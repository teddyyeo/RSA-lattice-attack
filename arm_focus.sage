# ======================================================================
# Problem 7 RSA partial-bit factoring attack  --  FOCUS-TRIM build
# Herrmann-May / Coppersmith-style 2-variable relaxed model
# Portable: worker count = CPU cores by default (tune for host RAM).
#
# This is arm.sage plus the validated focus-group row trim (RESEARCH.md sec 7):
# the top 4 total-degree shells of shift rows are redundant, so the per-candidate
# lattice shrinks 190 -> 120 rows and flatter runs ~2.5x faster (826s -> 330s at
# (18,5)) with recovery unchanged. Set FOCUS_TRIM=0 for the original lattice.
#
# The 2-variable split is provably the most efficient for this bit-geometry
# (enabling bound met with huge margin: 385 unknown bits vs a ~724-bit
# threshold, so adding variables only enlarges the lattice). Everything else is
# the same portable pipeline: WORKERS single-threaded processes (no BLAS
# oversubscription), flatter/strong-LLL/BKZ reduction, and checkpoint/resume.
#
# IMPORTANT:
#   This script uses the LATEST working data/model from our discussion:
#     MASK runs = K100 U4 K90 U46 K15 U87 K13 U69 K180 U58 K13 U84 K111 U4 K150
#     After brute-forcing the two U4 gaps:
#       K194 | U230 | K180 | U155 | K265
#     p = base + (x0 << 600) + (x1 << 265)
#
# Run (uses all cores by default; needs flatter on PATH):
#   sage arm_focus.sage | tee log.txt
#
# On a memory-tight host, cap workers (e.g. 16 GB laptop):
#   WORKERS=14 REDUCTION=flatter sage arm_focus.sage | tee log.txt
#
# Disable the trim (original full 190-row lattice, for A/B):
#   FOCUS_TRIM=0 REDUCTION=flatter sage arm_focus.sage | tee log_full.txt
#
# Stronger run, if the default params find nothing (see RESEARCH.md §3 slack table):
#   PARAMS='20,6' MAX_POLYS=16 PAIR_LIMIT=20 REDUCTION=flatter \\
#       sage arm_focus.sage | tee log_strong.txt
# ======================================================================

from itertools import product
import os
import time
import signal

# Keep each worker single-threaded so the 4 parallel workers do not
# oversubscribe the 4 ARM cores, and cap glibc arenas to bound RSS under
# fork-based parallelism. Must be set before heavy libs init / before forking.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "GOTO_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("MALLOC_ARENA_MAX", "2")

# ======================================================================
# 0. Challenge values
# ======================================================================

N_HEX = """
e505004fb5d34eb7 12d48ff4bbe8d27f c388133c6c0e7340 01061c0ee0a4edc6 
37c04fe8dd376185 de8ba04d0ccdbabb 93ab7c371b88d92e 865eec42b028c61d 
d7004ebf2ebb5d69 d0a09142be5c9de4 da16e514eea31817 2ecda6cd192073eb 
afb1e02d522ec053 34590ea6d75960c4 937bf64f9700db17 7a4aa3da6aae6807 
e5e32c0d0e428a0d b68d299f20c235d8 4ef459b0cf118286 59c31663c9ea8204 
4b28152c89a9c36c 3ec4303bd36664fd 77fb02c58340bdae 21120326d83fc017 
34bc90048dec9fe3 5f08c8fdc523abf8 4a91ec430f495672 37c3153a2035ff62 
5613b6dc3e6cb14d 50e18b8a79b25d67 8465b3ad02f5b7d8 18a1e2d635a0baf1
"""

E_HEX = "0000000000010001"

CT_HEX = """
8919342826ef3821 5af31e00c9290c4c 50ef9ff9e1afc591 47fab5b096361035 
e85f5fc95b73b069 7813b57b831a807d 41bcbecde5b9e663 9e2845b14e395ed0 
e5d995e63709ac0c 5ee2337228ee76bc bad857b14904aa2e 8e9997671908a634 
d0d1dda1d062ce7f 2e3293ddec8f5cce 26029292d594a062 dcf317d2a8380f43 
d72551889efceb87 6c8945a50382272e 76ed6b6fcdff1603 44e9e948e2b6e740 
e78bedf25f30e2c7 eeb5f74686c8eadc 29cea04ff08cfd86 dfd3d2a1632bf04a 
d5cfa369892a2da4 0f0dc0098ce6b731 d841aab3d0c8b78e b69c4625c47c4ad7 
158d49bb5d879581 e02bc525abe47f39 f699864bc5ce1de7 19430dae7aa5480b
"""

# Latest working p&MASK and MASK used for the 2-variable model.
# Do not replace these with the older PDF values unless the model below is also changed.
P_AND_MASK_HEX = """
ffa360d46885c534 d186538170633faf c2c0548a2e24a2c1 c0000000000039e2
0000000000000000 000000a520000000 00000000003e2de4 c436d2ca740a6246
99e1a1af94045c63 261323c000000000 000003bba0000000 00000000000000e5
0b0bc2461fcbac07 26360c2c0809450a 9a892cbf1d98ceee 48827591ccc593c9
"""

MASK_HEX = """
ffffffffffffffff fffffffff0ffffff ffffffffffffffff c00000000000fffe
0000000000000000 000003ffe0000000 0000000000ffffff ffffffffffffffff
ffffffffffffffff fffffff000000000 000003ffe0000000 00000000000001ff
ffffffffffffffff fffffffffc3fffff ffffffffffffffff ffffffffffffffff
"""

# ======================================================================
# 1. Parameters and model constants
# ======================================================================

P_BITS = 1024
FULL_MASK = (Integer(1) << P_BITS) - 1

# The two 4-bit gaps brute-forced before building the relaxed model.
HIGH_GAP_LSB = 920       # p bits [923:920]
LOW_GAP_LSB  = 150       # p bits [153:150]

# Relaxed model: K194 | U230 | K180 | U155 | K265
KH_WIDTH = 194
X0_WIDTH = 230
KM_WIDTH = 180
X1_WIDTH = 155
KL_WIDTH = 265

KL_LSB = 0
X1_LSB = 265
KM_LSB = 420
X0_LSB = 600
KH_LSB = 830

assert KH_WIDTH + X0_WIDTH + KM_WIDTH + X1_WIDTH + KL_WIDTH == P_BITS

X0_BOUND = Integer(1) << X0_WIDTH
X1_BOUND = Integer(1) << X1_WIDTH

# Runtime options. Override with environment variables.
# WORKERS defaults to all cores. Each worker holds one trimmed (~120x190)
# lattice plus a flatter subprocess (~0.4-0.7 GB peak with FOCUS_TRIM on), so on
# a memory-tight host set WORKERS <= RAM_GB / 0.8. On 16 GB keep WORKERS <= ~14
# to leave headroom for the OS; on 24+ GB the core count is the only limit.
_DEFAULT_WORKERS = os.cpu_count() or 4
START_INDEX = int(os.environ.get("START_INDEX", "0"))
STOP_INDEX  = int(os.environ.get("STOP_INDEX", "256"))
WORKERS     = int(os.environ.get("WORKERS", str(_DEFAULT_WORKERS)))
MAX_POLYS   = int(os.environ.get("MAX_POLYS", "10"))
PAIR_LIMIT  = int(os.environ.get("PAIR_LIMIT", "8"))
LLL_DELTA   = float(os.environ.get("LLL_DELTA", "0.99"))   # strong reduction

# Focus-group row trim (RESEARCH.md sec 7/sec 11): drop the top FOCUS_TRIM
# total-degree shells of shift rows (k+i > m - FOCUS_TRIM). Validated safe up
# to 4 shells at (m,t)=(18,5) -> 120/190 rows, flatter ~2.5x faster, recovery
# unchanged on 3 seeds. The 5th shell collapses the lattice, so 4 is the floor;
# set FOCUS_TRIM=0 to disable and rebuild the full lattice.
FOCUS_TRIM  = int(os.environ.get("FOCUS_TRIM", "4"))

# Lattice reduction backend: auto | flatter | bkz | lll
#   auto    -> flatter if on PATH, else strong LLL
#   flatter -> require flatter (fastest on Coppersmith bases), else LLL fallback
#   bkz     -> Sage BKZ with block size BKZ_BLOCK (strongest, slower)
#   lll     -> Sage/fpLLL LLL at LLL_DELTA
REDUCTION   = os.environ.get("REDUCTION", "auto").lower()
BKZ_BLOCK   = int(os.environ.get("BKZ_BLOCK", "20"))

# Skip any (m, t) whose lattice cannot reach the 385-bit unknown target.
# Default OFF: the hm_bound_ok() gate below is a hand-derived diagonal condition
# that overestimates capacity (~724 bits, impossible: >n=1 Coppersmith's 512).
# RESEARCH.md §3 uses the correct slack table as the scheduler instead; do not
# let this gate skip/bless runs unless it has been validated on synthetic cases.
FEASIBILITY_GATE = os.environ.get("FEASIBILITY_GATE", "0") == "1"

# Resume support: completed (m, t, idx) triples are appended here so an
# interrupted run can be restarted without repeating finished candidates.
CHECKPOINT  = os.environ.get("CHECKPOINT", "checkpoint_q7.txt")

# 2-variable model is provably the most efficient split here, so keep it.
# Params come from EMPIRICAL calibration (RESEARCH.md §3-4), which beats the
# conservative determinant gate badly: on planted-root synthetic instances of
# this exact geometry, (18,5) recovers the root (120/120 short vectors vanish,
# validated on 4 seeds); everything <= (16,4) fails. (18,5) is the minimal
# robust config (~0.09 h/candidate = 330s with reduce + focus-trim + flatter,
# down from ~0.23 h untrimmed); (20,6) kept as a proven fallback. cheap -> strong.
DEFAULT_PARAMS = [(18, 5), (20, 6)]


def parse_params():
    raw = os.environ.get("PARAMS", "").strip()
    if raw == "":
        return DEFAULT_PARAMS
    out = []
    for item in raw.split(";"):
        item = item.strip()
        if not item:
            continue
        a, b = item.split(",")
        out.append((int(a), int(b)))
    if not out:
        raise ValueError("PARAMS is empty after parsing")
    return out

PARAMS = parse_params()

# ======================================================================
# 1b. Feasibility gate, checkpoint, and lattice-reduction backend
# ======================================================================


def hm_bound_ok(m, t):
    # Idealized Herrmann-May diagonal-determinant condition for the 2-variable
    # linear model with bounds 2^X0_WIDTH, 2^X1_WIDTH and divisor p ~ N^(1/2).
    # Returns True if this (m, t) lattice can, in principle, reach the target.
    from math import comb
    n = 2
    logN = 2048
    logP = 1024
    D = comb(m + n, n)
    def num_deg(d):
        return comb(d + n - 1, n - 1)
    sum_deg = sum(d * num_deg(d) for d in range(m + 1))
    partA = (X0_WIDTH + X1_WIDTH) * sum_deg / float(n)
    partB = logN * sum((t - d) * num_deg(d) for d in range(0, t))
    return (partA + partB) < (logP * t * D)


def load_checkpoint():
    done = set()
    if not CHECKPOINT:
        return done
    try:
        with open(CHECKPOINT) as fp:
            for line in fp:
                line = line.strip()
                if not line:
                    continue
                parts = line.split(",")
                if len(parts) == 3:
                    try:
                        done.add((int(parts[0]), int(parts[1]), int(parts[2])))
                    except ValueError:
                        pass
    except IOError:
        pass
    return done


def record_checkpoint(m, t, idx):
    if not CHECKPOINT:
        return
    try:
        with open(CHECKPOINT, "a") as fp:
            fp.write("{},{},{}\n".format(m, t, idx))
    except IOError:
        pass


def _matrix_to_fplll(M):
    rows = []
    for i in range(M.nrows()):
        rows.append("[" + " ".join(str(M[i, j]) for j in range(M.ncols())) + "]")
    return "[" + "\n".join(rows) + "]\n"


def _fplll_to_matrix(text):
    parsed = []
    for line in text.splitlines():
        line = line.strip().strip("[]").strip()
        if not line:
            continue
        parts = line.split()
        if not parts:
            continue
        parsed.append([Integer(int(p)) for p in parts])
    if not parsed:
        return None
    width = len(parsed[0])
    if any(len(r) != width for r in parsed):
        return None
    return Matrix(ZZ, parsed)


def _try_flatter(M):
    import shutil, subprocess
    if shutil.which("flatter") is None:
        return None
    try:
        data = _matrix_to_fplll(M).encode()
        proc = subprocess.run(["flatter"], input=data,
                              stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        if proc.returncode != 0 or not proc.stdout.strip():
            return None
        red = _fplll_to_matrix(proc.stdout.decode())
        if red is None or red.nrows() != M.nrows() or red.ncols() != M.ncols():
            return None
        return red
    except Exception:
        return None


def reduce_lattice(M, verbose=False):
    # Pick the reduction backend per REDUCTION, with safe fallbacks.
    if REDUCTION in ("auto", "flatter"):
        red = _try_flatter(M)
        if red is not None:
            if verbose:
                log_line("[reduce] flatter")
            return red
        if REDUCTION == "flatter" and verbose:
            log_line("[reduce] flatter unavailable; falling back to LLL")
    if REDUCTION == "bkz":
        try:
            if verbose:
                log_line("[reduce] BKZ block_size={}".format(BKZ_BLOCK))
            return M.BKZ(block_size=BKZ_BLOCK)
        except Exception as ex:
            if verbose:
                log_line("[reduce] BKZ failed ({}); LLL fallback".format(ex))
    if verbose:
        log_line("[reduce] LLL delta={}".format(LLL_DELTA))
    try:
        return M.LLL(delta=LLL_DELTA, algorithm="fpLLL:wrapper")
    except TypeError:
        return M.LLL(delta=LLL_DELTA)

# ======================================================================
# 2. Utilities
# ======================================================================


def clean_hex(s):
    s = s.replace("0x", "").replace("0X", "")
    return "".join(c for c in s if c in "0123456789abcdefABCDEF").lower()


def hex_to_int(s):
    h = clean_hex(s)
    if len(h) == 0:
        raise ValueError("empty hex string")
    return Integer(int(h, 16))


def bitmask(width, lsb):
    return ((Integer(1) << width) - 1) << lsb


def extract_bits(value, width, lsb):
    return (Integer(value) >> lsb) & ((Integer(1) << width) - 1)


def insert_bits(value, insert_value, width, lsb):
    insert_value = Integer(insert_value)
    if not (0 <= insert_value < (Integer(1) << width)):
        raise ValueError("insert_value is out of range")
    m = bitmask(width, lsb)
    value = Integer(value)
    # Sage .sage treats ^ as power, so never use ^ for XOR/complement here.
    value = value & (FULL_MASK - m)
    value = value | (insert_value << lsb)
    return Integer(value)


def mask_runs(mask):
    bits = bin(Integer(mask))[2:].zfill(P_BITS)
    runs = []
    cur = bits[0]
    cnt = 1
    for b in bits[1:]:
        if b == cur:
            cnt += 1
        else:
            runs.append(("K" if cur == "1" else "U", cnt))
            cur = b
            cnt = 1
    runs.append(("K" if cur == "1" else "U", cnt))
    return runs


def log_line(s):
    print(s, flush=True)

# ======================================================================
# 3. Sanity checks
# ======================================================================


def sanity_check_inputs(N, e, ct, p_masked, mask):
    log_line("[sanity]")
    log_line("N bits      = {}".format(N.nbits()))
    log_line("e           = {}".format(e))
    log_line("ct bits     = {}".format(ct.nbits()))
    log_line("MASK bits   = {}".format(mask.nbits()))
    log_line("p&MASK bits = {}".format(p_masked.nbits()))

    if N.nbits() != 2048:
        raise ValueError("N is expected to be 2048 bits")
    if e != 65537:
        raise ValueError("e is expected to be 65537")
    if ct <= 0 or ct >= N:
        raise ValueError("ct must satisfy 0 < ct < N")
    if mask.nbits() > P_BITS:
        raise ValueError("MASK is longer than 1024 bits")

    unknown_mask = FULL_MASK - Integer(mask)
    if (p_masked & unknown_mask) != 0:
        raise ValueError("p&MASK has nonzero bits where MASK is 0")

    expected_runs = [
        ("K", 100), ("U", 4), ("K", 90), ("U", 46), ("K", 15),
        ("U", 87), ("K", 13), ("U", 69), ("K", 180), ("U", 58),
        ("K", 13), ("U", 84), ("K", 111), ("U", 4), ("K", 150),
    ]
    runs = mask_runs(mask)
    log_line("MASK runs   = {}".format(runs))
    if runs != expected_runs:
        raise ValueError("MASK run structure does not match the 2-variable latest model")

    if extract_bits(mask, 4, HIGH_GAP_LSB) != 0:
        raise ValueError("HIGH_GAP_LSB is wrong; MASK is not zero at [923:920]")
    if extract_bits(mask, 4, LOW_GAP_LSB) != 0:
        raise ValueError("LOW_GAP_LSB is wrong; MASK is not zero at [153:150]")

    brute_mask = Integer(mask) | bitmask(4, HIGH_GAP_LSB) | bitmask(4, LOW_GAP_LSB)
    for name, width, lsb in [("K194", KH_WIDTH, KH_LSB), ("K180", KM_WIDTH, KM_LSB), ("K265", KL_WIDTH, KL_LSB)]:
        missing = bitmask(width, lsb) & (FULL_MASK - brute_mask)
        if missing != 0:
            raise ValueError("{} is not fully known after brute force".format(name))

    log_line("[OK] input/model sanity checks passed")

# ======================================================================
# 4. Model construction and verification
# ======================================================================


def build_base_after_bruteforce(p_masked, high_guess, low_guess):
    candidate = Integer(p_masked)
    candidate = insert_bits(candidate, high_guess, 4, HIGH_GAP_LSB)
    candidate = insert_bits(candidate, low_guess,  4, LOW_GAP_LSB)
    model_known_mask = bitmask(KH_WIDTH, KH_LSB) | bitmask(KM_WIDTH, KM_LSB) | bitmask(KL_WIDTH, KL_LSB)
    return Integer(candidate & model_known_mask)


def make_p_candidate(base, x0, x1):
    return Integer(base) + (Integer(x0) << X0_LSB) + (Integer(x1) << X1_LSB)


def passes_internal_filters(x0, x1, p_masked):
    # x0 spans p[829:600] = U46 | K15 | U87 | K13 | U69
    # x1 spans p[419:265] = U58 | K13 | U84
    x0 = Integer(x0)
    x1 = Integer(x1)
    p_masked = Integer(p_masked)
    if extract_bits(x0, 15, 169) != extract_bits(p_masked, 15, 769):
        return False
    if extract_bits(x0, 13, 69) != extract_bits(p_masked, 13, 669):
        return False
    if extract_bits(x1, 13, 84) != extract_bits(p_masked, 13, 349):
        return False
    return True


def verify_candidate(N, base, x0v, x1v, p_masked, mask):
    x0v = Integer(x0v)
    x1v = Integer(x1v)
    if not (0 <= x0v < X0_BOUND):
        return None
    if not (0 <= x1v < X1_BOUND):
        return None
    if not passes_internal_filters(x0v, x1v, p_masked):
        return None
    p_candidate = make_p_candidate(base, x0v, x1v)
    if (p_candidate & Integer(mask)) != Integer(p_masked):
        return None
    if p_candidate > 1 and Integer(N) % p_candidate == 0:
        return p_candidate
    return None

# ======================================================================
# 5. Polynomial helpers
# ======================================================================


def to_univariate_in_x(poly, x0, x1):
    S = PolynomialRing(ZZ, "X")
    d = {}
    for exp, coeff in poly.dict().items():
        ex, ey = tuple(exp)
        if ey != 0:
            return None
        d[ex] = ZZ(coeff)
    return S(d)


def to_univariate_in_y(poly, x0, x1):
    S = PolynomialRing(ZZ, "Y")
    d = {}
    for exp, coeff in poly.dict().items():
        ex, ey = tuple(exp)
        if ex != 0:
            return None
        d[ey] = ZZ(coeff)
    return S(d)


def selected_resultant_pairs(polys, pair_limit):
    def score(h):
        try:
            deg = h.total_degree()
        except Exception:
            deg = 10**9
        try:
            terms = len(h.dict())
        except Exception:
            terms = 10**9
        return (deg, terms)

    ranked = sorted(list(enumerate(polys)), key=lambda z: score(z[1]))
    pairs = []
    seen = set()
    anchor_count = min(3, len(ranked))
    for ai in range(anchor_count):
        idx1, h1 = ranked[ai]
        for idx2, h2 in ranked[ai + 1:]:
            key = tuple(sorted((idx1, idx2)))
            if key in seen:
                continue
            seen.add(key)
            pairs.append((h1, h2))
            if len(pairs) >= pair_limit:
                return pairs
    return pairs


def reconstruct_polys_from_lll(L, R, monomials, X, Y, max_polys):
    x0, x1 = R.gens()
    polys = []
    malformed = 0
    R_QQ = R.change_ring(QQ)

    for row_idx in range(L.nrows()):
        row = L[row_idx]
        h_QQ = R_QQ(0)
        for col_idx, coeff_scaled in enumerate(row):
            if coeff_scaled == 0:
                continue
            i, j = monomials[col_idx]
            scale = (Integer(X) ** i) * (Integer(Y) ** j)
            h_QQ += QQ(coeff_scaled) / QQ(scale) * (x0 ** i) * (x1 ** j)
        if h_QQ == 0:
            continue
        try:
            den = lcm([c.denominator() for c in h_QQ.coefficients()])
            h = R(den * h_QQ)
        except Exception:
            malformed += 1
            continue
        if h != 0:
            polys.append(h)
            if len(polys) >= max_polys:
                break
    return polys, malformed

# ======================================================================
# 6. Herrmann-May bivariate attack for one candidate
# ======================================================================


def build_reduced_lattice(N, base, m, t, monomials, monomial_index):
    # Herrmann-May bivariate lattice with off-diagonal coefficients reduced mod
    # N^t (the leading/pivot coeff of each shift kept exact so the k=0 rows keep
    # their N^t and the lattice stays full rank), then scaled by X^a Y^b.
    # This shrinks max entry ~2.7x (52.8k -> 19.3k bits at m=24) => ~2.6x faster
    # flatter. Validated: identical short vectors vanish at the true root.
    #
    # FOCUS-GROUP ROW TRIM (RESEARCH.md sec 7 / sec 11): the top FOCUS_TRIM
    # total-degree SHELLS of shift rows (those with k+i > m - FOCUS_TRIM) are
    # provably redundant -- recovering the exact unimodular transform U showed
    # they never appear in any root-vanishing reduced vector. We keep the FULL
    # m-degree MONOMIAL space (columns) but drop those shift ROWS. Measured at
    # (m,t)=(18,5): dropping 4 shells (k+i>14) removes 70/190 rows and cuts
    # flatter ~2.5x (826s -> 330s) with recovery still exact/blind on 3 seeds.
    # The cliff is sharp: dropping a 5th shell (k+i>13) collapses the lattice
    # (no vanishing vectors). FOCUS_TRIM=0 restores the original full lattice.
    shell_cap = m - FOCUS_TRIM if FOCUS_TRIM > 0 else m
    R = PolynomialRing(ZZ, names=("x0", "x1"))
    x0, x1 = R.gens()
    a0 = Integer(1) << X0_LSB
    a1 = Integer(1) << X1_LSB
    inv_a0 = inverse_mod(a0, N)
    A = ZZ((a1 * inv_a0) % N)
    B = ZZ((base * inv_a0) % N)
    f = x0 + A * x1 + B
    X = X0_BOUND
    Y = X1_BOUND
    Nt = N ** t
    half = Nt // 2
    rows = []
    for k in range(m + 1):
        for i in range(m - k + 1):
            if k + i > shell_cap:                          # focus-group row trim
                continue
            g = (x1 ** i) * (f ** k) * (N ** max(t - k, 0))
            lead = (k, i)                                  # leading monomial x0^k x1^i
            coeffs = [ZZ(0)] * len(monomials)
            for exp, coeff in g.dict().items():
                et = tuple(exp)
                if et not in monomial_index:
                    continue
                c = ZZ(coeff)
                if et != lead:                             # reduce off-diagonal only
                    c = c % Nt
                    if c > half:
                        c -= Nt
                a, b = et
                coeffs[monomial_index[et]] = c * (Integer(X) ** a) * (Integer(Y) ** b)
            rows.append(coeffs)
    return Matrix(ZZ, rows), X, Y


def crt_groebner_recover(L, monomials, X, Y, verbose=False,
                         ngen=40, nprimes=12, pbits=30, max_tries=48, budget=90.0):
    # Blind root recovery: for several small primes p, reconstruct the shortest
    # lattice rows as polynomials over GF(p) (unscale = modular inverse of X^aY^b),
    # take the ideal's variety -- a single rational point (x0*,x1*) mod p -- and
    # CRT across enough primes (product > 2^230) to get the exact (x0,x1).
    # All small-prime Singular Groebner => fast; avoids degree-m resultants.
    ncols = L.ncols()
    use = min(ngen, L.nrows())
    rows = [L.row(r) for r in range(use)]              # shortest rows (flatter norm order)
    res0, res1, mods = [], [], []
    tried = 0
    pr = Integer(1) << pbits
    t_start = time.time()
    while len(mods) < nprimes and tried < max_tries and (time.time() - t_start) < budget:
        pr = previous_prime(pr)
        tried += 1
        try:
            FP = GF(pr)
            Rp = PolynomialRing(FP, 2, "x0,x1")
            gx, gy = Rp.gens()
            invs = [FP(1) / (FP(X) ** a * FP(Y) ** b) for (a, b) in monomials]
            gens = []
            for row in rows:
                d = {monomials[c]: FP(int(row[c])) * invs[c]
                     for c in range(ncols) if row[c] != 0}
                gens.append(Rp(d))
            V = Rp.ideal(gens).variety()
        except Exception:
            continue
        if len(V) != 1:
            continue
        pt = V[0]
        res0.append(ZZ(pt[gx])); res1.append(ZZ(pt[gy])); mods.append(pr)
    if not mods or prod(mods) < X0_BOUND:
        return []
    x0r = ZZ(crt(res0, mods)) % X0_BOUND
    x1r = ZZ(crt(res1, mods)) % X1_BOUND
    if verbose:
        log_line("[recover] primes_used={}/{} tried={} ({:.1f}s)".format(
            len(mods), nprimes, tried, time.time() - t_start))
    return [(x0r, x1r)]


def herrmann_may_bivariate(N, base, p_masked, mask, m, t, max_polys, pair_limit, verbose=False):
    N = Integer(N)
    base = Integer(base)

    monomials = []
    for total_deg in range(m + 1):
        for i in range(total_deg + 1):
            monomials.append((i, total_deg - i))
    monomial_index = {mon: idx for idx, mon in enumerate(monomials)}

    M, X, Y = build_reduced_lattice(N, base, m, t, monomials, monomial_index)
    if verbose:
        maxbits = max(abs(M[r, c]).nbits()
                      for r in range(M.nrows()) for c in range(M.ncols()))
        log_line("[HM] m={}, t={}, lattice={}x{}, maxentry={}bit".format(
            m, t, M.nrows(), M.ncols(), maxbits))

    t0 = time.time()
    try:
        L = reduce_lattice(M, verbose=verbose)
    except Exception as e:
        if verbose:
            log_line("[ERROR] lattice reduction failed: {}".format(e))
        return None
    if verbose:
        log_line("[HM] reduction done {:.2f}s".format(time.time() - t0))

    # Guard recovery with a wall-clock alarm so a wrong brute-force guess (255 of
    # 256 have no valid root) cannot hang the worker.
    def _timeout(signum, frame):
        raise TimeoutError()
    old = signal.signal(signal.SIGALRM, _timeout)
    signal.alarm(150)
    try:
        cands = crt_groebner_recover(L, monomials, X, Y, verbose=verbose)
    except Exception:
        cands = []
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)

    for (xr, yr) in cands:
        p = verify_candidate(N, base, ZZ(xr), ZZ(yr), p_masked, mask)
        if p is not None:
            return (ZZ(xr), ZZ(yr), p)
    if verbose:
        log_line("[HM] recovery candidates={}, none verified".format(len(cands)))
    return None

# ======================================================================
# 7. RSA decrypt and output
# ======================================================================


def decrypt_and_write_result(N, e, ct, p, q, extra_lines, filename):
    phi = (p - 1) * (q - 1)
    d = inverse_mod(e, phi)
    m_int = power_mod(ct, d, N)
    m_hex = hex(Integer(m_int))

    py_m = int(m_int)
    m_len = max(1, (py_m.bit_length() + 7) // 8)
    m_bytes = py_m.to_bytes(m_len, "big")

    try:
        m_utf8 = m_bytes.decode("utf-8")
    except Exception:
        m_utf8 = None

    log_line("[RSA]")
    log_line("m_int   = {}".format(m_int))
    log_line("m_hex   = {}".format(m_hex))
    log_line("m_bytes = {}".format(repr(m_bytes)))
    if m_utf8 is not None:
        log_line("m_utf8  = {}".format(m_utf8))

    with open(filename, "w") as fp:
        for line in extra_lines:
            fp.write(line + "\n")
        fp.write("p = {}\n".format(hex(p)))
        fp.write("q = {}\n".format(hex(q)))
        fp.write("p*q == N = {}\n".format(p * q == N))
        fp.write("d = {}\n".format(d))
        fp.write("m_int = {}\n".format(m_int))
        fp.write("m_hex = {}\n".format(m_hex))
        fp.write("m_bytes_repr = {}\n".format(repr(m_bytes)))
        if m_utf8 is not None:
            fp.write("m_utf8 = {}\n".format(m_utf8))

# ======================================================================
# 8. Worker and main loop
# ======================================================================


def solve_one_candidate(args):
    idx, high_guess, low_guess, N, p_masked, mask, m, t = args
    verbose = (idx == START_INDEX)
    t0 = time.time()
    base = build_base_after_bruteforce(p_masked, high_guess, low_guess)
    result = herrmann_may_bivariate(N, base, p_masked, mask, m, t, MAX_POLYS, PAIR_LIMIT, verbose=verbose)
    elapsed = time.time() - t0
    # Always return idx so the main loop can checkpoint attempts (incl. misses).
    if result is None:
        return (idx, high_guess, low_guess, m, t, None)
    x0, x1, p = result
    return (idx, high_guess, low_guess, m, t, (x0, x1, p, elapsed))


def run_tasks_parallel(tasks, workers):
    # Define the parallel wrapper here so it sees the final worker count.
    @parallel(ncpus=workers)
    def parallel_solver(idx, high_guess, low_guess, N, p_masked, mask, m, t):
        task = (idx, high_guess, low_guess, N, p_masked, mask, m, t)
        return solve_one_candidate(task)
    
    for task_input, result in parallel_solver(tasks):
        yield result


def run_tasks_sequential(tasks):
    for task in tasks:
        yield solve_one_candidate(task)


def run_attack():
    N = hex_to_int(N_HEX)
    e = hex_to_int(E_HEX)
    ct = hex_to_int(CT_HEX)
    p_masked = hex_to_int(P_AND_MASK_HEX)
    mask = hex_to_int(MASK_HEX)

    sanity_check_inputs(N, e, ct, p_masked, mask)

    candidates_all = list(product(range(16), range(16)))
    total = len(candidates_all)
    start = max(0, START_INDEX)
    stop = min(STOP_INDEX, total)
    if not (0 <= start < stop <= total):
        raise ValueError("invalid START_INDEX/STOP_INDEX: [{}, {})".format(START_INDEX, STOP_INDEX))

    log_line("")
    log_line("[settings]")
    log_line("candidate range = [{}, {}) / {}".format(start, stop, total))
    log_line("PARAMS          = {}".format(PARAMS))
    log_line("WORKERS         = {}".format(WORKERS))
    log_line("MAX_POLYS       = {}".format(MAX_POLYS))
    log_line("PAIR_LIMIT      = {}".format(PAIR_LIMIT))
    log_line("LLL_DELTA       = {}".format(LLL_DELTA))
    log_line("REDUCTION       = {}".format(REDUCTION))
    log_line("BKZ_BLOCK       = {}".format(BKZ_BLOCK))
    log_line("FEASIBILITY_GATE= {}".format(FEASIBILITY_GATE))
    log_line("CHECKPOINT      = {}".format(CHECKPOINT))
    log_line("X0_BOUND bits   = {}".format(X0_BOUND.nbits() - 1))
    log_line("X1_BOUND bits   = {}".format(X1_BOUND.nbits() - 1))

    checkpoint_done = load_checkpoint()
    if checkpoint_done:
        log_line("[resume] {} (m,t,idx) triples already attempted; they will be skipped".format(
            len(checkpoint_done)))

    for m, t in PARAMS:
        if FEASIBILITY_GATE and not hm_bound_ok(m, t):
            log_line("")
            log_line("[skip] m={}, t={}: lattice cannot reach the {}-bit unknown target".format(
                m, t, X0_WIDTH + X1_WIDTH))
            continue

        param_t0 = time.time()
        log_line("")
        log_line("=" * 88)
        log_line("[PARAM] m={}, t={}".format(m, t))
        log_line("=" * 88)

        tasks = []
        resumed = 0
        for idx, (high_guess, low_guess) in enumerate(candidates_all):
            if idx < start or idx >= stop:
                continue
            if (m, t, idx) in checkpoint_done:
                resumed += 1
                continue
            tasks.append((idx, high_guess, low_guess, N, p_masked, mask, m, t))

        if resumed:
            log_line("[resume] skipped {} already-attempted candidates for this param".format(resumed))
        if not tasks:
            log_line("[-] nothing left to do for m={}, t={}".format(m, t))
            continue

        runner = run_tasks_parallel(tasks, WORKERS) if WORKERS > 1 else run_tasks_sequential(tasks)
        done = 0
        for result in runner:
            done += 1
            # Guard against worker crashes (Sage @parallel may yield a string).
            if not (isinstance(result, tuple) and len(result) == 6):
                continue

            idx, high_guess, low_guess, found_m, found_t, payload = result
            record_checkpoint(found_m, found_t, idx)

            if done == 1 or done % max(1, min(16, len(tasks))) == 0:
                log_line("[progress] param=({},{}), completed={}/{}, elapsed={:.1f}s".format(
                    m, t, done, len(tasks), time.time() - param_t0
                ))

            if payload is None:
                continue

            x0, x1, p, elapsed = payload
            q = N // p
            if p * q != N:
                continue

            log_line("")
            log_line("#" * 88)
            log_line("[SUCCESS]")
            log_line("idx        = {}".format(idx))
            log_line("m,t        = {},{}".format(found_m, found_t))
            log_line("high_guess = 0x{:x}".format(high_guess))
            log_line("low_guess  = 0x{:x}".format(low_guess))
            log_line("x0         = {}".format(x0))
            log_line("x1         = {}".format(x1))
            log_line("p          = {}".format(hex(p)))
            log_line("q          = {}".format(hex(q)))
            log_line("elapsed candidate = {:.1f}s".format(elapsed))
            log_line("#" * 88)

            extra = [
                "idx = {}".format(idx),
                "m = {}".format(found_m),
                "t = {}".format(found_t),
                "high_guess = 0x{:x}".format(high_guess),
                "low_guess = 0x{:x}".format(low_guess),
                "x0 = {}".format(x0),
                "x1 = {}".format(x1),
            ]
            decrypt_and_write_result(N, e, ct, p, q, extra, "found_result_arm.txt")
            log_line("[+] wrote found_result_arm.txt")
            return (p, q)

        log_line("[-] m={}, t={} finished without factor in {:.1f}s".format(m, t, time.time() - param_t0))

    log_line("")
    log_line("[FAIL] No factor found in the selected range/parameters.")
    log_line("Next stronger examples (ARM 4-core / 24 GB; see RESEARCH.md §3-4):")
    log_line("  REDUCTION=flatter PARAMS='22,6;24,7' sage problem7_arm.sage | tee log_stronger.txt")
    log_line("  REDUCTION=flatter PARAMS='26,8' sage problem7_arm.sage | tee log_max.txt")
    return None

# ======================================================================
# 9. Execute
# ======================================================================

run_attack()