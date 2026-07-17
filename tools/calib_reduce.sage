# Calibration driver for the "reduce" path.
#
# Plants a KNOWN root with the real 2-variable geometry (230-bit x0 window at
# bit 600, 155-bit x1 window at bit 265), builds the HM lattice with off-diagonal
# coefficients reduced mod N^t (diagonal/pivot kept exact so the k=0 rows keep
# their N^t and the lattice stays full rank), reduces with flatter, and checks:
#   (1) flatter wall time,
#   (2) how many reduced rows actually VANISH at the planted root (ground truth),
#   (3) whether a resultant pair recovers the planted (x0*, x1*).
#
# REDUCE=0 -> arm.sage's original (un-reduced) construction, for A/B comparison.
# MT="24,7" chooses (m,t).  SEED chooses the planted instance.
import os, time, subprocess, signal

os.environ.setdefault("OMP_NUM_THREADS", "1")

X0_LSB = 600; X1_LSB = 265
X0_W = 230; X1_W = 155
X0_BOUND = Integer(1) << X0_W
X1_BOUND = Integer(1) << X1_W
REDUCE = os.environ.get("REDUCE", "1") == "1"
SEED   = int(os.environ.get("SEED", "12345"))

def log(s): print(s, flush=True)

def cmod(c, mod):
    r = c % mod
    if r > mod // 2: r -= mod
    return r

def plant_instance():
    set_random_seed(SEED)
    while True:
        p = random_prime(Integer(1) << 1024, lbound=Integer(1) << 1023)
        if p.nbits() == 1024: break
    while True:
        q = random_prime(Integer(1) << 1024, lbound=Integer(1) << 1023)
        if q.nbits() == 1024: break
    N = p * q
    x0s = (p >> X0_LSB) & (X0_BOUND - 1)
    x1s = (p >> X1_LSB) & (X1_BOUND - 1)
    base = p - (x0s << X0_LSB) - (x1s << X1_LSB)
    assert base + (x0s << X0_LSB) + (x1s << X1_LSB) == p
    return N, p, base, ZZ(x0s), ZZ(x1s)

def matrix_to_fplll(M):
    return "[" + "\n".join(
        "[" + " ".join(str(M[i, j]) for j in range(M.ncols())) + "]"
        for i in range(M.nrows())) + "]\n"

def build_matrix(N, base, m, t):
    R = PolynomialRing(ZZ, names=("x0", "x1")); x0, x1 = R.gens()
    a0 = Integer(1) << X0_LSB; a1 = Integer(1) << X1_LSB
    inv_a0 = inverse_mod(a0, N)
    A = ZZ((a1 * inv_a0) % N); B = ZZ((base * inv_a0) % N)
    f = x0 + A*x1 + B
    X = X0_BOUND; Y = X1_BOUND
    Nt = N ** t
    monomials = []
    for td in range(m + 1):
        for i in range(td + 1):
            monomials.append((i, td - i))
    midx = {mon: k for k, mon in enumerate(monomials)}
    rows = []
    for k in range(m + 1):
        for i in range(m - k + 1):
            g = (x1 ** i) * (f ** k) * (N ** max(t - k, 0))
            lead = (k, i)                         # leading monomial x0^k x1^i
            coeffs = [ZZ(0)] * len(monomials)
            for exp, c in g.dict().items():
                et = tuple(exp)
                if et not in midx: continue
                c = ZZ(c)
                if REDUCE and et != lead:
                    c = cmod(c, Nt)               # shrink off-diagonal only
                a, b = et
                coeffs[midx[et]] = c * (X ** a) * (Y ** b)
            rows.append(coeffs)
    return Matrix(ZZ, rows), monomials, R, X, Y

def reconstruct(L, monomials, R, X, Y, cap):
    x0, x1 = R.gens()
    RQ = R.change_ring(QQ)
    polys = []
    for r in range(min(cap, L.nrows())):
        h = RQ(0)
        for c in range(L.ncols()):
            v = L[r, c]
            if v == 0: continue
            i, j = monomials[c]
            h += QQ(v) / QQ((Integer(X) ** i) * (Integer(Y) ** j)) * (x0 ** i) * (x1 ** j)
        if h == 0: continue
        den = lcm([co.denominator() for co in h.coefficients()])
        polys.append(R(den * h))
    return polys

class _TO(Exception): pass
def _alarm(s, f): raise _TO()

def recover_root(polys, R, base, N, x0s, x1s):
    # BLIND recovery from the shortest reconstructed polys (mirrors the real run):
    # linear 2x2 if possible, else guarded resultant on the two lowest-degree polys.
    x0, x1 = R.gens()
    def height(h): return max((abs(c) for c in h.coefficients()), default=0)
    cand = sorted(polys, key=lambda h: (h.total_degree(), height(h)))[:12]
    # Strategy A: >=2 independent linear polys -> 2x2 solve
    lins = [h for h in cand if h.total_degree() == 1]
    for a in range(len(lins)):
        for b in range(a + 1, len(lins)):
            h1, h2 = lins[a], lins[b]
            a1, b1, c1 = h1.monomial_coefficient(x0), h1.monomial_coefficient(x1), h1.constant_coefficient()
            a2, b2, c2 = h2.monomial_coefficient(x0), h2.monomial_coefficient(x1), h2.constant_coefficient()
            det = a1 * b2 - a2 * b1
            if det == 0: continue
            xr = (-c1 * b2 + c2 * b1) / det
            yr = (-a1 * c2 + a2 * c1) / det
            if xr in ZZ and yr in ZZ:
                return (ZZ(xr), ZZ(yr), "linear")
    # Strategy B: guarded resultant on the two lowest-degree polys
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(120)
    try:
        for a in range(min(4, len(cand))):
            for b in range(a + 1, min(4, len(cand))):
                res = cand[a].resultant(cand[b], x1)
                if res == 0: continue
                d = {}; ok = True
                for e, c in res.dict().items():
                    ex, ey = tuple(e)
                    if ey != 0: ok = False; break
                    d[ex] = ZZ(c)
                if not (ok and d): continue
                for xr, _ in PolynomialRing(ZZ, "T")(d).roots():
                    xr = ZZ(xr)
                    hy = cand[a].subs({x0: xr}); dy = {}; okk = True
                    for e, c in hy.dict().items():
                        ex, ey = tuple(e)
                        if ex != 0: okk = False; break
                        dy[ey] = ZZ(c)
                    if okk and dy:
                        for yr, _ in PolynomialRing(ZZ, "T")(dy).roots():
                            return (xr, ZZ(yr), "resultant")
    except _TO:
        return ("timeout", None, "resultant")
    finally:
        signal.alarm(0)
    return None

def main():
    targets = []
    for it in os.environ.get("MT", "24,7").split(";"):
        a, b = it.split(","); targets.append((int(a), int(b)))
    N, p, base, x0s, x1s = plant_instance()
    log("[plant] seed={} N.nbits={} x0*.nbits={} x1*.nbits={}".format(
        SEED, N.nbits(), x0s.nbits(), x1s.nbits()))
    for (m, t) in targets:
        log("[cfg] m={} t={} REDUCE={}".format(m, t, REDUCE))
        tc = time.time()
        M, monomials, R, X, Y = build_matrix(N, base, m, t)
        maxbits = max(abs(M[i, j]).nbits() for i in range(M.nrows()) for j in range(M.ncols()))
        log("[build] dim={}x{} maxentry={}bit in {:.1f}s".format(
            M.nrows(), M.ncols(), maxbits, time.time() - tc))
        data = matrix_to_fplll(M).encode()
        tr = time.time()
        proc = subprocess.run(["flatter"], input=data,
                              stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        red_s = time.time() - tr
        if proc.returncode != 0 or not proc.stdout.strip():
            log("[flatter] FAILED after {:.1f}s".format(red_s)); continue
        L = Matrix(ZZ, [[Integer(int(z)) for z in ln.strip().strip("[]").split()]
                        for ln in proc.stdout.decode().splitlines()
                        if ln.strip().strip("[]").strip()])
        log("[flatter] {:.1f}s  ({:.2f} h)".format(red_s, red_s / 3600.0))
        # Save the reduced lattice so post-analysis never needs to re-reduce.
        lpath = "L_m{}_t{}_s{}_{}.txt".format(m, t, SEED, "red" if REDUCE else "raw")
        with open(lpath, "w") as fp:
            fp.write(matrix_to_fplll(L))
        log("[saved] {}".format(lpath))
        # Ground-truth check ONLY (no resultants -- those hung for 17h).
        # Success criterion (RESEARCH.md): >=2 INDEPENDENT rows vanish at the
        # true root. We know the root, so we just count + rank; no root solving.
        tg = time.time()
        x0g, x1g = R.gens()
        polys = reconstruct(L, monomials, R, X, Y, cap=120)
        vanish = [h for h in polys if h(x0s, x1s) == 0]
        rank = 0
        if vanish:
            vm = Matrix(ZZ, [[h.monomial_coefficient(x0g ** i * x1g ** j)
                              for (i, j) in monomials] for h in vanish])
            rank = vm.rank()
        log("[gt] reconstructed={} vanish_at_root={} indep_rank={} ({:.1f}s)".format(
            len(polys), len(vanish), rank, time.time() - tg))
        # End-to-end BLIND root recovery on the SAME in-process instance.
        tr2 = time.time()
        rec = recover_root(polys, R, base, N, x0s, x1s)
        recovered = False; method = "none"
        if rec is not None and rec[0] != "timeout":
            xr, yr, method = rec
            p_rec = base + (ZZ(xr) << X0_LSB) + (ZZ(yr) << X1_LSB)
            recovered = (xr == x0s and yr == x1s and N % p_rec == 0)
        elif rec is not None:
            method = "timeout"
        log("[recover] method={} recovered={} ({:.1f}s)".format(method, recovered, time.time() - tr2))
        log("[RESULT] m={} t={} REDUCE={} flatter={:.1f}s vanish={} indep={} RECOVERED={} SUCCESS={}".format(
            m, t, REDUCE, red_s, len(vanish), rank, recovered, recovered))

main()
