# Focus-group LEARNING (no flatter needed).
#
# We have B_red = saved flatter output for the SEED planted instance, and B_orig
# is deterministic. Recover the unimodular transform U with U*B_orig = B_red.
# A reduced row r is PRODUCTIVE if it vanishes at the planted root (x0*,x1*).
# An original row j is COLD if it is never used (U[r,j]==0) by ANY productive
# reduced row -> dropping it cannot break the existing vanishing relations.
# Print the cold set grouped by shift index (k,i) so we see the structure.
import os, time
X0_LSB = 600; X1_LSB = 265
X0_W = 230; X1_W = 155
X0_BOUND = Integer(1) << X0_W
X1_BOUND = Integer(1) << X1_W
SEED = int(os.environ.get("SEED", "12345"))
M    = int(os.environ.get("M", "18"))
T    = int(os.environ.get("T", "5"))
LFILE = os.environ.get("LFILE", "L_m18_t5_s12345_red.txt")
def log(s): print(s, flush=True)

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
    return N, p, base, ZZ(x0s), ZZ(x1s)

def cmod(c, mod):
    r = c % mod
    if r > mod // 2: r -= mod
    return r

def build_matrix(N, base, m, t, keep=None):
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
    shift_index = []          # row r -> (k,i)
    rows = []
    for k in range(m + 1):
        for i in range(m - k + 1):
            if keep is not None and (k, i) not in keep:
                continue
            g = (x1 ** i) * (f ** k) * (N ** max(t - k, 0))
            lead = (k, i)
            coeffs = [ZZ(0)] * len(monomials)
            for exp, c in g.dict().items():
                et = tuple(exp)
                if et not in midx: continue
                c = ZZ(c)
                if et != lead:
                    c = cmod(c, Nt)
                a, b = et
                coeffs[midx[et]] = c * (X ** a) * (Y ** b)
            rows.append(coeffs)
            shift_index.append((k, i))
    return Matrix(ZZ, rows), monomials, shift_index, X, Y

def load_matrix(path):
    rows = []
    for ln in open(path):
        ln = ln.strip().strip("[]").strip()
        if ln: rows.append([Integer(int(z)) for z in ln.split()])
    return Matrix(ZZ, rows)

def main():
    N, p, base, x0s, x1s = plant_instance()
    t0 = time.time()
    Bo, monomials, shift_index, X, Y = build_matrix(N, base, M, T)
    log("[orig] dim={}x{} ({:.1f}s)".format(Bo.nrows(), Bo.ncols(), time.time()-t0))
    Br = load_matrix(LFILE)
    log("[red ] dim={}x{} loaded".format(Br.nrows(), Br.ncols()))

    # which reduced rows vanish at the planted root
    pmon = [ZZ(x0s) ** a * ZZ(x1s) ** b for (a, b) in monomials]
    invden = []   # x0^i x1^j / (X^i Y^j) -> integer test: row . (pmon/scale)==0
    # reduced row R has R[c] = coeff * X^i Y^j ; poly value = sum R[c]*(x0s^i x1s^j)/(X^i Y^j)
    # multiply through by lcm not needed: test rational == 0 <=> sum R[c]*x0s^i*x1s^j / (X^i Y^j) == 0
    scale = [ (Integer(X) ** a) * (Integer(Y) ** b) for (a,b) in monomials ]
    def vanishes(row):
        acc = QQ(0)
        for c in range(len(row)):
            v = row[c]
            if v == 0: continue
            acc += QQ(v) * QQ(pmon[c]) / QQ(scale[c])
            # early bail impossible; sum small
        return acc == 0
    tv = time.time()
    productive = [r for r in range(Br.nrows()) if vanishes(Br.row(r))]
    log("[vanish] productive(reduced rows vanishing at root)={}/{} ({:.1f}s)".format(
        len(productive), Br.nrows(), time.time()-tv))

    # recover U: U * Bo = Br  -> U = Br * Bo^{-1}; solve_left
    tu = time.time()
    U = Bo.solve_left(Br)
    intg = all(U[i,j] in ZZ for i in range(U.nrows()) for j in range(U.ncols()))
    log("[U] {}x{} integral={} ({:.1f}s)".format(U.nrows(), U.ncols(), intg, time.time()-tu))

    # cold original row j: U[r,j]==0 for every productive reduced row r
    ncol = U.ncols()
    used = [False]*ncol
    contrib = [0]*ncol
    for r in productive:
        urow = U.row(r)
        for j in range(ncol):
            if urow[j] != 0:
                used[j] = True
                contrib[j] += 1
    cold = [j for j in range(ncol) if not used[j]]
    hot  = [j for j in range(ncol) if used[j]]
    log("[cold] original rows never used by any productive vector: {}/{}".format(len(cold), ncol))

    # group cold rows by shift index k
    from collections import defaultdict
    byk = defaultdict(list)
    for j in cold:
        k, i = shift_index[j]
        byk[k].append(i)
    for k in sorted(byk):
        log("   cold k={:2d}: i in {}".format(k, sorted(byk[k])))
    # contribution histogram: how many productive vectors use each hot row
    hotsorted = sorted(hot, key=lambda j: contrib[j])
    log("[hot ] least-used hot rows (row->(k,i)->#productive-uses):")
    for j in hotsorted[:15]:
        log("   row {:3d} (k={:2d},i={:2d}) used_by={}".format(j, shift_index[j][0], shift_index[j][1], contrib[j]))

    # write the conservative KEEP set (all hot rows) for the validation run
    keep = [shift_index[j] for j in hot]
    with open("focus_keep_m{}_t{}.txt".format(M,T), "w") as fp:
        fp.write(repr(sorted(keep)))
    log("[write] focus_keep_m{}_t{}.txt  keep={} drop={}".format(M, T, len(keep), len(cold)))

    # threshold-based (more aggressive) keep sets: drop rows used by <= THR
    # productive vectors.  THR=0 == conservative.  Saves each + a contrib dump.
    with open("focus_contrib_m{}_t{}.txt".format(M,T), "w") as fp:
        fp.write(repr([(shift_index[j], contrib[j]) for j in range(ncol)]))
    for thr in (1, 2, 3):
        keep_t = [shift_index[j] for j in range(ncol) if contrib[j] > thr]
        with open("focus_keep_m{}_t{}_thr{}.txt".format(M,T,thr), "w") as fp:
            fp.write(repr(sorted(keep_t)))
        log("[write] thr={} keep={} drop={}".format(thr, len(keep_t), ncol-len(keep_t)))

main()
