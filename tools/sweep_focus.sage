# Focus-group CLIFF sweep: peel the lattice by total-degree shell.
# For each Dcap in DCAPS, keep only shift rows with k+i <= Dcap (columns stay
# full m=18 monomial space), flatter, and test blind CRT recovery. Finds the
# smallest surviving shift-degree cap = the cliff.  SEED plants a known root.
import os, time, subprocess
X0_LSB = 600; X1_LSB = 265
X0_W = 230; X1_W = 155
X0_BOUND = Integer(1) << X0_W
X1_BOUND = Integer(1) << X1_W
SEED = int(os.environ.get("SEED", "12345"))
M    = int(os.environ.get("M", "18"))
T    = int(os.environ.get("T", "5"))
DCAPS = [int(x) for x in os.environ.get("DCAPS", "16,15,14,13").split(",")]
NPRIMES = 10; PBITS = 30; NGEN = 40
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

def build_matrix(N, base, m, t, dcap):
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
            if k + i > dcap:
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
    return Matrix(ZZ, rows), monomials, X, Y

def matrix_to_fplll(Mx):
    return "[" + "\n".join(
        "[" + " ".join(str(Mx[i, j]) for j in range(Mx.ncols())) + "]"
        for i in range(Mx.nrows())) + "]\n"

def crt_recover(L, monomials, X, Y, x0s, x1s):
    res0, res1, mods = [], [], []
    pr = Integer(1) << PBITS
    ncols = L.ncols()
    use = min(NGEN, L.nrows())
    rows = [L.row(r) for r in range(use)]
    tried = 0
    while len(mods) < NPRIMES and tried < 60:
        pr = previous_prime(pr); tried += 1
        try:
            FP = GF(pr)
            Rp = PolynomialRing(FP, 2, "x0,x1"); gx, gy = Rp.gens()
            invs = [FP(1)/(FP(X)**a * FP(Y)**b) for (a,b) in monomials]
            gens = []
            for row in rows:
                d = {monomials[c]: FP(int(row[c]))*invs[c] for c in range(ncols) if row[c]!=0}
                gens.append(Rp(d))
            V = Rp.ideal(gens).variety()
        except Exception:
            continue
        if len(V) != 1: continue
        pt = V[0]
        res0.append(ZZ(pt[gx])); res1.append(ZZ(pt[gy])); mods.append(pr)
    if not mods or prod(mods) < X0_BOUND:
        return None
    x0r = ZZ(crt(res0, mods)) % X0_BOUND
    x1r = ZZ(crt(res1, mods)) % X1_BOUND
    return (x0r, x1r, len(mods))

def main():
    N, p, base, x0s, x1s = plant_instance()
    log("[plant] seed={}".format(SEED))
    for dcap in DCAPS:
        tb = time.time()
        Mx, monomials, X, Y = build_matrix(N, base, M, T, dcap)
        log("[build Dcap={}] dim={}x{} ({:.1f}s)".format(dcap, Mx.nrows(), Mx.ncols(), time.time()-tb))
        data = matrix_to_fplll(Mx).encode()
        tr = time.time()
        proc = subprocess.run(["flatter"], input=data,
                              stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        red_s = time.time() - tr
        if proc.returncode != 0 or not proc.stdout.strip():
            log("[RESULT Dcap={}] dim={} flatter=FAIL SUCCESS=False".format(dcap, Mx.nrows())); continue
        L = Matrix(ZZ, [[Integer(int(z)) for z in ln.strip().strip("[]").split()]
                        for ln in proc.stdout.decode().splitlines()
                        if ln.strip().strip("[]").strip()])
        rec = crt_recover(L, monomials, X, Y, x0s, x1s)
        if rec is None:
            log("[RESULT Dcap={}] dim={} flatter={:.1f}s recover=NONE SUCCESS=False".format(
                dcap, Mx.nrows(), red_s)); continue
        x0r, x1r, np_ = rec
        p_rec = base + (x0r << X0_LSB) + (x1r << X1_LSB)
        ok = (x0r == x0s and x1r == x1s and N % p_rec == 0)
        log("[RESULT Dcap={}] dim={} flatter={:.1f}s ({:.3f}h) primes={} SUCCESS={}".format(
            dcap, Mx.nrows(), red_s, red_s/3600.0, np_, ok))

main()
