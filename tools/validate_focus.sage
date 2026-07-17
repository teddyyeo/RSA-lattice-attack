# Focus-group VALIDATION: build the TRIMMED lattice (only the kept shift rows),
# flatter it, and confirm blind CRT-Groebner recovery still returns (x0*,x1*).
# Compares flatter wall time against the full-lattice baseline.
#   KEEPFILE=focus_keep_m18_t5.txt  (repr of sorted [(k,i),...]);  empty -> full.
import os, time, subprocess
X0_LSB = 600; X1_LSB = 265
X0_W = 230; X1_W = 155
X0_BOUND = Integer(1) << X0_W
X1_BOUND = Integer(1) << X1_W
SEED = int(os.environ.get("SEED", "12345"))
M    = int(os.environ.get("M", "18"))
T    = int(os.environ.get("T", "5"))
KEEPFILE = os.environ.get("KEEPFILE", "")
NPRIMES = int(os.environ.get("NPRIMES", "10"))
PBITS   = int(os.environ.get("PBITS", "30"))
NGEN    = int(os.environ.get("NGEN", "40"))
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

def build_matrix(N, base, m, t, keep):
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
    return Matrix(ZZ, rows), monomials, X, Y

def matrix_to_fplll(M):
    return "[" + "\n".join(
        "[" + " ".join(str(M[i, j]) for j in range(M.ncols())) + "]"
        for i in range(M.nrows())) + "]\n"

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
    keep = None
    if KEEPFILE and os.path.exists(KEEPFILE):
        keep = set(tuple(t) for t in eval(open(KEEPFILE).read()))
    tag = "full" if keep is None else "trim{}".format(len(keep))
    tb = time.time()
    M_, monomials, X, Y = build_matrix(N, base, M, T, keep)
    maxbits = max(abs(M_[i,j]).nbits() for i in range(M_.nrows()) for j in range(M_.ncols()))
    log("[build:{}] dim={}x{} maxentry={}bit ({:.1f}s)".format(
        tag, M_.nrows(), M_.ncols(), maxbits, time.time()-tb))
    data = matrix_to_fplll(M_).encode()
    tr = time.time()
    proc = subprocess.run(["flatter"], input=data,
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    red_s = time.time() - tr
    if proc.returncode != 0 or not proc.stdout.strip():
        log("[flatter:{}] FAILED after {:.1f}s".format(tag, red_s)); return
    L = Matrix(ZZ, [[Integer(int(z)) for z in ln.strip().strip("[]").split()]
                    for ln in proc.stdout.decode().splitlines()
                    if ln.strip().strip("[]").strip()])
    log("[flatter:{}] {:.1f}s ({:.3f} h)".format(tag, red_s, red_s/3600.0))
    trc = time.time()
    rec = crt_recover(L, monomials, X, Y, x0s, x1s)
    if rec is None:
        log("[recover:{}] FAILED".format(tag))
        log("[RESULT] tag={} flatter={:.1f}s SUCCESS=False".format(tag, red_s)); return
    x0r, x1r, np_ = rec
    p_rec = base + (x0r << X0_LSB) + (x1r << X1_LSB)
    ok = (x0r == x0s and x1r == x1s and N % p_rec == 0)
    log("[recover:{}] x0_ok={} x1_ok={} N%p==0={} primes={} ({:.1f}s)".format(
        tag, x0r==x0s, x1r==x1s, N%p_rec==0, np_, time.time()-trc))
    log("[RESULT] tag={} dim={} flatter={:.1f}s ({:.3f}h) SUCCESS={}".format(
        tag, M_.nrows(), red_s, red_s/3600.0, ok))

main()
