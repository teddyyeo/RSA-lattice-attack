# RSA-lattice-attack

Optimized lattice attack against a partially revealed RSA prime — a
Herrmann–May / multivariate-Coppersmith attack that recovers a 1024-bit prime
`p` from 672 of its bits (352 unknown, in 7 discontinuous blocks) and decrypts
textbook RSA. Built for and solved on the 2026 Cryptanalysis Competition Q7.

**Result:** `N` factored, `ct` decrypted →
`FLAG{d1rty_b1t_l34k_c0pp3rsm1th_m33ts_str4t3gy}` (see
[`found_result_arm.txt`](found_result_arm.txt)). ~14 h on a 4-core ARM box.

## The attack in one paragraph

Brute-force the two 4-bit unknown gaps (256 combinations). For each, model the
remaining unknowns as a **2-variable** linear equation modulo the unknown
divisor `p ≈ √N` — the minimal viable split (the 565-bit post-brute-force span
exceeds the 512-bit single-block Coppersmith limit, so `n ≥ 2`; the big known
run in the middle is the natural seam). Build the Herrmann–May shift lattice,
reduce it with **flatter**, and recover the root. The correct brute-force pair
factors `N`.

## What makes it fast (see `RESEARCH.md` for the full evidence)

1. **Reduced lattice construction** — off-diagonal coefficients taken mod `Nᵗ`
   (pivots exact). Shrinks entries ~2.7× (52.8k → 19.3k bits) ⇒ **~2.6× faster
   flatter**. (`build_reduced_lattice` in `arm.sage`.)
2. **Empirically calibrated `(m,t) = (18,5)`** — not the determinant table's
   conservative `(24,7)`. Calibration on planted-root synthetic instances showed
   `(18,5)` recovers the root (validated on 4 seeds) while `≤(16,4)` fails; this
   alone is a ~3–4× speedup. **Don't trust the bound as a scheduler — calibrate.**
3. **Blind CRT small-prime Gröbner recovery** — the shortest lattice vectors are
   all degree-`m`, which defeats resultants (17 h hang) and the rank-degenerate
   kernel. Instead, take the ideal's variety mod several ~30-bit primes (a single
   rational point each, ~0.2 s via Singular) and **CRT** to the exact `(x0,x1)`.
   ~5–10 s, blind. (`crt_groebner_recover` in `arm.sage`.)

## Files

| Path | Purpose |
|---|---|
| `arm.sage`              | The solver: brute-force loop + reduced HM lattice + flatter + CRT recovery + RSA decrypt. Run it. |
| `RESEARCH.md`           | Decision record: model, feasibility, the reduce optimization, the empirical `(m,t)` floor, and the recovery method — with measured numbers and papers. |
| `PROGRESS.md`           | Earlier 5-variable investigation and why it was abandoned. |
| `context.md`            | Problem data (N, e, ct, MASK, p&MASK) and attack-design notes. |
| `q7.pdf`                | Original contest problem statement. |
| `found_result_arm.txt`  | Solution: p, q, d, and the decrypted plaintext (flag). |
| `tools/calib_reduce.sage` | Synthetic planted-root calibration driver: finds the minimal working `(m,t)` and validates the reduce construction. |
| `tools/crt_recover.sage`  | Standalone validation of the CRT small-prime Gröbner root recovery against a known root. |

## Running it

Requires SageMath (via `passagemath-standard`) and **flatter** on `PATH`
(build from source; on ARM use `-mcpu=native`). Then:

```bash
sage arm.sage | tee run.log      # defaults: (18,5)/(20,6), 4 workers, flatter, checkpointed
```

Environment knobs (all optional): `PARAMS='18,5;20,6'`, `WORKERS=4`,
`REDUCTION=flatter`, `START_INDEX`/`STOP_INDEX` (candidate range),
`CHECKPOINT=checkpoint_q7.txt` (resume). On success it writes
`found_result_arm.txt` and stops.

The run needs Python 3.12, `flatter`, and the Sage stack (fpylll, sympy, gmpy2).
