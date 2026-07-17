# Q7 — RSA Partial Key Exposure: Progress & Analysis

**Contest:** 2026 암호분석경진대회 — Problem 7
**Status: unsolved.** Infrastructure is solid and validated on toy cases; the
real-scale attack has not yet succeeded. See "Open question" at the bottom.

---

## 1. Problem summary

- `N = pq`, both 1024-bit, `N` is 2048-bit.
- 672 of `p`'s 1024 bits are known (`P_AND_MASK`), 352 unknown, in 7
  discontinuous blocks (parsed by `work/blocks.py` from `MASK`):

  | offset | length | treatment |
  |---|---|---|
  | 150 | 4  | brute force (16) |
  | 265 | 84 | lattice unknown x0 |
  | 362 | 58 | lattice unknown x1 |
  | 600 | 69 | lattice unknown x2 |
  | 682 | 87 | lattice unknown x3 |
  | 784 | 46 | lattice unknown x4 |
  | 920 | 4  | brute force (16) |

- Plan: brute-force the two 4-bit blocks (256 combinations), and for each,
  recover the remaining 5 unknowns (344 bits total, `β=0.5` since `p≈√N`) via
  the **Herrmann–May "Solving Linear Equations Modulo Divisors"** (ASIACRYPT
  2008) multivariate Coppersmith lattice attack.
- Once `p` is recovered: `q=N/p`, `d=e⁻¹ mod (p-1)(q-1)`, `pt=ct^d mod N`.

## 2. Infrastructure built

- **SageMath** installed via `pip install passagemath-standard` — Ubuntu
  24.04 dropped the `sagemath` apt package (only DB sub-packages remain in
  `universe`); the official pip-packaged Sage works fine and is much faster
  to install than building/conda.
- **`flatter`** (keeganryan/flatter) built from source. This was **critical**:
  Sage's default `Matrix.LLL()` and even direct `fpylll.LLL.reduction()`
  essentially never terminate on these lattices once `t>0` mixes in
  full-`N`-sized rows into an otherwise tiny-entry basis (highly skewed
  profile). flatter reduced a 462×462 lattice (m=6) with entries up to
  12,283 bits in ~8 minutes; the alternatives were killed after 100–300s with
  no progress.
- Vendored **`jvdsn/crypto-attacks`** (MIT) for its `PartialInteger` (models
  known/unknown bit-block structure) and `herrmann_may_multivariate.py`.
- Vendored **`kionactf/coppersmith`** — a more actively-maintained,
  CTF-community implementation with automated parameter tuning and native
  flatter support. See §5 — this looks like the more promising path.

## 3. Bugs found (both fixed locally, see README's patch step for kionactf's)

### 3.1 `jvdsn/crypto-attacks`: wrong Howgrave-Graham modulus (not patched upstream, we just stopped using this library)

`herrmann_may_multivariate.modular_multivariate` calls
`reconstruct_polynomials(L, f, N, ...)`, passing **plain `N`** as the
size-filter threshold, regardless of `t`/`β`. Compare
`howgrave_graham.modular_univariate` (the n=1 case), which correctly uses
`N**m`. The right threshold for this construction is `p^t ≥ N^(β·t)` — using
plain `N` is only correct by coincidence when `β·t=1`.

**Verified empirically**: reconstructed all 462 rows of a reduced m=6,t=1
lattice with **zero** filtering, and checked each against the known
(synthetic) true root directly — **none vanish**. This is independent of any
threshold reasoning: m=6 genuinely has no usable vectors, for any t. (We also
separately confirmed the bug causes false accepts: with the buggy `N`
threshold, some rows pass the filter without actually vanishing at the root,
corrupting downstream root-finding — either hanging in Sage's slow
Groebner-then-`solve()` fallback burning 30+GB RAM, or converging to a wrong
answer.)

Net effect: after finding this, we stopped trusting jvdsn's implementation
for this attack and moved to kionactf's (§5), which has a properly-derived
closed-form feasibility/parameter formula (§4) and, after its own bugfix,
solves the toy case correctly and fast.

### 3.2 `kionactf/coppersmith`: real Python bug, patched

`rootfind_ZZ.py`'s `solve_ZZ_symbolic_linear_near_bounds_internal` calls
`do_lattice_reduction(mattrans, lllopt_symbolic_linear)` — a positional arg
— but `do_lattice_reduction(mat, **kwds)` only accepts keyword args. Fixed
by changing the call to `do_lattice_reduction(mattrans, **lllopt_symbolic_linear)`
(see README §5 for the exact patch command). Also needed:
`context.rootfindZZopt['lllopt_symbolic_linear'] = {'algorithm': FLATTER}`
since the default (`FPLLL_BKZ`) shells out to an `fplll` CLI binary we never
built (only the library, via `libfplll-dev`).

## 4. Theoretical feasibility (this part is solid — corrects an earlier wrong analysis)

I initially hand-derived a closed-form "does a short-enough vector exist"
volumetric formula to avoid spending hours of lattice reduction on doomed
`(m,t)` pairs (`work/bound_check.py`). **This formula is buggy** — it
predicted failure for a toy case that we'd already empirically verified
succeeds (m=2,t=1 on 5×6-bit unknown blocks), and separately predicted the
real problem is infeasible even at m=200. Do not trust `bound_check.py`.

`kionactf/coppersmith`'s `coppersmith_linear()` has a properly-implemented
version of the same idea (matches the literature): a closed-form
`log_N_X_bound` (max tolerable fraction of unknown bits) and an `mestimate`
(required `m` for the formula's guarantee). Computed for our real numbers
(n=5, β=0.5, blocks 84/58/69/87/46 = 344 bits):

```
log_N_X (unknown fraction)  = 344/2048 = 0.16797
log_N_X_bound (max allowed) = 0.17638   →  361.2 bits allowed
feasible: YES (344 < 361, matches context.md's claim almost exactly)
mestimate ≈ 67   →  lattice dimension C(67+5,5) ≈ 14,000,000  (infeasible)
```

The margin (17 bits out of 1024, ~1.7%) is very tight, and convergence to the
asymptotic bound is slow near the boundary — the "guaranteed" `m` this
formula gives is astronomically large. **However**, this formula is known to
be conservative/pessimistic in practice (same story as `epsilon` in jvdsn's
`factorize_p`, which similarly gives absurd default `m` estimates) — real
CTF successes routinely beat the "guaranteed" `m` by a lot. context.md's own
estimate of `m≈4-7` is presumably based on this kind of optimism, not the
rigorous bound.

We also checked whether brute-forcing *more* bits (shrinking the lattice
unknowns further, at the cost of more combinations) helps enough to reach a
tractable `m`: it doesn't. Even sacrificing 100 *more* bits (2¹⁰⁰ extra
combinations — already absurd) only brings `mestimate` down to ~10.

## 5. Empirical status

- **Toy case (5 unknowns × 6 bits, tiny margin) validates the whole
  pipeline**: jvdsn's construction "succeeded" at m=2,t=1, but on closer
  inspection (§3.1) this was actually a lucky pass through the buggy filter,
  not a rigorous success. kionactf's construction (correct filter) required
  **m=3,t=1** to genuinely succeed on the same toy case — ~3 seconds total.
- **Real scale (5 unknowns, 84/58/69/87/46 bits) with jvdsn**: m=6 tested
  exhaustively across t=0,1,2,3 (each requiring a fresh multi-hour flatter
  reduction) — all fail (ground-truth-verified: 0 of 462 reconstructed
  polynomials vanish at the true root for any t at m=6). m=8 (1287-dim) and
  m=10 (3003-dim) were left running in the background for ~10 hours with no
  result before being killed (session paused for cleanup/reassessment) — see
  §4, this matches the "guaranteed m≈67" finding; m=8/m=10 were long shots.
- **Real scale with kionactf** (more promising, in progress when paused):
  - m=3,t=1: 0 polys found
  - m=3,t=2 / t=3: 1 poly found (not enough — need ≥5 independent for n=5)
  - m=4,t=1/2/3: 1 poly found, ~95-108s each
  - m=5,t=2: was running when the session was paused for cleanup. **This is
    the next thing to try** — see below.

## 6. Open question / next steps

The core uncertainty is unresolved: **what `m` does kionactf's more
efficient construction actually need for our real 344-bit/5-block/β=0.5
problem?** It's already looking like more than jvdsn's naive construction
needed for the toy case proportionally, but kionactf's lattice may still be
meaningfully more efficient per-dimension than jvdsn's (fewer wasted shift
polynomials), which could make a real difference at this scale. Concretely:

1. Resume calibration of `kionactf/coppersmith`'s `coppersmith_linear_core`
   on the real bit sizes (`work/kiona_real_test.py`), continuing the m=5,6,7…
   sweep that was interrupted. Use `context.rootfindZZopt` with
   `TRIANGULATE`-only and `FLATTER` for the symbolic-linear step (see the
   working config in `kiona_real_test.py`).
2. If kionactf's construction also requires impractically large `m`, the
   likely conclusion is that this exact 5-block/344-bit configuration is
   not practically solvable via a generic Herrmann-May-style implementation,
   and either (a) a further-improved construction from the literature is
   needed (Takayasu–Kunihiro 2013 "Better Lattice Constructions...", or the
   2015 "Solving Linear Equations Modulo Unknown Divisors: Revisited" —
   neither has been investigated for an existing implementation yet), or
   (b) the intended solution uses a different technique than plain
   Herrmann-May entirely, or (c) context.md's `m≈4-7` claim should be
   revisited/re-verified against its original source.
3. Once/if a working `(m,t)` is found and verified on synthetic data via
   `gt_calibrate.py`-style ground-truth checking, run the real 256-way
   brute force (`work/solve.py` / `work/brute_runner.sh`) in parallel across
   the 20 available cores.

## 7. Miscellaneous notes for whoever picks this up

- The `work/` directory was moved mid-session from `/root/cryptocontest` into
  this repo; `venv/` needed its baked-in absolute paths patched (see
  README's "Note on relocating"). Already fixed as of this commit.
- `PARI` needs a bigger stack for the real-scale kernel computation
  (`pari.allocatemem(8_000_000_000)` — the 1GB default overflows).
- Every long-running lattice reduction should be launched with the tool's
  proper background-execution flag, not manual shell `&`/`nohup` — the
  latter got reaped when its originating tool call returned, once, silently
  wasting ~10 minutes of compute.
