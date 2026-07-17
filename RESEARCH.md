# RSA Lattice Attack — Optimization Plan & Evidence

Forward-looking research notes for Q7 (RSA partial-key exposure). This is the
decision record: what to do, and the reasoning/evidence/papers behind each
choice. Dead-end detail from earlier attempts has been dropped; only what
changes the next run is kept.

---

## 0. Conclusion — VALIDATED END-TO-END (measured this session)

The full pipeline is proven on planted-root synthetic instances of the exact
problem geometry (230-bit x0 window at bit 600, 155-bit x1 window at bit 265,
β=0.5). `arm.sage` now implements all of it.

1. **2-variable relaxed model.** Geometry forces it (§2); 5-variable is a dead
   end (needs `m≈67`, dimension `C(m+5,5)`).

2. **Minimal robust config is `(m,t) = (18,5)`, NOT `(24,7)`.** The determinant
   slack table (§3) is *hugely conservative* — it forbade everything below
   `(24,7)`, but empirically **(18,5) recovers the root** (120/120 shortest
   vectors vanish at the true root, full rank; validated on 4 independent
   seeds), while everything **≤ (16,4) fails** (0 vanishing). Sharp cliff
   between (16,4) and (18,5). Use `(18,5)`; `(20,6)` as a proven fallback.
   → This is RESEARCH's thesis (bounds don't predict practical success) made
   concrete: ~6 units of `m` and a ~3-4x speedup were left on the table by
   trusting the gate.

3. **Reduce the lattice entries mod `Nᵗ`** (off-diagonal only, pivots exact).
   Measured: max entry 52.8k → 19.3k bits, flatter **3.7 h → 1.43 h at m=24
   (2.6×)**; at `(18,5)` one reduction is **~0.37 h**. Correctness preserved
   (identical short vectors vanish). Implemented in `build_reduced_lattice`.

4. **Root recovery = blind CRT small-prime Gröbner** (§11). The shortest vectors
   are all degree-`m`, so pairwise resultants hang (the old 17-h trap) and the
   naive kernel is rank-degenerate. Instead: reconstruct the shortest rows as
   polys mod several ~30-bit primes, take the ideal's variety (a *single*
   rational point per prime, ~0.2 s each via Singular), and CRT (product >
   2²³⁰) to the exact `(x0,x1)`. **~5–10 s, blind**, general across `m`.
   Implemented in `crt_groebner_recover`.

5. **Per-candidate ≈ 0.37 h** (flatter-dominated; recovery negligible). Real
   256-way run on 4 single-threaded workers: **~24 h worst, ~12 h expected**
   (early exit) — down from the original ~10 days.

6. **`flatter` mandatory** (built N1-native, `-mcpu=neoverse-n1`); `WORKERS=4`
   single-threaded (`OMP=1`, `MALLOC_ARENA_MAX=2`); `FEASIBILITY_GATE=0` — the
   `hm_bound_ok` gate is wrong (implies ~724-bit capacity, impossible) and is
   superseded by the empirical floor above.

7. **Do not brute-force extra bits to lower `m`** — on 4 cores the ×2 candidate
   blow-up per bit outweighs the cheaper lattice (§5). Stay at 256 candidates.

---

## 1. Problem shape

`N = p·q`, both primes 1024-bit, `N` 2048-bit. 672 bits of `p` known via
`MASK`; 352 unknown in 7 blocks. The two 4-bit blocks (offsets 150, 920) are
brute-forced (256 combinations). After that brute force the run structure is:

```
K194 | U230 | K180 | U155 | K265        (MSB → LSB, offsets below)
```

At the true root the recovery polynomial is linear in the unknown blocks and
vanishes modulo the unknown divisor `p ≈ N^0.5`. This is exactly the
**Herrmann–May "linear equations modulo an unknown divisor"** setting.

---

## 2. Why the 2-variable split (geometry, not preference)

After brute-forcing the two 4-bit gaps, the unknown region spans bits
**[265, 830) = 565 bits contiguous**, cleaved by the 180-bit known run at
[420, 600). This dictates the number of variables:

- **n = 1 is impossible.** Single-block Coppersmith ("factoring with high bits
  known") recovers an unknown block up to `N^(β²) = 2^512` bits of `p` at
  β = 0.5. Our span is 565 > 512. So ≥ 2 variables are required.
- **n = 2 is the minimum, and the K180 run is the natural seam:**

  ```
  p = base + x0·2^600 + x1·2^265
  x0 < 2^230   spans p[829:600]  = U46 | K15 | U87 | K13 | U69
  x1 < 2^155   spans p[419:265]  = U58 | K13 | U84
  ```

- **More variables are strictly worse.** Herrmann–May's cost is exponential in
  the number of blocks: lattice dimension `C(m+n, n)` and the recoverable-bit
  bound both degrade with `n`. The dead 5-variable model needed `m ≈ 67`
  (dimension `≈ 14M`); n = 2 keeps dimension `(m+1)(m+2)/2` (hundreds).

**Cost of the relaxation:** the lattice sees 385 root-bound bits (230 + 155),
not the true 344 unknown bits — the 41 interior known bits (K15, K13, K13) are
*not* used to shrink the lattice. They are recovered and validated afterward
(internal-bit filters, §7). This 41-bit waste is the price of dropping from 5
variables to 2, and it is worth it.

*Evidence:* Coppersmith 1996 (univariate/bivariate small roots, the 512-bit
n = 1 limit); Herrmann–May 2008 (the exact n-block model and its
exponential-in-n cost).

---

## 3. Feasibility & parameter schedule (the core evidence)

The binding question is not "is it feasible?" (it is) but "which `(m,t)` is
strong enough?" The determinant / Howgrave-Graham gate for this basis is
explicit because the shift set is triangular.

**Shift set** (current baseline, `x1`-only):

```
g_{i,k} = x1^i · f^k · N^max(t-k, 0),   0 ≤ k ≤ m,  0 ≤ i ≤ m-k
```

**Closed form** (monomials ordered by total degree):

```
w   = (m+1)(m+2)/2                     # dimension
C   = m(m+1)(m+2)/6
E_N = t(t+1)(3m - t + 4)/6
S   = log2(X0) + log2(X1) = 230 + 155 = 385

avg det bits per row = (C·S + E_N·log2 N) / w
gate (usable short vector exists): (C·S + E_N·log2 N) / w  <  t·log2(p) = t·1024
```

**Asymptotic capacity.** With `τ = t/m`, dropping finite-m terms:

```
S/log2(N) < τ³ − 3τ² + 3βτ
```

For β = 0.5 the RHS is maximized at `τ = 1 − √(1−β) = 0.29289`, giving capacity
`0.2071 · 2048 ≈ 424 bits`. Our root is 385 bits → **inside by only ~39 bits**.
Near-boundary instances predictably need large lattices (this is Miller–
Narayanan–Venkatesan's central warning: enabling conditions do not guarantee
practical success near the boundary).

**Finite-m slack (bits above the gate; negative = cannot work):**

| m | best t | τ | dim | slack |
|---|---|---|---|---|
| 18 | 6 | 0.28 | 190 | −89 ❌ |
| 20 | 6 | 0.30 | 231 | −22 ❌ |
| 22 | 7 | 0.27 | 276 | −18 ❌ |
| **24** | **7** | **0.29** | **325** | **+30 ✅** |
| 26 | 8 | 0.31 | 378 | +44 ✅ |
| 28 | 8 | 0.29 | 435 | +79 ✅ |
| 30 | 9 | 0.30 | 496 | +102 ✅ |

**Consequences:**
- Viable regime **starts at (24,7)**. Recommended ladder: `24,7 ; 26,8 ; 28,8`;
  `(30,9)` as the safe-but-expensive fallback.
- `arm.sage`'s `DEFAULT_PARAMS = [(12,4),(15,5),(18,6)]` and its "stronger"
  hints `(20,6);(22,7)` are all below the gate. **Delete them.**
- `arm.sage`'s `hm_bound_ok` implies a ~724-bit capacity — impossible (exceeds
  the n = 1 limit of 512). It is the same class of hand-derived formula that
  misled an earlier attempt. **Disable `FEASIBILITY_GATE`;** use this slack
  table as the scheduler instead.

---

## 4. Hardware reality (4-core Neoverse-N1, 24 GB, no swap)

| Resource | Value | Implication |
|---|---|---|
| CPU | 4× Neoverse-N1, 1 thread/core | 4-way candidate parallelism; no SMT |
| RAM | 24 GB, **0 swap** | fine at m ≤ 30 (§ below); OOM = kill, so cap per-worker RSS |
| Disk | ~99 GB free | ample for checkpoints/dumps |
| Reducer | `flatter` **not yet installed** | build it first — hard blocker |

**Memory at the real m.** Per-worker peak RSS ≈ 1–2 GB (m=24→~1 GB,
m=30→~2 GB); 4 workers ≈ 6–8 GB, well under 24 GB. **No need to reduce
`WORKERS` for memory.** Keep the single-thread env caps
(`OMP/OPENBLAS/MKL_NUM_THREADS=1`, `MALLOC_ARENA_MAX=2`) so the `flatter`
subprocess does not oversubscribe the 4 cores.

**Time is the binding constraint.** At the required m=24–30 these lattices are
comparable to or larger than the old n=5 462-dim/12k-bit lattice (~8 min each).
Realistic per-candidate cost is **minutes to tens of minutes**; 256 candidates
× tiers ÷ 4 cores ⇒ a **multi-hour-to-multi-day** job. Exact numbers to be
measured in calibration (§6) — do not start a blind real sweep.

---

## 5. Trap: brute-forcing extra bits does NOT help here

RESEARCH-worthy but rejected on this box. Brute-forcing `b` extra high bits of
`x0` reduces `S` by `b`, improving slack by `C/w = m/3` bits per bit — so a few
extra bits can flip a marginal `(m,t)` and let you drop `m`. **But each extra
bit doubles the candidate count.** Accounting:

- Dropping (28,8)→~(22,7) needs ~3 extra bits ⇒ ×8 candidates (256→2048).
- The (22,7) lattice is only ~3× cheaper per candidate than (28,8).
- Net: ~2.5× **more** total work.

On a 4-core box the candidate multiplier dominates. **Stay at 256 candidates,
pay the higher m.** (Extra-bit brute is a valid lever only on many-core
clusters.)

---

## 6. Calibrate before committing cores (highest-yield task)

Build a 2-variable ground-truth calibration driver (analogous to the old
`gt_calibrate.py`, but for this model). It must **not** touch the real unknown
root. On synthetic same-mask instances with a known `p`, it answers:

> For this exact 2-variable shape, which shift rows + reducer produce two
> independent true-root-vanishing relations at the lowest cost?

Log per trial: `(m, t, shift_policy, reducer, seed)`, lattice dimension, row
tags `(k, a, b)`, reduced-row index, polynomial weighted degree/height, whether
it vanishes at the true `(x0,x1)`, and whether a resultant pair recovers the
root. Two payoffs:

1. **Pins the minimal working `(m,t,policy)`** → the real sweep runs one
   validated configuration instead of a blind ladder.
2. **Benchmarks per-candidate wall time on this ARM box** → real ETA.

*Evidence:* ground-truth calibration is what rescued the earlier investigation;
Ajani–Bright reinforce the same operational idea — use fast algebraic checks as
pruning signals, not just final solvers.

---

## 7. Levers to lower required `m` (construction, not size)

The way to buy back the 4-core budget is smaller/smarter lattices at the same
strength, not bigger ones.

1. **Focus-group row selection** (Miller–Narayanan–Venkatesan). Coppersmith
   lattices carry many rows that satisfy the determinant story but never appear
   in a reduced vector that vanishes at the root. Learn the productive
   sublattice on smaller synthetic instances (heatmap over shift params
   `(k, a, b)`), keep the hot rows, drop cold ones → lower effective dimension.

2. **Asymmetric shift shape** (Bloemer–May, Jochemsz–May, Coron). The current
   `x1`-only family is a baseline, not sacred. Try a weighted-degree family
   that also mixes `x0`:

   ```
   g = x0^a · x1^b · f^k · N^max(t-k,0),   with   230·a + 155·b ≤ W_k
   ```

   Start with `a ∈ {0,1}` (then `{0,1,2}`); let focus-group testing decide
   which rows survive. Useful constructions follow the polynomial's monomial
   geometry, not a blind total-degree cap.

3. **Reduce once, mine many rows.** The right relation pair may be rows 5–40.
   **Increase `MAX_POLYS` before increasing `m`.** Rank pairs by weighted
   degree + unscaled height + nonzero resultant; score by synthetic true-root
   vanishing, not just degree/term-count.

4. **Cheaper post-LLL via internal-bit filters.** Apply the interior known-bit
   constraints as congruence filters on resultant roots *before* candidate
   reconstruction: `x0` must match its 15-bit and 13-bit islands, `x1` its
   13-bit island. Roots outside those residue classes are discarded early. Keep
   resultant-first extraction (Groebner blew memory on bad relations before).

5. **Reducer is a first-class parameter.** `flatter` for all real work; Sage LLL
   only for tiny calibration. Compare reducers by *usable vanishing relations
   per second*, not wall time. A slower reducer wins if it lets us lower `m`.

6. **Tiered passes.** Cheap (low delta, small `MAX_POLYS`/`PAIR_LIMIT`, broad
   `(m,t,policy)` scan) → medium (stronger, only policies that worked
   synthetically) → strong (`flatter`/BKZ, best 1–2 policies only).

---

## 8. Execution order

1. Build environment: `flatter` (+ deps) and Sage. **Blocker.**
2. Fix `arm.sage`: `PARAMS = 24,7;26,8;28,8`, `FEASIBILITY_GATE=0`,
   `WORKERS=4`, `flatter` reducer, per-worker RSS cap.
3. Write the synthetic calibration driver (§6); find minimal `(m,t,policy)` +
   benchmark. Iterate §7 levers here (cheap, single-core).
4. Real 256-way run at the validated config, tiered, checkpointed, launched via
   a proper background flag (not shell `&`).
5. On success: `q = N/p`, `d = e⁻¹ mod (p−1)(q−1)`, `pt = ct^d mod N`.

---

## 9. Fallback (only if both lattice paths fail)

SAT + Coppersmith hybrid (Ajani–Bright): encode `p·q = N` with the known mask,
let a SAT solver assign low bits of `p`, call univariate Coppersmith to test
whether a partial assignment can extend to a factor, learn blocking clauses.
Credible for arbitrary-position leaked bits, but a heavy new implementation and
their experiments are far below 2048-bit scale. Treat as last resort, not a
shortcut.

---

## 10. Papers (why each matters going forward)

- **Herrmann–May, ASIACRYPT 2008** — the exact model; cost exponential in
  block count ⇒ minimize variables. https://doi.org/10.1007/978-3-540-89255-7_25
- **Coppersmith, EUROCRYPT 1996** — base small-root tool; the n = 1 512-bit
  limit that forces n ≥ 2. https://doi.org/10.1007/3-540-68339-9_14 ,
  https://doi.org/10.1007/3-540-68339-9_16
- **Miller–Narayanan–Venkatesan (focus groups)** — bounds don't predict
  practical success near the boundary; learn productive sublattices; reducer
  choice matters. https://arxiv.org/abs/1708.09445
- **Bloemer–May / Jochemsz–May / Coron** — small-root success is monomial/shift
  selection, not one determinant inequality; keep it bivariate and solve with
  bivariate-specific methods. https://doi.org/10.1007/11426639_15 ,
  https://doi.org/10.1007/11935230_18 , https://doi.org/10.1007/978-3-540-74143-5_21
- **Ajani–Bright, ISSAC 2024** — SAT + Coppersmith fallback; algebraic checks as
  pruning signals. https://arxiv.org/abs/2406.20071

---

## 11. Measured results & the recovery method (this session)

Everything below is measured on planted-root synthetic instances with the exact
geometry, single-threaded on the 4-core Neoverse-N1. All of it is in `arm.sage`.

### Reduce optimization (`build_reduced_lattice`)
The raw construction expands `f^k` with un-reduced ~2048-bit coefficients
(`A`,`B` are full residues after making `f` monic), so entries reach `~2048·m`
bits (52,795 at m=24) and flatter took **3.7 h**. Reducing off-diagonal
coefficients mod `Nᵗ` (pivots kept exact → full rank) caps entries near
`log2(Nᵗ)+log2(X²⁴) ≈ 20k` bits (measured 19,326). flatter then: **1.43 h at
m=24 (2.6×), ~0.37 h at (18,5)**. Speedup is ~linear in bit-length (flatter is
precision-bound). Correctness verified: the same short vectors vanish at the
true root.

### Empirical (m,t) floor
| (m,t) | dim | flatter (reduced) | vanish@root |
|---|---|---|---|
| ≤ (16,4) | ≤153 | 0.1–223 s | **0** ❌ |
| **(18,5)** | 190 | **~0.37 h** | 120/120 ✅ |
| (20,6) | 231 | ~0.48 h | 120/120 ✅ |
| (22,6) | 276 | ~0.70 h | 120/120 ✅ |
| (24,7) | 325 | 1.43 h | 120/120 ✅ |

Sharp cliff at (16,4)→(18,5); `(18,5)` robust across seeds {12345,7,99,2024}.
The slack table (§3) is correct as a *lower bound on capacity* but wildly
pessimistic as a *scheduler* — do not use it to pick `m`; calibrate.

### Root recovery (`crt_groebner_recover`) — the hard-won piece
The shortest reduced vectors are all **degree-m** polynomials, which defeats the
obvious methods:
- pairwise **resultants** of degree-18/20 bivariates hang (the original 17-hour
  runaway; a single pair timed out at 120 s);
- the **kernel** of the vanishing rows is ~35-dimensional (196/231 rows vanish,
  can't pin a 231-monomial vector), and the multiplication operator on it is
  rank-deficient, so the eigenvalue method needs border-basis handling;
- **large-P Gröbner** (P > 2²³⁰ for an exact lift) uses Sage's slow generic
  path (fast Singular needs P < 2³¹).

The working method: the ideal is cleanly **zero-dimensional with a single
rational point** mod a small prime. So compute the variety mod ~10–12 distinct
~30-bit primes (Singular, ~0.2 s each) and **CRT** the residues (product > 2²³⁰)
to the exact `(x0,x1)`. **Blind** (feeds the shortest lattice rows; no knowledge
of the root), **~5–10 s**, general across `m`. Guarded by a 150 s SIGALRM per
candidate so wrong brute-force guesses (255/256) fail fast. Unscaling to poly
coefficients is a modular inverse of `X^a Y^b`, cheap mod p.

### Operational
- `flatter` built `-mcpu=neoverse-n1`, on PATH, visible to Sage.
- Per-candidate ≈ 0.37 h; real 256-way run ~12–24 h on 4 workers, checkpointed.
- Memory: fine (per-worker ~1–2 GB, no swap needed at m≤30).
- Note: concurrent flatter workers do **not** materially slow each other
  (measured: 3 parallel (18,5) all ~0.37 h), but naive Python rational
  reconstruction under that load is ~1000× slower — the CRT recovery avoids it
  by working mod small primes.
