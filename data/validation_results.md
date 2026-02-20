# Collatz LIS -- Validation Results

**Date:** 2026-02-20
**Script:** poc2_collatz_validate.py

---

## Summary: What's Real and What's Not

The validation suite reveals that the Collatz LIS finding is **real but more nuanced than
we initially claimed**. Several important corrections to the narrative:

### What IS real:
- The LIS ratio (~0.57 at N=1M) is genuinely lower than random permutations (~0.81)
- The ratio converges to a positive constant (estimated ~0.47 by log fit)
- The deviation from random IS about temporal structure, not just value distribution
- LIS code is verified correct (brute-force match on 1000 cases)

### What needs correction:

**1. The BDJ comparison is WRONG as framed.**
Comparing Collatz to 2*sqrt(k) (random permutation prediction) is not the right null model.
Collatz trajectories are multiplicative sequences with specific step distributions, not
permutations. The random permutation ratio is 0.81 (not 1.0) at these sequence lengths,
so our "43% deviation" is actually closer to a **27% deviation** (0.59 vs 0.81).

**2. The multiplicative random walk nearly matches (0.64 vs 0.59).**
A simple random walk with matched up/down fraction (33%) and step magnitudes (3x up, 0.5x
down) produces ratio 0.64 -- only 0.05 from Collatz's 0.59. The "deep structure" claim
is weaker than we thought. Most of the ratio is explained by the **step distribution alone**.

**3. The Bernoulli model also nearly matches (0.57 vs 0.59).**
An independent Bernoulli process with p=0.316 up-probability produces ratio 0.57.
This means the up/down fraction ALONE explains most of the LIS ratio.
The autocorrelation structure (lag-1 = -0.46, strongly alternating) does NOT make it worse --
the Markov model actually gives a LOWER ratio (0.48), not higher.

**4. Trajectory correlation is MASSIVE (97.4% shared suffixes).**
Almost all trajectories share most of their length as a common suffix to 1. This means
our "1 million independent samples" are really ~500 effective independent samples. Our
bootstrap CIs need to be widened by ~6x, making the CI roughly [0.55, 0.59] instead of
[0.565, 0.568].

**5. The ratio has NOT converged yet.**
The best-fit convergence model (r = 0.469 + 1.38/log(N)) predicts the asymptotic limit
at ~0.47, but the log-convergence is slow. At N=1B, the predicted ratio is still 0.54.
We don't know the true limit with high confidence.

---

## Detailed Findings

### V1: LIS Code Correctness -- PASS
- All 10 named test cases pass
- 0/1000 mismatches against brute-force on random sequences of length 10
- Algorithm is correct

### V2: Null Model Comparison

| Model | Ratio | Gap from Collatz |
|---|---|---|
| Collatz actual | 0.593 | baseline |
| Random permutation (BDJ) | 0.812 | +0.219 |
| Multiplicative RW (matched params) | 0.638 | +0.045 |
| Geometric RW (exact step dist) | 1.666 | +1.073 |
| Shuffled Collatz values | 0.813 | +0.219 |

**Key insight:** Shuffled Collatz (0.813) matches random permutation (0.812). This proves
that the LOW ratio is about **temporal ordering**, not about the value distribution. When
you destroy the ordering, the ratio goes back to the random baseline.

The multiplicative RW (0.638) is the closest match. The remaining gap (0.045) is the
"genuine Collatz structure" beyond simple step statistics.

The geometric RW with exact step distributions gives 1.67 (too high) because sampling
actual Collatz step ratios creates artificially high variance in value space.

### V3: Trajectory Correlation

- **97.4%** of each trajectory is shared suffix with previously seen trajectories
- 93% of trajectories share >90% of their length
- Effective independent sample size: ~500 out of 20,000
- **This is a real problem for our statistics.** Bootstrap CIs are 6x too narrow.
- However: the PREFIX (up to peak) is where LIS lives, and prefixes are more independent
  - Full trajectory ratio: 0.590
  - Prefix-only ratio: 0.853
  - **The low ratio comes from the SUFFIX, not the prefix!**

**Critical realization:** LIS of the full trajectory is LOW because the long shared
descending suffix dilutes the LIS. The LIS/2sqrt(k) ratio is low partly because k is
inflated by the long monotone descent, while LIS can't grow during descent. This is a
**mechanical artifact of trajectory structure**, not a deep number-theoretic signal.

### V4: Convergence Modeling

| Model | Asymptotic Limit | RMSE |
|---|---|---|
| r = a + b/log(N) | 0.469 | 0.0023 (best) |
| r = a + b/sqrt(N) | 0.572 | 0.0064 |
| r = a + b/N^0.25 | 0.549 | 0.0032 |

The log model fits best, predicting convergence to ~0.47. But the models disagree on the
limit (0.47 to 0.57 range), and all give decent fits. We cannot confidently determine the
asymptotic limit from N=1M data alone.

Predictions:
- N=10M: 0.555 (log), 0.573 (sqrt)
- N=100M: 0.544 (log), 0.573 (sqrt)
- N=1B: 0.536 (log), 0.572 (sqrt)

### V5: Autocorrelation Analysis

- Lag-1 autocorrelation of up/down indicator: **-0.458** (strongly ALTERNATING)
- Average up-run length: 1.00 (up-steps are almost always isolated)
- Average down-run length: 1.98 (down-steps come in pairs -- this is the 3n+1 -> /2 -> /2 pattern)

| Model | Ratio |
|---|---|
| Collatz actual | 0.594 |
| Independent Bernoulli (p=0.316) | 0.571 |
| Markov chain (matched autocorr) | 0.478 |
| sqrt(0.316/0.5) theoretical | 0.795 |

**The Bernoulli model nearly matches Collatz.** The autocorrelation does NOT explain the gap --
in fact, the strongly alternating pattern (negative autocorrelation) gives the Markov model
a LOWER ratio than Collatz, not higher. The up-fraction alone (31.6%) is the primary driver.

The sqrt(p/0.5) = 0.795 theoretical prediction is wrong because it assumes LIS growth is
proportional to sqrt(fraction of increasing steps). The actual relationship is more complex
because BDJ applies to the rank-permutation, not directly to the step-fraction.

---

## Revised Assessment of What's Publishable

### Strong claims (supported):
1. Collatz LIS/2sqrt(k) ratio is significantly lower than random baseline (~0.59 vs ~0.81)
2. The ratio is about temporal ordering (shuffling destroys it)
3. The up-step fraction (31.6%) is the primary driver
4. The ratio appears to converge to a positive constant (est. 0.47-0.57)

### Weak claims (need qualification):
1. ~~"Massive Tracy-Widom deviation"~~ -- comparing to TW is the wrong framework;
   Collatz trajectories are not random permutations
2. ~~"Ratio ~ 0.57 is a universal constant"~~ -- it's still drifting and model-dependent
3. ~~"Novel deep structure"~~ -- most of the effect is explained by the 31.6% up-fraction
   and the long descending suffix

### What IS genuinely interesting:
1. The **gap between multiplicative RW (0.64) and Collatz (0.59)** encodes structure
   beyond step-size distribution. It's small (0.05) but reproducible.
2. The **prefix-only ratio (0.85) vs full ratio (0.59)** reveals that Collatz's low LIS
   ratio is primarily a SUFFIX effect -- the long descent to 1 dominates.
3. The **LIS/LDS asymmetry (0.36)** is a clean quantitative characterization.
4. The **parameter sensitivity** (ratio varies with b in 3n+b) is genuinely novel.

### Honest paper framing:
Not "Collatz has deep non-random LIS structure" but rather "LIS as a pseudorandomness
measure for Collatz trajectories: what the ratio tells us about trajectory geometry."
The contribution is the FRAMEWORK (using LIS as a probe) and the DECOMPOSITION (showing
which properties drive the ratio), not a single dramatic deviation number.

---

## Status: VALIDATED with corrections. Still publishable, but with a more honest framing.
