# CollatzLIS

Longest increasing subsequences of Collatz trajectories: a null model decomposition.

We study the LIS of Collatz trajectories as a quantitative probe of pseudorandomness. The mean ratio LIS/2sqrt(k) is **0.567** at N=10^6, significantly below the **0.812** random baseline. A hierarchy of null models decomposes the deviation: the 31.6% up-step fraction and the long descending suffix account for most of the gap, with a residual of ~0.05 beyond matched stochastic models.

## Key results

| Finding | Value |
|---------|-------|
| Collatz LIS ratio at N=10^6 | 0.567 |
| Random permutation baseline | 0.812 |
| Shuffled Collatz baseline | 0.813 |
| Bernoulli (p=0.316) model | 0.571 |
| Multiplicative RW model | 0.638 |
| Residual beyond matched models | ~0.05 |
| Prefix-only ratio | 0.853 (near baseline) |
| LIS/LDS asymmetry | 0.36 |
| 3n+b family range | 0.49--0.61 |

## Prerequisites

- Python 3.11+
- No external dependencies (standard library only)

## Run

```bash
# Main experiment suite (scaling, parameter sensitivity, structural analysis)
python src/collatz_lis.py

# Validation suite (correctness, null models, trajectory correlation, convergence)
python src/validate.py
```

The main experiment runs scaling from N=1,000 to N=1,000,000 and takes several minutes at the largest scale. The validation suite runs five independent checks and takes approximately 5-10 minutes.

## Project structure

```
src/
  collatz_lis.py      # Main experiment suite (5 experiments)
                      #   1. Scaling (N=1k to 1M)
                      #   2. Parameter sensitivity (3n+b family)
                      #   3. Bootstrap statistics
                      #   4. Structural analysis (directionality, peak, runs)
                      #   5. Random baseline control
  validate.py         # Validation suite (5 validations)
                      #   V1. LIS algorithm correctness (brute-force check)
                      #   V2. Null model comparison (5 models)
                      #   V3. Trajectory correlation (suffix sharing)
                      #   V4. Convergence modeling (3 functional forms)
                      #   V5. Autocorrelation analysis (Bernoulli/Markov)

data/
  experiment_results.md     # Detailed findings from experiment suite
  validation_results.md     # Validation findings with corrections

paper/
  paper.tex                 # Manuscript (LaTeX, article class)
```

## The null model hierarchy

The central methodological contribution is a hierarchy of stochastic models, each capturing one additional structural property of Collatz trajectories:

| Model | What it captures | Ratio | Gap from Collatz |
|-------|-----------------|-------|-----------------|
| Random permutation | Nothing (BDJ baseline) | 0.812 | +0.219 |
| Shuffled Collatz | Value distribution only | 0.813 | +0.220 |
| Multiplicative RW | Step-size distribution | 0.638 | +0.045 |
| Bernoulli (p=0.316) | Up-fraction only | 0.571 | -0.022 |
| Markov (matched autocorr) | Up-fraction + alternation | 0.478 | -0.115 |
| **Collatz actual** | **All structure** | **0.593** | **baseline** |

Null model comparisons at N=50,000. Scaling experiments (Table 1 in the paper) extend to N=10^6.

## Validated claims

| Claim | Status |
|-------|--------|
| LIS ratio significantly below random baseline | Confirmed (0.567 vs 0.812 at N=10^6) |
| Deviation reflects temporal ordering | Confirmed (shuffling restores to 0.813) |
| Up-step fraction is the primary driver | Confirmed (Bernoulli p=0.316 matches within 0.02) |
| Ratio appears to converge | Likely, but limit uncertain (0.47--0.57 range) |
| LIS algorithm correct | Verified against brute-force on 1000 random sequences |

## Known limitations

- **Trajectory correlation:** 97.4% shared suffixes reduce effective sample size to ~500, widening bootstrap CIs by ~6x
- **Convergence not settled:** models disagree on asymptotic limit (0.47 to 0.57); N>10^8 needed
- **No analytical derivation:** all results are empirical
- **Prefix-only ratio is near-normal (0.853):** the low full ratio is primarily a suffix effect

## References

- Baik, Deift, Johansson (1999). On the distribution of the length of the longest increasing subsequence of random permutations. *J. Amer. Math. Soc.* 12, 1119-1178.
- Lagarias (1985). The 3x+1 problem and its generalizations. *Amer. Math. Monthly* 92, 3-23.
- Tao (2022). Almost all orbits of the Collatz map attain almost bounded values. *Forum of Mathematics, Pi* 10, e12.
- Aldous, Diaconis (1999). Longest increasing subsequences: from patience sorting to the Baik-Deift-Johansson theorem. *Bull. Amer. Math. Soc.* 36, 413-432.

## License

[MIT](LICENSE)
