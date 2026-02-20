"""
Collatz LIS -- Validation Suite
Are we actually seeing what we think we're seeing?

Key concerns:
  1. BDJ is for random permutations. Collatz trajectories aren't permutations.
     The right null model is random sequences with matched properties.
  2. Trajectories share suffixes (merge). Our samples are correlated.
  3. Ratio still drifting at N=1M. Is it converging or going to 0?
  4. Structural gap: 31.6% increasing predicts 0.79, we see 0.57.
  5. Sanity: is the LIS code even correct?
"""
import math
import time
import bisect
import random
from collections import defaultdict

# ============================================================
# Core (same as v2)
# ============================================================

def collatz_trajectory(n):
    traj = [n]
    while n != 1:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        traj.append(n)
    return traj

def lis_length(seq):
    if not seq:
        return 0
    tails = []
    for val in seq:
        pos = bisect.bisect_left(tails, val)
        if pos == len(tails):
            tails.append(val)
        else:
            tails[pos] = val
    return len(tails)

# ============================================================
# VALIDATION 1: Is the LIS code correct?
# ============================================================

def validate_lis_code():
    print("=" * 70)
    print("VALIDATION 1: LIS Algorithm Correctness")
    print("=" * 70)

    tests = [
        ([1, 2, 3, 4, 5], 5, "sorted ascending"),
        ([5, 4, 3, 2, 1], 1, "sorted descending"),
        ([3, 1, 4, 1, 5, 9, 2, 6], 4, "known sequence [3,1,4,1,5,9,2,6]"),
        ([10, 9, 2, 5, 3, 7, 101, 18], 4, "leetcode example"),
        ([1], 1, "single element"),
        ([1, 1, 1, 1], 1, "all equal (strict increase)"),
        ([2, 1], 1, "two elements descending"),
        ([1, 2], 2, "two elements ascending"),
        (list(range(100)), 100, "0..99 ascending"),
        (list(range(100, 0, -1)), 1, "100..1 descending"),
    ]

    all_pass = True
    for seq, expected, name in tests:
        result = lis_length(seq)
        ok = result == expected
        if not ok:
            all_pass = False
        print(f"  {'PASS' if ok else 'FAIL'}: {name}: got {result}, expected {expected}")

    # Brute-force verification on small random sequences
    print(f"\n  Brute-force verification on 1000 random sequences (n=10)...")

    def lis_brute(seq):
        """O(2^n) brute force for verification."""
        n = len(seq)
        best = 0
        for mask in range(1 << n):
            subseq = [seq[i] for i in range(n) if (mask >> i) & 1]
            if all(subseq[i] < subseq[i+1] for i in range(len(subseq)-1)):
                best = max(best, len(subseq))
        return best

    mismatches = 0
    for _ in range(1000):
        seq = [random.randint(1, 50) for _ in range(10)]
        fast = lis_length(seq)
        brute = lis_brute(seq)
        if fast != brute:
            mismatches += 1
            print(f"    MISMATCH: seq={seq}, fast={fast}, brute={brute}")

    print(f"  Brute-force mismatches: {mismatches}/1000")
    print(f"  LIS algorithm: {'VERIFIED' if mismatches == 0 and all_pass else 'FAILED'}")
    return all_pass and mismatches == 0


# ============================================================
# VALIDATION 2: Correct null model
# BDJ applies to random PERMUTATIONS. Collatz trajectories are
# NOT permutations -- they have:
#   - Non-uniform value ranges (can go much higher than starting n)
#   - Strong autocorrelation (each step depends on previous)
#   - Specific step-size distributions
# The right comparison is: sequences with matched statistical properties.
# ============================================================

def validate_null_model(n_max=50000):
    print(f"\n{'='*70}")
    print(f"VALIDATION 2: Correct Null Model (N={n_max:,})")
    print(f"{'='*70}\n")

    # Collect Collatz trajectory properties
    coll_data = []
    step_ratios_up = []    # ratio t[i+1]/t[i] when going up
    step_ratios_down = []  # ratio t[i+1]/t[i] when going down

    for n in range(2, min(n_max + 1, 20001)):
        traj = collatz_trajectory(n)
        k = len(traj)
        l = lis_length(traj)
        coll_data.append((k, l))

        for i in range(k - 1):
            if traj[i+1] > traj[i]:
                step_ratios_up.append(traj[i+1] / traj[i])
            else:
                step_ratios_down.append(traj[i+1] / traj[i])

    # Characterize step distributions
    avg_up = sum(step_ratios_up) / len(step_ratios_up) if step_ratios_up else 0
    avg_down = sum(step_ratios_down) / len(step_ratios_down) if step_ratios_down else 0
    frac_up = len(step_ratios_up) / (len(step_ratios_up) + len(step_ratios_down))

    print(f"  Collatz step statistics:")
    print(f"    Fraction of up-steps: {frac_up:.4f}")
    print(f"    Average up-step ratio (t[i+1]/t[i]): {avg_up:.4f}")
    print(f"    Average down-step ratio:             {avg_down:.4f}")
    print(f"    (Up steps are ~3x+1, then half. So typical up ratio ~ (3n+1)/(n) ~ 3)")
    print(f"    (Down steps are always /2, so ratio = 0.5)")

    # NULL MODEL A: Random permutation (BDJ baseline)
    print(f"\n  --- Null Model A: Random Permutations ---")
    rand_perm_ratios = []
    for k, _ in coll_data[:5000]:
        perm = list(range(k))
        random.shuffle(perm)
        l = lis_length(perm)
        rand_perm_ratios.append(l / (2 * math.sqrt(k)) if k > 0 else 0)
    rpm = sum(rand_perm_ratios) / len(rand_perm_ratios)
    print(f"  Mean ratio: {rpm:.4f} (should be ~1.0 for large k)")

    # NULL MODEL B: Random walk with matched step distribution
    # At each step: with prob frac_up, multiply by avg_up; else multiply by avg_down
    print(f"\n  --- Null Model B: Multiplicative Random Walk (matched step dist) ---")
    rw_ratios = []
    for k, _ in coll_data[:5000]:
        traj = [1000.0]  # start value
        for _ in range(k - 1):
            if random.random() < frac_up:
                traj.append(traj[-1] * avg_up)
            else:
                traj.append(traj[-1] * avg_down)
        l = lis_length(traj)
        rw_ratios.append(l / (2 * math.sqrt(k)) if k > 0 else 0)
    rwm = sum(rw_ratios) / len(rw_ratios)
    print(f"  Mean ratio: {rwm:.4f}")

    # NULL MODEL C: Geometric random walk (multiply by fixed ratio each step)
    # Expected behavior: should have same inc fraction as Collatz
    print(f"\n  --- Null Model C: Geometric RW (exact Collatz step distribution) ---")
    # Use the ACTUAL step-ratio distributions, sampled with replacement
    all_steps = [(r, True) for r in step_ratios_up[:10000]] + \
                [(r, False) for r in step_ratios_down[:10000]]
    grw_ratios = []
    for k, _ in coll_data[:5000]:
        traj = [random.randint(100, 10000)]
        for _ in range(k - 1):
            ratio, is_up = random.choice(all_steps)
            traj.append(traj[-1] * ratio)
        l = lis_length(traj)
        grw_ratios.append(l / (2 * math.sqrt(k)) if k > 0 else 0)
    grwm = sum(grw_ratios) / len(grw_ratios)
    print(f"  Mean ratio: {grwm:.4f}")

    # NULL MODEL D: Shuffled Collatz trajectories
    # Same VALUES as Collatz but in random order -- destroys temporal structure
    print(f"\n  --- Null Model D: Shuffled Collatz Trajectories ---")
    shuf_ratios = []
    for n in range(2, min(5001, n_max + 1)):
        traj = collatz_trajectory(n)
        k = len(traj)
        shuffled = list(traj)
        random.shuffle(shuffled)
        l = lis_length(shuffled)
        shuf_ratios.append(l / (2 * math.sqrt(k)) if k > 0 else 0)
    sm = sum(shuf_ratios) / len(shuf_ratios)
    print(f"  Mean ratio: {sm:.4f}")

    # Collatz actual
    coll_ratios = [l / (2 * math.sqrt(k)) for k, l in coll_data[:5000] if k > 1]
    cm = sum(coll_ratios) / len(coll_ratios)

    # Summary
    print(f"\n  === NULL MODEL COMPARISON ===")
    print(f"  {'Model':>40} | {'Ratio':>8} | {'vs Collatz':>10}")
    print(f"  {'-'*40}-+-{'-'*8}-+-{'-'*10}")
    print(f"  {'Collatz actual':>40} | {cm:8.4f} | {'baseline':>10}")
    print(f"  {'A: Random permutation (BDJ)':>40} | {rpm:8.4f} | {rpm-cm:+10.4f}")
    print(f"  {'B: Multiplicative RW (matched params)':>40} | {rwm:8.4f} | {rwm-cm:+10.4f}")
    print(f"  {'C: Geometric RW (exact step dist)':>40} | {grwm:8.4f} | {grwm-cm:+10.4f}")
    print(f"  {'D: Shuffled Collatz values':>40} | {sm:8.4f} | {sm-cm:+10.4f}")

    print(f"\n  Interpretation:")
    if abs(grwm - cm) < 0.05:
        print(f"  ** Model C (exact step dist) MATCHES Collatz -> the ratio is fully")
        print(f"     explained by step-size distribution. No additional mystery.")
    elif abs(rwm - cm) < 0.05:
        print(f"  ** Model B (matched params) matches -> ratio explained by up/down")
        print(f"     fractions and average magnitudes alone.")
    elif abs(sm - cm) < 0.05:
        print(f"  ** Model D (shuffled) matches -> ratio is about the VALUE distribution,")
        print(f"     not the temporal ordering. LIS is insensitive to Collatz dynamics.")
    else:
        print(f"  ** No null model matches Collatz. The ratio encodes genuine temporal")
        print(f"     structure of the Collatz map beyond step-size distribution.")
        if sm > cm + 0.05:
            print(f"     Shuffled > Collatz means ordering matters (temporal autocorrelation).")
        if grwm > cm + 0.05:
            print(f"     Geometric RW > Collatz means step DEPENDENCIES matter")
            print(f"     (not just marginal distribution).")

    return {'collatz': cm, 'perm': rpm, 'mult_rw': rwm, 'geo_rw': grwm, 'shuffled': sm}


# ============================================================
# VALIDATION 3: Trajectory correlation (shared suffixes)
# ============================================================

def validate_correlation(n_max=50000):
    print(f"\n{'='*70}")
    print(f"VALIDATION 3: Trajectory Correlation (N={n_max:,})")
    print(f"{'='*70}\n")

    # How many trajectories share the suffix after hitting value v?
    # Measure: at what point do trajectories merge?
    # For each n, find the first value that appeared in a previous trajectory.

    seen_values = {}  # value -> first n that produced it
    merge_points = []  # (n, fraction_of_trajectory_that_is_shared)

    for n in range(2, n_max + 1):
        traj = collatz_trajectory(n)
        k = len(traj)

        # Find earliest merge point
        merge_idx = k  # no merge
        for i, v in enumerate(traj):
            if v in seen_values and seen_values[v] != n:
                merge_idx = i
                break

        shared_frac = (k - merge_idx) / k if k > 0 else 0
        merge_points.append(shared_frac)

        # Record all values as seen by this trajectory
        for v in traj:
            if v not in seen_values:
                seen_values[v] = n

    avg_shared = sum(merge_points) / len(merge_points)

    print(f"  Average fraction of trajectory that is shared suffix: {avg_shared:.4f}")
    print(f"  (1.0 = entire trajectory is shared; 0.0 = completely unique)")

    # Distribution of shared fractions
    print(f"\n  Shared suffix fraction distribution:")
    hist = [0] * 10
    for f in merge_points:
        idx = min(int(f * 10), 9)
        hist[idx] += 1
    max_h = max(hist)
    for i in range(10):
        lo = i / 10
        hi = (i + 1) / 10
        bar = '#' * int(50 * hist[i] / max_h) if max_h > 0 else ''
        print(f"    [{lo:.1f}-{hi:.1f}) {bar} ({hist[i]})")

    # What does this mean for our LIS measurement?
    # If trajectories share most of their suffix, then LIS values are correlated.
    # But LIS is dominated by the UNIQUE prefix (where the trajectory goes up),
    # not the shared suffix (which is monotonically decreasing toward 1).

    # Test: compute LIS only on the unique prefix vs full trajectory
    print(f"\n  LIS of unique prefix vs full trajectory (first 10000):")
    prefix_ratios = []
    full_ratios = []
    for n in range(2, min(10001, n_max + 1)):
        traj = collatz_trajectory(n)
        k = len(traj)

        # Find merge point
        merge_idx = k
        for i, v in enumerate(traj):
            if v < n and v in range(2, n):  # crude: if value is less than n, likely shared
                # More precise: check if we've passed through a "small" value
                pass
        # Actually, let's just use the peak as the divider
        peak_val = max(traj)
        peak_idx = traj.index(peak_val)

        prefix = traj[:peak_idx + 1]  # up to peak
        suffix = traj[peak_idx:]       # from peak onward

        l_full = lis_length(traj)
        l_prefix = lis_length(prefix)
        l_suffix = lis_length(suffix)

        if k > 5:
            full_ratios.append(l_full / (2 * math.sqrt(k)))
            if len(prefix) > 1:
                prefix_ratios.append(l_prefix / (2 * math.sqrt(len(prefix))))

    avg_full = sum(full_ratios) / len(full_ratios)
    avg_prefix = sum(prefix_ratios) / len(prefix_ratios) if prefix_ratios else 0

    print(f"    Full trajectory ratio:   {avg_full:.4f}")
    print(f"    Prefix-only ratio:       {avg_prefix:.4f}")
    print(f"    (If similar: LIS is dominated by prefix, shared suffixes don't matter)")
    print(f"    (If different: shared suffixes significantly affect the LIS)")

    # Effective sample size
    # With ~{avg_shared:.0%} correlation, effective N is reduced
    effective_n = n_max * (1 - avg_shared)
    print(f"\n  Effective independent sample size: ~{effective_n:,.0f} (from {n_max:,})")
    print(f"  CI width should be multiplied by ~{math.sqrt(n_max / effective_n):.2f}x")

    return avg_shared


# ============================================================
# VALIDATION 4: Ratio convergence modeling
# ============================================================

def validate_convergence():
    print(f"\n{'='*70}")
    print(f"VALIDATION 4: Ratio Convergence Modeling")
    print(f"{'='*70}\n")

    # Recompute at many checkpoints and fit convergence model
    checkpoints = []
    for exp in range(3, 7):  # 10^3 to 10^6
        for mult in [1, 2, 5]:
            v = mult * (10 ** exp)
            if v <= 1000000:
                checkpoints.append(v)
    checkpoints = sorted(set(checkpoints))

    running_sum_lis = 0
    running_sum_k = 0
    running_count = 0
    cp_idx = 0
    cp_data = []

    for n in range(2, 1000001):
        traj = collatz_trajectory(n)
        k = len(traj)
        l = lis_length(traj)
        running_sum_lis += l
        running_sum_k += k
        running_count += 1

        if cp_idx < len(checkpoints) and n == checkpoints[cp_idx]:
            avg_lis = running_sum_lis / running_count
            avg_k = running_sum_k / running_count
            ratio = avg_lis / (2 * math.sqrt(avg_k))
            cp_data.append((n, ratio))
            cp_idx += 1

    # Fit: ratio = a + b/log(N) (logarithmic convergence model)
    # If ratio = a + b/log(N), then log(N) * (ratio - a) = b
    # Try several models:
    print(f"  Convergence data:")
    for n_val, ratio in cp_data:
        print(f"    N={n_val:>10,}: ratio={ratio:.5f}")

    # Model A: ratio -> constant (r = a + b/N^c)
    # Approximate: fit last 5 points to r = a + b/sqrt(N) via least squares
    print(f"\n  Fitting convergence models on last 8 data points...")
    fit_data = cp_data[-8:]

    # Model: r = a + b / log(N)
    # Solve via simple linear regression: r = a + b * (1/log(N))
    xs = [1.0 / math.log(n) for n, _ in fit_data]
    ys = [r for _, r in fit_data]
    n_fit = len(xs)
    sum_x = sum(xs)
    sum_y = sum(ys)
    sum_xy = sum(x * y for x, y in zip(xs, ys))
    sum_xx = sum(x * x for x in xs)

    b_log = (n_fit * sum_xy - sum_x * sum_y) / (n_fit * sum_xx - sum_x ** 2) if (n_fit * sum_xx - sum_x ** 2) != 0 else 0
    a_log = (sum_y - b_log * sum_x) / n_fit

    # Model: r = a + b / sqrt(N)
    xs2 = [1.0 / math.sqrt(n) for n, _ in fit_data]
    sum_x2 = sum(xs2)
    sum_xy2 = sum(x * y for x, y in zip(xs2, ys))
    sum_xx2 = sum(x * x for x in xs2)
    b_sqrt = (n_fit * sum_xy2 - sum_x2 * sum_y) / (n_fit * sum_xx2 - sum_x2 ** 2) if (n_fit * sum_xx2 - sum_x2 ** 2) != 0 else 0
    a_sqrt = (sum_y - b_sqrt * sum_x2) / n_fit

    # Model: r = a + b / N^0.25
    xs3 = [1.0 / (n ** 0.25) for n, _ in fit_data]
    sum_x3 = sum(xs3)
    sum_xy3 = sum(x * y for x, y in zip(xs3, ys))
    sum_xx3 = sum(x * x for x in xs3)
    b_q = (n_fit * sum_xy3 - sum_x3 * sum_y) / (n_fit * sum_xx3 - sum_x3 ** 2) if (n_fit * sum_xx3 - sum_x3 ** 2) != 0 else 0
    a_q = (sum_y - b_q * sum_x3) / n_fit

    # Compute residuals
    def rmse(xs_model, a, b):
        preds = [a + b * x for x in xs_model]
        return math.sqrt(sum((p - y) ** 2 for p, y in zip(preds, ys)) / n_fit)

    rmse_log = rmse(xs, a_log, b_log)
    rmse_sqrt = rmse(xs2, a_sqrt, b_sqrt)
    rmse_q = rmse(xs3, a_q, b_q)

    print(f"\n  {'Model':>30} | {'a (limit)':>10} | {'b':>10} | {'RMSE':>10}")
    print(f"  {'-'*30}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}")
    print(f"  {'r = a + b/log(N)':>30} | {a_log:10.5f} | {b_log:10.4f} | {rmse_log:10.6f}")
    print(f"  {'r = a + b/sqrt(N)':>30} | {a_sqrt:10.5f} | {b_sqrt:10.4f} | {rmse_sqrt:10.6f}")
    print(f"  {'r = a + b/N^0.25':>30} | {a_q:10.5f} | {b_q:10.4f} | {rmse_q:10.6f}")

    # Predictions at N=10M and N=100M
    print(f"\n  Predictions:")
    for n_pred in [10_000_000, 100_000_000, 1_000_000_000]:
        p_log = a_log + b_log / math.log(n_pred)
        p_sqrt = a_sqrt + b_sqrt / math.sqrt(n_pred)
        p_q = a_q + b_q / (n_pred ** 0.25)
        print(f"    N={n_pred:>13,}: log={p_log:.5f}  sqrt={p_sqrt:.5f}  q={p_q:.5f}")

    best_model = min([("log", a_log, rmse_log), ("sqrt", a_sqrt, rmse_sqrt), ("N^0.25", a_q, rmse_q)],
                     key=lambda x: x[2])
    print(f"\n  Best fit: {best_model[0]} model (RMSE={best_model[2]:.6f})")
    print(f"  Estimated asymptotic limit: {best_model[1]:.5f}")

    if best_model[1] > 0.4:
        print(f"  -> Ratio converges to a POSITIVE constant (~{best_model[1]:.3f})")
        print(f"  -> This IS a meaningful characteristic of Collatz trajectories")
    elif best_model[1] > 0.1:
        print(f"  -> Ratio may converge to a small positive constant (~{best_model[1]:.3f})")
        print(f"  -> Need larger N to confirm")
    else:
        print(f"  -> Ratio may be going to 0. The 'constant' could be an artifact of finite N.")
        print(f"  -> CAUTION: the finding may not survive at larger scales")

    return best_model


# ============================================================
# VALIDATION 5: Autocorrelation structure
# Why does 31.6% increasing predict 0.79 but we see 0.57?
# ============================================================

def validate_autocorrelation(n_max=20000):
    print(f"\n{'='*70}")
    print(f"VALIDATION 5: Step Autocorrelation (the 0.79 vs 0.57 gap)")
    print(f"{'='*70}\n")

    # Hypothesis: increasing steps are CLUSTERED (not independent).
    # If up-steps come in bursts followed by long down-runs,
    # the LIS can't grow as much as if ups were uniformly scattered.

    # Measure: autocorrelation of the up/down indicator sequence
    all_autocorr = []
    all_run_lengths_up = []
    all_run_lengths_down = []

    for n in range(2, n_max + 1):
        traj = collatz_trajectory(n)
        k = len(traj)
        if k < 10:
            continue

        # Binary indicator: 1 = up step, 0 = down step
        indicators = [1 if traj[i+1] > traj[i] else 0 for i in range(k-1)]

        # Autocorrelation at lag 1
        mean_ind = sum(indicators) / len(indicators)
        var_ind = sum((x - mean_ind)**2 for x in indicators) / len(indicators)
        if var_ind > 0 and len(indicators) > 2:
            cov1 = sum((indicators[i] - mean_ind) * (indicators[i+1] - mean_ind)
                       for i in range(len(indicators)-1)) / (len(indicators)-1)
            autocorr = cov1 / var_ind
            all_autocorr.append(autocorr)

        # Run lengths
        cur_run = 1
        cur_type = indicators[0]
        for i in range(1, len(indicators)):
            if indicators[i] == cur_type:
                cur_run += 1
            else:
                if cur_type == 1:
                    all_run_lengths_up.append(cur_run)
                else:
                    all_run_lengths_down.append(cur_run)
                cur_run = 1
                cur_type = indicators[i]
        if cur_type == 1:
            all_run_lengths_up.append(cur_run)
        else:
            all_run_lengths_down.append(cur_run)

    avg_autocorr = sum(all_autocorr) / len(all_autocorr) if all_autocorr else 0
    avg_up_run = sum(all_run_lengths_up) / len(all_run_lengths_up) if all_run_lengths_up else 0
    avg_down_run = sum(all_run_lengths_down) / len(all_run_lengths_down) if all_run_lengths_down else 0

    print(f"  Average lag-1 autocorrelation of up/down indicator: {avg_autocorr:.4f}")
    print(f"    (0 = independent steps; >0 = clustered; <0 = alternating)")
    print(f"  Average up-run length:   {avg_up_run:.2f}")
    print(f"  Average down-run length: {avg_down_run:.2f}")

    # Generate MATCHED null: same fraction of ups, but independent (Bernoulli)
    print(f"\n  Control: Independent Bernoulli with same up-fraction (0.316)...")
    bernoulli_ratios = []
    for n in range(2, min(5001, n_max + 1)):
        traj = collatz_trajectory(n)
        k = len(traj)
        if k < 5:
            continue

        # Generate random sequence with same up-fraction as Collatz
        # Use multiplicative model: up -> *3, down -> *0.5
        vals = [float(random.randint(100, 10000))]
        for _ in range(k - 1):
            if random.random() < 0.316:
                vals.append(vals[-1] * (2.5 + random.random()))  # up
            else:
                vals.append(vals[-1] * 0.5)  # down
        l = lis_length(vals)
        bernoulli_ratios.append(l / (2 * math.sqrt(k)))

    avg_bernoulli = sum(bernoulli_ratios) / len(bernoulli_ratios)

    # Generate MATCHED null with same autocorrelation (Markov chain)
    print(f"  Control: Markov chain matching autocorrelation...")
    # 2-state Markov: P(up|up) and P(up|down) chosen to match
    # overall fraction = 0.316 and autocorrelation
    p_up = 0.316
    # For Markov chain with states {up, down}:
    # P(up|up) = p_uu, P(up|down) = p_ud
    # Stationary: p_up = p_ud / (1 - p_uu + p_ud)
    # Autocorrelation at lag 1: rho = p_uu - p_ud
    # From: p_up = p_ud / (1 - (p_uu - p_ud))
    # Let rho = avg_autocorr, then p_uu = p_up + rho*(1-p_up), p_ud = p_up - rho*p_up
    # (approximately, for small rho corrections)
    rho = max(-0.99, min(0.99, avg_autocorr))
    p_uu = p_up + rho * (1 - p_up)
    p_ud = p_up * (1 - rho)
    p_uu = max(0.01, min(0.99, p_uu))
    p_ud = max(0.01, min(0.99, p_ud))

    print(f"  Markov params: P(up|up)={p_uu:.3f}, P(up|down)={p_ud:.3f}")

    markov_ratios = []
    for n in range(2, min(5001, n_max + 1)):
        traj = collatz_trajectory(n)
        k = len(traj)
        if k < 5:
            continue

        vals = [float(random.randint(100, 10000))]
        state = 1 if random.random() < p_up else 0

        for _ in range(k - 1):
            if state == 1:
                vals.append(vals[-1] * (2.5 + random.random()))
                state = 1 if random.random() < p_uu else 0
            else:
                vals.append(vals[-1] * 0.5)
                state = 1 if random.random() < p_ud else 0
        l = lis_length(vals)
        markov_ratios.append(l / (2 * math.sqrt(k)))

    avg_markov = sum(markov_ratios) / len(markov_ratios)

    # Collatz actual (small sample)
    coll_ratios = []
    for n in range(2, min(5001, n_max + 1)):
        traj = collatz_trajectory(n)
        k = len(traj)
        if k < 5:
            continue
        l = lis_length(traj)
        coll_ratios.append(l / (2 * math.sqrt(k)))
    avg_coll = sum(coll_ratios) / len(coll_ratios)

    print(f"\n  === AUTOCORRELATION GAP ANALYSIS ===")
    print(f"  {'Model':>35} | {'Ratio':>8}")
    print(f"  {'-'*35}-+-{'-'*8}")
    print(f"  {'Collatz actual':>35} | {avg_coll:8.4f}")
    print(f"  {'Independent Bernoulli (p=0.316)':>35} | {avg_bernoulli:8.4f}")
    print(f"  {'Markov chain (matched autocorr)':>35} | {avg_markov:8.4f}")
    print(f"  {'sqrt(0.316/0.5) prediction':>35} | {math.sqrt(0.316/0.5):8.4f}")

    print(f"\n  Diagnosis:")
    if abs(avg_markov - avg_coll) < 0.03:
        print(f"  ** Markov model MATCHES Collatz -> autocorrelation fully explains the gap")
        print(f"  ** The 0.57 ratio is a consequence of step CLUSTERING, not deep structure")
    elif abs(avg_bernoulli - avg_coll) < 0.03:
        print(f"  ** Bernoulli matches -> simple up-fraction explains everything")
    else:
        print(f"  ** Neither model fully matches -> Collatz has structure beyond")
        print(f"     simple step-correlation. Higher-order dependencies or specific")
        print(f"     step-size distributions contribute to the LIS ratio.")
        gap_from_bernoulli = avg_bernoulli - avg_coll
        gap_from_markov = avg_markov - avg_coll
        print(f"     Gap from Bernoulli: {gap_from_bernoulli:+.4f}")
        print(f"     Gap from Markov:    {gap_from_markov:+.4f}")

    return {
        'autocorr': avg_autocorr,
        'avg_up_run': avg_up_run,
        'avg_down_run': avg_down_run,
        'collatz_ratio': avg_coll,
        'bernoulli_ratio': avg_bernoulli,
        'markov_ratio': avg_markov,
    }


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    random.seed(42)

    print("=" * 70)
    print("Collatz LIS -- Validation Suite")
    print("Are we actually seeing what we think we're seeing?")
    print("=" * 70)
    print()

    v1 = validate_lis_code()

    v2 = validate_null_model(n_max=20000)

    v3 = validate_correlation(n_max=20000)

    v4 = validate_convergence()

    v5 = validate_autocorrelation(n_max=10000)

    # ============================================================
    print(f"\n{'='*70}")
    print(f"VALIDATION SUMMARY")
    print(f"{'='*70}")

    print(f"\n  1. LIS code: {'CORRECT' if v1 else 'BUGGY'}")
    print(f"  2. Null model comparison:")
    print(f"     - Collatz ratio:       {v2['collatz']:.4f}")
    print(f"     - Random permutation:  {v2['perm']:.4f}")
    print(f"     - Multiplicative RW:   {v2['mult_rw']:.4f}")
    print(f"     - Exact step dist RW:  {v2['geo_rw']:.4f}")
    print(f"     - Shuffled Collatz:    {v2['shuffled']:.4f}")
    print(f"  3. Trajectory correlation: {v3:.1%} shared suffix")
    print(f"  4. Convergence: {v4[0]} model, limit ~{v4[1]:.4f}")
    print(f"  5. Autocorrelation: lag-1 = {v5['autocorr']:.4f}")
    print(f"     - Collatz:   {v5['collatz_ratio']:.4f}")
    print(f"     - Bernoulli: {v5['bernoulli_ratio']:.4f}")
    print(f"     - Markov:    {v5['markov_ratio']:.4f}")

    print(f"\n  OVERALL VERDICT:")
    concerns = []
    if not v1:
        concerns.append("LIS code is buggy!")
    if abs(v2['geo_rw'] - v2['collatz']) < 0.03:
        concerns.append("Exact-step-dist random walk matches Collatz -- finding may be trivial")
    if abs(v2['shuffled'] - v2['collatz']) < 0.03:
        concerns.append("Shuffled Collatz matches -- temporal structure doesn't matter")
    if v4[1] < 0.1:
        concerns.append("Ratio may converge to 0 -- not a meaningful constant")
    if abs(v5['markov_ratio'] - v5['collatz_ratio']) < 0.03:
        concerns.append("Simple Markov model explains everything -- no deep structure")

    if not concerns:
        print(f"  ALL VALIDATIONS PASS. The finding is robust.")
        print(f"  The Collatz LIS ratio encodes genuine structure beyond")
        print(f"  simple step-distribution or first-order autocorrelation.")
    else:
        print(f"  CONCERNS IDENTIFIED:")
        for c in concerns:
            print(f"    - {c}")
