#-*- coding: utf-8 -*#
from pathlib import Path

code = r'''"""
Experiment 09 - Prime-gap half-sequence modulo 2
Erdos #251 research

Goal:
    g_n = p_{n+1} - p_n
    h_n = g_n / 2
    a_n = h_n mod 2 = (g_n / 2) mod 2

This is a deliberately weakened problem:
    Does the binary sequence a_n show eventual periodic structure?

We test:
1. Distribution of 0/1
2. Candidate periods
3. Best-period search
4. Subword complexity p(k)
5. Random binary baseline
6. A simple stability check on prefixes

This experiment is exploratory. It does NOT prove anything
about irrationality of Erdos #251.
"""

import math
import random
from collections import Counter


# ------------------------------------------------------------
# Prime generation
# ------------------------------------------------------------

def sieve(n):
    """Return all primes <= n."""
    is_prime = bytearray(b"\x01") * (n + 1)
    is_prime[0:2] = b"\x00\x00"

    limit = int(math.isqrt(n))
    for p in range(2, limit + 1):
        if is_prime[p]:
            start = p * p
            is_prime[start:n + 1:p] = b"\x00" * (
                ((n - start) // p) + 1
            )

    return [i for i in range(2, n + 1) if is_prime[i]]


def first_n_primes(n):
    """Generate at least n primes using a safe sieve bound."""
    if n < 1:
        return []

    # Rosser-style practical upper bound for n >= 6.
    if n < 6:
        bound = 20
    else:
        bound = int(n * (math.log(n) + math.log(math.log(n)))) + 100

    while True:
        primes = sieve(bound)
        if len(primes) >= n:
            return primes[:n]
        bound *= 2


# ------------------------------------------------------------
# Build a_n = (g_n / 2) mod 2
# ------------------------------------------------------------

def half_gap_binary(num_gaps):
    """
    Return:
        gaps = [g_1, ..., g_num_gaps]
        a    = [h_2 mod 2, ..., h_(num_gaps) mod 2]
             where h_n = g_n / 2.

    g_1 = 1 is exceptional, so we start at g_2 = 2.
    """
    primes = first_n_primes(num_gaps + 1)
    gaps = [
        primes[i + 1] - primes[i]
        for i in range(num_gaps)
    ]

    # Ignore g_1=1; all later gaps are even.
    usable_gaps = gaps[1:]

    a = [(g // 2) % 2 for g in usable_gaps]

    return usable_gaps, a


# ------------------------------------------------------------
# Basic statistics
# ------------------------------------------------------------

def binary_distribution(seq):
    c = Counter(seq)
    n = len(seq)
    return {
        0: c.get(0, 0),
        1: c.get(1, 0),
        "p1": c.get(1, 0) / n if n else float("nan"),
    }


# ------------------------------------------------------------
# Period matching
# ------------------------------------------------------------

def period_match(seq, p):
    """
    Compare a[i] with a[i-p].
    Returns:
        match_rate, compared
    """
    if p <= 0 or p >= len(seq):
        return float("nan"), 0

    matches = sum(
        1 for i in range(p, len(seq))
        if seq[i] == seq[i - p]
    )
    compared = len(seq) - p

    return matches / compared, compared


def best_period(seq, max_period=100):
    results = []

    max_period = min(max_period, len(seq) - 1)

    for p in range(1, max_period + 1):
        rate, compared = period_match(seq, p)
        results.append((rate, p, compared))

    results.sort(reverse=True)
    return results


# ------------------------------------------------------------
# Subword complexity
# ------------------------------------------------------------

def subword_complexity(seq, k):
    """
    Number of distinct binary blocks of length k.
    """
    if k <= 0 or k > len(seq):
        return 0

    return len({
        tuple(seq[i:i+k])
        for i in range(len(seq) - k + 1)
    })


def complexity_table(seq, max_k=20):
    return [
        (k, subword_complexity(seq, k))
        for k in range(1, max_k + 1)
    ]


# ------------------------------------------------------------
# Random binary baseline
# ------------------------------------------------------------

def random_binary(n):
    return [random.getrandbits(1) for _ in range(n)]


def random_period_baseline(n, max_period, trials=200):
    """
    For each random sequence:
    find its best period among 1..max_period.
    """
    best_rates = []

    for _ in range(trials):
        seq = random_binary(n)
        best = best_period(seq, max_period)[0]
        best_rates.append(best[0])

    best_rates.sort()

    mean = sum(best_rates) / len(best_rates)

    def percentile(values, q):
        idx = int(q * (len(values) - 1))
        return values[idx]

    return {
        "mean": mean,
        "median": percentile(best_rates, 0.50),
        "p95": percentile(best_rates, 0.95),
        "p99": percentile(best_rates, 0.99),
        "max": max(best_rates),
        "all": best_rates,
    }


# ------------------------------------------------------------
# Prefix stability
# ------------------------------------------------------------

def prefix_stability(seq, lengths, periods):
    rows = []

    for L in lengths:
        prefix = seq[:L]

        for p in periods:
            rate, compared = period_match(prefix, p)
            rows.append((L, p, rate, compared))

    return rows


# ------------------------------------------------------------
# Main experiment
# ------------------------------------------------------------

def main():
    NUM_GAPS = 100_000
    MAX_PERIOD = 200
    MAX_K = 20

    print("=" * 78)
    print("EXPERIMENT 09 - PRIME-GAP HALF-SEQUENCE MODULO 2")
    print("=" * 78)

    gaps, seq = half_gap_binary(NUM_GAPS)

    print(f"Prime gaps used: {len(gaps)}")
    print(f"Binary sequence length: {len(seq)}")

    # --------------------------------------------------------
    # A. Distribution
    # --------------------------------------------------------
    print("\n" + "-" * 78)
    print("A - DISTRIBUTION OF h_n mod 2")
    print("-" * 78)

    dist = binary_distribution(seq)

    print(f"0 count: {dist[0]:,}")
    print(f"1 count: {dist[1]:,}")
    print(f"P(1):    {dist['p1']:.9f}")
    print(f"P(0):    {1 - dist['p1']:.9f}")

    # --------------------------------------------------------
    # B. Best periods
    # --------------------------------------------------------
    print("\n" + "-" * 78)
    print("B - BEST PERIODS")
    print("-" * 78)

    best = best_period(seq, MAX_PERIOD)

    print("Top 20 candidate periods:")
    print(f"{'period':>8} {'match %':>12} {'compared':>12}")

    for rate, p, compared in best[:20]:
        print(f"{p:>8} {100*rate:>11.6f}% {compared:>12}")

    # --------------------------------------------------------
    # C. Selected periods
    # --------------------------------------------------------
    print("\n" + "-" * 78)
    print("C - SELECTED PERIOD CHECK")
    print("-" * 78)

    selected = [1, 2, 3, 4, 5, 6, 7, 8,
                10, 12, 16, 20, 24, 32, 40,
                50, 64, 80, 100, 128, 160, 200]

    print(f"{'period':>8} {'match %':>12} {'compared':>12}")

    for p in selected:
        rate, compared = period_match(seq, p)
        print(f"{p:>8} {100*rate:>11.6f}% {compared:>12}")

    # --------------------------------------------------------
    # D. Subword complexity
    # --------------------------------------------------------
    print("\n" + "-" * 78)
    print("D - SUBWORD COMPLEXITY")
    print("-" * 78)

    table = complexity_table(seq, MAX_K)

    print(f"{'k':>5} {'p(k)':>12} {'p(k)/k':>12}")

    for k, pk in table:
        print(f"{k:>5} {pk:>12} {pk/k:>12.4f}")

    # --------------------------------------------------------
    # E. Morse-Hedlund diagnostic
    # --------------------------------------------------------
    print("\n" + "-" * 78)
    print("E - MORSE-HEDLUND DIAGNOSTIC")
    print("-" * 78)

    print(
        "For an infinite binary sequence, eventual periodicity has "
        "strong constraints on subword complexity."
    )
    print(
        "We only use this as a diagnostic: finite data cannot prove "
        "non-periodicity."
    )

    for k, pk in table:
        flag = ">" if pk > k else "<="
        print(f"k={k:2d}: p(k)={pk:6d} {flag} k")

    # --------------------------------------------------------
    # F. Random baseline
    # --------------------------------------------------------
    print("\n" + "-" * 78)
    print("F - RANDOM BINARY PERIOD BASELINE")
    print("-" * 78)

    baseline = random_period_baseline(
        n=len(seq),
        max_period=MAX_PERIOD,
        trials=200,
    )

    real_best_rate, real_best_p, _ = best[0]

    print(f"Real best period: {real_best_p}")
    print(f"Real best match:  {100*real_best_rate:.6f}%")

    print(f"Random mean:      {100*baseline['mean']:.6f}%")
    print(f"Random median:    {100*baseline['median']:.6f}%")
    print(f"Random 95th pct:  {100*baseline['p95']:.6f}%")
    print(f"Random 99th pct:  {100*baseline['p99']:.6f}%")
    print(f"Random maximum:   {100*baseline['max']:.6f}%")

    exceed = sum(
        x >= real_best_rate
        for x in baseline["all"]
    )

    empirical_p = (exceed + 1) / (len(baseline["all"]) + 1)

    print(f"Random trials >= real: {exceed}/{len(baseline['all'])}")
    print(f"Empirical p-value:     {empirical_p:.6f}")

    # --------------------------------------------------------
    # G. Prefix stability
    # --------------------------------------------------------
    print("\n" + "-" * 78)
    print("G - PREFIX STABILITY")
    print("-" * 78)

    lengths = [1_000, 2_000, 5_000, 10_000,
               20_000, 50_000, len(seq)]

    candidate_periods = [2, 4, 8, 16, 32, 64, 128]

    rows = prefix_stability(
        seq,
        lengths,
        candidate_periods,
    )

    print(f"{'L':>8} {'p':>6} {'match %':>12} {'compared':>12}")

    for L, p, rate, compared in rows:
        if math.isnan(rate):
            continue
        print(
            f"{L:>8} {p:>6} "
            f"{100*rate:>11.6f}% {compared:>12}"
        )

    # --------------------------------------------------------
    # H. Final interpretation
    # --------------------------------------------------------
    print("\n" + "=" * 78)
    print("EXPERIMENT 09 SUMMARY")
    print("=" * 78)

    print("""
We studied the weakened sequence

    a_n = h_n mod 2
        = (g_n / 2) mod 2.

This is equivalent to observing prime gaps modulo 4.

The experiment asks whether this binary sequence shows
simple eventual-periodic structure.

IMPORTANT:
This experiment is deliberately weaker than Erdős #251.
Even if a_n is non-periodic, that alone does NOT prove
irrationality of T_1, because the full dyadic sum also depends
on higher bits of h_n and on carry propagation.

The useful outcome is therefore:
    data -> candidate structural lemma
rather than:
    data -> proof of irrationality.
""")

    print("Done.")


if __name__ == "__main__":
    main()
'''

path = Path("/mnt/data/experiment_09_half_gap_mod2.py")
path.write_text(code, encoding="utf-8")

print(f"已生成：{path}")
print("运行：python experiment_09_half_gap_mod2.py")
