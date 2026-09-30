"""
Erdos #251 — Experiment 08
2-adic structure of the integer-shift numerator

We study

    T_N^(K) = A_N / 2^K

and

    D_h^(K)(N) = T_{N+h}^(K) - T_N^(K)
                = (A_{N+h} - A_N) / 2^K.

Experiment 07 found that D_h(N) was rigorously non-integral for
all tested N <= 100000 and h <= 32, with many values extremely
close to half-integers.

This experiment asks:

    What is the 2-adic structure of
        Delta_A = A_{N+h} - A_N ?

In particular we measure:

    v2(Delta_A)
    Delta_A / 2^(K-1) mod 2
    residue Delta_A mod 2^r
    distance of D_h^(K)(N) to the nearest half-integer.

Important:
    This is an exploratory computational experiment.
    It does NOT prove irrationality of Erdős #251.
"""

from math import log
import time


# ============================================================
# PARAMETERS
# ============================================================

N_MAX = 100_000
H_MAX = 32
K = 60

print("=" * 80)
print("ERDOS #251 — EXPERIMENT 08")
print("2-ADIC STRUCTURE OF INTEGER-SHIFT NUMERATORS")
print("=" * 80)
print(f"N range : 1 .. {N_MAX}")
print(f"h range : 1 .. {H_MAX}")
print(f"K       : {K}")
print()


# ============================================================
# PRIME GENERATION
# ============================================================

def generate_primes(n):
    if n < 6:
        limit = 15
    else:
        limit = int(n * (log(n) + log(log(n)))) + 100

    while True:
        sieve = bytearray(b"\x01") * (limit + 1)
        sieve[0:2] = b"\x00\x00"

        for p in range(2, int(limit ** 0.5) + 1):
            if sieve[p]:
                sieve[p * p : limit + 1 : p] = b"\x00" * (
                    ((limit - p * p) // p) + 1
                )

        primes = [i for i in range(2, limit + 1) if sieve[i]]

        if len(primes) >= n:
            return primes[:n]

        limit *= 2


needed = N_MAX + H_MAX + K + 5

start = time.time()
primes = generate_primes(needed)
print(f"Generated {len(primes):,} primes")
print(f"Time: {time.time() - start:.3f} s")
print()

gaps = [
    primes[i + 1] - primes[i]
    for i in range(len(primes) - 1)
]


# ============================================================
# EXACT DYADIC NUMERATORS A_N
# ============================================================

DEN = 1 << K

A = [0] * (N_MAX + H_MAX + 1)

num = 0

for j in range(1, K + 1):
    num += gaps[j] * (1 << (K - j))

A[1] = num

for N in range(1, N_MAX + H_MAX):
    A[N + 1] = (
        2 * A[N]
        - DEN * gaps[N]
        + gaps[N + K]
    )


# ============================================================
# v2 FUNCTION
# ============================================================

def v2(x):
    """
    2-adic valuation:
        v2(x) = largest r such that 2^r | x.

    For x = 0, return infinity.
    """
    if x == 0:
        return float("inf")

    x = abs(x)
    return (x & -x).bit_length() - 1


def distance_to_half_numerator(x):
    """
    x is the numerator of a dyadic number x / 2^K.

    Distance from x/2^K to the nearest half-integer is

        | x/2^K - (m + 1/2) |

    minimized over integer m.

    In numerator units, this is distance from x to
        (2m+1) * 2^(K-1).

    Returns the numerator of that distance.
    """
    half = 1 << (K - 1)

    r = x % (1 << K)

    # Odd residue relative to the half-integer grid.
    d = abs(r - half)

    # Periodic modulo 2^K.
    d = min(d, (1 << K) - d)

    return d


# ============================================================
# BASIC SAMPLE
# ============================================================

print("=" * 80)
print("SAMPLE VALUES")
print("=" * 80)

for h in [1, 2, 3, 5, 10, 20, 27, 32]:

    N = 31233

    delta = A[N + h] - A[N]

    value = delta / DEN
    frac = value - int(value)

    valuation = v2(delta)
    half_dist = distance_to_half_numerator(delta) / DEN

    print(
        f"h={h:2d} | "
        f"D={value: .12f} | "
        f"v2(Delta_A)={valuation:2} | "
        f"dist_half={half_dist:.3e}"
    )

print()


# ============================================================
# EXPERIMENT A
# DISTRIBUTION OF v2(Delta_A)
# ============================================================

print("=" * 80)
print("EXPERIMENT A — DISTRIBUTION OF v2(Delta_A)")
print("=" * 80)

valuation_counts = {}

total = 0

for h in range(1, H_MAX + 1):

    for N in range(1, N_MAX - h + 2):

        delta = A[N + h] - A[N]
        val = v2(delta)

        valuation_counts[val] = valuation_counts.get(val, 0) + 1
        total += 1

print(f"Total samples: {total:,}")
print()

print(
    f"{'v2':>8} "
    f"{'count':>14} "
    f"{'percentage':>12}"
)

print("-" * 40)

for val in sorted(valuation_counts, key=lambda x: (float("inf") if x == float("inf") else x)):

    count = valuation_counts[val]

    print(
        f"{str(val):>8} "
        f"{count:>14,} "
        f"{100 * count / total:>11.6f}%"
    )

print()


# ============================================================
# EXPERIMENT B
# v2 BY SHIFT h
# ============================================================

print("=" * 80)
print("EXPERIMENT B — v2 BY SHIFT h")
print("=" * 80)

print(
    f"{'h':>4} "
    f"{'mean v2':>12} "
    f"{'median-ish':>12} "
    f"{'v2=K-1':>12} "
    f"{'v2>=K-1':>12}"
)

print("-" * 65)

for h in range(1, H_MAX + 1):

    vals = []

    exact_kminus1 = 0
    at_least_kminus1 = 0

    for N in range(1, N_MAX - h + 2):

        delta = A[N + h] - A[N]
        val = v2(delta)

        vals.append(val)

        if val == K - 1:
            exact_kminus1 += 1

        if val >= K - 1:
            at_least_kminus1 += 1

    vals_sorted = sorted(vals)

    mean_v = sum(vals) / len(vals)
    median_v = vals_sorted[len(vals_sorted) // 2]

    print(
        f"{h:>4} "
        f"{mean_v:>12.4f} "
        f"{median_v:>12} "
        f"{exact_kminus1:>12,} "
        f"{at_least_kminus1:>12,}"
    )

print()


# ============================================================
# EXPERIMENT C
# DISTANCE TO HALF-INTEGER
# ============================================================

print("=" * 80)
print("EXPERIMENT C — DISTANCE TO HALF-INTEGER")
print("=" * 80)

thresholds = [
    1e-1,
    1e-2,
    1e-3,
    1e-4,
    1e-5,
    1e-6,
    1e-8,
    1e-10,
]

half_counts = {t: 0 for t in thresholds}
sample_count = 0

min_dist = None
min_record = None

for h in range(1, H_MAX + 1):

    for N in range(1, N_MAX - h + 2):

        delta = A[N + h] - A[N]

        dist = distance_to_half_numerator(delta) / DEN

        sample_count += 1

        for t in thresholds:
            if dist < t:
                half_counts[t] += 1

        if min_dist is None or dist < min_dist:

            min_dist = dist

            min_record = (
                N,
                h,
                delta,
                dist,
                v2(delta)
            )

print(f"Total samples: {sample_count:,}")
print()

print(
    f"{'threshold':>14} "
    f"{'count':>14} "
    f"{'percentage':>14}"
)

print("-" * 48)

for t in thresholds:

    c = half_counts[t]

    print(
        f"{t:>14.1e} "
        f"{c:>14,} "
        f"{100*c/sample_count:>13.6f}%"
    )

print()

print("Closest half-integer found:")
print(
    f"N={min_record[0]}, "
    f"h={min_record[1]}, "
    f"D_h^(K)={min_record[2]/DEN:.15f}, "
    f"distance={min_record[3]:.3e}, "
    f"v2={min_record[4]}"
)

print()


# ============================================================
# EXPERIMENT D
# MODULAR RESIDUES
# ============================================================

print("=" * 80)
print("EXPERIMENT D — RESIDUES MODULO 2^r")
print("=" * 80)

r_values = [4, 8, 12, 16, 20]

for r in r_values:

    modulus = 1 << r

    residue_counts = {}

    # Use a representative h set so this remains readable.
    for h in [1, 2, 3, 5, 10, 20, 27, 32]:

        counts = {}

        for N in range(1, N_MAX - h + 2):

            delta = A[N + h] - A[N]

            residue = delta % modulus

            counts[residue] = counts.get(residue, 0) + 1

        residue_counts[h] = counts

    print()
    print(f"Modulo 2^{r}:")
    print()

    for h in [1, 2, 3, 5, 10, 20, 27, 32]:

        counts = residue_counts[h]

        unique = len(counts)

        top = sorted(
            counts.items(),
            key=lambda x: x[1],
            reverse=True
        )[:8]

        print(
            f"h={h:2d} | "
            f"unique residues={unique:5d} | "
            f"top={top}"
        )

print()


# ============================================================
# EXPERIMENT E
# NORMALIZED PARITY
# ============================================================
#
# If v2(Delta_A) = K-1, then
#
#     Delta_A / 2^(K-1)
#
# is odd.
#
# This corresponds to D_h^(K) lying exactly on the
# half-integer grid modulo 1.
# ============================================================

print("=" * 80)
print("EXPERIMENT E — NORMALIZED PARITY")
print("=" * 80)

for h in range(1, H_MAX + 1):

    count_v_kminus1 = 0
    count_v_ge_kminus1 = 0
    count_odd_normalized = 0

    for N in range(1, N_MAX - h + 2):

        delta = A[N + h] - A[N]

        val = v2(delta)

        if val == K - 1:
            count_v_kminus1 += 1

            normalized = delta >> (K - 1)

            if normalized & 1:
                count_odd_normalized += 1

        if val >= K - 1:
            count_v_ge_kminus1 += 1

    print(
        f"h={h:2d} | "
        f"v2=K-1: {count_v_kminus1:8,} | "
        f"v2>=K-1: {count_v_ge_kminus1:8,} | "
        f"odd normalized: {count_odd_normalized:8,}"
    )

print()


# ============================================================
# EXPERIMENT F
# SEARCH FOR h-DEPENDENT PATTERN
# ============================================================

print("=" * 80)
print("EXPERIMENT F — h-DEPENDENT 2-ADIC PATTERN")
print("=" * 80)

print(
    "For each h we compute the average valuation and the "
    "most common v2 value."
)
print()

for h in range(1, H_MAX + 1):

    counts = {}

    for N in range(1, N_MAX - h + 2):

        val = v2(A[N + h] - A[N])

        counts[val] = counts.get(val, 0) + 1

    mode_val, mode_count = max(
        counts.items(),
        key=lambda x: x[1]
    )

    print(
        f"h={h:2d} | "
        f"mode v2={mode_val:2} "
        f"({100*mode_count/(N_MAX-h+1):.3f}%)"
    )

print()


# ============================================================
# EXPERIMENT G
# DIRECT TEST OF HALF-INTEGER PHENOMENON
# ============================================================

print("=" * 80)
print("EXPERIMENT G — HALF-INTEGER PHENOMENON BY h")
print("=" * 80)

print(
    f"{'h':>4} "
    f"{'mean dist':>14} "
    f"{'median dist':>14} "
    f"{'<1e-3':>10} "
    f"{'<1e-6':>10}"
)

print("-" * 60)

for h in range(1, H_MAX + 1):

    dists = []

    c1 = 0
    c2 = 0

    for N in range(1, N_MAX - h + 2):

        delta = A[N + h] - A[N]

        dist = distance_to_half_numerator(delta) / DEN

        dists.append(dist)

        if dist < 1e-3:
            c1 += 1

        if dist < 1e-6:
            c2 += 1

    dists.sort()

    mean_d = sum(dists) / len(dists)
    median_d = dists[len(dists)//2]

    print(
        f"{h:>4} "
        f"{mean_d:>14.6e} "
        f"{median_d:>14.6e} "
        f"{c1:>10,} "
        f"{c2:>10,}"
    )

print()


# ============================================================
# FINAL INTERPRETATION
# ============================================================

print("=" * 80)
print("INTERPRETATION")
print("=" * 80)

print(
    """
The main questions to inspect are:

1. Does v2(A[N+h]-A[N]) concentrate at a small value?

2. Does the valuation depend systematically on h?

3. Is v2(A[N+h]-A[N]) often equal to K-1?
   If so, D_h^(K)(N) lies exactly on a half-integer
   modulo 1.

4. If not exactly K-1, does the residue still concentrate
   strongly near the half-integer residue modulo 2^r?

5. Does the pattern persist for K=40, 50, 60?
   A genuine arithmetic pattern should survive changes in K.

6. Most importantly:
   can the observed 2-adic pattern be rewritten as a theorem
   about the prime-gap sequence rather than merely a numerical
   observation?

Do NOT interpret this experiment as a proof of irrationality.
The goal is to discover a lemma that could eventually connect
prime-gap arithmetic to the rationality obstruction.
"""
)

print("=" * 80)
print("EXPERIMENT 08 COMPLETE")
print("=" * 80)
