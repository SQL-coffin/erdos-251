"""
Erdos #251 — Experiment 07
Integer-shift non-integrality certificates for
D_h(N) = T_{N+h} - T_N.

Goal:
    If T_1 is rational, then for some h the exact tail shift
    D_h(N) is eventually an integer.

This experiment does NOT prove irrationality. It asks whether finite,
rigorously certified computations can show D_h(N) is non-integral
for many N and h.

Mathematical definitions
------------------------
g_n = p_{n+1} - p_n

T_N = sum_{j>=1} g_{N+j} / 2^j

T_N^(K) = sum_{j=1}^K g_{N+j} / 2^j

D_h(N) = T_{N+h} - T_N

D_h^(K)(N) = T_{N+h}^(K) - T_N^(K)

Tail bound
----------
For m >= 6 we use the standard explicit bound

    p_m < m (log m + log log m) < m^2.

Since g_m = p_{m+1} - p_m < p_{m+1},

    g_{N+j} < (N+j+1)^2

for the ranges used here.

Therefore

    0 < T_N - T_N^(K)
      < sum_{j=K+1}^infinity (N+j+1)^2 / 2^j.

The last sum is evaluated exactly.

If

    dist(D_h^(K)(N), Z) > E_N(K) + E_{N+h}(K),

then D_h(N) cannot be an integer. This is a rigorous finite
certificate, assuming the stated prime upper bound.

Important:
    A finite set of certificates is NOT enough to prove that
    D_h(N) is non-integral infinitely often.
"""

from math import log, log2
from pathlib import Path
import time


# ============================================================
# PARAMETERS
# ============================================================

N_MAX = 100_000
H_MAX = 32
K = 60

print("=" * 80)
print("ERDOS #251 — EXPERIMENT 07")
print("INTEGER-SHIFT NON-INTEGRALITY CERTIFICATES")
print("=" * 80)
print(f"N range : 1 .. {N_MAX}")
print(f"h range : 1 .. {H_MAX}")
print(f"Tail K  : {K}")
print()


# ============================================================
# PRIME GENERATION
# ============================================================

def generate_primes(n):
    """Return the first n primes using a simple sieve."""
    if n < 1:
        return []

    if n < 6:
        limit = 15
    else:
        # Safe practical upper estimate for the requested size.
        limit = int(n * (log(n) + log(log(n)))) + 100

    while True:
        sieve = bytearray(b"\x01") * (limit + 1)
        sieve[0:2] = b"\x00\x00"

        for p in range(2, int(limit ** 0.5) + 1):
            if sieve[p]:
                sieve[p * p : limit + 1 : p] = b"\x00" * (
                    ((limit - p * p) // p) + 1
                )

        primes = [
            i for i in range(2, limit + 1)
            if sieve[i]
        ]

        if len(primes) >= n:
            return primes[:n]

        limit *= 2


needed_primes = N_MAX + H_MAX + K + 5

start = time.time()
primes = generate_primes(needed_primes)
print(f"Generated primes : {len(primes):,}")
print(f"Largest prime    : {primes[-1]:,}")
print(f"Generation time  : {time.time() - start:.3f} s")
print()


# ============================================================
# PRIME GAPS
# ============================================================

gaps = [
    primes[i + 1] - primes[i]
    for i in range(len(primes) - 1)
]


# ============================================================
# EXACT FINITE TAIL NUMERATOR
# ============================================================
#
# T_N^(K) = A_N / 2^K
#
# A_N = sum_{j=1}^K g_{N+j} * 2^(K-j)
#
# We build A_N using the exact recurrence
#
# A_(N+1) = 2 A_N - 2^K g_(N+1) + g_(N+K+1)
#
# so there is no floating-point arithmetic in the certificate.
# ============================================================

DEN = 1 << K
WEIGHT = DEN

A = [0] * (N_MAX + H_MAX + 1)

# Initial numerator for T_1^(K)
numerator = 0

for j in range(1, K + 1):
    numerator += gaps[j] * (1 << (K - j))

A[1] = numerator

for N in range(1, N_MAX + H_MAX):
    A[N + 1] = (
        2 * A[N]
        - DEN * gaps[N]
        + gaps[N + K]
    )


# ============================================================
# EXACT TAIL ERROR BOUND
# ============================================================
#
# E_N(K) =
# sum_{j=K+1}^infinity (N+j+1)^2 / 2^j
#
# Put q = K+1 and A0 = N+1.
#
# sum_{j=q}∞ 1/2^j
#     = 1 / 2^(q-1)
#
# sum_{j=q}∞ j/2^j
#     = (q+1) / 2^(q-1)
#
# sum_{j=q}∞ j^2/2^j
#     = (q^2 + 2q + 3) / 2^(q-1)
#
# Hence E_N(K) has denominator 2^K exactly.
# ============================================================

def tail_bound_numerator(N):
    """
    Return B_N such that

        0 < T_N - T_N^(K) < B_N / 2^K.

    Uses g_m < m^2 for the relevant m.
    """
    q = K + 1
    a0 = N + 1

    return (
        a0 * a0
        + 2 * a0 * (q + 1)
        + (q * q + 2 * q + 3)
    )


# ============================================================
# CERTIFICATE TEST
# ============================================================

def distance_numerator_to_integer(x_num):
    """
    x = x_num / DEN.

    Return numerator of distance(x, Z), i.e.

        min(r, DEN-r)

    where r = x_num mod DEN.
    """
    r = x_num % DEN
    return min(r, DEN - r)


total_tests = 0
certified = 0
failed = 0

# Per-h statistics
results = {}

global_best_margin = -1
global_best = None

start = time.time()

for h in range(1, H_MAX + 1):

    h_certified = 0
    h_failed = 0
    min_margin = None
    worst_case = None
    best_case = None

    max_N = N_MAX - h + 1

    for N in range(1, max_N + 1):

        # Exact finite approximation:
        #
        # D_h^(K)(N)
        #   = (A[N+h] - A[N]) / 2^K
        #
        d_num = A[N + h] - A[N]

        # Exact distance of finite D from nearest integer.
        dist_num = distance_numerator_to_integer(d_num)

        # Rigorous upper bound on the absolute truncation error:
        #
        # |D_h - D_h^(K)|
        #   <= E_N(K) + E_(N+h)(K)
        #
        err_num = (
            tail_bound_numerator(N)
            + tail_bound_numerator(N + h)
        )

        margin_num = dist_num - err_num

        total_tests += 1

        if margin_num > 0:
            certified += 1
            h_certified += 1

            if best_case is None or margin_num > best_case[0]:
                best_case = (
                    margin_num,
                    N,
                    d_num,
                    dist_num,
                    err_num,
                )

            if (
                global_best is None
                or margin_num > global_best_margin
            ):
                global_best_margin = margin_num
                global_best = (
                    h,
                    N,
                    d_num,
                    dist_num,
                    err_num,
                    margin_num,
                )

        else:
            failed += 1
            h_failed += 1

            if worst_case is None or margin_num < worst_case[0]:
                worst_case = (
                    margin_num,
                    N,
                    d_num,
                    dist_num,
                    err_num,
                )

    results[h] = {
        "tested": max_N,
        "certified": h_certified,
        "failed": h_failed,
        "best": best_case,
        "worst": worst_case,
    }

print(f"Certificate computation time: {time.time() - start:.3f} s")
print()


# ============================================================
# RESULTS
# ============================================================

print("=" * 80)
print("RESULTS BY SHIFT h")
print("=" * 80)

print(
    f"{'h':>4} "
    f"{'tested':>10} "
    f"{'certified':>12} "
    f"{'failed':>10} "
    f"{'rate':>10}"
)

print("-" * 55)

for h in range(1, H_MAX + 1):

    r = results[h]

    rate = r["certified"] / r["tested"]

    print(
        f"{h:>4} "
        f"{r['tested']:>10,} "
        f"{r['certified']:>12,} "
        f"{r['failed']:>10,} "
        f"{100 * rate:>9.4f}%"
    )

print()


# ============================================================
# GLOBAL SUMMARY
# ============================================================

print("=" * 80)
print("GLOBAL CERTIFICATE SUMMARY")
print("=" * 80)

print(f"Total tests      : {total_tests:,}")
print(f"Certified        : {certified:,}")
print(f"Not certified    : {failed:,}")

if total_tests:
    print(
        f"Certificate rate : "
        f"{100 * certified / total_tests:.6f}%"
    )

print()


# ============================================================
# STRONGEST CERTIFICATE
# ============================================================

print("=" * 80)
print("STRONGEST CERTIFICATE FOUND")
print("=" * 80)

if global_best is not None:

    h, N, d_num, dist_num, err_num, margin_num = global_best

    print(f"h                 = {h}")
    print(f"N                 = {N}")
    print(f"D_h^(K) numerator = {d_num}")
    print(f"D_h^(K)           = {d_num / DEN:.15f}")
    print(f"distance to Z     = {dist_num / DEN:.15e}")
    print(f"tail error bound  = {err_num / DEN:.15e}")
    print(f"certificate margin= {margin_num / DEN:.15e}")

    print()
    print(
        "Therefore, for this particular (N,h), "
        "D_h(N) is rigorously non-integral."
    )

print()


# ============================================================
# STABILITY CHECK ACROSS K
# ============================================================
#
# We do a smaller sample with several K values.
# This is a numerical stability check, not an additional theorem.
# ============================================================

print("=" * 80)
print("K-STABILITY CHECK")
print("=" * 80)

K_values = [30, 40, 50, 60]

sample_N = min(N_MAX, 20_000)
sample_H = H_MAX

for K_test in K_values:

    den = 1 << K_test

    # Build exact finite numerators for this K.
    B = [0] * (sample_N + sample_H + 1)

    num = 0
    for j in range(1, K_test + 1):
        num += gaps[j] * (1 << (K_test - j))

    B[1] = num

    for N in range(1, sample_N + sample_H):
        B[N + 1] = (
            2 * B[N]
            - den * gaps[N]
            + gaps[N + K_test]
        )

    def bound_num(N):
        q = K_test + 1
        a0 = N + 1
        return (
            a0 * a0
            + 2 * a0 * (q + 1)
            + (q * q + 2 * q + 3)
        )

    cert = 0
    tests = 0

    for h in range(1, sample_H + 1):
        for N in range(1, sample_N - h + 2):

            d_num = B[N + h] - B[N]

            r = d_num % den
            dist = min(r, den - r)

            err = bound_num(N) + bound_num(N + h)

            if dist > err:
                cert += 1

            tests += 1

    print(
        f"K={K_test:>3} | "
        f"certified={cert:,}/{tests:,} | "
        f"rate={100*cert/tests:.6f}%"
    )

print()


# ============================================================
# INTERPRETATION
# ============================================================

print("=" * 80)
print("INTERPRETATION")
print("=" * 80)

print(
    """
What this experiment establishes:

1. For every certified (N,h), the finite computation proves
       D_h(N) != integer
   for the TRUE infinite tail, not merely for T^(K).

2. This is stronger than a floating-point observation because
   the finite value and the error bound are handled using exact
   integer arithmetic.

3. The computation still has a finite range:
       N <= N_MAX,
       h <= H_MAX.
   Therefore it does NOT prove that D_h(N) is non-integral
   infinitely often.

4. To turn this into a genuine partial theorem toward
   Erdos #251, we need a mathematical argument extending the
   observed non-integrality beyond every finite computational
   range, or a theorem forcing infinitely many certificates.

5. The most useful next question is therefore:
       Why are these D_h(N) so consistently close to
       half-integers / non-integers?
   We should inspect the arithmetic structure of
       2^K D_h^(K)(N)
   modulo powers of 2 and modulo odd denominators.
"""
)

print("=" * 80)
print("EXPERIMENT 07 COMPLETE")
print("=" * 80)
