# -*- coding: utf-8 -*-
"""
Experiment 12 - Rationality finite-state grid obstruction

Erdos #251 research.

If T_1 is rational, then its binary fractional-part orbit
is eventually periodic. For a period P, this implies
(2^P - 1) * T_N is an integer eventually.

We therefore test the stronger directly observable condition:
dist((2^P - 1) * T_N, Z) > tail certificate.

No floating-point arithmetic is used.
"""

from fractions import Fraction
import math


def generate_primes(limit):
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, int(math.isqrt(limit)) + 1):
        if sieve[p]:
            sieve[p*p:limit+1:p] = b"\x00" * (
                ((limit - p*p) // p) + 1
            )
    return [n for n in range(2, limit + 1) if sieve[n]]


def generate_gaps(primes):
    return [primes[i+1] - primes[i] for i in range(len(primes)-1)]


def compute_T(gaps, N, K):
    total = Fraction(0, 1)
    for j in range(1, K + 1):
        total += Fraction(gaps[N + j - 1], 2**j)
    return total


def scan_T_fast(gaps, N_start, N_end, K):
    T = compute_T(gaps, N_start, K)
    results = [T]
    for N in range(N_start, N_end):
        T = (
            2 * T
            - gaps[N]
            + Fraction(gaps[N + K], 2**K)
        )
        results.append(T)
    return results


def distance_to_integer(x):
    lower = x.numerator // x.denominator
    return min(
        abs(x - lower),
        abs((lower + 1) - x)
    )


def tail_bound(N, K):
    r = K + 1
    a = N + 1

    sum_j2 = Fraction(
        2 * (r**2 + 2*r + 3),
        2**r
    )
    sum_j = Fraction(
        2 * (r + 1),
        2**r
    )
    sum_1 = Fraction(2, 2**r)

    return sum_j2 + 2*a*sum_j + a*a*sum_1


# ------------------------------------------------------------
# TDD TESTS
# ------------------------------------------------------------

def test_grid_equivalence():
    T = Fraction(17, 13)
    q = 7

    # Integer part does not matter:
    # q*T integer iff q*frac(T) integer.
    frac = T - (T.numerator // T.denominator)

    assert (q*T).denominator == 13
    assert q*frac == q*T - q*(T.numerator // T.denominator)

    print("TEST 12A: Grid equivalence")
    print("PASS")


def test_fast_recurrence():
    primes = generate_primes(5000)
    gaps = generate_gaps(primes)

    K = 20
    fast = scan_T_fast(gaps, 1, 50, K)

    for i, N in enumerate(range(1, 51)):
        assert fast[i] == compute_T(gaps, N, K)

    print("TEST 12B: Fast recurrence")
    print("PASS")


def test_tail_certificate():
    primes = generate_primes(5000)
    gaps = generate_gaps(primes)

    N = 1
    K = 20
    P = 3
    q = 2**P - 1

    T_K = compute_T(gaps, N, K)
    dist = distance_to_integer(q*T_K)
    error = q * tail_bound(N, K)

    assert dist > error

    print("TEST 12C: Rigorous grid certificate")
    print("PASS")


def scan_periods(gaps, T_values, P_values, N_start, N_end, K):
    rows = []
    for P in P_values:
        q = 2**P - 1
        certified = 0
        closest = None

        for i, N in enumerate(range(N_start, N_end + 1)):
            value = q * T_values[i]
            dist = distance_to_integer(value)
            error = q * tail_bound(N, K)

            is_certified = dist > error
            certified += is_certified

            if closest is None or dist < closest[0]:
                closest = (dist, N, error)

        rows.append({
            "P": P,
            "q": q,
            "certified_count": certified,
            "closest": closest,
        })

    return rows


def main():
    test_grid_equivalence()
    test_fast_recurrence()
    test_tail_certificate()

    N_start = 1
    N_end = 100_000
    K = 60
    P_values = range(1, 33)

    print("=" * 72)
    print("EXPERIMENT 12 - FINITE-STATE GRID OBSTRUCTION")
    print("=" * 72)
    print(f"N range: {N_start}..{N_end}")
    print(f"K = {K}")
    print(f"P range: {min(P_values)}..{max(P_values)}")

    primes = generate_primes(2_000_000)
    gaps = generate_gaps(primes)

    T_values = scan_T_fast(
        gaps,
        N_start,
        N_end,
        K
    )

    rows = scan_periods(
        gaps,
        T_values,
        P_values,
        N_start,
        N_end,
        K
    )

    print()
    print(f"{'P':>3} {'q=2^P-1':>12} {'certified':>12}"
          f" {'closest N':>12} {'distance':>18} {'tail':>18}")

    for row in rows:
        dist, N, error = row["closest"]
        print(
            f"{row['P']:>3} "
            f"{row['q']:>12} "
            f"{row['certified_count']:>12} "
            f"{N:>12} "
            f"{float(dist):>18.10e} "
            f"{float(error):>18.10e}"
        )

    print()
    print("Interpretation:")
    print("A 100% certified row means every tested N has")
    print("(2^P - 1) T_N rigorously non-integral.")
    print("This is a finite obstruction only; P is also finite.")


if __name__ == "__main__":
    main()
