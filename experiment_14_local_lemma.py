from fractions import Fraction
import math

"""
Experiment 14 - Local dyadic lemma

Define
    Delta_j = g_{N+h+j} - g_{N+j}
    c_j = Delta_j / 2
    D_h(N) = sum_{j>=1} Delta_j / 2^j
           = sum_{j>=1} c_j / 2^(j-1).

Target interval:
    1/2 < D_h(N) < 1

Sufficient prefix pattern:
    c_1,c_2,c_3 = (1,-1,1)
gives
    D^(3) = 3/4.

If the rigorously bounded remaining tail has absolute value < 1/4,
then the infinite D lies in (1/2,1).

The negative mirror pattern
    (-1,1,-1)
gives D^(3) = -3/4 and targets (-1,-1/2).
"""


def generate_primes(limit):
    sieve = bytearray([1]) * (limit + 1)
    sieve[0:2] = bytearray([0, 0])

    for p in range(2, int(math.isqrt(limit)) + 1):
        if sieve[p]:
            sieve[p*p:limit+1:p] = bytearray(
                ((limit - p*p)//p) + 1
            )

    return [n for n in range(2, limit + 1) if sieve[n]]


def generate_gaps(primes):
    return [
        primes[i+1] - primes[i]
        for i in range(len(primes) - 1)
    ]


def compute_T(gaps, N, K):
    total = Fraction(0, 1)

    for j in range(1, K + 1):
        total += Fraction(
            gaps[N + j - 1],
            2**j
        )

    return total


def scan_T_fast(gaps, N_start, N_end, K):
    T = compute_T(gaps, N_start, K)
    values = [T]

    for N in range(N_start, N_end):
        T = (
            2*T
            - gaps[N]
            + Fraction(gaps[N + K], 2**K)
        )
        values.append(T)

    return values


def tail_bound_difference(N, h, K):
    """
    Bound the infinite tail of

        D_h(N) - D_h^(K)(N)

    using
        g_m < (m+1)^2.
    """
    r = K + 1
    a = N + 1

    # Sum_{j=r}^inf (a+j)^2 / 2^j
    # = Sum j^2/2^j + 2a Sum j/2^j + a^2 Sum 1/2^j.
    def square_tail(x):
        sum_j2 = Fraction(
            2*(r*r + 2*r + 3),
            2**r
        )
        sum_j = Fraction(
            2*(r + 1),
            2**r
        )
        sum_1 = Fraction(2, 2**r)

        return (
            sum_j2
            + 2*x*sum_j
            + x*x*sum_1
        )

    return (
        square_tail(a)
        + square_tail(a + h)
    )


def finite_D(gaps, N, h, K):
    total = Fraction(0, 1)

    for j in range(1, K + 1):
        delta = gaps[N + h + j - 1] - gaps[N + j - 1]
        total += Fraction(delta, 2**j)

    return total


# ------------------------------------------------------------
# TDD TESTS
# ------------------------------------------------------------

def test_pattern_value():
    c = [1, -1, 1]

    D3 = sum(
        Fraction(c[j - 1], 2**(j - 1))
        for j in range(1, 4)
    )

    assert D3 == Fraction(3, 4)

    c = [-1, 1, -1]

    D3 = sum(
        Fraction(c[j - 1], 2**(j - 1))
        for j in range(1, 4)
    )

    assert D3 == Fraction(-3, 4)

    print("TEST 14A: Local pattern value")
    print("PASS")


def test_tail_bound():
    primes = generate_primes(5000)
    gaps = generate_gaps(primes)

    N = 10
    h = 2
    K = 20

    Dk = finite_D(gaps, N, h, K)
    err = tail_bound_difference(N, h, K)

    # The certificate must at least be non-negative.
    assert err > 0
    assert isinstance(Dk, Fraction)

    print("TEST 14B: Difference tail bound")
    print("PASS")


def test_mirror_interval():
    positive_prefix = Fraction(3, 4)
    negative_prefix = Fraction(-3, 4)

    assert (
        Fraction(1, 2)
        < positive_prefix
        < Fraction(1, 1)
    )

    assert (
        Fraction(-1, 1)
        < negative_prefix
        < Fraction(-1, 2)
    )

    print("TEST 14C: Target intervals")
    print("PASS")


def scan_pattern(gaps, T_values, N_start, N_end, h, K):
    results = []

    for N in range(N_start, N_end + 1):
        c = [
            (
                gaps[N + h + j - 1]
                - gaps[N + j - 1]
            ) // 2
            for j in range(1, 4)
        ]

        pattern = tuple(c)

        if pattern not in {
            (1, -1, 1),
            (-1, 1, -1),
        }:
            continue

        Dk = finite_D(gaps, N, h, K)
        error = tail_bound_difference(N, h, K)

        if pattern == (1, -1, 1):
            certified = (
                Fraction(1, 2)
                < Dk - error
                and Dk + error < Fraction(1, 1)
            )
            target = "positive"
        else:
            certified = (
                Fraction(-1, 1)
                < Dk - error
                and Dk + error < Fraction(-1, 2)
            )
            target = "negative"

        results.append({
            "N": N,
            "h": h,
            "pattern": pattern,
            "Dk": Dk,
            "error": error,
            "certified": certified,
            "target": target,
        })

    return results


def main():
    test_pattern_value()
    test_tail_bound()
    test_mirror_interval()

    N_start = 1
    N_end = 100_000
    H_max = 32
    K = 100

    T_end = N_end + H_max + 1

    print("=" * 78)
    print("EXPERIMENT 14 - LOCAL DYADIC LEMMA")
    print("=" * 78)
    print(f"N range: {N_start}..{N_end}")
    print(f"h range: 1..{H_max}")
    print(f"K = {K}")

    primes = generate_primes(2_000_000)
    gaps = generate_gaps(primes)

    T_values = scan_T_fast(
        gaps,
        N_start,
        T_end,
        K
    )

    total_candidates = 0
    total_certified = 0

    print()
    print(
        f"{'h':>3} {'pattern':>14} "
        f"{'candidates':>12} {'certified':>12}"
    )

    for h in range(1, H_max + 1):
        rows = scan_pattern(
            gaps,
            T_values,
            N_start,
            N_end,
            h,
            K
        )

        positive = [
            x for x in rows
            if x["pattern"] == (1, -1, 1)
        ]

        negative = [
            x for x in rows
            if x["pattern"] == (-1, 1, -1)
        ]

        for label, subset in [
            ("(+,-,+)", positive),
            ("(-,+,-)", negative),
        ]:
            certified = sum(
                x["certified"]
                for x in subset
            )

            print(
                f"{h:>3} {label:>14} "
                f"{len(subset):>12} "
                f"{certified:>12}"
            )

            total_candidates += len(subset)
            total_certified += certified

    print()
    print(f"Total pattern candidates: {total_candidates}")
    print(f"Total certified:          {total_certified}")

    print()
    print("LEMMA:")
    print("If c1,c2,c3 = (1,-1,1) and the")
    print("rigorous remaining tail is < 1/4 in absolute value,")
    print("then 1/2 < D_h(N) < 1.")
    print()
    print("Mirror:")
    print("If c1,c2,c3 = (-1,1,-1) and the")
    print("rigorous remaining tail is < 1/4 in absolute value,")
    print("then -1 < D_h(N) < -1/2.")
    print()
    print("These are sufficient local certificates, not")
    print("claims that such patterns occur infinitely often.")


if __name__ == "__main__":
    main()
