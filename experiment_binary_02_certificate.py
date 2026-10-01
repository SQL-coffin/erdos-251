from fractions import Fraction


def generate_primes(limit):
    sieve = [True] * (limit + 1)
    sieve[0:2] = [False, False]

    for p in range(2, int(limit ** 0.5) + 1):
        if sieve[p]:
            for multiple in range(p * p, limit + 1, p):
                sieve[multiple] = False

    return [n for n, is_prime in enumerate(sieve) if is_prime]


def generate_prime_gaps(primes):
    return [
        primes[i + 1] - primes[i]
        for i in range(len(primes) - 1)
    ]


def T_finite(gaps, N, K):
    total = Fraction(0, 1)

    for j in range(1, K + 1):
        total += Fraction(gaps[N + j - 1], 2 ** j)

    return total


def tail_bound(K):
    """
    Bound:

        R_K = sum_{j=K+1}^∞ g_(j+1) / 2^j

    using the elementary bound

        g_(j+1) < (j+2)^2.
    """

    bound = Fraction(0, 1)

    # Compute a very safe finite upper bound.
    # Remaining terms after this are negligible.
    for j in range(K + 1, K + 200):
        bound += Fraction((j + 2) ** 2, 2 ** j)

    return bound


def binary_prefix_integer(x, m):
    """
    Return the integer B corresponding to the first m
    binary digits of x in [0,1).
    """
    return (x * (2 ** m)).numerator // (
        x * (2 ** m)
    ).denominator


def certify_prefix(gaps, N, K, m):
    """
    Certify the first m binary digits of T_N.

    We construct an interval containing the true T_N.
    """

    T_K = T_finite(gaps, N, K)

    E = tail_bound(K)

    lower = T_K
    upper = T_K + E

    # We only need the fractional part.
    lower_frac = lower - lower.numerator // lower.denominator
    upper_frac = upper - upper.numerator // upper.denominator

    B = binary_prefix_integer(lower_frac, m)

    interval_lower = Fraction(B, 2 ** m)
    interval_upper = Fraction(B + 1, 2 ** m)

    certified = (
        lower_frac >= interval_lower
        and
        upper_frac < interval_upper
    )

    return {
        "certified": certified,
        "K": K,
        "m": m,
        "prefix_integer": B,
        "interval_lower": interval_lower,
        "interval_upper": interval_upper,
        "lower": lower_frac,
        "upper": upper_frac,
        "error": E,
    }


def binary_string(B, m):
    return format(B, f"0{m}b")


def main():

    print("=" * 60)
    print("BINARY WORKFLOW — EXPERIMENT 02")
    print("Part 1C: Rigorous Binary Prefix Certificate")
    print("=" * 60)

    PRIME_LIMIT = 2_000_000
    N = 1
    K = 100

    primes = generate_primes(PRIME_LIMIT)
    gaps = generate_prime_gaps(primes)

    print(f"Primes generated : {len(primes):,}")
    print(f"K                : {K}")
    print()

    for m in [10, 20, 30, 40, 50]:

        result = certify_prefix(
            gaps,
            N,
            K,
            m
        )

        prefix = binary_string(
            result["prefix_integer"],
            m
        )

        print("-" * 60)
        print(f"m = {m}")
        print(f"Candidate prefix : 0.{prefix}")
        print(f"Certified        : {result['certified']}")
        print()

        if result["certified"]:

            print("CERTIFICATE: PASS")

            print(
                f"Interval: "
                f"[{result['interval_lower']}, "
                f"{result['interval_upper']})"
            )

        else:

            print("CERTIFICATE: FAIL")

        print(
            f"Tail bound E_K = "
            f"{float(result['error']):.3e}"
        )

    print()
    print("=" * 60)
    print("EXPERIMENT COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()