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
    """
    Exact finite approximation:

        T_N^(K) = sum_{j=1}^K g_(N+j) / 2^j
    """

    total = Fraction(0, 1)

    for j in range(1, K + 1):
        total += Fraction(
            gaps[N + j - 1],
            2 ** j
        )

    return total


def tail_bound(N, K):
    """
    Rigorous upper bound for

        R_N(K)
        =
        T_N - T_N^(K)

    using

        g_(N+j) < (N+j+1)^2.
    """

    bound = Fraction(0, 1)

    for j in range(K + 1, K + 200):

        bound += Fraction(
            (N + j + 1) ** 2,
            2 ** j
        )

    return bound


def distance_to_integer(x):
    """
    Exact distance from rational x to nearest integer.
    """

    floor_x = x.numerator // x.denominator
    fractional = x - floor_x

    return min(
        fractional,
        1 - fractional
    )


def C_finite(gaps, N, P, K):
    """
    Finite approximation to

        C_P(N) = (2^P - 1) T_N.
    """

    T = T_finite(gaps, N, K)

    multiplier = 2 ** P - 1

    return multiplier * T


def certificate_for_C(gaps, N, P, K):
    """
    Test whether we can rigorously certify

        C_P(N) is NOT an integer.
    """

    C_K = C_finite(
        gaps,
        N,
        P,
        K
    )

    raw_tail = tail_bound(
        N,
        K
    )

    multiplier = 2 ** P - 1

    scaled_tail = multiplier * raw_tail

    distance = distance_to_integer(
        C_K
    )

    margin = distance - scaled_tail

    certified = margin > 0

    return {
        "certified": certified,
        "C_K": C_K,
        "distance": distance,
        "tail_error": scaled_tail,
        "margin": margin
    }


def main():

    print("=" * 70)
    print("ERDOS #251 — EXPERIMENT 10")
    print("Fixed-Period Integer-Shift Obstruction")
    print("=" * 70)

    PRIME_LIMIT = 2_000_000
    K = 80
    P_MAX = 32
    N_MAX = 100

    print(f"Prime limit : {PRIME_LIMIT}")
    print(f"K           : {K}")
    print(f"P range     : 1..{P_MAX}")
    print(f"N range     : 1..{N_MAX}")
    print()

    print("Generating primes...")

    primes = generate_primes(
        PRIME_LIMIT
    )

    gaps = generate_prime_gaps(
        primes
    )

    print(
        f"Primes generated : {len(primes):,}"
    )

    print(
        f"Gaps generated   : {len(gaps):,}"
    )

    print()

    print("-" * 70)
    print(
        "P | N | certified | distance | tail error | margin"
    )
    print("-" * 70)

    for P in range(1, P_MAX + 1):

        found = False

        for N in range(1, N_MAX + 1):

            result = certificate_for_C(
                gaps,
                N,
                P,
                K
            )

            if result["certified"]:

                print(
                    f"{P:2d} | "
                    f"{N:3d} | "
                    f"YES | "
                    f"{float(result['distance']):.6e} | "
                    f"{float(result['tail_error']):.6e} | "
                    f"{float(result['margin']):.6e}"
                )

                found = True
                break

        if not found:

            print(
                f"{P:2d} | "
                f"--- | "
                f"NO certificate found"
            )

    print()
    print("=" * 70)
    print("EXPERIMENT COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()