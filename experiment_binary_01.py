from fractions import Fraction


# ============================================================
# BINARY WORKFLOW — EXPERIMENT 01
# Part 1: Binary fractional dynamics
#
# Goal:
#   Verify the exact relation
#
#       f_{N+1} = {2 f_N}
#
# where
#
#       f_N = fractional_part(T_N)
#
#       T_N = sum_{j>=1} g_{N+j} / 2^j
#
# We use finite dyadic approximations only for computation.
# ============================================================


def generate_primes(limit):
    """Generate primes up to limit using the sieve of Eratosthenes."""
    sieve = [True] * (limit + 1)
    sieve[0:2] = [False, False]

    for p in range(2, int(limit ** 0.5) + 1):
        if sieve[p]:
            for multiple in range(p * p, limit + 1, p):
                sieve[multiple] = False

    return [n for n, is_prime in enumerate(sieve) if is_prime]


def generate_prime_gaps(primes):
    """Return g_n = p_(n+1) - p_n."""
    return [
        primes[i + 1] - primes[i]
        for i in range(len(primes) - 1)
    ]


def T_finite(gaps, N, K):
    """
    Exact finite approximation:

        T_N^(K) = sum_{j=1}^K g_{N+j} / 2^j

    N uses 1-based mathematical indexing.
    """
    total = Fraction(0, 1)

    for j in range(1, K + 1):
        total += Fraction(gaps[N + j - 1], 2 ** j)

    return total


def fractional_part(x):
    """Exact fractional part."""
    return x - x.numerator // x.denominator


def binary_fraction(x, digits=32):
    """
    Convert a Fraction in [0,1) into its first `digits`
    binary fractional digits.
    """
    bits = []

    for _ in range(digits):
        x *= 2

        if x >= 1:
            bits.append("1")
            x -= 1
        else:
            bits.append("0")

    return "".join(bits)

def test_doubling_map(gaps, N, K):
    """
    Test the exact finite-truncation identity:

        T_(N+1)^(K)
        =
        2 T_N^(K)
        - g_(N+1)
        + g_(N+K+1) / 2^K
    """

    T_N = T_finite(gaps, N, K)
    T_next = T_finite(gaps, N + 1, K)

    g_next = gaps[N]
    g_tail = gaps[N + K]

    # Exact finite-truncation correction term
    correction = Fraction(g_tail, 2 ** K)

    # Correct finite identity
    predicted_T_next = (
        2 * T_N
        - g_next
        + correction
    )

    identity_passed = (predicted_T_next == T_next)

    # Fractional parts
    f_N = fractional_part(T_N)
    f_next = fractional_part(T_next)

    predicted_fractional = fractional_part(
        predicted_T_next
    )

    return (
        identity_passed,
        f_N,
        f_next,
        predicted_fractional,
        correction
    )


def test_binary_shift(gaps, N, K, digits=32):
    """
    Test the binary left-shift relation.

        binary(f_(N+1))
        =
        binary(f_N)[1:]

    """
    T_N = T_finite(gaps, N, K)
    T_next = T_finite(gaps, N + 1, K)

    f_N = fractional_part(T_N)
    f_next = fractional_part(T_next)

    bits_N = binary_fraction(f_N, digits)
    bits_next = binary_fraction(f_next, digits)

    expected = bits_N[1:] + "?"

    # Compare only the first digits that are meaningful
    # after the shift.
    passed = bits_next[:-1] == bits_N[1:]

    return passed, bits_N, bits_next, expected


def main():

    print("=" * 60)
    print("BINARY WORKFLOW — EXPERIMENT 01")
    print("Part 1: Binary fractional dynamics")
    print("=" * 60)

    # --------------------------------------------------------
    # Parameters
    # --------------------------------------------------------

    PRIME_LIMIT = 2_000_000
    K = 60
    TEST_COUNT = 100
    BINARY_DIGITS = 32

    print(f"Prime limit : {PRIME_LIMIT}")
    print(f"K           : {K}")
    print(f"Tests       : {TEST_COUNT}")
    print()

    # --------------------------------------------------------
    # Generate prime gaps
    # --------------------------------------------------------

    print("Generating primes...")

    primes = generate_primes(PRIME_LIMIT)
    gaps = generate_prime_gaps(primes)

    print(f"Primes generated : {len(primes):,}")
    print(f"Gaps generated   : {len(gaps):,}")
    print()

    # --------------------------------------------------------
    # Test 1 — Doubling map
    # --------------------------------------------------------

    print("-" * 60)
    print("TEST 1 — DOUBLING MAP")
    print("-" * 60)

    passed_count = 0

    for N in range(1, TEST_COUNT + 1):

        passed, f_N, f_next, predicted, correction = test_doubling_map(
            gaps,
            N,
            K
        )

        if passed:
            passed_count += 1
        else:
            print(f"FAIL at N = {N}")
            print(f"f_N       = {f_N}")
            print(f"f_(N+1)   = {f_next}")
            print(f"predicted = {predicted}")

    print(f"Passed: {passed_count}/{TEST_COUNT}")

    if passed_count == TEST_COUNT:
        print("STATUS: PASS")
    else:
        print("STATUS: FAIL")

    # --------------------------------------------------------
    # Test 2 — Binary shift
    # --------------------------------------------------------

    print()
    print("-" * 60)
    print("TEST 2 — BINARY LEFT SHIFT")
    print("-" * 60)

    passed_count = 0

    for N in range(1, TEST_COUNT + 1):

        passed, bits_N, bits_next, expected = test_binary_shift(
            gaps,
            N,
            K,
            BINARY_DIGITS
        )

        if passed:
            passed_count += 1
        else:
            print(f"FAIL at N = {N}")
            print(f"f_N     = 0.{bits_N}")
            print(f"f_(N+1) = 0.{bits_next}")

    print(f"Passed: {passed_count}/{TEST_COUNT}")

    if passed:
        passed_count += 1
    else:
        print(f"FAIL at N = {N}")
        print(f"f_N             = {f_N}")
        print(f"f_(N+1)         = {f_next}")
        print(f"predicted       = {predicted}")
        print(f"tail correction = {correction}")

    # --------------------------------------------------------
    # Test 3 — Binary prefix stability
    # --------------------------------------------------------

    print()
    print("-" * 60)
    print("TEST 3 — BINARY PREFIX STABILITY")
    print("-" * 60)

    N = 1
    K_VALUES = [20, 40, 60, 80, 100]
    PREFIX_LENGTHS = [10, 20, 30, 40]

    binary_results = {}

    for K_test in K_VALUES:
        T = T_finite(gaps, N, K_test)
        f = fractional_part(T)

        binary_results[K_test] = binary_fraction(
            f,
            max(PREFIX_LENGTHS)
        )

    for prefix_length in PREFIX_LENGTHS:

        print()
        print(f"Prefix length: {prefix_length} bits")

        reference = binary_results[K_VALUES[-1]][:prefix_length]

        stable = True

        for K_test in K_VALUES:

            current = binary_results[K_test][:prefix_length]

            print(
                f"K={K_test:3d} : "
                f"0.{current}"
            )

            if current != reference:
                stable = False

        if stable:
            print("STATUS: STABLE")
        else:
            print("STATUS: NOT STABLE")    
    # --------------------------------------------------------
    # Example
    # --------------------------------------------------------

    print()
    print("-" * 60)
    print("EXAMPLE")
    print("-" * 60)

    N = 1

    T_N = T_finite(gaps, N, K)
    T_next = T_finite(gaps, N + 1, K)

    f_N = fractional_part(T_N)
    f_next = fractional_part(T_next)

    bits_N = binary_fraction(f_N, BINARY_DIGITS)
    bits_next = binary_fraction(f_next, BINARY_DIGITS)

    print(f"T_{N}^(K)     = {T_N}")
    print(f"f_{N}         = {f_N}")
    print(f"binary       = 0.{bits_N}")
    print()

    print(f"T_{N+1}^(K)   = {T_next}")
    print(f"f_{N+1}       = {f_next}")
    print(f"binary       = 0.{bits_next}")

    print()
    print("=" * 60)
    print("EXPERIMENT COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()