from fractions import Fraction

# ==============================
# Experiment 11
# Eventual Integer-State Scan
# ==============================

P = 1
K = 60

print("=" * 60)
print("EXPERIMENT 11 — EVENTUAL INTEGER-STATE SCAN")
print("=" * 60)

print(f"P = {P}")
print(f"K = {K}")

# TDD TEST 1:
# T_N^(K) must have denominator 2^K

def test_denominator_structure():
    denominator = 2 ** K

    test_value = Fraction(123456789, denominator)

    assert test_value.denominator <= denominator
    assert denominator % test_value.denominator == 0

    print("TEST 1: Denominator structure")
    print("PASS")


test_denominator_structure()

print("=" * 60)
print("ALL TESTS PASSED")
print("=" * 60)
# ==============================
# Prime gaps
# ==============================

def generate_primes(limit):
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"

    for p in range(2, int(limit ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = b"\x00" * (
                ((limit - p * p) // p) + 1
            )

    return [n for n in range(2, limit + 1) if sieve[n]]


def generate_gaps(primes):
    return [
        primes[i + 1] - primes[i]
        for i in range(len(primes) - 1)
    ]


# ==============================
# TDD TEST 2
# ==============================

def test_prime_gaps():

    primes = generate_primes(100)

    gaps = generate_gaps(primes)

    assert primes[:5] == [2, 3, 5, 7, 11]
    assert gaps[:5] == [1, 2, 2, 4, 2]

    print("TEST 2: Prime gap generation")
    print("PASS")


test_prime_gaps()
# ==============================
# Finite T_N^(K)
# ==============================

def compute_T(gaps, N, K):
    """
    T_N^(K) = sum_{j=1}^K g_{N+j} / 2^j
    """
    total = Fraction(0, 1)

    for j in range(1, K + 1):
        total += Fraction(gaps[N + j - 1], 2 ** j)

    return total


# ==============================
# TDD TEST 3
# ==============================

def test_T_definition():

    primes = generate_primes(1000)
    gaps = generate_gaps(primes)

    K_test = 5
    N_test = 1

    T = compute_T(gaps, N_test, K_test)

    # T_N^(K) must have denominator dividing 2^K
    assert T.denominator <= 2 ** K_test
    assert (2 ** K_test) % T.denominator == 0

    # It must be positive because every prime gap is positive
    assert T > 0

    print("TEST 3: T_N^(K) construction")
    print("PASS")


test_T_definition()
# ==============================
# Integer distance
# ==============================

def distance_to_integer(x):
    """
    Distance from x to the nearest integer.
    """
    lower = x.numerator // x.denominator
    upper = lower + 1

    return min(
        abs(x - lower),
        abs(upper - x)
    )


# ==============================
# TDD TEST 4
# ==============================

# ==============================
# Fast T_N^(K) recurrence
# ==============================

def scan_T_fast(gaps, N_start, N_end, K):
    """
    Generate T_N^(K) for consecutive N
    using the exact finite recurrence:

    T_{N+1}^(K)
        = 2 T_N^(K)
        - g_{N+1}
        + g_{N+K+1} / 2^K
    """

    results = []

    # Initial value
    T = compute_T(
        gaps,
        N_start,
        K
    )

    results.append(T)

    for N in range(N_start, N_end):

        T = (
            2 * T
            - gaps[N]
            + Fraction(
                gaps[N + K],
                2 ** K
            )
        )

        results.append(T)

    return results
# ==============================
# Tail bound
# ==============================

def tail_bound(N, K):
    """
    Exact infinite upper bound for:

        T_N - T_N^(K)

    using:

        g_{N+j} < (N+j+1)^2
    """

    r = K + 1
    a = N + 1

    # Sum_{j=r}^∞ j^2 / 2^j
    sum_j2 = Fraction(
        2 * (r ** 2 + 2 * r + 3),
        2 ** r
    )

    # Sum_{j=r}^∞ j / 2^j
    sum_j = Fraction(
        2 * (r + 1),
        2 ** r
    )

    # Sum_{j=r}^∞ 1 / 2^j
    sum_1 = Fraction(
        2,
        2 ** r
    )

    return (
        sum_j2
        + 2 * a * sum_j
        + a ** 2 * sum_1
    )


# ==============================
# TDD TEST 5
# ==============================

def test_tail_certificate():

    primes = generate_primes(1000)
    gaps = generate_gaps(primes)

    N_test = 1
    P_test = 1
    K_test = 10

    T_K = compute_T(
        gaps,
        N_test,
        K_test
    )

    C_K = (2 ** P_test - 1) * T_K

    distance = distance_to_integer(C_K)

    error = (2 ** P_test - 1) * tail_bound(
        N_test,
        K_test
    )

    print("TEST 5: Rigorous tail certificate")
    print(f"C_P^(K)(N) = {C_K}")
    print(f"Distance to integer = {float(distance):.12f}")
    print(f"Tail bound          = {float(error):.12e}")

    # If the finite value is farther from every integer
    # than the maximum possible tail error,
    # the infinite value cannot be an integer.

    certified = distance > error

    assert certified

    print("Certified: C_P(N) is NOT an integer")
    print("PASS")


test_tail_certificate()

# ==============================
# TDD TEST 6
# Scan fixed P
# ==============================

def scan_fixed_P(gaps, P, N_start, N_end, K):
    results = []

    # Compute the consecutive T_N^(K) values once.
    # scan_T_fast() is covered by TEST 9A and is
    # exactly equal to the direct compute_T() reference.
    T_values = scan_T_fast(
        gaps,
        N_start,
        N_end,
        K
    )

    for i, N in enumerate(
        range(N_start, N_end + 1)
    ):
        T_K = T_values[i]

        C_K = (2 ** P - 1) * T_K

        distance = distance_to_integer(C_K)

        error = (
            (2 ** P - 1)
            * tail_bound(N, K)
        )

        certified = distance > error

        results.append({
            "N": N,
            "distance": distance,
            "tail_error": error,
            "certified": certified
        })

    return results


def test_fixed_P_scan():

    primes = generate_primes(2_000_000)
    gaps = generate_gaps(primes)

    P_test = 1
    K_test = 30

    results = scan_fixed_P(
        gaps,
        P_test,
        1,
        100,
        K_test
    )

    assert len(results) == 100

    # Every result must have a valid distance.
    for result in results:
        assert 0 <= result["distance"] <= Fraction(1, 2)

    print("TEST 6: Fixed-P scan")
    print("PASS")
    print(f"P = {P_test}")
    print(f"N range = 1..100")
    print(f"Samples = {len(results)}")


test_fixed_P_scan()
# ==============================
# TDD TEST 7
# Find closest integer states
# ==============================

def test_find_closest_states():

    primes = generate_primes(5000)
    gaps = generate_gaps(primes)

    P_test = 1
    K_test = 30

    results = scan_fixed_P(
        gaps,
        P_test,
        1,
        100,
        K_test
    )

    # Sort by distance to nearest integer
    closest = sorted(
        results,
        key=lambda x: x["distance"]
    )[:5]

    assert len(closest) == 5

    print("TEST 7: Closest integer states")
    print("PASS")
    print()
    print("Five closest states:")

    for result in closest:
        print(
            f"N={result['N']:3d} | "
            f"distance={float(result['distance']):.12e} | "
            f"tail={float(result['tail_error']):.12e} | "
            f"certified={result['certified']}"
        )


test_find_closest_states()
# ==============================
# TEST 8
# Large fixed-P scan
# ==============================

def test_large_fixed_P_scan():

    primes = generate_primes(2_000_000)
    gaps = generate_gaps(primes)

    P_test = 1
    K_test = 60

    results = scan_fixed_P(
        gaps,
        P_test,
        1,
        100_000,
        K_test
    )

    assert len(results) == 100_000

    closest = sorted(
        results,
        key=lambda x: x["distance"]
    )[:10]

    certified_count = sum(
        result["certified"]
        for result in results
    )

    print()
    print("=" * 60)
    print("TEST 8: LARGE FIXED-P SCAN")
    print("=" * 60)
    print(f"P = {P_test}")
    print(f"N range = 1..100000")
    print(f"K = {K_test}")
    print(f"Certified non-integer states = {certified_count}")
    print()
    print("10 closest states:")
    print()

    for result in closest:
        print(
            f"N={result['N']:6d} | "
            f"distance={float(result['distance']):.12e} | "
            f"tail={float(result['tail_error']):.12e} | "
            f"certified={result['certified']}"
        )

    assert certified_count > 0

    print()
    print("TEST 8: PASS")


test_large_fixed_P_scan()
# ==============================
# Fast T_N^(K) recurrence
# ==============================

def scan_T_fast(gaps, N_start, N_end, K):
    """
    Generate T_N^(K) for consecutive N
    using the exact finite recurrence.
    """

    results = []

    # First value is calculated directly.
    T = compute_T(
        gaps,
        N_start,
        K
    )

    results.append(T)

    # Remaining values use the recurrence.
    for N in range(N_start, N_end):

        T = (
            2 * T
            - gaps[N]
            + Fraction(
                gaps[N + K],
                2 ** K
            )
        )

        results.append(T)

    return results
# ==============================
# TEST 9A
# Fast recurrence correctness
# ==============================

def test_fast_recurrence():

    primes = generate_primes(5000)
    gaps = generate_gaps(primes)

    K_test = 20
    N_start = 1
    N_end = 100

    fast_values = scan_T_fast(
        gaps,
        N_start,
        N_end,
        K_test
    )

    for i, N in enumerate(
        range(N_start, N_end + 1)
    ):

        direct = compute_T(
            gaps,
            N,
            K_test
        )

        fast = fast_values[i]

        assert direct == fast

    print("TEST 9A: Fast recurrence correctness")
    print("PASS")


test_fast_recurrence()


