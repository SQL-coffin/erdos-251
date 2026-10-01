from fractions import Fraction
import math

"""
Experiment 13 - Small tail differences + gap mismatch

D_h(N) = T_{N+h} - T_N.

Exact identity:
D_h(N+1) = 2 D_h(N) - g_{N+h+1} + g_{N+1}.

We search for finite witnesses where both neighboring
tail differences have absolute value < 1, while the
corresponding gaps are different.
"""

def generate_primes(limit):
    sieve = bytearray([1]) * (limit + 1)
    sieve[0:2] = bytearray([0, 0])

    for p in range(2, int(math.isqrt(limit)) + 1):
        if sieve[p]:
            sieve[p*p:limit+1:p] = bytearray(
                ((limit - p*p) // p) + 1
            )

    return [n for n in range(2, limit + 1) if sieve[n]]


def generate_gaps(primes):
    return [primes[i+1] - primes[i]
            for i in range(len(primes) - 1)]


def compute_T(gaps, N, K):
    total = Fraction(0, 1)
    for j in range(1, K + 1):
        total += Fraction(gaps[N + j - 1], 2**j)
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


def tail_bound(N, K):
    r = K + 1
    a = N + 1

    sum_j2 = Fraction(
        2*(r**2 + 2*r + 3),
        2**r
    )
    sum_j = Fraction(
        2*(r + 1),
        2**r
    )
    sum_1 = Fraction(2, 2**r)

    return (
        sum_j2
        + 2*a*sum_j
        + a*a*sum_1
    )


def test_recurrence_identity():
    primes = generate_primes(5000)
    gaps = generate_gaps(primes)

    K = 20
    T = scan_T_fast(gaps, 1, 100, K)

    for h in range(1, 6):
        for N in range(1, 50):
            d = T[N + h - 1] - T[N - 1]
            d_next = T[N + h] - T[N]

            boundary = Fraction(
                gaps[N + h + K] - gaps[N + K],
                2**K
            )

            rhs = (
                2*d
                - gaps[N + h]
                + gaps[N]
                + boundary
            )

            assert d_next == rhs

    print("TEST 13A: D_h recurrence")
    print("PASS")


def test_small_tail_certificate():
    primes = generate_primes(5000)
    gaps = generate_gaps(primes)

    K = 20
    T = scan_T_fast(gaps, 1, 20, K)

    h = 1
    N = 1

    d_k = T[N + h - 1] - T[N - 1]
    error = tail_bound(N, K) + tail_bound(N + h, K)

    assert abs(d_k) + error < 1

    print("TEST 13B: Small-tail certificate")
    print("PASS")


def certify_small(value, error):
    return abs(value) + error < 1


def scan_witnesses(
    gaps,
    T_values,
    tail_errors,
    h,
    N_start,
    N_end,
):
    witnesses = []

    for N in range(N_start, N_end + 1):
        d_k = (
            T_values[N + h - 1]
            - T_values[N - 1]
        )
        d_next_k = (
            T_values[N + h]
            - T_values[N]
        )

        error = tail_errors[N] + tail_errors[N + h]
        error_next = (
            tail_errors[N + 1]
            + tail_errors[N + h + 1]
        )

        small_now = certify_small(d_k, error)
        small_next = certify_small(d_next_k, error_next)
        gap_mismatch = gaps[N + h] != gaps[N]

        if small_now and small_next and gap_mismatch:
            witnesses.append({
                "N": N,
                "h": h,
                "d": d_k,
                "d_next": d_next_k,
                "error": error,
                "error_next": error_next,
                "gap_now": gaps[N],
                "gap_shift": gaps[N + h],
            })

    return witnesses


def main():
    test_recurrence_identity()
    test_small_tail_certificate()

    K = 100
    N_start = 1
    N_end = 100_000
    H_max = 32
    T_end = N_end + H_max + 1

    print("=" * 78)
    print("EXPERIMENT 13 - SMALL TAIL DIFFERENCE + GAP MISMATCH")
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

    tail_errors = {
        N: tail_bound(N, K)
        for N in range(1, T_end + 1)
    }

    all_witnesses = {}

    print()
    print(
        f"{'h':>3} {'witnesses':>12} "
        f"{'first N':>12} {'last N':>12}"
    )

    for h in range(1, H_max + 1):
        witnesses = scan_witnesses(
            gaps,
            T_values,
            tail_errors,
            h,
            N_start,
            N_end
        )

        all_witnesses[h] = witnesses

        if witnesses:
            print(
                f"{h:>3} {len(witnesses):>12} "
                f"{witnesses[0]['N']:>12} "
                f"{witnesses[-1]['N']:>12}"
            )
        else:
            print(
                f"{h:>3} {0:>12} "
                f"{'-':>12} {'-':>12}"
            )

    print()
    print("Example late witnesses:")

    shown = 0

    for h, witnesses in all_witnesses.items():
        late = [
            w for w in witnesses
            if w["N"] >= 90_000
        ]

        if not late:
            continue

        w = late[0]

        print(
            f"h={h:2d} N={w['N']:6d} "
            f"| gaps=({w['gap_now']},{w['gap_shift']}) "
            f"| |D|~{float(abs(w['d'])):.3e} "
            f"| |D_next|~{float(abs(w['d_next'])):.3e}"
        )

        shown += 1
        if shown >= 5:
            break

    print()
    print("Interpretation:")
    print(
        "A witness has both neighboring shifts rigorously "
        "inside (-1,1), while the two gaps differ."
    )
    print(
        "For eventual integer-shift period h, such a sufficiently "
        "late witness contradicts eventual integrality."
    )
    print("Finite witnesses alone do not prove irrationality.")


if __name__ == "__main__":
    main()
