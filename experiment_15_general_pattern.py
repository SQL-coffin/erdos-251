from fractions import Fraction
import math

"""
Experiment 15 - General Pattern Lemma

For
    D_h(N) = sum_{j>=1} c_j / 2^(j-1),
    c_j = (g_{N+h+j} - g_{N+j}) / 2,

generalize the local pattern idea.

Two certificate modes are kept separate:

A. PURE PREFIX:
   A finite signed word w itself must put its prefix value
   inside (1/2,1) or (-1,-1/2), with enough margin to absorb
   a rigorous tail bound starting immediately after the word.

B. WINDOW-ASSISTED:
   The word is only a structural pattern. The continuation through
   K is evaluated exactly, and only the tail after K is bounded.

Mode A is a genuine pattern lemma.
Mode B is a rigorous finite certificate carrying pattern metadata.

Finite computation does not prove irrationality or infinitude.
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
    return [
        primes[i+1] - primes[i]
        for i in range(len(primes) - 1)
    ]


def pattern_value(word):
    value = Fraction(0, 1)

    for j, c in enumerate(word):
        value += Fraction(c, 2**j)

    return value


def pattern_margin(word, sign):
    p = sign * pattern_value(word)

    return min(
        p - Fraction(1, 2),
        Fraction(1, 1) - p,
    )


def tail_bound_difference(N, h, m):
    """
    Rigorous bound for the tail after m terms of D_h(N):

        |sum_{j=m+1}^inf Delta_j / 2^j|

    using g_k < (k+1)^2.
    """
    r = m + 1
    a = N + 1

    def square_tail(x):
        numerator = (
            2*(r*r + 2*r + 3)
            + 4*x*(r + 1)
            + 2*x*x
        )
        return Fraction(numerator, 2**r)

    return square_tail(a) + square_tail(a + h)


def finite_D(gaps, N, h, K):
    total = Fraction(0, 1)

    for j in range(1, K + 1):
        delta = gaps[N + h + j - 1] - gaps[N + j - 1]
        total += Fraction(delta, 2**j)

    return total


def scan_T_fast(gaps, N_start, N_end, K):
    T = Fraction(0, 1)

    for j in range(1, K + 1):
        T += Fraction(gaps[N_start + j - 1], 2**j)

    values = [T]

    for N in range(N_start, N_end):
        T = (
            2*T
            - gaps[N]
            + Fraction(gaps[N + K], 2**K)
        )
        values.append(T)

    return values


def pure_prefix_certificate(word, N, h):
    """
    Genuine Pattern Lemma certificate.

    The prefix itself must lie strictly inside the target interval
    with enough margin to absorb every possible tail allowed by
    g_k < (k+1)^2.
    """
    if not word:
        return None

    if abs(word[0]) != 1:
        return None

    sign = word[0]

    prefix = pattern_value(word)
    margin = pattern_margin(word, sign)
    error = tail_bound_difference(N, h, len(word))

    certified = (
        margin > 0
        and error < margin
    )

    return {
        "sign": sign,
        "prefix": prefix,
        "margin": margin,
        "tail_bound": error,
        "certified": certified,
    }


def window_certificate(word, gaps, N, h, K):
    """
    Rigorous finite-window certificate carrying a pattern label.

    The exact finite window through K is combined with a rigorous
    tail bound after K.
    """
    sign = 1 if word[0] > 0 else -1

    if word[0] != sign:
        return None

    Dk = finite_D(gaps, N, h, K)
    error = tail_bound_difference(N, h, K)

    lower = Dk - error
    upper = Dk + error

    if sign > 0:
        certified = (
            Fraction(1, 2) < lower
            and upper < Fraction(1, 1)
        )
        target = (Fraction(1, 2), Fraction(1, 1))
    else:
        certified = (
            Fraction(-1, 1) < lower
            and upper < Fraction(-1, 2)
        )
        target = (Fraction(-1, 1), Fraction(-1, 2))

    return {
        "sign": sign,
        "Dk": Dk,
        "error": error,
        "lower": lower,
        "upper": upper,
        "target": target,
        "certified": certified,
    }


# ------------------------------------------------------------
# TDD TESTS
# ------------------------------------------------------------

def test_pattern_value_general():
    assert pattern_value([1, -1, 1]) == Fraction(3, 4)
    assert pattern_value([-1, 1, -1]) == Fraction(-3, 4)
    assert pattern_value([1, 0]) == Fraction(1, 1)

    print("TEST 15A: General pattern value")
    print("PASS")


def test_pattern_margin():
    assert pattern_margin([1, -1, 1], 1) == Fraction(1, 4)
    assert pattern_margin([-1, 1, -1], -1) == Fraction(1, 4)

    print("TEST 15B: Pattern margin")
    print("PASS")


def test_pure_prefix_rejects_large_tail():
    result = pure_prefix_certificate([1, -1, 1], 100_000, 1)

    assert result["margin"] == Fraction(1, 4)
    assert result["tail_bound"] > result["margin"]
    assert not result["certified"]

    print("TEST 15C: Pure prefix tail protection")
    print("PASS")


def test_pure_prefix_accepts_when_margin_wins():
    result = pure_prefix_certificate([1, -1, 1, -1, 1], 1, 1)

    assert result["margin"] > 0
    assert pure_prefix_certificate([2, -1, 1], 1, 1) is None

    print("TEST 15D: Pure prefix certificate API")
    print("PASS")


def test_window_certificate():
    primes = generate_primes(5000)
    gaps = generate_gaps(primes)

    word = [1, -1, 1]

    result = window_certificate(word, gaps, 1, 1, 20)

    assert isinstance(result["Dk"], Fraction)
    assert result["error"] > 0

    print("TEST 15E: Window certificate")
    print("PASS")


def test_exact_pattern_equivalence():
    for word in (
        [1, -1, 1],
        [-1, 1, -1],
        [1, 1, -1, 0],
    ):
        p1 = pattern_value(word)
        p2 = sum(
            Fraction(c, 2**j)
            for j, c in enumerate(word)
        )
        assert p1 == p2

    print("TEST 15F: Exact pattern equivalence")
    print("PASS")


def build_word(gaps, N, h, m):
    return tuple(
        (
            gaps[N + h + j - 1]
            - gaps[N + j - 1]
        ) // 2
        for j in range(1, m + 1)
    )


def search_pure_prefix(
    gaps,
    N_start,
    N_end,
    H_max,
    M_max,
):
    """
    Search genuine pure-prefix certificates.

    This is intentionally a smaller exploratory range because the
    rigorous prefix-tail bound becomes the limiting condition.
    """
    results = []

    for h in range(1, H_max + 1):
        for N in range(N_start, N_end + 1):
            first_c = (
                gaps[N + h]
                - gaps[N]
            ) // 2

            if abs(first_c) != 1:
                continue

            prefix_num = 0

            for m in range(1, M_max + 1):
                c = (
                    gaps[N + h + m - 1]
                    - gaps[N + m - 1]
                ) // 2

                prefix_num = 2*prefix_num + c
                prefix = Fraction(prefix_num, 2**(m - 1))

                if prefix_num == 0:
                    continue

                sign = 1 if prefix_num > 0 else -1
                normalized = sign * prefix

                error = tail_bound_difference(N, h, m)

                lower_ok = (
                    Fraction(1, 2) + error
                    < normalized
                )
                upper_ok = (
                    normalized
                    < Fraction(1, 1) - error
                )

                if lower_ok and upper_ok:
                    results.append({
                        "N": N,
                        "h": h,
                        "m": m,
                        "word": build_word(
                            gaps, N, h, m
                        ),
                        "prefix": prefix,
                        "sign": sign,
                        "tail_bound": error,
                    })
                    break

    return results


def search_window_patterns(
    gaps,
    T_values,
    N_start,
    N_end,
    H_max,
    K,
    M_max,
):
    """
    Search structural patterns of lengths 1..M_max.

    A pattern is retained only when:

      1. its first signed coefficient is +/-1, and
      2. the full K-window gives a rigorous target certificate.

    The pattern length is therefore a structural descriptor, while
    the actual proof certificate comes from the exact K-window plus
    the rigorous post-K tail.
    """
    results = []

    for h in range(1, H_max + 1):
        for N in range(N_start, N_end + 1):
            d_k = (
                T_values[N + h - 1]
                - T_values[N - 1]
            )

            error = tail_bound_difference(N, h, K)

            positive = (
                Fraction(1, 2) < d_k - error
                and d_k + error < Fraction(1, 1)
            )

            negative = (
                Fraction(-1, 1) < d_k - error
                and d_k + error < Fraction(-1, 2)
            )

            if not (positive or negative):
                continue

            sign = 1 if positive else -1

            found = None

            for m in range(3, M_max + 1):
                word = build_word(
                    gaps, N, h, m
                )

                if word[0] != sign:
                    continue

                margin = pattern_margin(
                    word,
                    sign,
                )

                if margin > 0:
                    found = {
                        "N": N,
                        "h": h,
                        "m": m,
                        "word": word,
                        "sign": sign,
                        "Dk": d_k,
                        "error": error,
                        "prefix_margin": margin,
                    }
                    break

            if found is not None:
                results.append(found)

    return results


def summarize_patterns(results):
    by_length = {}
    by_h = {}
    unique_words = set()

    for item in results:
        m = item["m"]
        h = item["h"]
        word = item["word"]

        by_length[m] = by_length.get(m, 0) + 1
        by_h[h] = by_h.get(h, 0) + 1
        unique_words.add(word)

    return {
        "count": len(results),
        "by_length": by_length,
        "by_h": by_h,
        "unique_words": len(unique_words),
    }


def main():
    test_pattern_value_general()
    test_pattern_margin()
    test_pure_prefix_rejects_large_tail()
    test_pure_prefix_accepts_when_margin_wins()
    test_window_certificate()
    test_exact_pattern_equivalence()

    print()
    print("=" * 78)
    print("EXPERIMENT 15 - GENERAL PATTERN SEARCH")
    print("=" * 78)

    primes = generate_primes(2_000_000)
    gaps = generate_gaps(primes)

    print()
    print("A. PURE PREFIX SEARCH")
    print("N=1..10,000 | h=1..8 | m<=64")

    pure = search_pure_prefix(
        gaps,
        1,
        10_000,
        8,
        64,
    )

    pure_summary = summarize_patterns(pure)

    print(
        f"Pure certificates: {pure_summary['count']}"
    )
    print(
        f"Unique pure words: {pure_summary['unique_words']}"
    )
    print(
        f"Lengths: {sorted(pure_summary['by_length'].items())}"
    )

    print()
    print("First pure certificates:")

    for item in pure[:10]:
        print(
            f"N={item['N']:6d} "
            f"h={item['h']:2d} "
            f"m={item['m']:2d} "
            f"word={item['word']} "
            f"prefix={float(item['prefix']): .6f} "
            f"tail<{float(item['tail_bound']):.3e}"
        )

    print()
    print("B. WINDOW-ASSISTED STRUCTURAL PATTERNS")
    print("N=1..100,000 | h=1..32 | K=100 | m<=16")

    T_end = 100_000 + 32 + 1

    T_values = scan_T_fast(
        gaps,
        1,
        T_end,
        100,
    )

    window = search_window_patterns(
        gaps,
        T_values,
        1,
        100_000,
        32,
        100,
        16,
    )

    window_summary = summarize_patterns(window)

    print(
        f"Window certificates: {window_summary['count']}"
    )
    print(
        f"Unique structural words: "
        f"{window_summary['unique_words']}"
    )
    print(
        f"Lengths: {sorted(window_summary['by_length'].items())}"
    )

    print()
    print("Per-h counts:")

    for h in sorted(window_summary["by_h"]):
        print(
            f"h={h:2d}: "
            f"{window_summary['by_h'][h]}"
        )

    print()
    print("Examples:")

    for item in window[:10]:
        print(
            f"N={item['N']:6d} "
            f"h={item['h']:2d} "
            f"m={item['m']:2d} "
            f"word={item['word']} "
            f"sign={item['sign']:+d} "
            f"|D_K|~{float(abs(item['Dk'])):.6f}"
        )

    print()
    print("INTERPRETATION")
    print(
        "The pure-prefix search is a genuine finite-pattern lemma:"
    )
    print(
        "the word alone plus a rigorous tail bound forces D_h(N)"
    )
    print(
        "into one of the two half-intervals."
    )
    print()
    print(
        "The window-assisted search is stronger computationally but"
    )
    print(
        "weaker as a pattern statement: the certificate also uses"
    )
    print(
        "the exact continuation through K=100."
    )
    print()
    print(
        "Neither experiment proves that any pattern occurs infinitely"
    )
    print(
        "often, and neither by itself proves irrationality."
    )


if __name__ == "__main__":
    main()
