# -*- coding: utf-8 -*-
import math
import time
import sympy as sp

# ============================================================
# Experiment 02
# Erdős Problem #251
#
# Study:
#     T_N = sum_{j>=1} g_{N+j} / 2^j
#
# where
#     g_n = p_{n+1} - p_n
#
# We want to test:
#
#     {T_N} ?= {T_{N+1}} / 2
#
# ============================================================


# ------------------------------------------------------------
# 1. Generate primes using Sieve of Eratosthenes
# ------------------------------------------------------------

def generate_primes(n):
    """
    Generate the first n prime numbers.
    """
    # Rough upper bound for nth prime
    if n < 6:
        limit = 15
    else:
        limit = int(n * (math.log(n) + math.log(math.log(n)))) + 100

    while True:
        sieve = bytearray(b"\x01") * (limit + 1)
        sieve[0:2] = b"\x00\x00"

        for p in range(2, int(math.sqrt(limit)) + 1):
            if sieve[p]:
                sieve[p * p : limit + 1 : p] = b"\x00" * (
                    ((limit - p * p) // p) + 1
                )

        primes = [i for i in range(2, limit + 1) if sieve[i]]

        if len(primes) >= n:
            return primes[:n]

        limit *= 2


# ------------------------------------------------------------
# 2. Parameters
# ------------------------------------------------------------

N_MAX = 100_000
K = 40

print("=" * 60)
print("Erdos #251 - Experiment 02")
print("=" * 60)

print(f"Number of N values : {N_MAX}")
print(f"Tail depth K       : {K}")
print()


# ------------------------------------------------------------
# 3. Generate primes
#
# We need N_MAX + K + 2 primes because T_N uses
# gaps beyond N.
# ------------------------------------------------------------

start = time.time()

primes = generate_primes(N_MAX + K + 2)

print(f"Generated {len(primes):,} primes")
print(f"Time: {time.time() - start:.3f} seconds")
print()


# ------------------------------------------------------------
# 4. Prime gaps
# ------------------------------------------------------------

gaps = [
    primes[i + 1] - primes[i]
    for i in range(len(primes) - 1)
]


# ------------------------------------------------------------
# 5. Compute T_N^(K)
#
# T_N^(K) = sum_{j=1}^K g_{N+j}/2^j
#
# Python index:
#
# g_n in mathematical notation corresponds to
# gaps[n-1] in Python.
# ------------------------------------------------------------

def compute_T(N):
    total = 0.0

    for j in range(1, K + 1):
        total += gaps[N + j - 1] / (2 ** j)

    return total


# ------------------------------------------------------------
# 6. Fractional part
# ------------------------------------------------------------

def fractional_part(x):
    return x - math.floor(x)


def distance_to_integer(x):
    f = fractional_part(x)
    return min(f, 1.0 - f)


# ------------------------------------------------------------
# 7. Experiment
# ------------------------------------------------------------

print("Testing the fractional-part relation:")
print()

max_error = 0.0

interesting = []

for N in range(1, N_MAX):

    T_N = compute_T(N)
    T_next = compute_T(N + 1)

    frac_N = fractional_part(T_N)
    frac_next = fractional_part(T_next)

    predicted = frac_next / 2

    error = abs(frac_N - predicted)

    max_error = max(max_error, error)

    if error > 1e-10:
        print("POSSIBLE FAILURE")
        print(f"N        = {N}")
        print(f"T_N      = {T_N}")
        print(f"T_N+1    = {T_next}")
        print(f"frac(T_N)      = {frac_N}")
        print(f"frac(T_N+1)/2  = {predicted}")
        print(f"error           = {error}")
        print()

    d = distance_to_integer(T_N)

    if d < 1e-4:
        interesting.append(
            (N, T_N, frac_N, d)
        )


# ------------------------------------------------------------
# 8. Results
# ------------------------------------------------------------

print("=" * 60)
print("RESULT")
print("=" * 60)

print(f"Maximum error:")
print(max_error)

print()

if max_error < 1e-10:
    print("The tested relation appears numerically valid:")
    print()
    print("    {T_N} = {T_(N+1)} / 2")
else:
    print("The relation is NOT numerically exact.")
    
print()


# ------------------------------------------------------------
# 9. Interesting near-integer values
# ------------------------------------------------------------

print("=" * 60)
print("T_N values very close to integers")
print("=" * 60)

for N, T, frac, d in interesting[:30]:

    nearest = round(T)

    print(
        f"N={N:6d}   "
        f"T={T:.12f}   "
        f"nearest={nearest:6d}   "
        f"distance={d:.3e}"
    )


# ------------------------------------------------------------
# 10. Examine the first few values manually
# ------------------------------------------------------------

print()
print("=" * 60)
print("First 20 fractional parts")
print("=" * 60)

for N in range(1, 21):

    T = compute_T(N)
    f = fractional_part(T)

    print(
        f"N={N:3d}   "
        f"T_N={T:.10f}   "
        f"frac={f:.10f}"
    )


print()
print("=" * 60)
print("Experiment complete")
print("=" * 60)
total = 0.0

for j in range(1, 11):
    gap = gaps[1 + j - 1]
    term = gap / (2 ** j)

    total += term

    print(
        f"j={j:2d}   "
        f"gap={gap:2d}   "
        f"term={term:.10f}   "
        f"sum={total:.10f}"
    )

print()
print("T_1 using first 10 terms =", total)
print("T_1 using K=40 =", compute_T(1))
print()
print("=" * 60)
print("CHECK RECURSION")
print("=" * 60)

T1 = compute_T(1)
T2 = compute_T(2)

left = T1
right = (gaps[1] + T2) / 2

print("T1 =", T1)
print("T2 =", T2)
print("g2 =", gaps[1])
print()
print("(g2 + T2) / 2 =", right)
print()
print("Difference =", abs(left - right))
print()
print("=" * 60)
print("EXACT INTEGER CHECK")
print("=" * 60)

K_test = 40

numerator_1 = 0

for j in range(1, K_test + 1):
    numerator_1 += gaps[j] * (2 ** (K_test - j))

T1_exact = numerator_1 / (2 ** K_test)

print("Numerator =", numerator_1)
print("Denominator =", 2 ** K_test)
print("T1 exact rational =", T1_exact)

print()
print("Float T1 =", compute_T(1))
print("Difference =", abs(T1_exact - compute_T(1)))
from fractions import Fraction

print()
print("=" * 60)
print("EXACT RECURSION CHECK")
print("=" * 60)


def compute_T_exact(N, K=40):

    numerator = 0

    for j in range(1, K + 1):
        numerator += gaps[N + j - 1] * (2 ** (K - j))

    denominator = 2 ** K

    return Fraction(numerator, denominator)


T1_exact = compute_T_exact(1)
T2_exact = compute_T_exact(2)

right_exact = (gaps[1] + T2_exact) / 2

print("T1 exact =", T1_exact)
print("T2 exact =", T2_exact)
print()

print("(g2 + T2) / 2 =", right_exact)
print()

print("Are they exactly equal?")
print(T1_exact == right_exact)
print("g42 =", gaps[41])
print()
print("=" * 60)
print("RAW PRIME / GAP CHECK")
print("=" * 60)

for n in range(1, 44):
    print(
        f"p_{n} = {primes[n-1]:4d}"
        f"    g_{n} = {gaps[n-1]:2d}"
    )

from fractions import Fraction

print()
print("=" * 60)
print("CROSS-CHECK T1")
print("=" * 60)

# 方法 1：使用我们原来的 compute_T_exact()
T1_method_1 = compute_T_exact(1)

# 方法 2：完全逐项相加
T1_method_2 = Fraction(0, 1)

for j in range(1, 41):
    term = Fraction(gaps[j], 2 ** j)
    T1_method_2 += term

print("Method 1:")
print(T1_method_1)

print()

print("Method 2:")
print(T1_method_2)

print()

print("Are they exactly equal?")
print(T1_method_1 == T1_method_2)

print()

print("Difference:")
print(T1_method_1 - T1_method_2)
print()
print("=" * 60)
print("MANUAL CHECK: T_2")
print("=" * 60)

total = 0.0

for j in range(1, 11):
    gap = gaps[2 + j - 1]
    term = gap / (2 ** j)

    total += term

    print(
        f"j={j:2d}   "
        f"gap={gap:2d}   "
        f"term={term:.10f}   "
        f"sum={total:.10f}"
    )

print()
print("T_2 using first 10 terms =", total)
print("T_2 using K=40 =", compute_T(2))
from fractions import Fraction

print()
print("=" * 60)
print("EXACT INTEGER CHECK: T_2")
print("=" * 60)

K_test = 40

numerator_2 = 0

for j in range(1, K_test + 1):
    numerator_2 += gaps[2 + j - 1] * (2 ** (K_test - j))

denominator_2 = 2 ** K_test

T2_exact_integer = Fraction(
    numerator_2,
    denominator_2
)

print("Numerator =", numerator_2)
print("Denominator =", denominator_2)
print("T2 exact rational =", T2_exact_integer)

print()
print("Float T2 =", compute_T(2))
print(
    "Difference =",
    float(T2_exact_integer) - compute_T(2)
)
print()
print("=" * 60)
print("CHECK T_2 BOUNDARY GAP")
print("=" * 60)

print("p_42 =", primes[41])
print("p_43 =", primes[42])
print("p_44 =", primes[43])

print()
print("g_42 =", gaps[41])
print("g_43 =", gaps[42])
from fractions import Fraction

T2_exact = compute_T_exact(2)
T3_exact = compute_T_exact(3)

lhs = (gaps[2] + T3_exact) / 2
rhs = T2_exact + Fraction(gaps[42], 2**41)

print()
print("=" * 60)
print("EXACT T_2 RECURSION CHECK")
print("=" * 60)

print("T2 =", T2_exact)
print("T3 =", T3_exact)

print()
print("(g3 + T3) / 2 =", lhs)
print("T2 + g43/2^41 =", rhs)

print()
print("Exactly equal:", lhs == rhs)
print("Difference:", lhs - T2_exact)
print("Expected boundary:", Fraction(gaps[42], 2**41))
print()
print("=" * 60)
print("LOCAL RECURSION STRUCTURE")
print("=" * 60)

for N in range(1, 11):
    current = compute_T(N)
    next_value = compute_T(N + 1)
    gap = gaps[N]

    predicted = (gap + next_value) / 2
    error = predicted - current

    print(
        f"N={N:2d} | "
        f"g_(N+1)={gap:2d} | "
        f"T_N={current:.10f} | "
        f"T_(N+1)={next_value:.10f} | "
        f"error={error:.3e}"
    )
    print()
print("=" * 60)
print("REMOVE THE FIRST GAP CONTRIBUTION")
print("=" * 60)

for N in range(1, 11):
    T = compute_T(N)
    next_T = compute_T(N + 1)
    gap = gaps[N]

    R = T - gap / 2
    expected = next_T / 2

    print(
        f"N={N:2d} | "
        f"R_N={R:.10f} | "
        f"T_(N+1)/2={expected:.10f} | "
        f"difference={R-expected:.3e}"
    )
from fractions import Fraction

print()
print("=" * 60)
print("EXPANDED RECURSION CHECK")
print("=" * 60)

T1_exact = compute_T_exact(1)

for N in range(2, 11):

    TN_exact = compute_T_exact(N)

    # 2^(N-1) * T1
    main_term = (2 ** (N - 1)) * T1_exact

    # Sum_{k=2}^{N} 2^(N-k) * g_k
    gap_sum = Fraction(0, 1)

    for k in range(2, N + 1):
        gap_sum += (2 ** (N - k)) * gaps[k - 1]

    expanded = main_term - gap_sum

    difference = TN_exact - expanded

    print(
        f"N={N:2d} | "
        f"T_N={TN_exact} | "
        f"expanded={expanded} | "
        f"difference={difference}"
    )
print()
print("=" * 60)
print("BOUNDARY TERM PATTERN")
print("=" * 60)

for N in range(2, 11):

    TN = compute_T_exact(N)
    T1 = compute_T_exact(1)

    gap_sum = Fraction(0, 1)

    for k in range(2, N + 1):
        gap_sum += (2 ** (N - k)) * gaps[k - 1]

    expanded = (2 ** (N - 1)) * T1 - gap_sum

    difference = TN - expanded

    print(
        f"N={N:2d} | "
        f"difference={difference} | "
        f"g_(N+41)={gaps[N + 40]}"
    )
print()
print("=" * 60)
print("FRACTIONAL PART RECURSION CHECK")
print("=" * 60)

for N in range(1, 11):

    TN = compute_T(N)
    direct_frac = TN % 1

    predicted = (2 ** (N - 1)) * compute_T(1)
    predicted_frac = predicted % 1

    difference = direct_frac - predicted_frac

    print(
        f"N={N:2d} | "
        f"direct={direct_frac:.12f} | "
        f"predicted={predicted_frac:.12f} | "
        f"difference={difference:.3e}"
    )
print()
print("=" * 60)
print("BINARY REPRESENTATION OF T_1^(40)")
print("=" * 60)

T1_exact = compute_T_exact(1)

numerator = T1_exact.numerator
denominator = T1_exact.denominator

print("T1 exact =", T1_exact)
print("T1 decimal =", float(T1_exact))

print()
print("Numerator binary:")
print(bin(numerator))

print()
print("Denominator:")
print(denominator)
print("Denominator binary:")
print(bin(denominator))

print()
print("T1 binary expansion:")
print(format(float(T1_exact), ".40f"))
print()
print("=" * 60)
print("PRIME GAP CONTRIBUTIONS IN BINARY")
print("=" * 60)

for j in range(1, 11):
    gap = gaps[j]
    contribution = Fraction(gap, 2**j)

    print(
        f"j={j:2d} | "
        f"gap={gap:2d} | "
        f"contribution={contribution} | "
        f"decimal={float(contribution):.10f}"
    )
print()
print("=" * 60)
print("EXACT BINARY CHECK: T_1^(10)")
print("=" * 60)

T10 = Fraction(0, 1)

for j in range(1, 11):
    T10 += Fraction(gaps[j], 2**j)

print("T1^(10) =", T10)
print("Decimal  =", float(T10))

# Binary expansion
integer_part = T10.numerator // T10.denominator
fraction_part = T10 - integer_part

bits = ""

for _ in range(20):
    fraction_part *= 2

    if fraction_part >= 1:
        bits += "1"
        fraction_part -= 1
    else:
        bits += "0"

print()
print("Binary:")
print(f"{integer_part}.{bits}")
print()
print("=" * 60)
print("HALF-GAP TRANSFORMATION CHECK")
print("=" * 60)

T1_gap = Fraction(0, 1)
T1_half_gap = Fraction(0, 1)

for j in range(1, K + 1):
    g = gaps[j]
    h = g // 2

    # Original:
    T1_gap += Fraction(g, 2**j)

    # Since g = 2h:
    T1_half_gap += Fraction(h, 2**(j - 1))

print("Original T1^(40):")
print(T1_gap)

print()
print("Half-gap representation:")
print(T1_half_gap)

print()
print("Exactly equal:", T1_gap == T1_half_gap)
print("Difference:", T1_gap - T1_half_gap)
print()
print("=" * 60)
print("HALF-GAP SEQUENCE")
print("=" * 60)

T10_half = Fraction(0, 1)

for m in range(10):
    n = m + 2

    g = gaps[n - 1]
    h = g // 2

    contribution = Fraction(h, 2**m)
    T10_half += contribution

    print(
        f"m={m:2d} | "
        f"n={n:2d} | "
        f"g_n={g:2d} | "
        f"h_n={h:2d} | "
        f"contribution={contribution}"
    )

print()
print("T1^(10) from half-gaps:")
print(T10_half)

print()
print("Decimal:")
print(float(T10_half))

print()
print("Expected:")
print(Fraction(1201, 512))

print()
print("Exactly equal:", T10_half == Fraction(1201, 512))
print()
print("=" * 60)
print("AUTOMATIC BINARY CARRY CHECK")
print("=" * 60)

# --------------------------------------------------
# 1. Build raw dyadic coefficients
#
# T1^(40) = sum h_(m+2) / 2^m
#
# h_n = g_n / 2
# --------------------------------------------------

raw = []

for m in range(K):
    h = gaps[m + 1] // 2
    raw.append(h)

print("Half-gap sequence:")
print(raw)


# --------------------------------------------------
# 2. Carry propagation
#
# raw[m] = number of units at 2^(-m)
#
# Two units at position m -> one unit at position m-1
# --------------------------------------------------

digits = raw.copy()

carry = 0

# Start from the smallest binary position
for m in range(K - 1, -1, -1):

    total = digits[m] + carry

    digits[m] = total % 2
    carry = total // 2


# carry is now the integer part
integer_part = carry


# --------------------------------------------------
# 3. Construct binary representation
# --------------------------------------------------


integer_part = 2 * carry + digits[0]

binary_fraction = "".join(str(bit) for bit in digits[1:])

binary_string = f"{integer_part:b}.{binary_fraction}"

print()
print("Integer part:")
print(integer_part)

print()
print("Binary fraction:")
print(binary_fraction)

print()
print("Full binary representation:")
print(binary_string)


# --------------------------------------------------
# 4. Convert binary representation back to exact Fraction
# --------------------------------------------------

binary_value = Fraction(integer_part, 1)

for m, bit in enumerate(digits):
    if bit == 1:
        binary_value += Fraction(1, 2 ** m )


# --------------------------------------------------
# 5. Compare with original exact T1
# --------------------------------------------------

original = compute_T_exact(1)

print()
print("=" * 60)
print("EXACT VERIFICATION")
print("=" * 60)

print("Original:")
print(original)

print()
print("From binary:")
print(binary_value)

print()
print("Exactly equal:", original == binary_value)

print()
print("Difference:")
print(original - binary_value)
print()
print("=" * 80)
print("CARRY TABLE")
print("=" * 80)

carry = 0

print(
    f"{'m':>3} "
    f"{'h':>4} "
    f"{'carry in':>10} "
    f"{'total':>8} "
    f"{'bit':>5} "
    f"{'carry out':>10}"
)

print("-" * 80)

for m in range(K - 1, -1, -1):

    h = raw[m]
    carry_in = carry

    total = h + carry_in

    bit = total % 2
    carry = total // 2

    print(
        f"{m:>3} "
        f"{h:>4} "
        f"{carry_in:>10} "
        f"{total:>8} "
        f"{bit:>5} "
        f"{carry:>10}"
    )
print()
print("=" * 70)
print("CARRY GROWTH EXPERIMENT")
print("=" * 70)

for K_test in [10, 20, 40, 80, 160, 320]:

    raw_test = []

    for m in range(K_test):
        h = gaps[m + 1] // 2
        raw_test.append(h)

    carry = 0
    carry_values = []

    for m in range(K_test - 1, -1, -1):

        total = raw_test[m] + carry

        bit = total % 2
        carry = total // 2

        carry_values.append(carry)

    print(
        f"K = {K_test:>3} | "
        f"max carry = {max(carry_values):>3} | "
        f"final carry = {carry:>3}"
    )
print()
print("=" * 85)
print("CARRY vs HALF-GAP GROWTH")
print("=" * 85)

print(
    f"{'K':>5} "
    f"{'max h':>8} "
    f"{'max carry':>12} "
    f"{'avg carry':>12} "
    f"{'carry/max_h':>15}"
)

print("-" * 85)

for K_test in [10, 20, 40, 80, 160, 320]:

    raw_test = []

    for m in range(K_test):
        h = gaps[m + 1] // 2
        raw_test.append(h)

    carry = 0
    carry_values = []

    for m in range(K_test - 1, -1, -1):

        total = raw_test[m] + carry

        bit = total % 2
        carry = total // 2

        carry_values.append(carry)

    max_h = max(raw_test)
    max_carry = max(carry_values)
    avg_carry = sum(carry_values) / len(carry_values)

    ratio = max_carry / max_h

    print(
        f"{K_test:>5} "
        f"{max_h:>8} "
        f"{max_carry:>12} "
        f"{avg_carry:>12.3f} "
        f"{ratio:>15.3f}"
    )
print()
print("=" * 80)
print("BINARY PERIOD DETECTION")
print("=" * 80)

# --------------------------------------------------
# Generate binary digits again
# --------------------------------------------------

K_period = 320

raw_period = []

for m in range(K_period):
    h = gaps[m + 1] // 2
    raw_period.append(h)

digits_period = raw_period.copy()

carry = 0

for m in range(K_period - 1, -1, -1):

    total = digits_period[m] + carry

    digits_period[m] = total % 2
    carry = total // 2


# digits[0] belongs to the integer part
# We only study the fractional binary digits.

bits = digits_period[1:]

print("Number of binary fractional bits:", len(bits))


# --------------------------------------------------
# Test possible periods
# --------------------------------------------------

print()
print(
    f"{'Period':>8} "
    f"{'Matches':>10} "
    f"{'Compared':>10} "
    f"{'Match %':>10}"
)

print("-" * 45)

# Only test relatively short periods
for p in range(1, 41):

    matches = 0
    compared = 0

    # Ignore the first 50 bits as a possible transient.
    start = 50

    for i in range(start, len(bits) - p):

        compared += 1

        if bits[i] == bits[i + p]:
            matches += 1

    percentage = 100 * matches / compared

    print(
        f"{p:>8} "
        f"{matches:>10} "
        f"{compared:>10} "
        f"{percentage:>9.2f}%"
    )
print()
print("=" * 80)
print("BINARY PREFIX STABILITY TEST")
print("=" * 80)


def get_binary_bits(K_test):
    """
    Calculate the fractional binary digits of T1^(K).
    """

    raw = []

    for m in range(K_test):
        h = gaps[m + 1] // 2
        raw.append(h)

    digits = raw.copy()

    carry = 0

    # Carry from right to left
    for m in range(K_test - 1, -1, -1):

        total = digits[m] + carry

        digits[m] = total % 2
        carry = total // 2

    # digits[0] is the lowest integer bit.
    # Fractional bits start at digits[1].
    return digits[1:]


# --------------------------------------------------
# Generate binary expansions
# --------------------------------------------------

bits_80 = get_binary_bits(80)
bits_160 = get_binary_bits(160)
bits_320 = get_binary_bits(320)


# --------------------------------------------------
# Compare two sequences
# --------------------------------------------------

def common_prefix(a, b):

    n = min(len(a), len(b))

    for i in range(n):

        if a[i] != b[i]:
            return i

    return n


# --------------------------------------------------
# Compare
# --------------------------------------------------

prefix_80_160 = common_prefix(bits_80, bits_160)
prefix_160_320 = common_prefix(bits_160, bits_320)
prefix_80_320 = common_prefix(bits_80, bits_320)


print()
print("Common prefix lengths:")
print()

print(
    "K=80  vs K=160 :",
    prefix_80_160,
    "bits"
)

print(
    "K=160 vs K=320 :",
    prefix_160_320,
    "bits"
)

print(
    "K=80  vs K=320 :",
    prefix_80_320,
    "bits"
)


# --------------------------------------------------
# Show first 100 bits
# --------------------------------------------------

print()
print("First 100 bits:")
print()

print("K=80 :")
print("".join(map(str, bits_80[:100])))

print()

print("K=160:")
print("".join(map(str, bits_160[:100])))

print()

print("K=320:")
print("".join(map(str, bits_320[:100])))
print()
print("=" * 90)
print("FULL BINARY PERIOD SCAN")
print("=" * 90)

# --------------------------------------------------
# Use the most reliable prefix
# K = 320
# --------------------------------------------------

K_scan = 320

bits_scan = get_binary_bits(K_scan)

# Ignore the initial transient
start = 50

print(f"Total fractional bits : {len(bits_scan)}")
print(f"Analysis starts at    : bit {start}")
print()


# --------------------------------------------------
# Scan periods 1 ... 150
# --------------------------------------------------

results = []

for p in range(1, 151):

    mismatches = 0
    compared = 0

    for i in range(start, len(bits_scan) - p):

        compared += 1

        if bits_scan[i] != bits_scan[i + p]:
            mismatches += 1

    matches = compared - mismatches
    match_percent = 100 * matches / compared

    results.append(
        (p, mismatches, compared, match_percent)
    )


# --------------------------------------------------
# Print all results
# --------------------------------------------------

print(
    f"{'Period':>8} "
    f"{'Mismatch':>10} "
    f"{'Compared':>10} "
    f"{'Match %':>10}"
)

print("-" * 45)

for p, mismatches, compared, match_percent in results:

    print(
        f"{p:>8} "
        f"{mismatches:>10} "
        f"{compared:>10} "
        f"{match_percent:>9.2f}%"
    )


# --------------------------------------------------
# Find best candidates
# --------------------------------------------------

print()
print("=" * 90)
print("BEST PERIOD CANDIDATES")
print("=" * 90)

best = sorted(
    results,
    key=lambda x: x[3],
    reverse=True
)

for p, mismatches, compared, match_percent in best[:10]:

    print(
        f"Period {p:>3} | "
        f"Match = {match_percent:>6.2f}% | "
        f"Mismatch = {mismatches:>3}"
    )
print()
print("=" * 90)
print("PERIOD CANDIDATE STABILITY TEST")
print("=" * 90)


# --------------------------------------------------
# 1. Generate enough primes
# --------------------------------------------------


K_values = [160, 320, 640, 1280]

max_K = max(K_values)

primes_large = list(
    sp.primerange(1, 200000)
)

gaps_large = [
    primes_large[i + 1] - primes_large[i]
    for i in range(max_K + 5)
]


# --------------------------------------------------
# 2. Generate binary digits
# --------------------------------------------------

def get_binary_bits_large(K_test):

    raw = []

    for m in range(K_test):

        h = gaps_large[m + 1] // 2

        raw.append(h)

    digits = raw.copy()

    carry = 0

    for m in range(K_test - 1, -1, -1):

        total = digits[m] + carry

        digits[m] = total % 2

        carry = total // 2

    return digits[1:]


# --------------------------------------------------
# 3. Period measurement
# --------------------------------------------------

candidate_periods = [
    141,
    44,
    113,
    80,
    103
]


def measure_period(bits, p, start=50):

    mismatches = 0
    compared = 0

    for i in range(start, len(bits) - p):

        compared += 1

        if bits[i] != bits[i + p]:

            mismatches += 1

    matches = compared - mismatches
    
    if compared == 0:
        return (
            matches,
            mismatches,
            compared,
            float("nan")
        )

    return (
        matches,
        mismatches,
        compared,
        100 * matches / compared
    )


# --------------------------------------------------
# 4. Test each K
# --------------------------------------------------

print()

print(
    f"{'K':>6} "
    f"{'Period':>8} "
    f"{'Match %':>10} "
    f"{'Mismatch':>10} "
    f"{'Compared':>10}"
)

print("-" * 55)


for K_test in K_values:

    bits = get_binary_bits_large(K_test)

    for p in candidate_periods:

        matches, mismatches, compared, percent = \
            measure_period(bits, p)

        print(
            f"{K_test:>6} "
            f"{p:>8} "
            f"{percent:>9.2f}% "
            f"{mismatches:>10} "
            f"{compared:>10}"
        )
print()
print("=" * 90)
print("PERIOD CANDIDATE STABILITY TEST")
print("=" * 90)




# --------------------------------------------------
# 2. Generate prime gaps
# --------------------------------------------------

gaps_large = [
    primes_large[i + 1] - primes_large[i]
    for i in range(max_K + 5)
]


# --------------------------------------------------
# 3. Generate binary digits
# --------------------------------------------------

def get_binary_bits_large(K_test):

    raw = []

    for m in range(K_test):

        h = gaps_large[m + 1] // 2

        raw.append(h)

    digits = raw.copy()

    carry = 0

    for m in range(K_test - 1, -1, -1):

        total = digits[m] + carry

        digits[m] = total % 2
        carry = total // 2

    return digits[1:]


# --------------------------------------------------
# 4. Period measurement
# --------------------------------------------------

candidate_periods = [
    141,
    44,
    113,
    80,
    103
]


def measure_period(bits, p, start=50):

    mismatches = 0
    compared = 0

    for i in range(start, len(bits) - p):

        compared += 1

        if bits[i] != bits[i + p]:
            mismatches += 1

    matches = compared - mismatches

    if compared == 0:
        return (
            matches,
            mismatches,
            compared,
            float("nan")
        )

    return (
        matches,
        mismatches,
        compared,
        100 * matches / compared
    )


# --------------------------------------------------
# 5. Test each K
# --------------------------------------------------

print()

print(
    f"{'K':>6} "
    f"{'Period':>8} "
    f"{'Match %':>10} "
    f"{'Mismatch':>10} "
    f"{'Compared':>10}"
)

print("-" * 55)


for K_test in K_values:

    bits = get_binary_bits_large(K_test)

    for p in candidate_periods:

        matches, mismatches, compared, percent = \
            measure_period(bits, p)

        print(
            f"{K_test:>6} "
            f"{p:>8} "
            f"{percent:>9.2f}% "
            f"{mismatches:>10} "
            f"{compared:>10}"
        )
print("\n" + "=" * 70)
print("FRACTIONAL PART / PRIME GAP TEST")
print("=" * 70)

T = []

for N in range(1, 22):
    value = 0.0

    for j in range(1, 100):
        value += gaps[N + j - 1] / (2 ** j)

    T.append(value)

for N in range(1, 10):
    M = int(T[N - 1])
    f = T[N - 1] - M

    next_T = T[N]

    predicted_gap = 2 * M - int(next_T)

    print(
        f"N={N:2d} | "
        f"T_N={T[N-1]:.10f} | "
        f"floor={M:2d} | "
        f"frac={f:.10f} | "
        f"g_(N+1)={gaps[N]:2d} | "
        f"predicted={predicted_gap:2d}"
    )
print("\n" + "=" * 95)
print("PRIME GAP / BINARY CARRY TEST")
print("=" * 95)

print(
    f"{'N':>3} "
    f"{'T_N':>14} "
    f"{'frac':>12} "
    f"{'bit':>5} "
    f"{'carry':>6} "
    f"{'g_(N+1)':>9} "
    f"{'floor(T_N)':>11} "
    f"{'floor(T_N+1)':>13}"
)

print("-" * 95)

for N in range(1, 21):

    M = int(T[N - 1])
    f = T[N - 1] - M

    # First binary digit of the fractional part
    bit = 1 if f >= 0.5 else 0

    # Carry produced by doubling the fractional part
    carry = bit

    next_floor = int(T[N])

    gap = gaps[N]

    print(
        f"{N:3d} "
        f"{T[N-1]:14.10f} "
        f"{f:12.10f} "
        f"{bit:5d} "
        f"{carry:6d} "
        f"{gap:9d} "
        f"{M:11d} "
        f"{next_floor:13d}"
    )
print("\n" + "=" * 80)
print("BINARY BIT / PRIME GAP STATISTICS")
print("=" * 80)

# ---------------------------------------------------------
# Generate T_N for N = 1 ... 1000
# ---------------------------------------------------------

T_stat = []

for N in range(1, 1001):
    value = 0.0

    for j in range(1, 100):
        value += gaps[N + j - 1] / (2 ** j)

    T_stat.append(value)


# ---------------------------------------------------------
# Extract binary first bits
# ---------------------------------------------------------

bits = []

for value in T_stat:
    frac = value - int(value)

    bit = 1 if frac >= 0.5 else 0

    bits.append(bit)


# ---------------------------------------------------------
# A. Count 0 / 1
# ---------------------------------------------------------

count_0 = bits.count(0)
count_1 = bits.count(1)

print("\nA. BINARY BIT DISTRIBUTION")
print("-" * 40)

print(f"0: {count_0} ({100 * count_0 / len(bits):.2f}%)")
print(f"1: {count_1} ({100 * count_1 / len(bits):.2f}%)")


# ---------------------------------------------------------
# B. Consecutive bit pairs
# ---------------------------------------------------------

pairs = {
    "00": 0,
    "01": 0,
    "10": 0,
    "11": 0
}

for i in range(len(bits) - 1):
    pair = str(bits[i]) + str(bits[i + 1])
    pairs[pair] += 1


print("\nB. CONSECUTIVE BINARY PAIRS")
print("-" * 40)

for pair in ["00", "01", "10", "11"]:
    print(f"{pair}: {pairs[pair]}")


# ---------------------------------------------------------
# C. Prime gap distribution conditioned on bit
# ---------------------------------------------------------

gap_by_bit = {
    0: {},
    1: {}
}

for N in range(1000):

    bit = bits[N]

    # g_(N+1)
    gap = gaps[N + 1]

    if gap not in gap_by_bit[bit]:
        gap_by_bit[bit][gap] = 0

    gap_by_bit[bit][gap] += 1


print("\nC. PRIME GAP DISTRIBUTION CONDITIONED ON BINARY BIT")
print("-" * 60)

for bit in [0, 1]:

    total = sum(gap_by_bit[bit].values())

    print(f"\nBit = {bit}   Total = {total}")

    for gap in sorted(gap_by_bit[bit]):
        count = gap_by_bit[bit][gap]

        print(
            f"gap = {gap:2d} : "
            f"{count:4d} "
            f"({100 * count / total:6.2f}%)"
        )


# ---------------------------------------------------------
# D. Specific gaps: 2, 4, 6, 8, 10
# ---------------------------------------------------------

print("\nD. CONDITIONAL PROBABILITIES")
print("-" * 60)

for gap in [2, 4, 6, 8, 10]:

    print(f"\nGap = {gap}")

    for bit in [0, 1]:

        total = sum(gap_by_bit[bit].values())

        count = gap_by_bit[bit].get(gap, 0)

        probability = 100 * count / total

        print(
            f"P(gap={gap} | bit={bit}) "
            f"= {probability:.2f}%"
        )
print("\n" + "=" * 80)
print("BINARY PREFIX PERIOD TEST")
print("=" * 80)


def longest_periodic_suffix(bits, L, max_period=None):

    sequence = bits[:L]

    if max_period is None:
        max_period = L // 2

    best_period = None
    best_matches = -1

    for p in range(1, max_period + 1):

        matches = 0
        compared = 0

        # Compare the suffix against itself shifted by p
        for i in range(p, L):
            compared += 1

            if sequence[i] == sequence[i - p]:
                matches += 1

        if compared > 0:
            match_rate = matches / compared

            if match_rate > best_matches:
                best_matches = match_rate
                best_period = p

    return best_period, best_matches


print(
    f"{'Length':>8} "
    f"{'Best Period':>14} "
    f"{'Match %':>12}"
)

print("-" * 40)


for L in [100, 200, 400, 800]:

    period, match_rate = longest_periodic_suffix(
        bits,
        L,
        L // 2
    )

    print(
        f"{L:8d} "
        f"{period:14d} "
        f"{100 * match_rate:11.2f}%"
    )
print("\n" + "=" * 80)
print("MONTE CARLO RANDOM BINARY BASELINE")
print("=" * 80)

import random


def best_period_match(sequence, max_period):

    L = len(sequence)

    best_period = None
    best_match = -1.0

    for p in range(1, max_period + 1):

        matches = 0
        compared = L - p

        for i in range(p, L):
            if sequence[i] == sequence[i - p]:
                matches += 1

        if compared > 0:
            match_rate = matches / compared

            if match_rate > best_match:
                best_match = match_rate
                best_period = p

    return best_period, best_match


# ---------------------------------------------------------
# Test only short periods
# ---------------------------------------------------------

L = 800
MAX_PERIOD = 50
TRIALS = 1000

real_sequence = bits[:L]

real_period, real_match = best_period_match(
    real_sequence,
    MAX_PERIOD
)

print("\nREAL T1 SEQUENCE")
print("-" * 40)

print(f"Length:       {L}")
print(f"Max period:   {MAX_PERIOD}")
print(f"Best period:  {real_period}")
print(f"Best match:   {100 * real_match:.2f}%")


# ---------------------------------------------------------
# Monte Carlo
# ---------------------------------------------------------

random_matches = []

for trial in range(TRIALS):

    random_sequence = [
        random.randint(0, 1)
        for _ in range(L)
    ]

    _, match = best_period_match(
        random_sequence,
        MAX_PERIOD
    )

    random_matches.append(match)


# ---------------------------------------------------------
# Statistics
# ---------------------------------------------------------

random_matches.sort()

average_random = sum(random_matches) / TRIALS

median_random = random_matches[TRIALS // 2]

percentile_95 = random_matches[int(0.95 * TRIALS)]

percentile_99 = random_matches[int(0.99 * TRIALS)]

maximum_random = random_matches[-1]


print("\nRANDOM BASELINE")
print("-" * 40)

print(f"Trials:       {TRIALS}")
print(f"Average:      {100 * average_random:.2f}%")
print(f"Median:       {100 * median_random:.2f}%")
print(f"95th pct:     {100 * percentile_95:.2f}%")
print(f"99th pct:     {100 * percentile_99:.2f}%")
print(f"Maximum:      {100 * maximum_random:.2f}%")


# ---------------------------------------------------------
# Compare real sequence with random baseline
# ---------------------------------------------------------

count_exceeded = sum(
    1
    for value in random_matches
    if value >= real_match
)

empirical_p = (
    (count_exceeded + 1)
    / (TRIALS + 1)
)


print("\nCOMPARISON")
print("-" * 40)

print(
    f"Random sequences >= real: "
    f"{count_exceeded}/{TRIALS}"
)

print(
    f"Empirical p-value: "
    f"{empirical_p:.4f}"
)
print("\n" + "=" * 80)
print("EXPERIMENT 05 - CONDITIONAL ENTROPY")
print("=" * 80)

from collections import defaultdict, Counter
import math
import random


# ---------------------------------------------------------
# 1. Count transitions:
#
#     g_n  ->  g_(n+1)
#
# ---------------------------------------------------------

transition_counts = defaultdict(Counter)
gap_counts = Counter()

for i in range(len(gaps) - 1):

    current_gap = gaps[i]
    next_gap = gaps[i + 1]

    transition_counts[current_gap][next_gap] += 1
    gap_counts[current_gap] += 1


# ---------------------------------------------------------
# 2. Calculate conditional entropy
#
# H(X | Y=y)
#
#     = - sum_x P(x|y) log2 P(x|y)
# ---------------------------------------------------------

def conditional_entropy_for_gap(current_gap):

    counts = transition_counts[current_gap]

    total = sum(counts.values())

    entropy = 0.0

    for count in counts.values():

        probability = count / total

        entropy -= probability * math.log2(probability)

    return entropy


# ---------------------------------------------------------
# 3. Overall conditional entropy
#
# H(g_(n+1) | g_n)
#
# = sum_y P(y) H(g_(n+1) | g_n=y)
# ---------------------------------------------------------

total_transitions = sum(gap_counts.values())

conditional_entropy = 0.0

for current_gap, count in gap_counts.items():

    probability_current = count / total_transitions

    entropy_current = conditional_entropy_for_gap(current_gap)

    conditional_entropy += (
        probability_current * entropy_current
    )


print()
print("REAL PRIME-GAP SEQUENCE")
print("-" * 50)

print(
    f"H(g_(n+1) | g_n) = "
    f"{conditional_entropy:.6f} bits"
)


# ---------------------------------------------------------
# 4. Compare with shuffled gap sequence
# ---------------------------------------------------------

shuffled_gaps = gaps.copy()

random.seed(251)

random.shuffle(shuffled_gaps)


# Count transitions in shuffled sequence

shuffled_transition_counts = defaultdict(Counter)
shuffled_gap_counts = Counter()

for i in range(len(shuffled_gaps) - 1):

    current_gap = shuffled_gaps[i]
    next_gap = shuffled_gaps[i + 1]

    shuffled_transition_counts[current_gap][next_gap] += 1
    shuffled_gap_counts[current_gap] += 1


# Calculate shuffled conditional entropy

shuffled_total = sum(shuffled_gap_counts.values())

shuffled_entropy = 0.0

for current_gap, count in shuffled_gap_counts.items():

    probability_current = count / shuffled_total

    counts = shuffled_transition_counts[current_gap]

    total = sum(counts.values())

    entropy_current = 0.0

    for next_count in counts.values():

        probability = next_count / total

        entropy_current -= (
            probability * math.log2(probability)
        )

    shuffled_entropy += (
        probability_current * entropy_current
    )


print()
print("SHUFFLED PRIME-GAP SEQUENCE")
print("-" * 50)

print(
    f"H(g_(n+1) | g_n) = "
    f"{shuffled_entropy:.6f} bits"
)


# ---------------------------------------------------------
# 5. Difference
# ---------------------------------------------------------

difference = shuffled_entropy - conditional_entropy

print()
print("COMPARISON")
print("-" * 50)

print(
    f"Real entropy     : "
    f"{conditional_entropy:.6f} bits"
)

print(
    f"Shuffled entropy : "
    f"{shuffled_entropy:.6f} bits"
)

print(
    f"Difference       : "
    f"{difference:.6f} bits"
)

print()
print("=" * 80)
print("EXPERIMENT 05 COMPLETE")
print("=" * 80)
# ============================================================
# EXPERIMENT 05B - SHUFFLED BASELINE (100 RUNS)
# ============================================================

print("\n" + "=" * 80)
print("EXPERIMENT 05B - SHUFFLED BASELINE")
print("=" * 80)

SHUFFLE_RUNS = 100

shuffle_entropies = []

for run in range(SHUFFLE_RUNS):

    test_gaps = gaps.copy()
    random.shuffle(test_gaps)

    transition_counts_test = defaultdict(Counter)
    gap_counts_test = Counter()

    for i in range(len(test_gaps) - 1):

        current_gap = test_gaps[i]
        next_gap = test_gaps[i + 1]

        transition_counts_test[current_gap][next_gap] += 1
        gap_counts_test[current_gap] += 1

    total_test = sum(gap_counts_test.values())

    entropy_test = 0.0

    for current_gap, count in gap_counts_test.items():

        probability_current = count / total_test

        counts = transition_counts_test[current_gap]
        total_current = sum(counts.values())

        entropy_current = 0.0

        for next_count in counts.values():

            probability = next_count / total_current

            entropy_current -= (
                probability * math.log2(probability)
            )

        entropy_test += (
            probability_current * entropy_current
        )

    shuffle_entropies.append(entropy_test)


# ------------------------------------------------------------
# Statistics
# ------------------------------------------------------------

shuffle_mean = sum(shuffle_entropies) / len(shuffle_entropies)

shuffle_variance = sum(
    (x - shuffle_mean) ** 2
    for x in shuffle_entropies
) / len(shuffle_entropies)

shuffle_std = math.sqrt(shuffle_variance)

shuffle_min = min(shuffle_entropies)
shuffle_max = max(shuffle_entropies)

real_entropy = conditional_entropy


print()
print(f"Number of shuffles : {SHUFFLE_RUNS}")
print(f"Real entropy       : {real_entropy:.6f}")
print(f"Shuffle mean       : {shuffle_mean:.6f}")
print(f"Shuffle std        : {shuffle_std:.6f}")
print(f"Shuffle minimum    : {shuffle_min:.6f}")
print(f"Shuffle maximum    : {shuffle_max:.6f}")

print()
print(
    f"Real - Shuffle mean = "
    f"{real_entropy - shuffle_mean:.6f} bits"
)

print(
    f"Z-score = "
    f"{(real_entropy - shuffle_mean) / shuffle_std:.3f}"
)

print()
print("=" * 80)
print("EXPERIMENT 05B COMPLETE")
print("=" * 80)
# ============================================================
# EXPERIMENT 06 - SECOND-ORDER CONDITIONAL ENTROPY
#
# H(g_(n+2) | g_n, g_(n+1))
# ============================================================

print("\n" + "=" * 80)
print("EXPERIMENT 06 - SECOND-ORDER CONDITIONAL ENTROPY")
print("=" * 80)


# ------------------------------------------------------------
# 1. Count:
#
# (g_n, g_(n+1)) -> g_(n+2)
# ------------------------------------------------------------

pair_transition_counts = defaultdict(Counter)
pair_counts = Counter()

for i in range(len(gaps) - 2):

    gap_a = gaps[i]
    gap_b = gaps[i + 1]
    gap_c = gaps[i + 2]

    pair = (gap_a, gap_b)

    pair_transition_counts[pair][gap_c] += 1
    pair_counts[pair] += 1


# ------------------------------------------------------------
# 2. Calculate:
#
# H(g_(n+2) | g_n, g_(n+1))
# ------------------------------------------------------------

total_pairs = sum(pair_counts.values())

second_order_entropy = 0.0

for pair, count in pair_counts.items():

    probability_pair = count / total_pairs

    next_counts = pair_transition_counts[pair]
    total_next = sum(next_counts.values())

    entropy_pair = 0.0

    for next_count in next_counts.values():

        probability = next_count / total_next

        entropy_pair -= (
            probability * math.log2(probability)
        )

    second_order_entropy += (
        probability_pair * entropy_pair
    )


# ------------------------------------------------------------
# 3. Compare with first-order entropy
# ------------------------------------------------------------

print()
print("REAL PRIME-GAP SEQUENCE")
print("-" * 50)

print(
    f"H(g_(n+1) | g_n)               = "
    f"{conditional_entropy:.6f} bits"
)

print(
    f"H(g_(n+2) | g_n, g_(n+1))     = "
    f"{second_order_entropy:.6f} bits"
)


information_gain = (
    conditional_entropy - second_order_entropy
)

print()
print(
    f"Additional information from "
    f"second gap = {information_gain:.6f} bits"
)


# ------------------------------------------------------------
# 4. Shuffled baseline
# ------------------------------------------------------------

SHUFFLE_RUNS_2 = 100

second_shuffle_entropies = []

for run in range(SHUFFLE_RUNS_2):

    test_gaps = gaps.copy()
    random.shuffle(test_gaps)

    test_pair_transition_counts = defaultdict(Counter)
    test_pair_counts = Counter()

    for i in range(len(test_gaps) - 2):

        gap_a = test_gaps[i]
        gap_b = test_gaps[i + 1]
        gap_c = test_gaps[i + 2]

        pair = (gap_a, gap_b)

        test_pair_transition_counts[pair][gap_c] += 1
        test_pair_counts[pair] += 1

    total_test_pairs = sum(test_pair_counts.values())

    entropy_test = 0.0

    for pair, count in test_pair_counts.items():

        probability_pair = count / total_test_pairs

        next_counts = test_pair_transition_counts[pair]
        total_next = sum(next_counts.values())

        entropy_pair = 0.0

        for next_count in next_counts.values():

            probability = next_count / total_next

            entropy_pair -= (
                probability * math.log2(probability)
            )

        entropy_test += (
            probability_pair * entropy_pair
        )

    second_shuffle_entropies.append(entropy_test)


# ------------------------------------------------------------
# 5. Statistics
# ------------------------------------------------------------

second_shuffle_mean = (
    sum(second_shuffle_entropies)
    / len(second_shuffle_entropies)
)

second_shuffle_variance = sum(
    (x - second_shuffle_mean) ** 2
    for x in second_shuffle_entropies
) / len(second_shuffle_entropies)

second_shuffle_std = math.sqrt(
    second_shuffle_variance
)

second_shuffle_min = min(
    second_shuffle_entropies
)

second_shuffle_max = max(
    second_shuffle_entropies
)


print()
print("SHUFFLED BASELINE")
print("-" * 50)

print(
    f"Number of shuffles              = "
    f"{SHUFFLE_RUNS_2}"
)

print(
    f"Real second-order entropy       = "
    f"{second_order_entropy:.6f}"
)

print(
    f"Shuffle mean                    = "
    f"{second_shuffle_mean:.6f}"
)

print(
    f"Shuffle std                     = "
    f"{second_shuffle_std:.6f}"
)

print(
    f"Shuffle minimum                 = "
    f"{second_shuffle_min:.6f}"
)

print(
    f"Shuffle maximum                 = "
    f"{second_shuffle_max:.6f}"
)


# ------------------------------------------------------------
# 6. Z-score
# ------------------------------------------------------------

second_z_score = (
    second_order_entropy
    - second_shuffle_mean
) / second_shuffle_std


print()
print(
    f"Real - Shuffle mean             = "
    f"{second_order_entropy - second_shuffle_mean:.6f}"
)

print(
    f"Z-score                         = "
    f"{second_z_score:.3f}"
)


print()
print("=" * 80)
print("EXPERIMENT 06 COMPLETE")
print("=" * 80)