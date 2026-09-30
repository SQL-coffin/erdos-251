import math
import time
import sympy as sp


# ============================================================
# Experiment 02
# Erd?s Problem #251
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

required_primes = N_MAX + K + 2

# Generate enough primes with SymPy.
primes = list(sp.primerange(1, 2_000_000))

if len(primes) < required_primes:
    raise RuntimeError(
        f"Not enough primes generated: {len(primes)} < {required_primes}"
    )

primes = primes[:required_primes]

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

# ?¹æ? 1ï¼šä½¿?¨æ?ä»¬å??¥ç? compute_T_exact()
T1_method_1 = compute_T_exact(1)

# ?¹æ? 2ï¼šå??¨é€é¡¹?¸å?
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

    # No valid comparisons if the period reaches beyond the
    # available suffix.
    if start >= len(bits) or start + p >= len(bits):
        return (0, 0, 0, float("nan"))

    mismatches = 0
    compared = 0

    for i in range(start, len(bits) - p):

        compared += 1

        if bits[i] != bits[i + p]:
            mismatches += 1

    matches = compared - mismatches
    match_percent = 100 * matches / compared

    return (
        matches,
        mismatches,
        compared,
        match_percent
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

        percent_text = "N/A" if compared == 0 else f"{percent:>9.2f}%"

        print(
            f"{K_test:>6} "
            f"{p:>8} "
            f"{percent_text:>10} "
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
