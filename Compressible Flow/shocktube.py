"""
SHOCK TUBE CALCULATOR for TI-84 Plus CE Python (MicroPython).
Calculates incident and reflected shock parameters from diaphragm ratio,
shock pressure ratio, or incident shock Mach number.
"""
import math
import sys


def solve(func, target, x0, x1):
    """Secant solver to invert monotonic relations."""
    for _ in range(40):
        f0, f1 = func(x0) - target, func(x1) - target
        if abs(f1) < 1e-6 or abs(x1 - x0) < 1e-6:
            return x1
        x0, x1 = x1, x1 - f1 * (x1 - x0) / (f1 - f0)
    return x1


def ask(prompt, default):
    """Helper for numeric input with a default on [ENTER]."""
    val = input("%s [%s]: " % (prompt, default))
    return float(val) if val else float(default)


def p4_p1_from_p2_p1(a1, a4, g1, g4, p2_p1):
    """Shock tube relation: diaphragm ratio P4/P1 from shock ratio P2/P1."""
    denom = math.sqrt(2 * g1 * (2 * g1 + (g1 + 1) * (p2_p1 - 1)))
    factor = 1 - ((g4 - 1) * (a1 / a4) * (p2_p1 - 1)) / denom
    if factor <= 0:
        return 1e9
    return p2_p1 * (factor ** ((-2 * g4) / (g4 - 1)))


# --- User Interface ---
print("--- SHOCK TUBE ---")
print("1:P2/P1  2:P4/P1  3:Ms")
opt = int(input("Known [1-3]: "))
val = float(input("Value: "))

print("\n-- Driven Gas (1) --")
g1 = ask("gamma1", 1.4)
T1 = ask("T1 (K)", 300)
R1 = ask("R1", 287)

print("-- Driver Gas (4) --")
g4 = ask("gamma4", 1.4)
T4 = ask("T4 (K)", 300)
R4 = ask("R4", 287)

# Speed of sound in driver (4) and driven (1) sections
a1 = math.sqrt(g1 * R1 * T1)
a4 = math.sqrt(g4 * R4 * T4)

# Resolve P2/P1 and P4/P1
if opt == 1:
    p2_p1 = val
    if p2_p1 < 1.0:
        print("Error: P2/P1 >= 1")
        sys.exit()
    p4_p1 = p4_p1_from_p2_p1(a1, a4, g1, g4, p2_p1)
elif opt == 2:
    p4_p1 = val
    if p4_p1 < 1.0:
        print("Error: P4/P1 >= 1")
        sys.exit()
    p2_p1 = solve(lambda p: p4_p1_from_p2_p1(a1, a4, g1, g4, p), p4_p1, 1.5, 4.0)
elif opt == 3:
    ms = val
    if ms < 1.0:
        print("Error: Ms >= 1")
        sys.exit()
    p2_p1 = 1 + ((2 * g1) / (g1 + 1)) * (ms**2 - 1)
    p4_p1 = p4_p1_from_p2_p1(a1, a4, g1, g4, p2_p1)
else:
    print("error")
    sys.exit()

# Incident shock Mach number Ms, shock speed W, induced flow velocity up
ms = math.sqrt((p2_p1 - 1) * (g1 + 1) / (2 * g1) + 1)
w = ms * a1
up = (2 * a1 / (g1 + 1)) * (ms - 1 / ms)

# Reflected Shock Mach number Mr and pressure ratios
if ms <= 1.0001:
    mr = 1.0
    p5_p2 = 1.0
else:
    term = math.sqrt(
        1 + (2 * (g1 - 1) / (g1 + 1) ** 2) * (ms**2 - 1) * (g1 + 1 / ms**2)
    )
    a_quad = (ms / (ms**2 - 1)) * term
    mr = (1 + math.sqrt(1 + 4 * a_quad**2)) / (2 * a_quad)
    p5_p2 = 1 + ((2 * g1) / (g1 + 1)) * (mr**2 - 1)

p5_p1 = p5_p2 * p2_p1

# Display results formatted for TI-84 screen (26x10 chars)
print("\nP4/P1: %.4f" % p4_p1)
print("P2/P1: %.4f" % p2_p1)
print("P5/P1: %.4f" % p5_p1)
print("P5/P2: %.4f" % p5_p2)
print("Ms   : %.4f" % ms)
print("Mr   : %.4f" % mr)
print("W    : %.1f m/s" % w)
print("up   : %.1f m/s" % up)
