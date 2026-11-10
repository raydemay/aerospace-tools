"""
NORMAL SHOCK CALCULATOR for TI-84 Plus CE Python (MicroPython).
Calculates flow properties across a normal shock wave given any one ratio.
"""
import math
import sys


def solve(func, target, x0=1.2, x1=2.0):
    """Secant solver to invert monotonic compressible flow relations."""
    for _ in range(40):
        f0, f1 = func(x0) - target, func(x1) - target
        if abs(f1) < 1e-6 or abs(x1 - x0) < 1e-6:
            return x1
        x0, x1 = x1, x1 - f1 * (x1 - x0) / (f1 - f0)
    return x1


# --- Normal Shock Forward Relations (given M1, gamma) ---
def p_ratio(g, m):
    """Static pressure ratio P2/P1."""
    return 1 + (2 * g / (g + 1)) * (m**2 - 1)


def r_ratio(g, m):
    """Density ratio rho2/rho1."""
    return (g + 1) * m**2 / (2 + (g - 1) * m**2)


def t_ratio(g, m):
    """Temperature ratio T2/T1."""
    return p_ratio(g, m) / r_ratio(g, m)


def m2_from_m1(g, m):
    """Downstream Mach number M2."""
    return math.sqrt((1 + 0.5 * (g - 1) * m**2) / (g * m**2 - 0.5 * (g - 1)))


def p0_ratio(g, m):
    """Total pressure ratio p02/p01 from entropy increase."""
    ds = (g * 287 / (g - 1)) * math.log(t_ratio(g, m)) - 287 * math.log(
        p_ratio(g, m)
    )
    return math.exp(-ds / 287)


def p02_p1_ratio(g, m):
    """Downstream total to upstream static pressure ratio p02/p1."""
    p01_p1 = (1 + 0.5 * (g - 1) * m**2) ** (g / (g - 1))
    return p0_ratio(g, m) * p01_p1


# --- User Interface ---
print("--- NORMAL SHOCK ---")
gin = input("gamma [1.4]: ")
g = float(gin) if gin else 1.4

print("1:M1    2:M2    3:P2/P1\n4:r2/r1 5:T2/T1 6:p02/p01\n7:p02/p1")
opt = int(input("Known [1-7]: "))
val = float(input("Value: "))

# Resolve upstream Mach number M1
if opt == 1:
    m1 = val
elif opt == 2:
    if val >= 1.0:
        print("M2 must be < 1.0")
        sys.exit()
    m1 = solve(lambda m: m2_from_m1(g, m), val, 1.2, 2.5)
else:
    lookup = {3: p_ratio, 4: r_ratio, 5: t_ratio, 6: p0_ratio, 7: p02_p1_ratio}
    m1 = solve(lambda m: lookup[opt](g, m), val, 1.2, 2.5)

# Calculate downstream properties
m2 = m2_from_m1(g, m1)
p2p1 = p_ratio(g, m1)
r2r1 = r_ratio(g, m1)
t2t1 = p2p1 / r2r1
p02p01 = p0_ratio(g, m1)
p02p1 = p02_p1_ratio(g, m1)

# Display results formatted for TI-84 screen (26x10 chars)
print("\nM1: %.4f  M2: %.4f" % (m1, m2))
print("P2/P1   : %.4f" % p2p1)
print("rho2/rho: %.4f" % r2r1)
print("T2/T1   : %.4f" % t2t1)
print("p02/p01 : %.4f" % p02p01)
print("p02/p1  : %.4f" % p02p1)
