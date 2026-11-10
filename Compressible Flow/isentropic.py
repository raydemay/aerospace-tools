"""
ISENTROPIC FLOW CALCULATOR for TI-84 Plus CE Python (MicroPython).
Calculates 1D isentropic flow properties given any one flow parameter.
"""
import math
import sys


def solve(func, target, x0, x1):
    """Secant solver to invert monotonic compressible flow relations."""
    for _ in range(40):
        f0, f1 = func(x0) - target, func(x1) - target
        if abs(f1) < 1e-6 or abs(x1 - x0) < 1e-6:
            return x1
        x0, x1 = x1, x1 - f1 * (x1 - x0) / (f1 - f0)
    return x1


# --- Forward Isentropic Relations (given M, gamma) ---
def t0_t(g, m):
    """Stagnation-to-static temperature ratio T0/T."""
    return 1 + 0.5 * (g - 1) * m**2


def p0_p(g, m):
    """Stagnation-to-static pressure ratio P0/P."""
    return t0_t(g, m) ** (g / (g - 1))


def r0_r(g, m):
    """Stagnation-to-static density ratio rho0/rho."""
    return t0_t(g, m) ** (1 / (g - 1))


def a_astar(g, m):
    """Area ratio A/A* from Mach number."""
    if m <= 0:
        return 1e9
    return math.sqrt(
        (1 / m**2) * ((2 / (g + 1)) * t0_t(g, m)) ** ((g + 1) / (g - 1))
    )


def prandtl_meyer(g, m):
    """Prandtl-Meyer angle nu (deg) for supersonic flow."""
    if m <= 1.0:
        return 0.0
    term1 = math.sqrt((g + 1) / (g - 1)) * math.atan(
        math.sqrt((g - 1) / (g + 1) * (m**2 - 1))
    )
    term2 = math.atan(math.sqrt(m**2 - 1))
    return (term1 - term2) * (180 / math.pi)


# --- User Interface ---
print("--- ISENTROPIC FLOW ---")
gin = input("gamma [1.4]: ")
g = float(gin) if gin else 1.4

print("1:Mach  2:T0/T  3:P0/P\n4:r0/r  5:A/A*  6:mu\n7:Prandtl-Meyer")
opt = int(input("Known [1-7]: "))
val = float(input("Value: "))

# Resolve Mach number M
if opt == 1:
    m = val
elif opt == 2:
    m = math.sqrt(max(0.0, (val - 1) * 2 / (g - 1)))
elif opt == 3:
    m = math.sqrt(max(0.0, ((val ** ((g - 1) / g)) - 1) * 2 / (g - 1)))
elif opt == 4:
    m = math.sqrt(max(0.0, ((val ** (g - 1)) - 1) * 2 / (g - 1)))
elif opt == 5:
    c_in = input("1:Sub 2:Sup [2]: ")
    if c_in == "1":
        m = solve(lambda x: a_astar(g, x), val, 0.2, 0.4)
    else:
        m = solve(lambda x: a_astar(g, x), val, 1.5, 2.5)
elif opt == 6:
    m = 1 / math.sin(val * math.pi / 180)
elif opt == 7:
    m = solve(lambda x: prandtl_meyer(g, x), val, 1.5, 2.5)
else:
    print("error")
    sys.exit()

# Compute all flow ratios
t0t = t0_t(g, m)
p0p = p0_p(g, m)
r0r = r0_r(g, m)
a_ast = a_astar(g, m)
mu = (math.asin(1 / m) * 180 / math.pi) if m >= 1.0 else 90.0
nu = prandtl_meyer(g, m)

# Display results formatted for TI-84 screen (26x10 chars)
print("\nMach    : %.4f" % m)
print("T0/T    : %.4f" % t0t)
print("P0/P    : %.4f" % p0p)
print("rho0/rho: %.4f" % r0r)
print("A/A*    : %.4f" % a_ast)
print("Mach ang: %.2f deg" % mu)
print("PM angle: %.2f deg" % nu)
