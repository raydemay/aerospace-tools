"""
OBLIQUE SHOCK CALCULATOR for TI-84 Plus CE Python (MicroPython).
Calculates flow properties across an oblique shock wave given any one parameter.
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


# --- Oblique Shock Relations ---
# Normal shock functions parameterized by normal Mach number Mn1
def p_ratio(g, mn):
    """Static pressure ratio P2/P1."""
    return 1 + (2 * g / (g + 1)) * (mn**2 - 1)


def r_ratio(g, mn):
    """Density ratio rho2/rho1."""
    return (g + 1) * mn**2 / (2 + (g - 1) * mn**2)


def t_ratio(g, mn):
    """Temperature ratio T2/T1."""
    return p_ratio(g, mn) / r_ratio(g, mn)


def mn2_from_mn1(g, mn):
    """Downstream normal Mach number Mn2."""
    return math.sqrt((1 + 0.5 * (g - 1) * mn**2) / (g * mn**2 - 0.5 * (g - 1)))


def p0_ratio(g, mn):
    """Total pressure ratio across shock p02/p01 from entropy change."""
    ds = (g * 287 / (g - 1)) * math.log(t_ratio(g, mn)) - 287 * math.log(
        p_ratio(g, mn)
    )
    return math.exp(-ds / 287)


def theta_from_beta(m1, g, b):
    """Deflection angle theta (rad) from M1 and wave angle beta (rad)."""
    if b >= math.pi / 2:
        return 0.0
    num = 2 * (1 / math.tan(b)) * (m1**2 * math.sin(b) ** 2 - 1)
    den = m1**2 * (g + math.cos(2 * b)) + 2
    return math.atan(num / den)


def m2_from_mn1(m1, g, mn1):
    """Downstream Mach number M2 from normal component Mn1."""
    sin_b = min(1.0, mn1 / m1)
    b = math.asin(sin_b)
    th = theta_from_beta(m1, g, b)
    mn2 = mn2_from_mn1(g, mn1)
    return mn2 / math.sin(b - th) if (b - th) > 1e-5 else mn2


# --- User Interface ---
print("--- OBLIQUE SHOCK ---")
m1 = float(input("M1: "))
gin = input("gamma [1.4]: ")
g = float(gin) if gin else 1.4

print("1:Beta  2:Theta 3:P2/P1\n4:r2/r1 5:T2/T1 6:p02/p01\n7:p02/p1 8:M2")
opt = int(input("Known [1-8]: "))
val = float(input("Value (deg if angle): "))

# Resolve normal Mach number Mn1 and wave angle beta (rad)
mu = math.asin(1 / m1)  # Mach angle (minimum physical wave angle)

if opt == 1:
    beta = val * math.pi / 180
    mn1 = m1 * math.sin(beta)
elif opt == 2:
    th = val * math.pi / 180
    c_in = input("1:Weak 2:Strong [1]: ")
    is_strong = c_in == "2"
    if is_strong:
        beta = solve(
            lambda b: theta_from_beta(m1, g, b),
            th,
            math.pi / 2 - 0.05,
            math.pi / 2 - 0.2,
        )
    else:
        beta = solve(lambda b: theta_from_beta(m1, g, b), th, mu + 0.05, mu + 0.2)
    mn1 = m1 * math.sin(beta)
elif opt == 8:
    mn1 = solve(lambda mn: m2_from_mn1(m1, g, mn), val, 1.1, 1.6)
    beta = math.asin(min(1.0, mn1 / m1))
else:
    lookup = {
        3: p_ratio,
        4: r_ratio,
        5: t_ratio,
        6: p0_ratio,
        7: lambda g, mn: p0_ratio(g, mn)
        * ((1 + 0.5 * (g - 1) * m1**2) ** (g / (g - 1))),
    }
    mn1 = solve(lambda mn: lookup[opt](g, mn), val, 1.1, 2.0)
    beta = math.asin(min(1.0, mn1 / m1))

# Downstream state properties
mn2 = mn2_from_mn1(g, mn1)
p2p1 = p_ratio(g, mn1)
r2r1 = r_ratio(g, mn1)
t2t1 = p2p1 / r2r1
p02p01 = p0_ratio(g, mn1)
p02p1 = p02p01 * ((1 + 0.5 * (g - 1) * m1**2) ** (g / (g - 1)))

theta = theta_from_beta(m1, g, beta) if opt != 2 else (val * math.pi / 180)
m2 = mn2 / math.sin(beta - theta) if (beta - theta) > 1e-5 else mn2

# Display results formatted for TI-84 screen (26x10 chars)
print("\nM1: %.4f  M2: %.4f" % (m1, m2))
print("Wave beta : %.2f deg" % (beta * 180 / math.pi))
print("Defl theta: %.2f deg" % (theta * 180 / math.pi))
print("P2/P1   : %.4f" % p2p1)
print("rho2/rho: %.4f" % r2r1)
print("T2/T1   : %.4f" % t2t1)
print("p02/p01 : %.4f" % p02p01)
print("p02/p1  : %.4f" % p02p1)
