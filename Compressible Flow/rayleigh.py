"""
RAYLEIGH FLOW CALCULATOR for TI-84 Plus CE Python (MicroPython).
Calculates 1D frictionless flow with heat addition/rejection in a constant-area duct given any one parameter.
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


# --- Forward Rayleigh Flow Relations (given M, gamma) ---
def t0_t0star(g, m):
    """Stagnation temperature ratio T0/T0*."""
    return ((g + 1) * m**2 * (2 + (g - 1) * m**2)) / ((1 + g * m**2) ** 2)


def t_tstar(g, m):
    """Static temperature ratio T/T*."""
    return ((g + 1) * m / (1 + g * m**2)) ** 2


def p_pstar(g, m):
    """Static pressure ratio P/P*."""
    return (g + 1) / (1 + g * m**2)


def p0_p0star(g, m):
    """Stagnation pressure ratio P0/P0*."""
    return ((g + 1) / (1 + g * m**2)) * (
        ((2 + (g - 1) * m**2) / (g + 1)) ** (g / (g - 1))
    )


def rho_rhostar(g, m):
    """Density ratio rho/rho*."""
    return (1 + g * m**2) / ((g + 1) * m**2)


def u_ustar(g, m):
    """Velocity ratio u/u*."""
    return ((g + 1) * m**2) / (1 + g * m**2)


# --- User Interface ---
print("--- RAYLEIGH FLOW ---")
print("1:M      2:T0/T0* 3:T/T*")
print("4:P/P*   5:P0/P0* 6:u/u*")
opt = int(input("Known [1-6]: "))
val = float(input("Value: "))
g_in = input("gamma [1.4]: ")
g = float(g_in) if g_in else 1.4

# Solve for Mach number M
if opt == 1:
    m = val
elif opt == 2:
    # T0/T0* has sub and super branches (max T0/T0* = 1 at M=1)
    if val > 1.0:
        print("Error: T0/T0* <= 1")
        sys.exit()
    sub = int(input("1:Sub 2:Super: "))
    x0, x1 = (0.2, 0.6) if sub == 1 else (1.5, 2.5)
    m = solve(lambda m: t0_t0star(g, m), val, x0, x1)
elif opt == 3:
    # T/T* max is at M = 1/sqrt(g)
    sub = int(input("1:Sub 2:Super: "))
    x0, x1 = (0.2, 0.5) if sub == 1 else (1.3, 2.0)
    m = solve(lambda m: t_tstar(g, m), val, x0, x1)
elif opt == 4:
    # P/P* is monotonic: P/P* = (g+1)/(1+g*M^2)
    arg = ((g + 1) / val - 1) / g
    if arg < 0:
        print("error")
        sys.exit()
    m = math.sqrt(arg)
elif opt == 5:
    # P0/P0* has sub and super branches
    sub = int(input("1:Sub 2:Super: "))
    x0, x1 = (0.2, 0.6) if sub == 1 else (1.5, 2.5)
    m = solve(lambda m: p0_p0star(g, m), val, x0, x1)
elif opt == 6:
    # u/u* = (g+1)M^2 / (1+g*M^2)
    denom = (g + 1) - g * val
    if denom <= 0 or val < 0:
        print("error")
        sys.exit()
    m = math.sqrt(val / denom)
else:
    print("error")
    sys.exit()

if m <= 0:
    print("error")
    sys.exit()

# Display results formatted for TI-84 screen (26x10 chars)
print("\nM     : %.4f" % m)
print("T0/T0*: %.4f" % t0_t0star(g, m))
print("T/T*  : %.4f" % t_tstar(g, m))
print("P/P*  : %.4f" % p_pstar(g, m))
print("P0/P0*: %.4f" % p0_p0star(g, m))
print("r/r*  : %.4f" % rho_rhostar(g, m))
print("u/u*  : %.4f" % u_ustar(g, m))
