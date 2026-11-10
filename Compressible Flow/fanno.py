"""
FANNO FLOW CALCULATOR for TI-84 Plus CE Python (MicroPython).
Calculates 1D adiabatic flow with friction in a constant-area duct given any one parameter.
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


# --- Forward Fanno Flow Relations (given M, gamma) ---
def t_tstar(g, m):
    """Static temperature ratio T/T*."""
    return (g + 1) / (2 + (g - 1) * m**2)


def p_pstar(g, m):
    """Static pressure ratio P/P*."""
    return (1 / m) * math.sqrt(t_tstar(g, m))


def p0_p0star(g, m):
    """Stagnation pressure ratio P0/P0*."""
    return (1 / m) * (((2 + (g - 1) * m**2) / (g + 1)) ** ((g + 1) / (2 * (g - 1))))


def rho_rhostar(g, m):
    """Density ratio rho/rho*."""
    return p_pstar(g, m) / t_tstar(g, m)


def u_ustar(g, m):
    """Velocity ratio u/u*."""
    return m * math.sqrt(t_tstar(g, m))


def fL_d(g, m):
    """Friction parameter 4fL*/D to reach sonic choke point."""
    return ((1 - m**2) / (g * m**2)) + ((g + 1) / (2 * g)) * math.log(
        ((g + 1) * m**2) / (2 + (g - 1) * m**2)
    )


# --- User Interface ---
print("--- FANNO FLOW ---")
print("1:M      2:T/T*   3:P/P*")
print("4:P0/P0* 5:4fL*/D 6:u/u*")
opt = int(input("Known [1-6]: "))
val = float(input("Value: "))
g_in = input("gamma [1.4]: ")
g = float(g_in) if g_in else 1.4

# Solve for Mach number M
if opt == 1:
    m = val
elif opt == 2:
    # T/T* = (g+1) / (2 + (g-1)M^2)
    arg = ((g + 1) / val - 2) / (g - 1)
    if arg < 0:
        print("error")
        sys.exit()
    m = math.sqrt(arg)
elif opt == 3:
    # P/P* is monotonic
    m = solve(lambda m: p_pstar(g, m), val, 0.5, 1.5)
elif opt == 4:
    # P0/P0* has sub and super branches
    sub = int(input("1:Sub 2:Super: "))
    x0, x1 = (0.2, 0.6) if sub == 1 else (1.5, 2.5)
    m = solve(lambda m: p0_p0star(g, m), val, x0, x1)
elif opt == 5:
    # 4fL*/D has sub and super branches
    sub = int(input("1:Sub 2:Super: "))
    x0, x1 = (0.3, 0.7) if sub == 1 else (1.5, 2.5)
    m = solve(lambda m: fL_d(g, m), val, x0, x1)
elif opt == 6:
    # u/u*
    denom = (g + 1) - (g - 1) * val**2
    if denom <= 0:
        print("error")
        sys.exit()
    m = math.sqrt((2 * val**2) / denom)
else:
    print("error")
    sys.exit()

if m <= 0:
    print("error")
    sys.exit()

# Display results formatted for TI-84 screen (26x10 chars)
print("\nM     : %.4f" % m)
print("T/T*  : %.4f" % t_tstar(g, m))
print("P/P*  : %.4f" % p_pstar(g, m))
print("P0/P0*: %.4f" % p0_p0star(g, m))
print("r/r*  : %.4f" % rho_rhostar(g, m))
print("4fL*/D: %.4f" % fL_d(g, m))
print("u/u*  : %.4f" % u_ustar(g, m))
