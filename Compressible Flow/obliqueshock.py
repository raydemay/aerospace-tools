import math
import sys

print("OBLIQUE SHOCK CALCULATOR\n")
M1 = float(input("Enter the upstream Mach number (M1):"))
selection = int(
    input(
        "What is the second known value?\n1: Wave Angle (Beta)\n2: Deflection Angle (theta)\n3: P2/P1\n4: rho2/rho1\n5: T2/T1\n6: p02/p01\n7: p02/p1\n8: M2\n\nEnter a number to continue: "
    )
)
knownValue = float(input("Enter the known value. Please enter any angles in degrees: "))
gamma = float(input("Enter value for gamma: "))

# Rayleigh Pitot Formula
def rayleighPitot(gamma, M):
    p = (
        ((gamma + 1) ** 2 * M**2 / (4 * gamma * M**2 - 2 * (gamma - 1)))
        ** (gamma / (gamma - 1))
        * (1 - gamma + 2 * gamma * M**2)
        / (gamma + 1)
    )
    return p


if selection == 1:
    # Wave angle is known
    beta = knownValue
    Mn1 = M1 * math.sin(beta * math.pi / 180)
elif selection == 2:
    # Deflection angle is known
    theta = knownValue
    # TODO Implement Eq. 4.19 from Anderson
    L = math.sqrt(
        (M1**2 - 1) ** 3
        - 3
        * (1 + (gamma - 1) * M1**2 / 2)
        * (1 + (gamma + 1) * M1**2 / 2)
        * math.tan(theta) ** 2
    )
    X = (
        (M1**2 - 1) ** 3
        - 9
        * (1 + (gamma - 1) * M1**2 / 2)
        * (1 + (gamma + 1) * M1**2 / 2 + (gamma + 1) * M1**4 / 4)
        * math.tan(theta) ** 2
    ) / L**3
    # check weak or strong solution
    A_weak = math.cos((4 * math.pi + math.acos(X)) / 3)
    A_strong = math.cos(math.acos(X) / 3)
    beta_weak = math.atan(M1**2 - 1 + 2 * L * A_weak) / (
        3 * (1 + (gamma - 1) * M1**2 / 2) * math.tan(theta)
    )
    beta_strong = math.atan(M1**2 - 1 + 2 * L * A_weak) / (
        3 * (1 + (gamma - 1) * M1**2 / 2) * math.tan(theta)
    )
elif selection == 3:
    # Static pressure ratio known
    P2_P1 = knownValue
    Mn1 = math.sqrt((P2_P1 - 1) * (gamma + 1) / (2 * gamma) + 1)
    beta = math.asin(Mn1 / M1)
elif selection == 4:
    # density ratio known
    rho2_rho1 = knownValue
    Mn1 = math.sqrt((2 * rho2_rho1) / ((gamma + 1) - rho2_rho1 * (gamma - 1)))
    beta = math.asin(Mn1 / M1)
elif selection == 5:
    # temperature ratio known
    T2_T1 = knownValue
    a = 2 * (gamma**2 - gamma)
    b = 4 * gamma - (gamma - 1) ** 2 - T2_T1 * (gamma + 1) ** 2
    c = -2 * (gamma - 1)
    Mn1 = math.sqrt((-b + math.sqrt(b**2 - 4 * a * c)) / (2 * a))
    beta = math.asin(Mn1 / M1)
elif selection == 6:
    # total pressure ratio known
    p02_p01 = knownValue
    M_guess = 2
    Mn1 = 1
    # Newton-Raphson Method
    while abs(M_guess - Mn1) > 1e-7:
        Mn1 = M_guess
        rho2_rho1 = (gamma + 1) * Mn1**2 / (2 + (gamma - 1) * Mn1**2)
        P1_P2 = (gamma + 1) / (2 * gamma * Mn1**2 - gamma + 1)
        f = (rho2_rho1 ** (gamma / (gamma - 1))) * (
            P1_P2 ** (1 / (gamma - 1))
        ) - p02_p01
        dP1P2 = -4 * gamma * Mn1 * (gamma + 1) / (2 * gamma * Mn1**2 - gamma + 1) ** 2
        drho2rho1 = 4 * (gamma + 1) * Mn1 / ((gamma - 1) * Mn1**2 + 2) ** 2
        fprime = (gamma / (gamma - 1)) * (
            rho2_rho1 ** (1 / (gamma - 1))
        ) * drho2rho1 * (P1_P2 ** (1 / (gamma - 1))) + (
            rho2_rho1 ** (gamma / (gamma - 1))
        ) * (
            1 / (gamma - 1)
        ) * (
            P1_P2 ** ((2 - gamma) / (gamma - 1))
        ) * dP1P2
        M_guess = Mn1 - (f / fprime)
    Mn1 = M_guess
elif selection == 7:
    # Pitot tube ratio known
    p02_p1 = knownValue
    M_guess = 2
    Mn1 = 1
    while abs(M_guess - Mn1) > 1e-7:
        Mn1 = M_guess
        f = rayleighPitot(gamma, Mn1) - p02_p1
        # Approximate derivative using 4th order central difference
        h = 1e-3  # This should be the step size with the lowest error
        fprime = (
            rayleighPitot(gamma, (Mn1 - 2 * h))
            - 8 * rayleighPitot(gamma, (Mn1 - h))
            + 8 * rayleighPitot(gamma, (Mn1 + h))
            - rayleighPitot(gamma, (Mn1 + 2 * h))
        ) / (12 * h)
        M_guess = Mn1 - (f / fprime)
    Mn1 = M_guess
elif selection == 8:
    # M2 known
    # TODO Implement NACA 1135 Equation 132
    M2 = knownValue
else:
    print("error")
    sys.exit()


# Normal shock components
# TODO

# Final values
print("M1:", M1)
print("M2:", M2)
print("Wave angle: ", beta, "degrees")
print("Deflection angle: ", theta, "degrees")
print("T2/T1 is ", T2_T1)
print("P2/P1 is ", P2_P1)
print("rho2/rho1 is", rho2_rho1)
print("p02/p01 is", p02_p01)
print("p02/p1 is", p02_p1)
