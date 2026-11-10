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
Mn1 = 0
beta = 0
theta = 0
M2 = 0
P2_P1 = 0
rho2_rho1 = 0
T2_T1 = 0
p02_p01 = 0
p02_p1 = 0

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
    if theta == 0:
        case = int(input("Weak shock (1) or Strong shock (2)? "))
        if case == 2:
            beta = 90.0
        else:
            beta = math.asin(1 / M1) * (180 / math.pi)
    else:
        theta_rad = theta * math.pi / 180
        tan_th = math.tan(theta_rad)
        L_sq = (
            (M1**2 - 1) ** 2
            - 3
            * (1 + (gamma - 1) * M1**2 / 2)
            * (1 + (gamma + 1) * M1**2 / 2)
            * tan_th**2
        )
        if L_sq < 0:
            print("Deflection angle too large; shock is detached.")
            sys.exit()
        L = math.sqrt(L_sq)
        X = (
            (M1**2 - 1) ** 3
            - 9
            * (1 + (gamma - 1) * M1**2 / 2)
            * (1 + (gamma - 1) * M1**2 / 2 + (gamma + 1) * M1**4 / 4)
            * tan_th**2
        ) / L**3
        if X > 1:
            X = 1.0
        elif X < -1:
            X = -1.0

        # check weak or strong solution
        A_weak = math.cos((4 * math.pi + math.acos(X)) / 3)
        A_strong = math.cos(math.acos(X) / 3)
        beta_weak = (
            math.atan(
                (M1**2 - 1 + 2 * L * A_weak)
                / (3 * (1 + (gamma - 1) * M1**2 / 2) * tan_th)
            )
            * (180 / math.pi)
        )
        beta_strong = (
            math.atan(
                (M1**2 - 1 + 2 * L * A_strong)
                / (3 * (1 + (gamma - 1) * M1**2 / 2) * tan_th)
            )
            * (180 / math.pi)
        )

        case = int(input("Weak shock (1) or Strong shock (2)? "))
        if case == 2:
            beta = beta_strong
        else:
            beta = beta_weak

    Mn1 = M1 * math.sin(beta * math.pi / 180)
elif selection == 3:
    # Static pressure ratio known
    P2_P1 = knownValue
    Mn1 = math.sqrt((P2_P1 - 1) * (gamma + 1) / (2 * gamma) + 1)
    beta = math.asin(Mn1 / M1) * (180 / math.pi)
elif selection == 4:
    # density ratio known
    rho2_rho1 = knownValue
    Mn1 = math.sqrt((2 * rho2_rho1) / ((gamma + 1) - rho2_rho1 * (gamma - 1)))
    beta = math.asin(Mn1 / M1) * (180 / math.pi)
elif selection == 5:
    # temperature ratio known
    T2_T1 = knownValue
    a = 2 * (gamma**2 - gamma)
    b = 4 * gamma - (gamma - 1) ** 2 - T2_T1 * (gamma + 1) ** 2
    c = -2 * (gamma - 1)
    Mn1 = math.sqrt((-b + math.sqrt(b**2 - 4 * a * c)) / (2 * a))
    beta = math.asin(Mn1 / M1) * (180 / math.pi)
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
    beta = math.asin(Mn1 / M1) * (180 / math.pi)
elif selection == 7:
    # Stagnation to upstream static pressure ratio known (p02/p1)
    p02_p1 = knownValue
    T0_T = 1 + 0.5 * (gamma - 1) * M1**2
    P0_P = T0_T ** (gamma / (gamma - 1))
    p02_p01 = p02_p1 / P0_P
    M_guess = 2
    Mn1 = 1
    # Newton-Raphson Method
    while abs(M_guess - Mn1) > 1e-7:
        Mn1 = M_guess
        P1_P2 = 1 / (1 + ((2 * gamma) / (gamma + 1)) * (Mn1**2 - 1))
        rho2_rho1 = (gamma + 1) * Mn1**2 / (2 + (gamma - 1) * Mn1**2)
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
    beta = math.asin(Mn1 / M1) * (180 / math.pi)
elif selection == 8:
    # M2 known
    # NACA 1135 Equation 132 solved as a quadratic in Mn1^2
    M2 = knownValue
    a = 4 * gamma + 2 * gamma * (gamma - 1) * M2**2
    b = (
        (4 * gamma - (gamma - 1) ** 2) * M2**2
        - (gamma + 1) ** 2 * M1**2
        - 4 * (gamma - 1)
    )
    c = -4 - 2 * (gamma - 1) * M2**2
    Mn1 = math.sqrt((-b + math.sqrt(b**2 - 4 * a * c)) / (2 * a))
    sin_beta = Mn1 / M1
    if sin_beta > 1:
        sin_beta = 1.0
    beta = math.asin(sin_beta) * (180 / math.pi)
else:
    print("error")
    sys.exit()


# Normal shock components
cp = gamma * 287 / (gamma - 1)
T0_T = 1 + 0.5 * (gamma - 1) * M1**2
P0_P = T0_T ** (gamma / (gamma - 1))

Mn2 = math.sqrt(
    (1 + (gamma - 1) * (Mn1**2) / 2) / (gamma * Mn1**2 - (gamma - 1) / 2)
)

if selection != 3:
    P2_P1 = 1 + ((2 * gamma) / (gamma + 1)) * (Mn1**2 - 1)
if selection != 4:
    rho2_rho1 = (gamma + 1) * Mn1**2 / (2 + (gamma - 1) * Mn1**2)
if selection != 5:
    T2_T1 = P2_P1 / rho2_rho1

entropyChange = cp * math.log(T2_T1) - 287 * math.log(P2_P1)

if selection != 6:
    p02_p01 = math.exp(-entropyChange / 287)
if selection != 7:
    p02_p1 = p02_p01 * P0_P

if selection != 1 and selection != 2:
    sin_beta = Mn1 / M1
    if sin_beta > 1:
        sin_beta = 1.0
    beta = math.asin(sin_beta) * (180 / math.pi)

if selection != 2:
    if beta == 90:
        theta = 0.0
    else:
        beta_rad = beta * math.pi / 180
        theta_rad = math.atan(
            2
            * (1 / math.tan(beta_rad))
            * (Mn1**2 - 1)
            / (M1**2 * (gamma + math.cos(2 * beta_rad)) + 2)
        )
        theta = theta_rad * 180 / math.pi

if selection != 8:
    M2 = Mn2 / math.sin((beta - theta) * math.pi / 180)

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
