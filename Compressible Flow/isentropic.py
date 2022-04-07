import sys
import math

print("ISENTROPIC FLOW CALCULATOR\n")
selection = int(
    input(
        "What is the known value?\n1: Mach number\n2: T0/T\n3: P0/P\n4: rho0/rho\n5: A/A*\n\nEnter a number to continue: "
    )
)
knownValue = float(input("Enter the known value: "))
gamma = float(input("Enter value for gamma: "))

M = 0
# Determine Mach number from inputs
if selection == 1:
    # Mach number is known
    M = knownValue
elif selection == 2:
    # T0/T is known
    M = math.sqrt((knownValue - 1) * 2 / (gamma - 1))
elif selection == 3:
    # P0/P is known
    M = math.sqrt(((knownValue ** ((gamma - 1) / gamma)) - 1) * 2 / (gamma - 1))
elif selection == 4:
    # rho0/rho is known
    M = math.sqrt(((knownValue ** (gamma - 1)) - 1) * 2 / (gamma - 1))
elif selection == 5:
    # A/A* known
    case = int(input("Subsonic (1) or Supersonic (2)?"))
    if case == 1:
        A_AStar = knownValue
        M_guess = 0.001
    elif case == 2:
        A_AStar = knownValue
        M_guess = 2
    R = A_AStar**2
    X_new = M_guess**2
    X = 0
    E = (gamma + 1) / (gamma - 1)
    P = (E - 1) / E
    Q = 1 / E
    while abs(X_new - X) > 0.0000001:
        # Newton-Raphson Method
        X = X_new
        f = (P + Q * X) ** (1 / Q) - R * X
        fprime = (P + Q * X) ** ((1 / Q) - 1) - R
        X_new = X - (f / fprime)
    M = math.sqrt(X)
else:
    print("error")
    sys.exit()

# Calculate ratios and Mach angle from given/calculated Mach Number
if selection != 2:
    T0_T = 1 + 0.5 * (gamma - 1) * M**2
if selection != 3:
    P0_P = T0_T ** (gamma / (gamma - 1))
if selection != 4:
    R0_R = T0_T ** (1 / (gamma - 1))
if selection != 5:
    A_AStar = math.sqrt(
        (1 / M**2) * ((2 / (gamma + 1)) * T0_T) ** ((gamma + 1) / (gamma - 1))
    )
if M > 1:
    Mach_angle = math.asin(1 / M) * (180 / math.pi)
    PrandtlMeyerAngle = (
        math.sqrt((gamma + 1) / (gamma - 1))
        * math.atan((math.sqrt(((gamma - 1) / (gamma + 1)) * (M**2 - 1))))
        - math.atan(math.sqrt(M**2 - 1))
    ) * (180 / math.pi)
else:
    Mach_angle = 90
    PrandtlMeyerAngle = 0
print("Mach number is ", "{:.6f}".format(M))
print("T0/T is ", "{:.6f}".format(T0_T))
print("P0/P is ", "{:.6f}".format(P0_P))
print("rho0/rho is", "{:.6f}".format(R0_R))
print("A/A*", "{:.6f}".format(A_AStar))
print("Mach angle is", "{:.6f}".format(Mach_angle), "degrees")
print("Prandtl-Meyer angle is ", "{:.6f}".format(PrandtlMeyerAngle), "degrees")
