import math
import sys

print("NORMAL SHOCK FLOW CALCULATOR\n")
selection = int(
    input(
        "What is the known value?\n1: M1\n2: M2\n3: P2/P1\n4: rho2/rho1\n5: T2/T1\n6: p02/p01\n7: p02/p1\n\nEnter a number to continue: "
    )
)
knownValue = float(input("Enter the known value: "))
gamma = float(input("Enter value for gamma: "))

M1 = 0
if selection == 1:
    # M1 known
    M1 = knownValue
elif selection == 2:
    # M2 known
    M2 = knownValue
    if knownValue <= 1:
        M1 = math.sqrt(
            (1 + (gamma - 1) * (M2**2) / 2) / (1 * gamma * M2**2 - (gamma - 1) / 2)
        )
    else:
        print("M2 must be subsonic")
        sys.exit()
elif selection == 3:
    # P2/P1 known
    p2_p1 = knownValue
    M1 = math.sqrt((p2_p1 - 1) * (gamma + 1) / (2 * gamma) + 1)
elif selection == 4:
    # rho2/rho1 known
    rho2_rho1 = knownValue
    M1 = math.sqrt((2 * rho2_rho1) / ((gamma + 1) - rho2_rho1 * (gamma - 1)))
elif selection == 5:
    # T2/T1 is known
    # Equation is a quadratic with repeated roots
    T2_T1 = knownValue
    a = 2 * (gamma**2 - gamma)
    b = 4 * gamma - (gamma - 1) ** 2 - T2_T1 * (gamma + 1) ** 2
    c = -2 * (gamma - 1)
    M1 = math.sqrt((-b + math.sqrt(b**2 - 4 * a * c)) / (2 * a))
elif selection == 6:
    # TODO: Make this math work
    print("This requires an iterative process and has not yet been programmed")
    sys.exit()
elif selection == 7:
    # TODO: Make this math work
    print("This requires an iterative process and has not yet been programmed")
    sys.exit()
else:
    print("error")
    sys.exit()
5
# Calculate ratios from given/calculated M1
cp = gamma * 287 / (gamma - 1)
T0_T = 1 + 0.5 * (gamma - 1) * M1**2
P0_P = T0_T ** (gamma / (gamma - 1))
if selection != 2:
    M2 = math.sqrt(
        (1 + (gamma - 1) * (M1**2) / 2) / (1 * gamma * M1**2 - (gamma - 1) / 2)
    )
if selection != 3:
    p2_p1 = 1 + ((2 * gamma) / (gamma + 1)) * (M1**2 - 1)
if selection != 4:
    rho2_rho1 = (gamma + 1) * M1**2 / (2 + (gamma - 1) * M1**2)
if selection != 5:
    T2_T1 = p2_p1 / rho2_rho1
entropyChange = cp * math.log(T2_T1) - 287 * math.log(p2_p1)
p02_p01 = math.exp(-entropyChange / 287)
p02_p1 = p02_p01 * P0_P
print("M1:", M1)
print("M2:", M2)
print("T2/T1 is ", T2_T1)
print("P2/P1 is ", p2_p1)
print("rho2/rho1 is", rho2_rho1)
print("p02/p01 is", p02_p01)
print("p02/p1 is", p02_p1)
