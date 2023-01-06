import math
import sys

print(
    "SHOCK TUBE CALCULATOR\nThis will calculate the pressure ratio and Mach number of shock (including reflected shock) based on diaphragm or shock pressure ratio"
)

selection = int(
    input(
        "What is the known value?\n1: P2/P1\n2: P4/P1\n\nEnter a number to continue: "
    )
)
knownValue = float(input("Enter the known value: "))
gammaDriver = float(input("Enter gamma for driver gas: "))
Tdriver = float(input("Enter temperature of driver gas in Kelvin: "))
Rdriver = float(input("Enter gas constant for driver gas: "))
gammaDriven = float(input("Enter gamma for driven gas: "))
Tdriven = float(input("Enter temperature of driven gas in Kelvin: "))
Rdriven = float(input("Enter gas constant for driven gas: "))


def pressureRatio(a1, a4, gamma1, gamma4, P2_P1):
    p = P2_P1 * (
        1
        - ((gamma4 - 1) * (a1 / a4) * (P2_P1 - 1))
        / math.sqrt(2 * gamma1 * (2 * gamma1 + (gamma1 + 1) * (P2_P1 - 1)))
    ) ** ((-2 * gamma4) / (gamma4 - 1))
    return p


a1 = math.sqrt(gammaDriven * Rdriven * Tdriven)
a4 = math.sqrt(gammaDriver * Rdriver * Tdriver)

if selection == 1:
    P2_P1 = knownValue
    P4_P1 = pressureRatio(a1, a4, gammaDriven, gammaDriver, P2_P1)
elif selection == 2:
    P4_P1 = knownValue
    P2_P1 = 0
    guess = 2
    while abs(guess - P2_P1) > 1e-8:
        P2_P1 = guess
        f = pressureRatio(a1, a4, gammaDriven, gammaDriver, P2_P1) - P4_P1
        # Approximate derivative using 4th order central difference
        h = 1e-3  # This should be the step size with the lowest error
        fprime = (
            pressureRatio(a1, a4, gammaDriven, gammaDriver, (P2_P1 - 2 * h))
            - 8 * pressureRatio(a1, a4, gammaDriven, gammaDriver, (P2_P1 - h))
            + 8 * pressureRatio(a1, a4, gammaDriven, gammaDriver, (P2_P1 + h))
            - pressureRatio(a1, a4, gammaDriven, gammaDriver, (P2_P1 + 2 * h))
        ) / (12 * h)
        guess = P2_P1 - (f / fprime)
else:
    print("error")
    sys.exit()

Ms = math.sqrt((P2_P1 - 1) * (gammaDriven + 1) / (2 * gammaDriven) + 1)

# Reflected Shock
a = (Ms / (Ms**2 - 1)) * math.sqrt(
    1
    + (
        (2 * (gammaDriven - 1) / (gammaDriven + 1) ** 2)
        * (Ms**2 - 1)
        * (gammaDriven + (1 / Ms**2))
    )
)
b = -1
c = -a
Mr = (-b + math.sqrt(b**2 - 4 * a * c)) / (2 * a)
P5_P2 = 1 + ((2 * gammaDriven) / (gammaDriven + 1)) * (Mr**2 - 1)

# Output values
print("P4/P1 = ", P4_P1)
print("P2_P1 = ", P2_P1)
print("Ms = ", Ms)
print("Mr = ", Mr)
print("P5/P2 = ", P5_P2)
