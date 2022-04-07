import sys
import math

print("OBLIQUE WAVE FLOW CALCULATOR\n")
selection = int(
    input(
        "What is the known value?\n1: M1\n2: M2\n3: P2/P1\n4: rho2/rho1\n5: T2/T1\n6: p02/p01\n7: Deflection Angle\n8: Wave Angle\n9: Mn1\n10: Mn2\n\nEnter a number to continue: "
    )
)
knownValue = float(input("Enter the known value: "))
gamma = float(input("Enter value for gamma: "))
