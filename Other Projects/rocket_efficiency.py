import numpy as np
from matplotlib import pyplot as plt


def prop_efficiency(c, u):
    # Propulsive efficiency equation
    eta = 2 * (u / c) / (1 + u**2 / c**2)
    return eta


def main():
    # Velocity of jet (m/s)
    c_1 = 2e3
    c_2 = 3.5e3
    c_3 = 5e3
    u = np.arange(0, 15e3)  # vehicle velocity

    # Get efficiency to plot
    eta_1 = prop_efficiency(c_1, u)
    eta_2 = prop_efficiency(c_2, u)
    eta_3 = prop_efficiency(c_3, u)
    fig, ax = plt.subplots()
    ax.plot(u, eta_1, u, eta_2, u, eta_3)
    ax.set(
        xlabel="Vehicle velocity (m/s)",
        ylabel=r"Efficiency $\eta_{prop}$",
        title="Propulsive Efficiency",
    )
    ax.legend((f"c = {c_1} m/s", f"c = {c_2} m/s", f"c = {c_3} m/s"))
    # fig.savefig("PropulsiveEfficiency.png")
    plt.show()


if __name__ == "__main__":
    main()
