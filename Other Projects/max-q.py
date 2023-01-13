import numpy as np
import matplotlib.pyplot as plt


def density(alt: float) -> float:
    """
    Returns the air density in kg/m^3 based on altitude.
    Equations from
    https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
    :param alt: Altitude in feet
    :return: Density in kg/m^3
    """
    if alt < 11000.0:
        temp = 15.04 - 0.00649 * alt
        p = 101.29 * ((temp + 273.1) / 288.08) ** 5.256
        rho = p / (0.2869 * (temp + 273.1))
    elif 11000.0 <= alt < 25000:
        temp = -56.46
        p = 22.65 * np.exp(1.73 - 0.000157 * alt)
        rho = p / (0.2869 * (temp + 273.1))
    else:
        temp = -131.21 + 0.00299 * alt
        p = 2.488 * ((temp + 273.1) / 216.6) ** -11.388
        rho = p / (0.2869 * (temp + 273.1))
    return rho


def velocity(time: float, acceleration: float) -> float:
    return acceleration * time


def altitude(time: float, acceleration: float) -> float:
    return 0.5 * acceleration * time**2


def main():
    thrust_sl = 7.607e6  # Newtons
    thrust_vac = 934e3  # Newtons
    cd = 0.5  # No idea what a typical rocket cd is
    m_0 = 549054  # starting mass (Falcon 9; kg)
    m_prop_first_stage = 395700  # first stage propellant mass (kg)
    m_prop_second_stage = 92670  # second stage propellant mass (kg)
    burn_time_first_stage = 162  # seconds
    pitch_final = 65 * np.pi / 180  # radians
    pitch_rate = pitch_final / (burn_time_first_stage / 4)
    S = 5  # frontal area (meters)
    m_dot_first_stage = (
        m_prop_first_stage / burn_time_first_stage
    )  # assuming constant thrust and mass flow rate

    a_0 = thrust_sl / (9.81 * m_0)  # initial acceleration

    plt.style.use("bmh")
    q_values = []
    time = np.arange(0.0, burn_time_first_stage + 0.01, 0.01)
    accel = a_0
    m = m_0
    pitch = 90 * np.pi / 180
    for elapsed_time in time:
        alt = altitude(elapsed_time, accel)
        # Dynamic pressure q = 0.5 * density * velocity^2
        q = 0.5 * density(alt) * velocity(elapsed_time, accel) ** 2
        drag = q * S * cd
        F_net = thrust_sl - drag
        q_values.append(q)
        m = m - (m_dot_first_stage * 0.01)
        if pitch > pitch_final:
            pitch = pitch - pitch_rate * elapsed_time
        else:
            pitch = pitch_final
        accel = (F_net / (m)) * np.sin(
            pitch
        )  # only want vertical component of acceleration

    # Plot w/ 1g accel
    plt.plot(time, q_values, "b-")
    max_val = max(q_values)
    ind = q_values.index(max_val)
    # Plot an arrow and text with the max value
    plt.annotate(
        f"{max_val:.2E} Pa @ {time[ind]:.2f} s",
        xy=(time[ind], max_val),
        xytext=(time[ind] + 10, max_val + 2),
        xycoords="data",
        arrowprops=dict(facecolor="black", shrink=0.05),
    )
    # plot the point of Max Q
    plt.plot(time[ind], max_val, "rx")

    # Configure plot and show
    plt.xlabel("Time (s)")
    plt.ylabel("Pressure (Pa)")
    plt.title("Dynamic pressure as a function of time")
    plt.show()


if __name__ == "__main__":
    main()
