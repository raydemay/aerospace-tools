import math
import matplotlib.pyplot as plt
import numpy as np
import scipy

# Initialize globals
mu = 3.986004e14  # standard gravitational parapeter of Earth
R_Earth = 6378.137e3  # Earth's radius
J2 = 1.08263e-3  # 2nd zonal harmonic


def inertial_to_orbital(r, v):
    # Calculate orbital parameters from given inputs

    h = np.cross(r, v)
    i_h = h / np.linalg.norm(h)
    c = np.cross(v, h) - (mu / np.linalg.norm(r)) * r
    e = np.linalg.norm(c) / mu
    i_e = c / np.linalg.norm(c)
    i_y = np.cross(i_h, i_e)
    energy = (np.linalg.norm(v) ** 2) / 2 - mu / np.linalg.norm(r)
    a = -mu / (2 * energy)
    i = math.acos(i_h[2])
    cos_LAN = -i_h[1] / math.sin(i)
    sin_LAN = i_h[0] / math.sin(i)
    LAN = np.arctan2(sin_LAN, cos_LAN)
    cos_omega = i_y[2] / math.sin(i)
    sin_omega = i_e[2] / math.sin(i)
    omega = np.arctan2(sin_omega, cos_omega)
    Roi = [np.transpose(i_e), np.transpose(i_y), np.transpose(i_h)]
    r_o = np.matmul(Roi, r)
    cos_f = r_o[0] / np.linalg.norm(r_o)
    sin_f = r_o[1] / np.linalg.norm(r_o)
    f = np.arctan2(sin_f, cos_f)
    params = [a, e, i, LAN, omega, f]
    return params


def orbit_eqn_J2_Earth(t, input):
    # Two-body orbit equation with J2 perturbation term
    # See Eqn(10.30) in Orbital Mechanics for Engineering Students 4th Edition by Curtis
    rx, ry, rz, vx, vy, vz = input
    r = np.linalg.norm(np.array([rx, ry, rz]))
    ax = -mu * rx / r**3 + (
        (1.5)
        * (J2 * mu * R_Earth**2 / r**4)
        * (rx / r)
        * (5 * (rz**2 / r**2) - 1)
    )
    ay = -mu * ry / r**3 + (
        (1.5)
        * (J2 * mu * R_Earth**2 / r**4)
        * (ry / r)
        * (5 * (rz**2 / r**2) - 1)
    )
    az = -mu * rz / r**3 + (
        (1.5)
        * (J2 * mu * R_Earth**2 / r**4)
        * (rz / r)
        * (5 * (rz**2 / r**2) - 3)
    )
    derivs = [vx, vy, vz, ax, ay, az]
    return derivs


def orbitalframe_rv(sma, ecc, inc, LAN, omega, f, steps):
    f_inc = np.linspace(0, 2 * math.pi, num=steps)  # true anomaly from 0 to 2pi
    r_inertialframe = np.zeros((3, np.size(f_inc)))
    p = sma * (1 - ecc**2)
    r = p / (1 + ecc * math.cos(f))
    B = np.transpose(ROI(omega, inc, LAN))
    r_orbitalframe = np.squeeze(
        np.array(
            [
                [np.multiply(r, np.cos(f_inc))],
                [np.multiply(r, np.sin(f_inc))],
                [np.zeros((np.size(f_inc)))],
            ]
        )
    )
    for i in range(0, steps):
        r_inertialframe[:, i] = (B @ r_orbitalframe[:, i]).T
    return r_inertialframe


def ROI(om, i, OM):
    ie1 = math.cos(om) * math.cos(OM) - math.sin(om) * math.cos(i) * math.sin(OM)
    ie2 = math.cos(om) * math.sin(OM) + math.sin(om) * math.cos(i) * math.cos(OM)
    ie3 = math.sin(om) * math.sin(i)
    iy1 = -(math.sin(om) * math.cos(OM) + math.cos(om) * math.cos(i) * math.sin(OM))
    iy2 = -math.sin(om) * math.sin(OM) + math.cos(om) * math.cos(i) * math.cos(OM)
    iy3 = math.cos(om) * math.sin(i)
    ih1 = math.sin(i) * math.sin(OM)
    ih2 = -math.sin(i) * math.cos(OM)
    ih3 = math.cos(i)
    rotations = np.array([[ie1, ie2, ie3], [iy1, iy2, iy3], [ih1, ih2, ih3]])
    return rotations


def plot_Earth_and_Orbits(vector):
    x = vector[:, 0, :]
    y = vector[:, 1, :]
    z = vector[:, 2, :]

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    # Make data
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x_Earth = R_Earth * np.outer(np.cos(u), np.sin(v))
    y_Earth = R_Earth * np.outer(np.sin(u), np.sin(v))
    z_Earth = R_Earth * np.outer(np.ones(np.size(u)), np.cos(v))

    # Plot the surface
    ax.plot_surface(x_Earth, y_Earth, z_Earth, alpha=0.5)

    # Plot orbits
    for i in range(0, len(x[0])):
        ax.plot(x[i, :], y[i, :], z[i, :], alpha=1)

    # Set an equal aspect ratio
    ax.set_aspect("equal")

    plt.show()


def main():
    # orbital params in inertial frame
    r0 = np.array(
        [-2.491984247928895e06, 4.793000892455519e05, -6.824701828788767e06]
    )  # m
    v0 = np.array(
        [6.708662359765611e03, 2.089444378549832e03, -2.302871311065931e03]
    )  # m/s
    steps = 100
    times = np.linspace(0, 365.25 * 24, num=steps)
    times_sec = times * 3600
    orbital_elements = np.zeros((6, steps))
    g = np.zeros((steps, 3, steps))

    # Get initial orbit conditions
    # a_0, e_0, i_0, LAN_0, omega_0, f_0 = inertial_to_orbital(r0, v0)
    # rI = orbitalframe_rv(a_0, e_0, i_0, LAN_0, omega_0, f_0, steps)

    # Numerically integrate
    init_vals = np.block([r0, v0])
    t_span = (0, times_sec[-1])

    sol = scipy.integrate.solve_ivp(
        orbit_eqn_J2_Earth,
        t_span,
        init_vals,
        method="DOP853",
        t_eval=times_sec,
        rtol=1e-10,
        atol=1e-12,
    )
    # split inertial frame r,v components into r and v arrays
    rx_i, ry_i, rz_i, vx_i, vy_i, vz_i = sol.y
    ri_vals = np.array([rx_i, ry_i, rz_i])
    vi_vals = np.array([vx_i, vy_i, vz_i])

    for j in range(0, steps):
        orbital_elements[:, j] = inertial_to_orbital(ri_vals[:, j], vi_vals[:, j])
    a, e, i, LAN, omega, f = orbital_elements
    for k in range(0, steps):
        g[k, :, :] = orbitalframe_rv(a[k], e[k], i[k], LAN[k], omega[k], f[k], steps)

    plot_Earth_and_Orbits(g)


if __name__ == "__main__":
    main()
