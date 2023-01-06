import math
import matplotlib
import numpy as np
import scipy
import time
import tqdm

# Initialize globals
mu = 3.986004e14
R_Earth = 6378.137e3
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


def main():
    # orbital params in inertial frame
    r0 = np.array([-2.491984247928895e06, 4.793000892455519e05, -6.824701828788767e06])
    v0 = np.array([6.708662359765611e03, 2.089444378549832e03, -2.302871311065931e03])
    f_inc = np.linspace(0, 2 * math.pi, num=360)

    times = np.linspace(0, 365.25 * 24, num=100)
    times_sec = times * 3600
    rI = np.zeros((3, np.size(f_inc)))
    # r_OF = []

    # Get initial orbit conditions
    orbit_params = inertial_to_orbital(r0, v0)
    a_0 = orbit_params[0]
    e_0 = orbit_params[1]
    i_0 = orbit_params[2]
    LAN_0 = orbit_params[3]
    omega_0 = orbit_params[4]
    f_0 = orbit_params[5]
    p_0 = a_0 * (1 - e_0**2)
    r1 = np.divide(p_0, (1 + e_0 * np.cos(f_inc)))
    ro = np.squeeze(
        np.array(
            [
                [np.multiply(r1, np.cos(f_inc))],
                [np.multiply(r1, np.sin(f_inc))],
                [np.zeros((np.size(f_inc)))],
            ]
        )
    )
    ROI1 = ROI(omega_0, i_0, LAN_0)
    RIO = ROI1.T

    for i in range(1, np.size(f_inc)):
        rI[:, i - 1] = (RIO @ ro[:, i - 1]).T

    # Figuring out avg orbit propagation from initial values
    T_0 = 2 * math.pi * math.sqrt(a_0**3 / mu)
    T_hr = T_0 / 3600
    K2_0 = (
        (-1.5)
        * (math.sqrt(mu) * J2 * R_Earth**2)
        / ((1 - e_0**2) ** 2 * a_0 ** (7 / 2))
    )
    omega_dot_avg = K2_0 * (2.5 * math.sin(i_0) ** 2 - 2)
    LAN_dot_avg = K2_0 * math.cos(i_0)
    LAN_dot = LAN_dot_avg * T_0

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
    rx, ry, rz, vx, vy, vz = sol.y
    print(rx.shape)


if __name__ == "__main__":
    main()
