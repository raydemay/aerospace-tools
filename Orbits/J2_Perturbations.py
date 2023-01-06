import math
import matplotlib
import numpy


def inertial_to_orbital(mu, r, v):
    # Calculate orbital parameters from given inputs

    h = numpy.cross(r, v)
    i_h = h / numpy.linalg.norm(h)
    c = numpy.cross(v, h) - (mu / numpy.linalg.norm(r)) * r
    e = numpy.linalg.norm(c) / mu
    i_e = c / numpy.linalg.norm(c)
    i_y = numpy.cross(i_h, i_e)
    energy = (numpy.linalg.norm(v) ** 2) / 2 - mu / numpy.linalg.norm(r)
    a = -mu / (2 * energy)
    i = math.acos(i_h[2])
    cos_LAN = -i_h[1] / math.sin(i)
    sin_LAN = i_h[0] / math.sin(i)
    LAN = numpy.arctan2(sin_LAN, cos_LAN)
    cos_omega = i_y[2] / math.sin(i)
    sin_omega = i_e[2] / math.sin(i)
    omega = numpy.arctan2(sin_omega, cos_omega)
    Roi = [numpy.transpose(i_e), numpy.transpose(i_y), numpy.transpose(i_h)]
    r_o = numpy.matmul(Roi, r)
    cos_f = r_o[0] / numpy.linalg.norm(r_o)
    sin_f = r_o[1] / numpy.linalg.norm(r_o)
    f = numpy.arctan2(sin_f, cos_f)
    params = [a, e, i, LAN, omega, f]
    return params


def main():
    # Initialize
    mu = 3.986004e14
    R_Earth = 6378.137e3
    # orbital params in inertial frame
    r0 = numpy.array(
        [-2.491984247928895e06, 4.793000892455519e05, -6.824701828788767e06]
    )
    v0 = numpy.array(
        [6.708662359765611e03, 2.089444378549832e03, -2.302871311065931e03]
    )
    f_inc = numpy.linspace(0, 2 * math.pi, num=360)
    J2 = 1.08263e-3  # 2nd zonal harmonic
    # times = [0, 1, 2, 4, 6, 10, 15, 25, 50, 75, 100, 150]
    times = numpy.linspace(0, 365.25 * 24)
    times_sec = times * 3600
    rI = numpy.zeros((3, numpy.size(f_inc)))
    # r_OF = []

    # Get initial orbit conditions
    orbit_params = inertial_to_orbital(mu, r0, v0)
    print(orbit_params)


if __name__ == "__main__":
    main()
