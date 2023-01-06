import math
import matplotlib
import numpy

def inertial_to_orbital(mu,r,v):
    h = numpy.cross(r,v)
    i_h = h/numpy.linalg.norm(h)
    c = numpy.cross(v,h) - (mu/numpy.linalg.norm(r))*r
    e = numpy.linalg.norm(c)/mu
    i_e = c/numpy.linalg.norm(c)
    i_y = numpy.cross(i_h,i_e)
    energy = (numpy.linalg.norm(v)^2)/2 - mu/numpy.linalg.norm(r)
    a = -mu/(2*energy)
    i = math.acos(i_h(3))
    cos_LAN = -i_h(2)/math.sin(i)
    sin_LAN = i_h(1)/math.sin(i)
    LAN = numpy.linalg.norm(sin_LAN, cos_LAN)
    cos_omega = i_y(3)/math.sin(i)
    sin_omega = i_e(3)/math.sin(i)
    omega = numpy.arctan2(sin_omega, cos_omega)
    r_o = Roi*r
    cos_f = r_o(1)/numpy.linalg.norm(r_o)
    sin_f = r_o(2)/numpy.linalg.norm(r_o)
    f = numpy.arctan2(sin_f, cos_f)
    return params
