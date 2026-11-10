"""
CONICAL FLOW CALCULATOR for TI-84 Plus CE Python (MicroPython).
Solves Taylor-Maccoll supersonic conical flow over a circular cone at zero angle of attack.
"""
import math
import sys


def solve_cone_from_shock(m_inf, th_s, g):
    """
    Integrates Taylor-Maccoll equation inward from conical shock wave
    angle th_s to the cone surface where polar velocity Vt = 0.
    Returns (th_c, vr_c) or None.
    """
    v_inf = m_inf * math.sqrt((g - 1) / (2 + (g - 1) * m_inf**2))
    vr = v_inf * math.cos(th_s)
    mn1 = m_inf * math.sin(th_s)
    if mn1 < 1.0:
        return None

    ratio = (2.0 + (g - 1) * mn1**2) / ((g + 1) * mn1**2)
    vt = -v_inf * math.sin(th_s) * ratio

    def tm(th, vr, vt):
        gm1_2 = (g - 1.0) / 2.0
        a2 = gm1_2 * (1.0 - vr**2 - vt**2)
        denom = a2 - vt**2
        if abs(denom) < 1e-12:
            return 0.0
        return (vt**2 * vr - a2 * (2.0 * vr + vt / math.tan(th))) / denom

    dtheta = -math.radians(0.25)
    theta = th_s
    while vt < 0 and theta > math.radians(0.2):
        # 4th-order Runge-Kutta
        k1_r, k1_t = vt, tm(theta, vr, vt)
        k2_r, k2_t = vt + 0.5 * dtheta * k1_t, tm(
            theta + 0.5 * dtheta, vr + 0.5 * dtheta * k1_r, vt + 0.5 * dtheta * k1_t
        )
        k3_r, k3_t = vt + 0.5 * dtheta * k2_t, tm(
            theta + 0.5 * dtheta, vr + 0.5 * dtheta * k2_r, vt + 0.5 * dtheta * k2_t
        )
        k4_r, k4_t = vt + dtheta * k3_t, tm(
            theta + dtheta, vr + dtheta * k3_r, vt + dtheta * k3_t
        )

        vr_next = vr + (dtheta / 6.0) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r)
        vt_next = vt + (dtheta / 6.0) * (k1_t + 2 * k2_t + 2 * k3_t + k4_t)
        theta_next = theta + dtheta

        if vt_next >= 0:
            frac = -vt / (vt_next - vt)
            return theta + frac * dtheta, vr + frac * (vr_next - vr)

        vr, vt, theta = vr_next, vt_next, theta_next

    return None


# --- User Interface ---
print("--- CONICAL FLOW ---")
print("1: th_c (Cone angle)")
print("2: th_s (Shock angle)")
opt = int(input("Known [1-2]: "))
val = float(input("Angle (deg): "))
m_in = input("M_inf [2.0]: ")
m_inf = float(m_in) if m_in else 2.0
g_in = input("gamma [1.4]: ")
g = float(g_in) if g_in else 1.4

if m_inf <= 1.0:
    print("Error: M_inf > 1")
    sys.exit()

mu = math.asin(1.0 / m_inf)

if opt == 1:
    th_c_target = math.radians(val)
    if th_c_target <= 0 or th_c_target >= math.radians(60):
        print("Error: Invalid th_c")
        sys.exit()

    # Secant shooting to find shock angle th_s matching th_c
    s0 = mu + 0.25 * th_c_target
    s1 = mu + 0.65 * th_c_target
    found = False
    for _ in range(30):
        r0 = solve_cone_from_shock(m_inf, s0, g)
        r1 = solve_cone_from_shock(m_inf, s1, g)
        if not r0 or not r1:
            break
        f0 = r0[0] - th_c_target
        f1 = r1[0] - th_c_target
        if abs(f1) < 1e-6 or abs(s1 - s0) < 1e-6:
            th_s = s1
            th_c, vr_c = r1
            found = True
            break
        s0, s1 = s1, s1 - f1 * (s1 - s0) / (f1 - f0)

    if not found:
        print("Shock detached or no sol")
        sys.exit()

elif opt == 2:
    th_s = math.radians(val)
    if th_s <= mu or th_s >= math.pi / 2:
        print("Error: th_s > mu (Mach angle)")
        sys.exit()
    res = solve_cone_from_shock(m_inf, th_s, g)
    if not res:
        print("Error: Cannot integrate")
        sys.exit()
    th_c, vr_c = res
else:
    print("error")
    sys.exit()

# Compute flow properties
mn1 = m_inf * math.sin(th_s)
a_c = math.sqrt(0.5 * (g - 1) * (1.0 - vr_c**2))
mc = vr_c / a_c

# Static pressure ratio across shock P2/P1
p2_p1 = 1.0 + (2.0 * g / (g + 1.0)) * (mn1**2 - 1.0)

# Stagnation pressure ratio across shock P02/P01
term1 = ((g + 1.0) * mn1**2) / (2.0 + (g - 1.0) * mn1**2)
term2 = (g + 1.0) / (2.0 * g * mn1**2 - (g - 1.0))
p02_p01 = (term1 ** (g / (g - 1.0))) * (term2 ** (1.0 / (g - 1.0)))

# Surface static pressure ratio Pc/Pinf
pc_pinf = (
    ((1.0 + 0.5 * (g - 1.0) * m_inf**2) / (1.0 + 0.5 * (g - 1.0) * mc**2))
    ** (g / (g - 1.0))
) * p02_p01

# Pressure coefficient Cp
cp = (2.0 / (g * m_inf**2)) * (pc_pinf - 1.0)

# Display results formatted for TI-84 screen (26x10 chars)
print("\nth_c  : %.2f deg" % math.degrees(th_c))
print("th_s  : %.2f deg" % math.degrees(th_s))
print("Mc    : %.4f" % mc)
print("Mn1   : %.4f" % mn1)
print("P2/P1 : %.4f (shock)" % p2_p1)
print("Pc/P1 : %.4f (cone)" % pc_pinf)
print("Cp    : %.4f" % cp)
print("P0c/P0: %.4f" % p02_p01)
