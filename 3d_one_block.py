"""
Main developer: Sang Inn Woo @ Incheon National University
The main goal of this code is to calculate the reliable soil thrust under 3d conditions
of a single track system (plate + grouser) for tracked vehicles over saturated clay ground.
Following the limit analysis concept (the upper bound theorem is used here),
this code assumes two geometrically acceptable failure modes (the block and triangular wedge failure modes).
The proposed soil thrust is the least value of soil thrusts from two failre modes.
"""

# import depending libraries
import numpy as np
from scipy.optimize import bisect

# function to solve
def fun_t(x, c1, c2) :
    return x ** 3 + c1 * x ** 2 - c2

def fun_p(x, c1, c2, c3, c4) :
    return x ** 2 + c1 * (c2 * x + c3) * np.sqrt(1 - x ** 2) + c4*x - 2

Fx_b = np.nan
Fx_t = np.nan
Fx_p = np.nan

# set height (H), length (L), and width (D) of soil block
H = 0.0357
L = 0.1239
D = 0.1239

# set vertical load acting the soil block
W_g = 0.

# set the undrained shear strength of the clay
c_u = 7.7

# set the angle at the tip of the plate
theta_tip = np.arctan(L/H)
sin_theta_tip = np.sin(theta_tip)


"""
Block failure mode
"""

# calculate the soil thrust for the block failure
Fx_b = c_u * L * (2 * H + D)


"""
Triangular wedge failure mode
"""

# check the start and end values for the bisection method
bisect_start_value = fun_t(0, W_g / c_u / H ** 2 + 2 * D / H, W_g / c_u / H ** 2 + D / H)
bisect_end_value = fun_t(sin_theta_tip, W_g / c_u / H ** 2 + 2 * D / H, W_g / c_u / H ** 2 + D / H)

# if solution exists
if bisect_start_value * bisect_end_value < 0 :

    # find angle of the failure surface with respect to the vertica line: triangle
    sin_theta_t = bisect(fun_t, 0, sin_theta_tip, args=(W_g / c_u / H ** 2 + 2 * D / H, W_g / c_u / H ** 2 + D / H))
    theta_rad_t = np.arcsin(sin_theta_t)
    theta_deg_t = theta_rad_t * 180 / np.pi

    # calculate the soil thrust for the triangular wedge failure
    Fx_t = W_g * np.cos(theta_rad_t) / np.sin(theta_rad_t) \
           + c_u * H * (D / np.sin(theta_rad_t) / np.cos(theta_rad_t) + H / np.cos(theta_rad_t))


"""
Trapezoidal wedge failure mode
"""

# check the start and end values for the bisection method
bisect_start_value = fun_p(sin_theta_tip, 2 / L, H, D, W_g / c_u / L ** 2)
bisect_end_value = fun_p(1, 2 / L, H, D, W_g / c_u / L ** 2)

# if solution exists
if bisect_start_value * bisect_end_value < 0 :

    # find angle of the failure surface with respect to the vertica line: trapezoidal
    sin_theta_p = bisect(fun_p, sin_theta_tip, 1, args=(2 / L, H, D, W_g / c_u / L ** 2))
    theta_rad_p = np.arcsin(sin_theta_p)
    theta_deg_p = theta_rad_p * 180 / np.pi
    Fx_p = W_g * np.cos(theta_rad_p) / np.sin(theta_rad_p) \
           + c_u * L * (D / np.sin(theta_rad_p) / np.sin(theta_rad_p) \
                        + 2 * H / np.sin(theta_rad_p) \
                        - L * np.cos(theta_rad_p) / np.sin(theta_rad_p) ** 2)

# Print min soil thrust as a solution
Fx = np.nanmin([Fx_b, Fx_t, Fx_p])

# print soil thrust values
print("Fxb = " + str(round(Fx_b,4)) + " kN")
print("Fxt = " + str(round(Fx_t,4)) + " kN")
print("Fxp = " + str(round(Fx_p,4)) + " kN")
print("Fx = " + str(round(Fx,4)) + " kN")