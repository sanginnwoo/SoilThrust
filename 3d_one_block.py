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

# declare a function to solve in the triangular wedge failure mode
def fun(x, c1, c2) :
    return x ** 3 + c1 * x ** 2 - c2

# set height (H), length (L), and width (D) of soil block
H = 0.0357
L = 0.1239
D = 0.1239

# set vertical load acting the soil block
W_g = 0.3454

# set the undrained shear strength of the clay
c_u = 7.7

# calculate the soil thrust for the block failure
Fx_b = c_u * L * (2 * H + D)

# find angle of the failure surface with respect to the vertica line
sin_theta = bisect(fun, 0, 1, args=(W_g / c_u / H ** 2 + 2 * D / H, W_g / c_u / H ** 2 + D / H))
theta_rad = numpy.arcsin(sin_theta)
theta_deg = theta_rad * 180 / numpy.pi

# calculate the soil thrust for the triangular wedge failure
Fx_t = W_g * np.cos(theta_rad) / np.sin(theta_rad) \
       + c_u * H * (D / np.sin(theta_rad) / np.cos(theta_rad) + H / np.cos(theta_rad))

# get soil thrust values with respect to undrained shear strength
Fx = np.minimum(Fx_t, Fx_b)

# print soil thrust values
print("Fxb = " + str(round(Fx_b,4)) + " kN")
print("Fxt = " + str(round(Fx_t,4)) + " kN")
print("Fx = " + str(round(Fx,4)) + " kN")