"""
Main developer: Sang Inn Woo @ Incheon National University
Pur

"""

# import depending libraries
import numpy as np
from scipy.optimize import bisect

# declare a function to solve
def fun(x, c1, c2) :
    return x ** 3 + c1 * x ** 2 - c2


# height (H), length (L), and width (D) of soil block
H = 0.0357
L = 0.1239
D = 0.1239
W_g = 0.3454
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

print("Fxb = " + str(round(Fx_b,4)) + " kN")
print("Fxt = " + str(round(Fx_t,4)) + " kN")
print("Fx = " + str(round(Fx,4)) + " kN")