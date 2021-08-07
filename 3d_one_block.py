import numpy
import numpy as np
from scipy.optimize import bisect
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def fun(x, c1, c2):
    return x ** 3 + c1 * x ** 2 - c2


# height (H), length (L), and width (D) of soil block
H = 0.03
L = 0.1
D = 0.3

# vertical force (kN) transmitted to the track system
sig_g = 5.6
W_g = sig_g * L * D

# undrained strength of clay (in kPa)
c_u_list = []
Fx_b = []
Fx_t = []

for i in range(1, 401):

    # increase c_u
    c_u = i*0.05
    c_u_list.append(c_u)

    # calculate the soil thrust for the block failure
    Fx_b.append(c_u * L * (2 * H + D))

    # find angle of the failure surface with respect to the vertica line
    sin_theta = bisect(fun, 0, 1, args=(W_g / c_u / H ** 2 + 2 * D / H, W_g / c_u / H ** 2 + D / H))
    theta_rad = numpy.arcsin(sin_theta)
    theta_deg = theta_rad * 180 / numpy.pi

    # calculate the soil thrust for the triangular wedge failure
    Fx_t.append(W_g * np.cos(theta_rad) / np.sin(theta_rad)
                + c_u * H * (D / np.sin(theta_rad) / np.cos(theta_rad) + H / np.cos(theta_rad)))

Fx = np.minimum(Fx_t, Fx_b)

plt.plot(c_u_list, Fx, color='silver', linewidth=8)
plt.plot(c_u_list, Fx_t, '-', color='black')
plt.plot(c_u_list, Fx_b, '--', color='black')
plt.show()
