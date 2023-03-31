import numpy as np
from scipy.optimize import bisect
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d


# function to solve
def fun_t(x, c1, c2) :
    return x ** 3 + c1 * x ** 2 - c2

def fun_p(x, c1, c2, c3, c4) :
    return x ** 2 + c1 * (c2 * x + c3) * np.sqrt(1 - x ** 2) + c4*x - 2


# height (H), length (L), and width (D) of soil block
H = 0.04
L = 0.15
W_g = 2.0

# angle of ground surface with respect to horizontal surface
psi_deg = 30
psi_rad = psi_deg * np.pi / 180
sin_psi = np.sin(psi_rad)
cos_psi = np.cos(psi_rad)

# undrained strength of clay (in kPa)
c_u_list = []
Fx_b = []
Fx_t = []
Fx_p = []
Fx = []

# set the angle at the tip of the plate
theta_tip = np.arctan(L/H)
sin_theta_tip = np.sin(theta_tip)

for i in range(1, 601) :
    # increase c_u
    c_u = i * 0.05
    c_u_list.append(c_u)

    """
    Block failure mode
    """

    # calculate the soil thrust for the block failure
    Fx_b.append(c_u * L - W_g * sin_psi)

    """
    Triangular wedge failure mode
    """

    # find angle of the failure surface with respect to the vertica line: triangle
    theta_rad = np.arctan(np.sqrt(1 + W_g * cos_psi / c_u / H))
    sin_theta = np.sin(theta_rad)
    cos_theta = np.cos(theta_rad)

    # calculate the soil thrust for the triangular wedge failure
    Fx_t.append(W_g * (cos_psi * cos_theta / sin_theta - sin_psi)
                + c_u * H / sin_theta / cos_theta )

    # Print min soil thrust as a solution
    Fx.append(np.nanmin([Fx_b[i-1], Fx_t[i-1]]))


# plot soil thrust values
plt.plot(c_u_list, Fx, color='silver', linewidth=8, label=r'$F_x$')
plt.plot(c_u_list, Fx_t, '--', color='grey', label=r'$F_{xt}$')
plt.plot(c_u_list, Fx_b, '-', color='grey', label=r'$F_{xb}$')

# set axes title
plt.xlabel("Undrained Shear Strength (kPa)")
plt.ylabel("Soil Thrust (kN/m)")

# set min values of x and y axes
plt.ylim(bottom=-1.0)
plt.xlim(left=0.)
plt.ylim(top=5.0)
plt.xlim(right=30.)

plt.tick_params(direction='in', bottom=True, top=True, left=True, right=True)
plt.grid(axis='both', color='lightgray', ls='-', lw=0.5)

plt.text(20, 0.02,
         r'$H$ = ' + str(round(H, 4)) + " m \n" +
         r'$L$ = ' + str(round(L, 4)) + " m \n" +
         r'$W_g$ = ' + str(round(W_g, 4)) + " kN/m",
         ha="left", va="bottom", size=12,
         )

# set legend
plt.legend(bbox_to_anchor=(0., 1.02, 1., .02),
           loc=4,
           ncol=4,
           mode="expand",
           borderaxespad=0.,
           frameon=False,
           fontsize=12)

plt.savefig("test.svg", format="svg")
plt.show()