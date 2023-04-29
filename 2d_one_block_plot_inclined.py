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
H = 0.036
L = 0.124
W_g = 1.36

# construct test data
c_u_test = [10.3, 11.7, 14.3, 28.6]
Fx_test = [1.21, 1.54, 1.85, 2.77]

# angle of ground surface with respect to horizontal surface
psi_deg = 0
psi_rad = psi_deg * np.pi / 180
sin_psi = np.sin(psi_rad)
cos_psi = np.cos(psi_rad)

# undrained strength of clay (in kPa)
c_u_list = []
Fx_b = []
Fx_t = []
Fx_p = []
Fx = []
theta = []

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

    theta.append(np.rad2deg(theta_rad))

    # calculate the soil thrust for the triangular wedge failure
    Fx_t.append(W_g * (cos_psi * cos_theta / sin_theta - sin_psi)
                + c_u * H / sin_theta / cos_theta )

    # Print min soil thrust as a solution
    Fx.append(np.nanmin([Fx_b[i-1], Fx_t[i-1]]))

# set the size of plot
fig, axs = plt.subplots(2, 1,
                        gridspec_kw={'height_ratios': [1, 2]})

# plot theta with respec to c_u
axs[0].plot(c_u_list, theta, color='black', label=r'$theta$')
axs[0].axhline(np.rad2deg(theta_tip), color='grey', linestyle='--')

# set min values of x and y axes
axs[0].set_ylim(45, 90)
axs[0].set_xlim(0, 30)

# set axes title
axs[0].set_ylabel("theta (deg)")

# Set grid
axs[0].set_yticks(np.arange(45, 91, 15))
axs[0].tick_params(direction='in', bottom=True, top=True, left=True, right=True)
axs[0].grid(axis='both', color='lightgray', ls='-', lw=0.5)

# plot soil thrust with respec to c_u
axs[1].plot(c_u_list, Fx, color='silver', linewidth=8, label=r'$F_x$')
axs[1].plot(c_u_list, Fx_t, '--', color='grey', label=r'$F_{xt}$')
axs[1].plot(c_u_list, Fx_b, '-', color='grey', label=r'$F_{xb}$')

# show markers
axs[1].plot(c_u_test, Fx_test, linewidth=0, marker='o', ms=6, mfc='white', mec='black', label='_nolegend_')

# set axes title
axs[1].set_xlabel("Undrained Shear Strength (kPa)")
axs[1].set_ylabel("Soil Thrust (kN/m)")

# set min values of x and y axes
axs[1].set_ylim(0, 4)
axs[1].set_xlim(0, 30)

# set grid
axs[1].tick_params(direction='in', bottom=True, top=True, left=True, right=True)
axs[1].grid(axis='both', color='lightgray', ls='-', lw=0.5)

# set annotations
axs[1].text(20, 0.05,
         r'$H$ = ' + str(round(H, 4)) + " m \n" +
         r'$L$ = ' + str(round(L, 4)) + " m \n" +
         r'$W_g$ = ' + str(round(W_g, 4)) + " kN/m",
         ha="left", va="bottom", size=10,
         )

# set legend
axs[1].legend(fontsize=10)

plt.savefig("test.svg", format="svg")
plt.show()