import numpy as np
from scipy.optimize import bisect
import matplotlib.pyplot as plt
import pandas as pd


# function to solve
def fun_t(x, c1, c2):
    return x ** 3 + c1 * x ** 2 - c2


def fun_p(x, c1, c2, c3, c4):
    return x ** 2 + c1 * (c2 * x + c3) * np.sqrt(1 - x ** 2) + c4 * x - 2


# height (H), length (L), and width (D) of soil block
H = 0.035
L = 0.095
D = 0.15
beta = np.radians(30)
W_g = 0.2

# set test number
test_no = 0

# vertical force (kN) transmitted to the track system
W_g_test = [0.0860, 0.1719, 0.2579, 0.3454]

beta_test = [np.radians(0), np.radians(0), np.radians(0), np.radians(0)]

# construct test data
c_u_test = [[7.5, 8.8, 21.7, 27.2],
            [7.5, 8.8, 21.7, 27.2],
            [7.7],
            [7.7]]
Fx_test = [[0.1463, 0.1643, 0.2963, 0.3513],
           [0.1788, 0.2008, 0.3288, 0.4275],
           [0.1924],
           [0.1982]]

# if test data is used, override variables
if test_no > 0:

    # set dimensions of the track system
    H = 0.0357
    L = 0.1239
    D = 0.1239

    # determine Wg from test number
    W_g = W_g_test[test_no - 1]

    beta = beta_test[test_no - 1]

    # set test data shown in plot
    c_u_show = c_u_test[test_no - 1]
    Fx_show = Fx_test[test_no - 1]

for b in range(-45, 60, 15):

    # undrained strength of clay (in kPa)
    c_u_list = []
    Fx_b = []
    Fx_t = []
    Fx_p = []
    Fx = []

    beta = np.radians(b)

    # set the angle at the tip of the plate
    theta_tip = np.arctan(L / H)
    sin_theta_tip = np.sin(theta_tip)

    for i in range(1, 601):
        # increase c_u
        c_u = i * 0.05
        c_u_list.append(c_u)

        """
        Block failure mode
        """

        # calculate the soil thrust for the block failure
        Fx_b.append(c_u * L * (2 * H + D) - W_g * np.sin(beta))

        """
        Triangular wedge failure mode
        """

        # check the start and end values for the bisection method
        bisect_start_value = fun_t(0,
                                   W_g * np.cos(beta) / c_u / H ** 2 + 2 * D / H,
                                   W_g * np.cos(beta) / c_u / H ** 2 + D / H)
        bisect_end_value = fun_t(sin_theta_tip,
                                 W_g * np.cos(beta) / c_u / H ** 2 + 2 * D / H,
                                 W_g * np.cos(beta) / c_u / H ** 2 + D / H)

        # if solution exists
        if bisect_start_value * bisect_end_value < 0:
            # find angle of the failure surface with respect to the vertica line: triangle
            sin_theta_t = bisect(fun_t, 0, sin_theta_tip, args=(W_g * np.cos(beta) / c_u / H ** 2 + 2 * D / H,
                                                                W_g * np.cos(beta) / c_u / H ** 2 + D / H))
            theta_rad_t = np.arcsin(sin_theta_t)
            theta_deg_t = theta_rad_t * 180 / np.pi

            # calculate the soil thrust for the triangular wedge failure
            Fx_t.append(W_g * (np.cos(theta_rad_t) * np.cos(beta) / np.sin(theta_rad_t) - np.sin(beta))
                        + c_u * H * (D / np.sin(theta_rad_t) / np.cos(theta_rad_t) + H / np.cos(theta_rad_t)))
        else:
            Fx_t.append(np.nan)

        # Print min soil thrust as a solution
        Fx.append(np.nanmin([Fx_b[i - 1], Fx_t[i - 1]]))

    # plot soil thrust values
    plt.plot(c_u_list, Fx, color='grey', linewidth=2, label=r'$F_x$')
    plt.plot(c_u_list, Fx_t, '--', color='black', linewidth=0.5, label=r'$F_{xt}$')
    plt.plot(c_u_list, Fx_b, '-', color='black', linewidth=0.5, label=r'$F_{xb}$')

# set axes title
plt.xlabel("Undrained Shear Strength (kPa)")
plt.ylabel("Soil Thrust (kN)")

# set min values of x and y axes
ymin = -0.2
ymax = 0.6
xmin = 0
xmax = 30

plt.ylim(bottom=ymin)
plt.ylim(top=ymax)
plt.xlim(left=xmin)
plt.xlim(right=xmax)

plt.tick_params(direction='in', bottom=True, top=True, left=True, right=True)
plt.grid(axis='both', color='lightgray', ls='-', lw=0.5)

txt_loc_x = xmin + 0.05 * (xmax - xmin)
txt_loc_y = ymin + 0.95 * (ymax - ymin)

plt.text(txt_loc_x, txt_loc_y,
         r'$H$ = ' + str(round(H, 4)) + " m \n" +
         r'$D$ = ' + str(round(D, 4)) + " m \n" +
         r'$L$ = ' + str(round(L, 4)) + " m \n" +
         r'$W_g$ = ' + str(round(W_g, 4)) + " kN",
         ha="left", va="top", size=12)

# show test results
if test_no > 0:

    fig_num = ''

    if test_no == 1:
        fig_num = '(a)'
    if test_no == 2:
        fig_num = '(b)'
    if test_no == 3:
        fig_num = '(c)'
    if test_no == 4:
        fig_num = '(d)'

    # show fig. number
    plt.text(1, 0.77, fig_num, ha="left", va="top", size=12)

    # show markers
    plt.plot(c_u_show, Fx_show, linewidth=0,
             marker='o', ms=8, mfc='white', mec='black',
             label='_nolegend_')

    for i in range(0, len(c_u_show)):
        if len(c_u_show) > 1:
            plt.text(c_u_show[i], Fx_show[i] - 0.02, str(test_no) + '-' + str(i + 1),
                     ha="left", va="top", size=10)
        else:
            plt.text(c_u_show[i], Fx_show[i] - 0.02, str(test_no),
                     ha="left", va="top", size=10)

# set legend
#plt.legend(bbox_to_anchor=(0., 1.02, 1., .02), loc=4, ncol=4,
#           mode="expand", borderaxespad=0., frameon=False, fontsize=12)

plt.savefig("test.svg", format="svg")
plt.show()

df = pd.DataFrame(list(zip(*[c_u_list, Fx_b, Fx_t, Fx_p, Fx]))).add_prefix('Col')
df.to_csv('file.csv', index=False)
