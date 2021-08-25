import numpy as np
from scipy.optimize import bisect
import matplotlib.pyplot as plt
import pandas as pd


# function to solve
def fun_t(x, c1, c2) :
    return x ** 3 + c1 * x ** 2 - c2

def fun_p(x, c1, c2, c3, c4) :
    return x ** 2 + c1 * (c2 * x + c3) * np.sqrt(1 - x ** 2) + c4*x - 2


# height (H), length (L), and width (D) of soil block
H = 0.04
L = 0.15
D = 0.15
W_g = 0.2

# set test number
test_no = 2

"""
Experimental Results (Case 1, 2, 3, 4)
"""

# vertical force (kN) transmitted to the track system
W_g_test = [0.0860, 0.1719, 0.2579, 0.3454]

# construct test data
c_u_test = [[7.5, 8.8, 21.7, 27.2],
            [7.5, 8.8, 21.7, 27.2],
            [7.7],
            [7.7]]
Fx_test = [[0.1463, 0.1643, 0.2963, 0.3513],
           [0.1788, 0.2008, 0.3288, 0.4275],
           [0.1924],
           [0.1982]]

"""
Numerical Results (Simulation A, B, C in Case 2)
"""

c_u_nume = [4.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0]
Fx_nume_A = [0.1350, 0.1868, 0.2356, 0.2798, 0.3196, 0.3737, 0.4581, 0.5383]
Fx_nume_B = [0.1040, 0.1509, 0.1927, 0.2323, 0.2637, 0.3010, 0.3692, 0.4326]
Fx_nume_C = [0.1008, 0.1447, 0.1811, 0.2172, 0.2550, 0.2918, 0.3529, 0.4061]

# if test data is used, override variables
if test_no > 0 :

    # set dimensions of the track system
    H = 0.0357
    L = 0.1239
    D = 0.1239

    # determine Wg from test number
    W_g = W_g_test[test_no - 1]

    # set test data shown in plot
    c_u_show = c_u_test[test_no - 1]
    Fx_show = Fx_test[test_no - 1]

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
    Fx_b.append(c_u * L * (2 * H + D))

    """
    Triangular wedge failure mode
    """

    # check the start and end values for the bisection method
    bisect_start_value = fun_t(0, W_g / c_u / H ** 2 + 2 * D / H, W_g / c_u / H ** 2 + D / H)
    bisect_end_value = fun_t(sin_theta_tip, W_g / c_u / H ** 2 + 2 * D / H, W_g / c_u / H ** 2 + D / H)

    # if solution exists
    if bisect_start_value * bisect_end_value < 0:
        # find angle of the failure surface with respect to the vertica line: triangle
        sin_theta_t = bisect(fun_t, 0, sin_theta_tip, args=(W_g / c_u / H ** 2 + 2 * D / H, W_g / c_u / H ** 2 + D / H))
        theta_rad_t = np.arcsin(sin_theta_t)
        theta_deg_t = theta_rad_t * 180 / np.pi

        # calculate the soil thrust for the triangular wedge failure
        Fx_t.append(W_g * np.cos(theta_rad_t) / np.sin(theta_rad_t)
                    + c_u * H * (D / np.sin(theta_rad_t) / np.cos(theta_rad_t) + H / np.cos(theta_rad_t)))
    else:
        Fx_t.append(np.nan)


    """
    Trapezoidal wedge failure mode
    """

    # check the start and end values for the bisection method
    bisect_start_value = fun_p(sin_theta_tip, 2 / L, H, D, W_g / c_u / L ** 2)
    bisect_end_value = fun_p(1, 2 / L, H, D, W_g / c_u / L ** 2)

    # if solution exists
    if bisect_start_value * bisect_end_value < 0:
        # find angle of the failure surface with respect to the vertica line: trapezoidal
        sin_theta_p = bisect(fun_p, sin_theta_tip, 1, args=(2 / L, H, D, W_g / c_u / L ** 2))
        theta_rad_p = np.arcsin(sin_theta_p)
        theta_deg_p = theta_rad_p * 180 / np.pi
        Fx_p.append(W_g * np.cos(theta_rad_p) / np.sin(theta_rad_p)
                    + c_u * L * (D / np.sin(theta_rad_p) / np.sin(theta_rad_p)
                                 + 2 * H / np.sin(theta_rad_p)
                                 - L * np.cos(theta_rad_p) / np.sin(theta_rad_p) ** 2))
    else:
        Fx_p.append(np.nan)



    # Print min soil thrust as a solution
    Fx.append(np.nanmin([Fx_b[i-1], Fx_t[i-1], Fx_p[i-1]]))


# plot soil thrust values
plt.plot(c_u_list, Fx, color='silver', linewidth=8, label=r'$F_x$')
plt.plot(c_u_list, Fx_t, '--', color='black', label=r'$F_{xt}$')
plt.plot(c_u_list, Fx_b, '-', color='black', label=r'$F_{xb}$')
plt.plot(c_u_list, Fx_p, ':', color='black', label=r'$F_{xp}$')

# set axes title
plt.xlabel("Undrained Shear Strength (kPa)")
plt.ylabel("Soil Thrust (kN)")

# set min values of x and y axes
plt.ylim(bottom=0.)
plt.xlim(left=0.)
plt.ylim(top=0.6)
plt.xlim(right=30.)

plt.tick_params(direction='in', bottom=True, top=True, left=True, right=True)
plt.grid(axis='both', color='lightgray', ls='-', lw=0.5)

plt.text(20, 0.02,
         r'$H$ = ' + str(round(H, 4)) + " m \n" +
         r'$D$ = ' + str(round(D, 4)) + " m \n" +
         r'$L$ = ' + str(round(L, 4)) + " m \n" +
         r'$W_g$ = ' + str(round(W_g, 4)) + " kN",
         ha="left", va="bottom", size=12,
         )


# show numerical results
if test_no == 2 :

    plt.plot(c_u_nume, Fx_nume_A, '-', linewidth=1, color='blue')
    plt.plot(c_u_nume, Fx_nume_B, '-', linewidth=1, color='blue')
    plt.plot(c_u_nume, Fx_nume_C, '-', linewidth=1, color='blue')


# show test results
if test_no > 0 :

    # show markers
    plt.plot(c_u_show, Fx_show, linewidth=0,
             marker='o', ms=8, mfc='white', mec='black',
             label='_nolegend_')

    for i in range(0, len(c_u_show)) :
        if len(c_u_show) > 1 :
            plt.text(c_u_show[i], Fx_show[i] - 0.02, str(test_no) + '-' + str(i + 1),
                     ha="left", va="top", size=10)
        else :
            plt.text(c_u_show[i], Fx_show[i] - 0.02, str(test_no),
                     ha="left", va="top", size=10)



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

df = pd.DataFrame(list(zip(*[c_u_list, Fx_b, Fx_t, Fx_p, Fx]))).add_prefix('Col')
df.to_csv('file.csv', index=False)