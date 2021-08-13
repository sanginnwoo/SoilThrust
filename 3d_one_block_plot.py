"""


"""

import numpy
import numpy as np
from scipy.optimize import bisect
import matplotlib.pyplot as plt


# function to solve
def fun(x, c1, c2) :
    return x ** 3 + c1 * x ** 2 - c2


# height (H), length (L), and width (D) of soil block
H = 0.02
L = 0.1
D = 0.2
W_g = 0.35

# set test number
test_no = 2

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

c_u_nume_f1 = [[],
            [4.0, 6.0, 8.0, 10.0,
             12.0, 15.0, 20.0, 25.0],
            [],
            []]
Fx_nume_f1 = [[],
           [0.1350, 0.1868, 0.2356, 0.2798,
            0.3196, 0.3737, 0.4581, 0.5383],
           [],
           []]

c_u_nume_f2 = [[],
            [4.0, 6.0, 8.0, 10.0,
             12.0, 15.0, 20.0, 25.0],
            [],
            []]
Fx_nume_f2 = [[],
           [0.1040, 0.1509, 0.1927, 0.2323,
            0.2637, 0.3010, 0.3692, 0.4326],
           [],
           []]

c_u_nume_f3 = [[],
            [4.0, 6.0, 8.0, 10.0,
             12.0, 15.0, 20.0, 25.0],
            [],
            []]
Fx_nume_f3 = [[],
           [0.0000, 0.1447, 0.1811, 0.2158,
            0.0000, 0.0000, 0.0000, 0.0000],
           [],
           []]

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

    c_u_nume_show_f1 = c_u_nume_f1[test_no - 1]
    Fx_nume_show_f1 = Fx_nume_f1[test_no - 1]

    c_u_nume_show_f2 = c_u_nume_f2[test_no - 1]
    Fx_nume_show_f2 = Fx_nume_f2[test_no - 1]

    c_u_nume_show_f3 = c_u_nume_f3[test_no - 1]
    Fx_nume_show_f3 = Fx_nume_f3[test_no - 1]


# undrained strength of clay (in kPa)
c_u_list = [0]
Fx_b = [0]
Fx_t = [0]

for i in range(1, 601) :
    # increase c_u
    c_u = i * 0.05
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

# get soil thrust values with respect to undrained shear strength
Fx = np.minimum(Fx_t, Fx_b)

# plot soil thrust values
plt.plot(c_u_list, Fx, color='silver', linewidth=8, label=r'$F_x$')
plt.plot(c_u_list, Fx_t, '-', color='black', label=r'$F_{xt}$')
plt.plot(c_u_list, Fx_b, '--', color='black', label=r'$F_{xb}$')

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

# show test results
if test_no > 0 :

    fig_num = ''

    if test_no == 1 :
        fig_num = '(a)'
    if test_no == 2 :
        fig_num = '(b)'
    if test_no == 3 :
        fig_num = '(c)'
    if test_no == 4 :
        fig_num = '(d)'

    # show fig. number
    plt.text(1, 0.57, fig_num, ha="left", va="top", size=12)

    # show markers
    plt.plot(c_u_show, Fx_show, linewidth=0, marker='o', ms=12, mfc='white', mec='black',
             label='Experimental Results')

    plt.plot(c_u_nume_show_f1, Fx_nume_show_f1, linewidth=1, marker='o', ms=6, mfc='red', mec='black',
             label='Numerical Results')

    plt.plot(c_u_nume_show_f2, Fx_nume_show_f2, linewidth=1, marker='o', ms=6, mfc='yellow', mec='black',
             label='Numerical Results')

    plt.plot(c_u_nume_show_f3, Fx_nume_show_f3, linewidth=1, marker='o', ms=6, mfc='blue', mec='black',
             label='Numerical Results')


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
           ncol=3,
           mode="expand",
           borderaxespad=0.,
           frameon=False,
           fontsize=12)

plt.savefig("test.svg", format="svg")
plt.show()
