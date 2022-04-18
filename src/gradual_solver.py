import numpy as np
import math
import sympy
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import fsolve


def func(x):
    r = 2
    v01 = x[0:3]
    v02 = x[3:6]
    v03 = x[6:9]
    v04 = x[9:12]
    v05 = x[12:15]
    v06 = x[15:18]
    v07 = x[18:21]
    v08 = x[21:24]
    v09 = x[24:27]
    v10 = x[27:30]
    v11 = x[30:33]
    v12 = x[33:36]

    eqa = np.array((v01[0] - v02[0]) ** 2 + (v01[1] - v02[1]) ** 2 + (v01[2] - v02[2]) ** 2) - r ** 2
    eqb = np.array((v01[0] - v03[0]) ** 2 + (v01[1] - v03[1]) ** 2 + (v01[2] - v03[2]) ** 2) - r ** 2
    eqc = np.array((v02[0] - v03[0]) ** 2 + (v02[1] - v03[1]) ** 2 + (v02[2] - v03[2]) ** 2) - r ** 2

    eq1 = np.array((v03[0] - v04[0]) ** 2 + (v03[1] - v04[1]) ** 2 + (v03[2] - v04[2]) ** 2) - r ** 2
    eq2 = np.array((v04[0] - v05[0]) ** 2 + (v04[1] - v05[1]) ** 2 + (v04[2] - v05[2]) ** 2) - r ** 2
    eq3 = np.array((v05[0] - v06[0]) ** 2 + (v05[1] - v06[1]) ** 2 + (v05[2] - v06[2]) ** 2) - r ** 2
    eq4 = np.array((v06[0] - v02[0]) ** 2 + (v06[1] - v02[1]) ** 2 + (v06[2] - v02[2]) ** 2) - r ** 2
    # cnxn to top rim
    eq5 = np.array((v04[0] - v01[0]) ** 2 + (v04[1] - v01[1]) ** 2 + (v04[2] - v01[2]) ** 2) - r ** 2
    eq6 = np.array((v05[0] - v01[0]) ** 2 + (v05[1] - v01[1]) ** 2 + (v05[2] - v01[2]) ** 2) - r ** 2
    eq7 = np.array((v06[0] - v01[0]) ** 2 + (v06[1] - v01[1]) ** 2 + (v06[2] - v01[2]) ** 2) - r ** 2

    ## body cnxn to top
    eq8 = np.array((v07[0] - v02[0]) ** 2 + (v07[1] - v02[1]) ** 2 + (v07[2] - v02[2]) ** 2) - r ** 2
    eq9 = np.array((v08[0] - v02[0]) ** 2 + (v08[1] - v02[1]) ** 2 + (v08[2] - v02[2]) ** 2) - r ** 2
    eq10 = np.array((v08[0] - v03[0]) ** 2 + (v08[1] - v03[1]) ** 2 + (v08[2] - v03[2]) ** 2) - r ** 2
    eq11 = np.array((v09[0] - v03[0]) ** 2 + (v09[1] - v03[1]) ** 2 + (v09[2] - v03[2]) ** 2) - r ** 2
    eq12 = np.array((v09[0] - v04[0]) ** 2 + (v09[1] - v04[1]) ** 2 + (v09[2] - v04[2]) ** 2) - r ** 2
    eq13 = np.array((v10[0] - v04[0]) ** 2 + (v10[1] - v04[1]) ** 2 + (v10[2] - v04[2]) ** 2) - r ** 2
    eq14 = np.array((v10[0] - v05[0]) ** 2 + (v10[1] - v05[1]) ** 2 + (v10[2] - v05[2]) ** 2) - r ** 2
    eq15 = np.array((v11[0] - v05[0]) ** 2 + (v11[1] - v05[1]) ** 2 + (v11[2] - v05[2]) ** 2) - r ** 2
    eq16 = np.array((v11[0] - v06[0]) ** 2 + (v11[1] - v06[1]) ** 2 + (v11[2] - v06[2]) ** 2) - r ** 2
    eq17 = np.array((v07[0] - v06[0]) ** 2 + (v07[1] - v06[1]) ** 2 + (v07[2] - v06[2]) ** 2) - r ** 2

    # bottom rim
    eq18 = np.array((v07[0] - v08[0]) ** 2 + (v07[1] - v08[1]) ** 2 + (v07[2] - v08[2]) ** 2) - r ** 2
    eq19 = np.array((v08[0] - v09[0]) ** 2 + (v08[1] - v09[1]) ** 2 + (v08[2] - v09[2]) ** 2) - r ** 2
    eq20 = np.array((v09[0] - v10[0]) ** 2 + (v09[1] - v10[1]) ** 2 + (v09[2] - v10[2]) ** 2) - r ** 2
    eq21 = np.array((v10[0] - v11[0]) ** 2 + (v10[1] - v11[1]) ** 2 + (v10[2] - v11[2]) ** 2) - r ** 2
    eq22 = np.array((v11[0] - v07[0]) ** 2 + (v11[1] - v07[1]) ** 2 + (v11[2] - v07[2]) ** 2) - r ** 2

    # Cnxn to bottom rim
    eq23 = np.array((v07[0] - v12[0]) ** 2 + (v07[1] - v12[1]) ** 2 + (v07[2] - v12[2]) ** 2) - r ** 2
    eq24 = np.array((v08[0] - v12[0]) ** 2 + (v08[1] - v12[1]) ** 2 + (v08[2] - v12[2]) ** 2) - r ** 2
    eq25 = np.array((v09[0] - v12[0]) ** 2 + (v09[1] - v12[1]) ** 2 + (v09[2] - v12[2]) ** 2) - r ** 2
    eq26 = np.array((v10[0] - v12[0]) ** 2 + (v10[1] - v12[1]) ** 2 + (v10[2] - v12[2]) ** 2) - r ** 2
    eq27 = np.array((v11[0] - v12[0]) ** 2 + (v11[1] - v12[1]) ** 2 + (v11[2] - v12[2]) ** 2) - r ** 2

    return [0, 0, 0, 0, 0, 0, eqa, eqb, eqc, eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11,
            eq12, eq13, eq14, eq15, eq16, eq17, eq18, eq19, eq20, eq21,
            eq22, eq23, eq24, eq25, eq26, eq27]


def fun_test(x, r):
    v01 = x[0:3]
    v02 = x[3:6]
    v03 = x[6:9]
    v04 = x[9:12]
    v05 = x[12:15]
    v06 = x[15:18]
    v07 = x[18:21]
    v08 = x[21:24]
    v09 = x[24:27]
    v10 = x[27:30]
    v11 = x[30:33]
    v12 = x[33:36]

    eq0 = np.array((v01[0] - v02[0]) ** 2 + (v01[1] - v02[1]) ** 2 + (v01[2] - v02[2]) ** 2) - r ** 2
    eq1 = np.array((v01[0] - v03[0]) ** 2 + (v01[1] - v03[1]) ** 2 + (v01[2] - v03[2]) ** 2) - r ** 2
    eq2 = np.array((v02[0] - v03[0]) ** 2 + (v02[1] - v03[1]) ** 2 + (v02[2] - v03[2]) ** 2) - r ** 2

    eq3 = np.array((v03[0] - v04[0]) ** 2 + (v03[1] - v04[1]) ** 2 + (v03[2] - v04[2]) ** 2) - r ** 2


    return [0, 0, 0, 0, 0, 0, eq0, eq1, eq2, eq3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0]


def graph_vec(x):
    v01 = x[0:3]
    v02 = x[3:6]
    v03 = x[6:9]
    v04 = x[9:12]
    v05 = x[12:15]
    v06 = x[15:18]
    v07 = x[18:21]
    v08 = x[21:24]
    v09 = x[24:27]
    v10 = x[27:30]
    v11 = x[30:33]
    v12 = x[33:36]

    # placing all verts together to form a dataframe which can be plotted
    DATA = pd.DataFrame([v01, v02, v03, v04, v05, v06, v07, v08, v09, v10, v11, v12], columns=["X", "Y", "Z"])

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(DATA['X'], DATA['Y'], DATA['Z'])

    # labeling the verts plotted
    n = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    for i, txt in enumerate(n):
        ax.text(DATA['X'][i], DATA['Y'][i], DATA['Z'][i], '%s' % (str(i)), size=8, zorder=1,
                color='k')

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    ax.set_box_aspect((np.ptp(DATA['X']), np.ptp(DATA['Y']), np.ptp(DATA['Z'])))

    plt.show()


phi = (1+math.sqrt(5))/2
init_verts = np.array([0, 1, phi, # v01
                        1, phi, 0, # v02
                        phi, 0, 1, # v03
                        0.25, -1, phi,  # v04
                       -phi, 0, 1,  # v05
                       -1, phi, 0,  # v06
                       0, 1, -phi,  # v07
                       phi, 0, -1,  # v08
                       1, -phi, 0,  # v09
                       -1, -phi, 0,  # v10
                       -phi, 0, -1,  # v11
                       0, -1, -phi,  # v12
               ])
init_verts = np.array([0, 1, phi, # v01
                        1, phi, 0, # v02
                        phi, 0, 1, # v03
                       (1+phi)/2, (0+phi)/2, (1+0)/2,  # v04
                       0, 0, 0,  # v05
                       0, 0, 0,  # v06
                       0, 0, 0,  # v07
                       0, 0, 0,  # v08
                       0, 0, 0,  # v09
                       0, 0, 0,  # v10
                       0, 0, 0,  # v11
                       0, 0, 0,  # v12
               ])

roots = fsolve(lambda x: fun_test(x, 2), init_verts)
print("inits")
print(init_verts)
print("roots")
print(roots)
print(fun_test(roots, 2))
graph_vec(roots)
