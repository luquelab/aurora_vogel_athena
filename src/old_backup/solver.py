import numpy as np
import math
import sympy
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import fsolve
from scipy.optimize import root
from scipy.optimize import basinhopping
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull
# init 3 vertices for equilateral triangle with side length 1
# golden ratio
phi=(1+math.sqrt(5))/2
# init verts
v01 = np.array([0, 0, 0])
v02 = np.array([math.sqrt(3)/2, 1/2, 0])
v03 = np.array([0, 1, 0])

v04x, v04y, v04z, v05x, v05y, v05z, v06x, v06y, v06z, v07x, v07y, v07z, v08x, v08y, v08z, v09x, v09y, v09z, v10x, v10y, v10z, v11x, v11y, v11z, v12x, v12y, v12z = sympy.symbols("v04x v04y v04z v05x v05y v05z v06x v06y v06z v07x v07y v07z v08x v08y v08z v09x v09y v09z v10x v10y v10z v11x v11y v11z v12x v12y v12z", real=True, positive=True)
wx, wy, wz = sympy.symbols("wx wy wz", real=True)


def grab_valid_soln(solutions, solv_type):
    if solv_type == "tetra":
        # sol via sympy solve
        # for sol in solutions:
        #     if sol[list(sol.keys())[-1]] >= 0:
        #         valid_soln = sol
        #         break

        # valid_soln_keys = list(valid_soln.keys())
        # valid_soln = np.array([valid_soln[valid_soln_keys[0]],
        #                        valid_soln[valid_soln_keys[1]],
        #                        valid_soln[valid_soln_keys[2]]])

        # temp fix, sympy nsolve
        valid_soln = solutions
        valid_soln = np.array([valid_soln[0], valid_soln[1], valid_soln[2]])
        return valid_soln

    if solv_type == "octa":
        for sol in solutions:
            i = 2
            val = True
            for k in range(int(len(sol)/3)):
                key = list(sol.keys())[i + 3 * k]
                if sol[key] < 0:
                    val = False
            if val:
                valid_soln = sol
                break
        solns = []
        i = 0
        for k in range(int(len(valid_soln) / 3)):
            v_keys = list(valid_soln.keys())
            solns.append(np.array([valid_soln[v_keys[0 + i]],
                                         valid_soln[v_keys[1 + i]],
                                         valid_soln[v_keys[2 + i]]]))
            i += 3
        return solns

########################################################################
# Successful Tetrahedron # new approach?
########################################################################
def tetrahedron_solver_new(dat):
    r1 = dat['r1']
    r2 = dat['r2']
    r3 = dat['r3']
    v01 = dat['v01']
    v02 = dat['v02']
    v03 = dat['v03']

    # eq1 = sympy.Eq((v01[0] - v04x)**2 + (v01[1] - v04y)**2 + (v01[2] - v04z)**2, r1 ** 2)
    # eq2 = sympy.Eq((v02[0] - v04x)**2 + (v02[1] - v04y)**2 + (v02[2] - v04z)**2, r2 ** 2)
    # eq3 = sympy.Eq((v03[0] - v04x)**2 + (v03[1] - v04y)**2 + (v03[2] - v04z)**2, r3 ** 2)
    eq1 = (v01[0] - v04x)**2 + (v01[1] - v04y)**2 + (v01[2] - v04z)**2 - r1 ** 2
    eq2 = (v02[0] - v04x) ** 2 + (v02[1] - v04y) ** 2 + (v02[2] - v04z) ** 2 - r2 ** 2
    eq3 = (v03[0] - v04x)**2 + (v03[1] - v04y)**2 + (v03[2] - v04z)**2 - r3 ** 2

    # return grab_valid_soln(sympy.solve([eq1, eq2, eq3]), "tetra")
    return grab_valid_soln(sympy.nsolve((eq1, eq2, eq3), (v04x, v04y, v04z), (1, 1, 1)), "tetra")

########################################################################
# Unsuccessful Tetrahedron
########################################################################
def tetrahedron_solver():
    # base
    eq1 = sympy.Eq((v02[0] - v04x)**2 + (v02[1] - v04y)**2 + (v02[2] - v04z)**2, 1.0**2)
    eq2 = sympy.Eq((v03[0] - v04x)**2 + (v03[1] - v04y)**2 + (v03[2] - v04z)**2, 1.0**2)
    # cnx to top
    eq3 = sympy.Eq((v04x - v01[0])**2 + (v04y - v01[1])**2 + (v04z - v01[2])**2, 1.0**2)

    return sympy.solve([eq1, eq2, eq3])
#    sympy.nsolve((eq1, eq2, eq3),
#                 (v04x, v04y, v04z),
#                 (0,0,0))

########################################################################
# Successful octahedron
########################################################################
def my_func(x, eqns):
    eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9 = eqns
    # v04x v04y v04z
    #            0
    # v05x v05y v05z
    #            1
    # v06x v06y v06z
    #  2    3    4


    eq1 = eq1.subs(v04z, x[0])
    eq2 = eq2.subs(v04z, x[0])
    eq3 = eq3.subs(v05z, x[1])
    eq4 = eq4.subs(v05z, x[1])
    eq5 = eq5.subs((v06x, v06y, v06z), (x[2], x[3], x[4]))
    eq6 = eq6.subs((v06x, v06y, v06z), (x[2], x[3], x[4]))
    eq7 = eq7.subs((v04z, v06x, v06y, v06z), (x[0], x[2], x[3], x[4]))
    eq8 = eq8.subs((v05z, v06x, v06y, v06z), (x[1], x[2], x[3], x[4]))
    eq9 = eq9.subs((v04z, v05z), (x[0], x[1]))

    return float(eq1.lhs + eq2.lhs + eq3.lhs + eq4.lhs + eq5.lhs + eq6.lhs + eq7.lhs + eq8.lhs + eq9.lhs)


def my_func_n(x):
    # v04x v04y v04z
    #  0    1    2
    # v05x v05y v05z
    #  3    4    5
    # v06x v06y v06z
    #  6    7    8
    eq1 = (v01[0] - x[0]) ** 2 + (v01[1] - x[1]) ** 2 + (v01[2] - x[2]) ** 2 - 1.0 ** 2
    eq2 = (v02[0] - x[0]) ** 2 + (v02[1] - x[1]) ** 2 + (v02[2] - x[2]) ** 2 - 1.0 ** 2
    eq3 = (v02[0] - x[3]) ** 2 + (v02[1] - x[4]) ** 2 + (v02[2] - x[5]) ** 2 - 1.0 ** 2
    eq4 = (v03[0] - x[3]) ** 2 + (v03[1] - x[4]) ** 2 + (v03[2] - x[5]) ** 2 - 1.0 ** 2
    eq5 = (v01[0] - x[6]) ** 2 + (v01[1] - x[7]) ** 2 + (v01[2] - x[8]) ** 2 - 1.0 ** 2
    eq6 = (v03[0] - x[6]) ** 2 + (v03[1] - x[7]) ** 2 + (v03[2] - x[8]) ** 2 - 1.0 ** 2
    eq7 = (x[0] - x[6]) ** 2 + (x[1] - x[7]) ** 2 + (x[2] - x[8]) ** 2 - 1.0 ** 2
    eq8 = (x[3] - x[6]) ** 2 + (x[4] - x[7]) ** 2 + (x[5] - x[8]) ** 2 - 1.0 ** 2
    eq9 = (x[0] - x[3]) ** 2 + (x[1] - x[4]) ** 2 + (x[2] - x[5]) ** 2 - 1.0 ** 2

    return eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9


def penalty_function(v):
    alpha = 10**-5
    penalty = 0
    my_pts = []
    for i in range(int(len(v)/3)):
        my_pts.append(np.array([v[i*3], v[i*3+1], v[i*3+2]]))
#    print("my pts are: " + str(my_pts))
    for i in range(len(my_pts)):
        for j in range(len(my_pts)):
            if i < j:
                if np.linalg.norm(my_pts[i] - my_pts[j]) < alpha:
#                    print("pts " + str(my_pts[i]) + " and " + str(my_pts[j]))
#                    print("have norm: " + str(np.linalg.norm(my_pts[i] - my_pts[j])))
#                    print("and first portion: " + str(1/np.linalg.norm(my_pts[i] - my_pts[j])))
#                    print("and penalty: " + str(abs(1/np.linalg.norm(my_pts[i] - my_pts[j]) - 1/alpha)))
                    penalty += abs(1/np.linalg.norm(my_pts[i] - my_pts[j]) - 1/alpha)
    return penalty


def my_func_n2(x):
    # v04x v04y v04z
    #  0    1    2
    # v05x v05y v05z
    #  3    4    5
    # v06x v06y v06z
    #  6    7    8
    v_pts = [v01[0], v01[1], v01[2],
             v02[0], v02[1], v02[2],
             v03[0], v03[1], v03[2]]
    penalty = penalty_function(np.array(v_pts + list(x)))
    eq1 = abs((v01[0] - x[0]) ** 2 + (v01[1] - x[1]) ** 2 + (v01[2] - x[2]) ** 2 - 1.0 ** 2) + penalty
    eq2 = abs((v02[0] - x[0]) ** 2 + (v02[1] - x[1]) ** 2 + (v02[2] - x[2]) ** 2 - 1.0 ** 2) + penalty
    eq3 = abs((v02[0] - x[3]) ** 2 + (v02[1] - x[4]) ** 2 + (v02[2] - x[5]) ** 2 - 1.0 ** 2) + penalty
    eq4 = abs((v03[0] - x[3]) ** 2 + (v03[1] - x[4]) ** 2 + (v03[2] - x[5]) ** 2 - 1.0 ** 2) + penalty
    eq5 = abs((v01[0] - x[6]) ** 2 + (v01[1] - x[7]) ** 2 + (v01[2] - x[8]) ** 2 - 1.0 ** 2) + penalty
    eq6 = abs((v03[0] - x[6]) ** 2 + (v03[1] - x[7]) ** 2 + (v03[2] - x[8]) ** 2 - 1.0 ** 2) + penalty
    eq7 = abs((x[0] - x[6]) ** 2 + (x[1] - x[7]) ** 2 + (x[2] - x[8]) ** 2 - 1.0 ** 2) + penalty
    eq8 = abs((x[3] - x[6]) ** 2 + (x[4] - x[7]) ** 2 + (x[5] - x[8]) ** 2 - 1.0 ** 2) + penalty
    eq9 = abs((x[0] - x[3]) ** 2 + (x[1] - x[4]) ** 2 + (x[2] - x[5]) ** 2 - 1.0 ** 2) + penalty
    print(penalty)
    print(str(eq1 + eq2 + eq3 + eq4 + eq5 + eq6 + eq7 + eq8 + eq9))
    return eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9
    # return eq1 + eq2 + eq3 + eq4 + eq5 + eq6 + eq7 + eq8 + eq9


def octahedron_nsolver_new():
    # round 1
    # example that works, very sensitive to input
    # x = root(my_func_n, x0=[1.0, 2.0, 3.0, 1.5, 2.5, 3.5, 4, 5, 8], method='lm').x
    x = root(my_func_n, x0=[1.0, 2.0, 3.0, 1.5, 2.5, 3.5, 4, 5, 8], method='lm').x
    # how to best punish pts that are the same....
    # x0_ran = 3*np.random.random_sample(9)
    # results = root(my_func_n2, x0=x0_ran, method='broyden2', tol=1e-7, options={'maxiter': 10000})
    # minimizer_kwargs = {"method": "L-BFGS-B", "bounds": [(low, high) for low, high in zip([0, 0, 0, 0, 0, 0, 0, 0, 0],
    #                                                                                [4, 4, 4, 4, 4, 4, 4, 4, 4])]}
    # results = basinhopping(my_func_n2, x0=x0_ran, minimizer_kwargs=minimizer_kwargs, niter=200)
    # print(results)
    # x = results.x
    print(my_func_n2(x))
    X = v01[0], v02[0], v03[0], x[0], x[3], x[6]
    Y = v01[1], v02[1], v03[1], x[1], x[4], x[7]
    Z = v01[2], v02[2], v03[2], x[2], x[5], x[8]
    fig = plt.figure()
    face_set = np.array([v01, v02, v03, np.array([x[0], x[1], x[2]]),
                         np.array([x[3], x[4], x[5]]),
                         np.array([x[6], x[7], x[8]])])
    print(face_set)
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.scatter3D(X, Y, Z, c='g')
    hull = ConvexHull(face_set)
    # draw the polygons of the convex hull
    color = ['b', 'r', 'r', 'r']
    c = 0
    for s in hull.simplices:
        tri = Poly3DCollection(face_set[s])
        tri.set_color(color[0])
        tri.set_edgecolor('k')
        tri.set_alpha(0.6)
        ax.add_collection3d(tri)
        c += 1
    # draw the vertices
    ax.scatter(face_set[:, 0], face_set[:, 1], face_set[:, 2], marker='o', color='purple')
    set_axes_equal(ax)
    plt.show()
    return 0


def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    # Took this from
    # https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


def octahedron_solver_new():
    eq1 = sympy.Eq((v01[0] - v04x)**2 + (v01[1] - v04y)**2 + (v01[2] - v04z)**2, 1.0**2)
    eq2 = sympy.Eq((v02[0] - v04x)**2 + (v02[1] - v04y)**2 + (v02[2] - v04z)**2, 1.0**2)
    eq3 = sympy.Eq((v02[0] - v05x) ** 2 + (v02[1] - v05y) ** 2 + (v02[2] - v05z) ** 2, 1.0 ** 2)
    eq4 = sympy.Eq((v03[0] - v05x) ** 2 + (v03[1] - v05y) ** 2 + (v03[2] - v05z) ** 2, 1.0 ** 2)
    eq5 = sympy.Eq((v01[0] - v06x) ** 2 + (v01[1] - v06y) ** 2 + (v01[2] - v06z) ** 2, 1.0 ** 2)
    eq6 = sympy.Eq((v03[0] - v06x) ** 2 + (v03[1] - v06y) ** 2 + (v03[2] - v06z) ** 2, 1.0 ** 2)
    eq7 = sympy.Eq((v04x - v06x) ** 2 + (v04y - v06y) ** 2 + (v04z - v06z) ** 2, 1.0 ** 2)
    eq8 = sympy.Eq((v05x - v06x) ** 2 + (v05y - v06y) ** 2 + (v05z - v06z) ** 2, 1.0 ** 2)
    eq9 = sympy.Eq((v04x - v05x) ** 2 + (v04y - v05y) ** 2 + (v04z - v05z) ** 2, 1.0 ** 2)
    # solns 1
    solns = sympy.solve([eq1, eq2, eq3, eq4])
    for i in range(len(solns)):
        soln = solns[i]
        v04x_n = soln[v04x]
        v04y_n = soln[v04y]
        v05x_n = soln[v05x]
        v05y_n = soln[v05y]
        # round 2
        eq1 = sympy.Eq((v01[0] - v04x_n) ** 2 + (v01[1] - v04y_n) ** 2 + (v01[2] - v04z) ** 2 - 1.0 ** 2, 0**2)
        eq2 = sympy.Eq((v02[0] - v04x_n) ** 2 + (v02[1] - v04y_n) ** 2 + (v02[2] - v04z) ** 2 - 1.0 ** 2, 0**2)
        eq3 = sympy.Eq((v02[0] - v05x_n) ** 2 + (v02[1] - v05y_n) ** 2 + (v02[2] - v05z) ** 2 - 1.0 ** 2, 0**2)
        eq4 = sympy.Eq((v03[0] - v05x_n) ** 2 + (v03[1] - v05y_n) ** 2 + (v03[2] - v05z) ** 2 - 1.0 ** 2, 0**2)
        eq5 = sympy.Eq((v01[0] - v06x) ** 2 + (v01[1] - v06y) ** 2 + (v01[2] - v06z) ** 2 - 1.0 ** 2, 0**2)
        eq6 = sympy.Eq((v03[0] - v06x) ** 2 + (v03[1] - v06y) ** 2 + (v03[2] - v06z) ** 2 - 1.0 ** 2, 0**2)
        eq7 = sympy.Eq((v04x_n - v06x) ** 2 + (v04y_n - v06y) ** 2 + (v04z - v06z) ** 2 - 1.0 ** 2, 0**2)
        eq8 = sympy.Eq((v05x_n - v06x) ** 2 + (v05y_n - v06y) ** 2 + (v05z - v06z) ** 2 - 1.0 ** 2, 0**2)
        eq9 = sympy.Eq((v04x_n - v05x_n) ** 2 + (v04y_n - v05y_n) ** 2 + (v04z - v05z) ** 2 - 1.0 ** 2, 0**2)

        # eq1 = (v01[0] - v04x_n) ** 2 + (v01[1] - v04y_n) ** 2 + (v01[2] - v04z) ** 2 - 1.0 ** 2
        # eq2 = (v02[0] - v04x_n) ** 2 + (v02[1] - v04y_n) ** 2 + (v02[2] - v04z) ** 2 - 1.0 ** 2
        # eq3 = (v02[0] - v05x_n) ** 2 + (v02[1] - v05y_n) ** 2 + (v02[2] - v05z) ** 2 - 1.0 ** 2
        # eq4 = (v03[0] - v05x_n) ** 2 + (v03[1] - v05y_n) ** 2 + (v03[2] - v05z) ** 2 - 1.0 ** 2
        # eq5 = (v01[0] - v06x) ** 2 + (v01[1] - v06y) ** 2 + (v01[2] - v06z) ** 2 - 1.0 ** 2
        # eq6 = (v03[0] - v06x) ** 2 + (v03[1] - v06y) ** 2 + (v03[2] - v06z) ** 2 - 1.0 ** 2
        # eq7 = (v04x_n - v06x) ** 2 + (v04y_n - v06y) ** 2 + (v04z - v06z) ** 2 - 1.0 ** 2
        # eq8 = (v05x_n - v06x) ** 2 + (v05y_n - v06y) ** 2 + (v05z - v06z) ** 2 - 1.0 ** 2
        # eq9 = (v04x_n - v05x) ** 2 + (v04y_n - v05y) ** 2 + (v04z - v05z) ** 2 - 1.0 ** 2

        print(sympy.solve([eq1, eq2, eq3, eq4, eq9]))

        # solns 2
        # print(sympy.nsolve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9,), (v04z, v05z,
        #                                                             v06x, v06y, v06z,), (0, 0, 0, 0, 0,)))
        # fsolve(my_func, args=([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9]), x0=[1.0, 1.0, 1.0, 1.0, 1.0])

    grab_valid_soln(sympy.solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9]), "octa")
    # This one runs forever
    return sympy.solve([eq1, eq2, eq3, eq4, eq5, eq6])

    # This one quits after some time
    # "RecursionError: maximum recursion depth exceeded in comparison"
    # sympy.nonlinsolve([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9],
    #          [v04x, v04y, v04z, v05x, v05y, v05z, v06x, v06y, v06z])

########################################################################
# Unsuccessful octahedron
########################################################################
def octahedron_solver():
    # hardcoded triangle side lengths 1
    # Top rim
    eq1 = sympy.Eq((v01[0] - v04x)**2 + (v01[1] - v04y)**2 + (v01[2] - v04z)**2, 1.0**2)
    eq2 = sympy.Eq((v01[0] - v05x)**2 + (v01[1] - v05y)**2 + (v01[2] - v05z)**2, 1.0**2)

    # Mid cnxn
    eq3 = sympy.Eq((v03[0] - v04x)**2 + (v03[1] - v04y)**2 + (v03[2] - v04z)**2, 1.0**2)
    eq4 = sympy.Eq((v04x - v05x)**2 + (v04y - v05y)**2 + (v04z - v05z)**2, 1.0**2)
    eq5 = sympy.Eq((v05x - v02[0])**2 + (v05y - v02[1])**2 + (v05z - v02[2])**2, 1.0**2)

    # Bottom rim
    eq6 = sympy.Eq((v02[0] - v06x)**2 + (v02[1] - v06y)**2 + (v02[2] - v06z)**2, 1.0**2)
    eq7 = sympy.Eq((v03[0] - v06x)**2 + (v03[1] - v06y)**2 + (v03[2] - v06z)**2, 1.0**2)
    eq8 = sympy.Eq((v04x - v06x)**2 + (v04y - v06y)**2 + (v04z - v06z)**2, 1.0**2)
    eq9 = sympy.Eq((v05x - v06x)**2 + (v05y - v06y)**2 + (v05z - v06z)**2, 1.0**2)

    # This one runs forever
    return sympy.solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9])

    # This one quits after some time
    # "RecursionError: maximum recursion depth exceeded in comparison"
    # sympy.nonlinsolve([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9],
    #          [v04x, v04y, v04z, v05x, v05y, v05z, v06x, v06y, v06z])

########################################################################
# Unsuccessful icosahedron
########################################################################


def icosahedron_solver1():
    # Top rim
    eq1 = sympy.Eq((v03[0] - v04x)**2 + (v03[1] - v04y)**2 + (v03[2] - v04z)**2, 1.0**2)
    eq2 = sympy.Eq((v04x - v05x)**2 + (v04y - v05y)**2 + (v04z - v05z)**2, 1.0**2)
    eq3 = sympy.Eq((v05x - v06x)**2 + (v05y - v06y)**2 + (v05z - v06z)**2, 1.0**2)
    eq4 = sympy.Eq((v06x - v02[0])**2 + (v06y - v02[1])**2 + (v06z - v02[2])**2, 1.0**2)
    # cnxn to top rim
    eq5 = sympy.Eq((v04x - v01[0])**2 + (v04y - v01[1])**2 + (v04z - v01[2])**2, 1.0**2)
    eq6 = sympy.Eq((v05x - v01[0])**2 + (v05y - v01[1])**2 + (v05z - v01[2])**2, 1.0**2)
    eq7 = sympy.Eq((v06x - v01[0])**2 + (v06y - v01[1])**2 + (v06z - v01[2])**2, 1.0**2)

    ## body cnxn to top
    eq8 = sympy.Eq((v07x - v02[0])**2 + (v07y - v02[1])**2 + (v07z - v02[2])**2, 1.0**2)
    eq9 = sympy.Eq((v08x - v02[0])**2 + (v08y - v02[1])**2 + (v08z - v02[2])**2, 1.0**2)
    eq10 = sympy.Eq((v08x - v03[0])**2 + (v08y - v03[1])**2 + (v08z - v03[2])**2, 1.0**2)
    eq11 = sympy.Eq((v09x - v03[0])**2 + (v09y - v03[1])**2 + (v09z - v03[2])**2, 1.0**2)
    eq12 = sympy.Eq((v09x - v04x)**2 + (v09y - v04y)**2 + (v09z - v04z)**2, 1.0**2)
    eq13 = sympy.Eq((v10x - v04x)**2 + (v10y - v04y)**2 + (v10z - v04z)**2, 1.0**2)
    eq14 = sympy.Eq((v10x - v05x)**2 + (v10y - v05y)**2 + (v10z - v05z)**2, 1.0**2)
    eq15 = sympy.Eq((v11x - v05x)**2 + (v11y - v05y)**2 + (v11z - v05z)**2, 1.0**2)
    eq16 = sympy.Eq((v11x - v06x)**2 + (v11y - v06y)**2 + (v11z - v06z)**2, 1.0**2)
    eq17 = sympy.Eq((v07x - v06x)**2 + (v07y - v06y)**2 + (v07z - v06z)**2, 1.0**2)

    # bottom rim
    eq18 = sympy.Eq((v07x - v08x)**2 + (v07y - v08y)**2 + (v07z - v08z)**2, 1.0**2)
    eq19 = sympy.Eq((v08x - v09x)**2 + (v08y - v09y)**2 + (v08z - v09z)**2, 1.0**2)
    eq20 = sympy.Eq((v09x - v10x)**2 + (v09y - v10y)**2 + (v09z - v10z)**2, 1.0**2)
    eq21 = sympy.Eq((v10x - v11x)**2 + (v10y - v11y)**2 + (v10z - v11z)**2, 1.0**2)
    eq22 = sympy.Eq((v11x - v07x)**2 + (v11y - v07y)**2 + (v11z - v07z)**2, 1.0**2)

    # Cnxn to bottom rim
    eq23 = sympy.Eq((v07x - v12x)**2 + (v07y - v12y)**2 + (v07z - v12z)**2, 1.0**2)
    eq24 = sympy.Eq((v08x - v12x)**2 + (v08y - v12y)**2 + (v08z - v12z)**2, 1.0**2)
    eq25 = sympy.Eq((v09x - v12x)**2 + (v09y - v12y)**2 + (v09z - v12z)**2, 1.0**2)
    eq26 = sympy.Eq((v10x - v12x)**2 + (v10y - v12y)**2 + (v10z - v12z)**2, 1.0**2)
    eq27 = sympy.Eq((v11x - v12x)**2 + (v11y - v12y)**2 + (v11z - v12z)**2, 1.0**2)

    # "RecursionError: maximum recursion depth exceeded in comparison"
    sympy.nonlinsolve([eq1, eq2, eq3, eq4, eq5, eq6, eq7,
                 eq8, eq9, eq10, eq11, eq12, eq13, eq14,
                 eq15, eq16, eq17, eq18, eq19, eq20, eq21,
                 eq22, eq23, eq24, eq25, eq26, eq27],
                      [v04x, v04y, v04z, v05x, v05y, v05z, v06x, v06y, v06z, v07x, v07y, v07z, v08x, v08y, v08z, v09x, v09y, v09z, v10x, v10y, v10z, v11x, v11y, v11z, v12x, v12y, v12z])

    # gets stuck
    #sympy.solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7,
    #             eq8, eq9, eq10, eq11, eq12, eq13, eq14,
    #             eq15, eq16, eq17, eq18, eq19, eq20, eq21,
    #             eq22, eq23, eq24, eq25, eq26, eq27])

    # gets stuck
    sympy.solve_poly_system([eq1, eq2, eq3, eq4, eq5, eq6, eq7,
                 eq8, eq9, eq10, eq11, eq12, eq13, eq14,
                 eq15, eq16, eq17, eq18, eq19, eq20, eq21,
                 eq22, eq23, eq24, eq25, eq26, eq27],v04x, v04y, v04z, v05x, v05y, v05z, v06x, v06y, v06z, v07x, v07y, v07z, v08x, v08y, v08z, v09x, v09y, v09z, v10x, v10y, v10z, v11x, v11y, v11z, v12x, v12y, v12z)

###################################################################


def icosahedron_solver2():
    # Top rim
    eq1 = (v03[0] - v04x)**2 + (v03[1] - v04y)**2 + (v03[2] - v04z)**2
    eq2 = (v04x - v05x)**2 + (v04y - v05y)**2 + (v04z - v05z)**2
    eq3 = (v05x - v06x)**2 + (v05y - v06y)**2 + (v05z - v06z)**2
    eq4 = (v06x - v02[0])**2 + (v06y - v02[1])**2 + (v06z - v02[2])**2
    # cnxn to top rim
    eq5 = (v04x - v01[0])**2 + (v04y - v01[1])**2 + (v04z - v01[2])**2
    eq6 = (v05x - v01[0])**2 + (v05y - v01[1])**2 + (v05z - v01[2])**2
    eq7 = (v06x - v01[0])**2 + (v06y - v01[1])**2 + (v06z - v01[2])**2

    ## body cnxn to top
    eq8 = (v07x - v02[0])**2 + (v07y - v02[1])**2 + (v07z - v02[2])**2
    eq9 = (v08x - v02[0])**2 + (v08y - v02[1])**2 + (v08z - v02[2])**2
    eq10 = (v08x - v03[0])**2 + (v08y - v03[1])**2 + (v08z - v03[2])**2
    eq11 = (v09x - v03[0])**2 + (v09y - v03[1])**2 + (v09z - v03[2])**2
    eq12 = (v09x - v04x)**2 + (v09y - v04y)**2 + (v09z - v04z)**2
    eq13 = (v10x - v04x)**2 + (v10y - v04y)**2 + (v10z - v04z)**2
    eq14 = (v10x - v05x)**2 + (v10y - v05y)**2 + (v10z - v05z)**2
    eq15 = (v11x - v05x)**2 + (v11y - v05y)**2 + (v11z - v05z)**2
    eq16 = (v11x - v06x)**2 + (v11y - v06y)**2 + (v11z - v06z)**2
    eq17 = (v07x - v06x)**2 + (v07y - v06y)**2 + (v07z - v06z)**2

    # bottom rim
    eq18 = (v07x - v08x)**2 + (v07y - v08y)**2 + (v07z - v08z)**2
    eq19 = (v08x - v09x)**2 + (v08y - v09y)**2 + (v08z - v09z)**2
    eq20 = (v09x - v10x)**2 + (v09y - v10y)**2 + (v09z - v10z)**2
    eq21 = (v10x - v11x)**2 + (v10y - v11y)**2 + (v10z - v11z)**2
    eq22 = (v11x - v07x)**2 + (v11y - v07y)**2 + (v11z - v07z)**2

    # Cnxn to bottom rim
    eq23 = (v07x - v12x)**2 + (v07y - v12y)**2 + (v07z - v12z)**2
    eq24 = (v08x - v12x)**2 + (v08y - v12y)**2 + (v08z - v12z)**2
    eq25 = (v09x - v12x)**2 + (v09y - v12y)**2 + (v09z - v12z)**2
    eq26 = (v10x - v12x)**2 + (v10y - v12y)**2 + (v10z - v12z)**2
    eq27 = (v11x - v12x)**2 + (v11y - v12y)**2 + (v11z - v12z)**2

    # "ZeroDivisionError: matrix is numerically singular"
    return sympy.nsolve((eq1, eq2, eq3, eq4, eq5, eq6, eq7,
            eq8, eq9, eq10, eq11, eq12, eq13, eq14,
            eq15, eq16, eq17, eq18, eq19, eq20, eq21,
            eq22, eq23, eq24, eq25, eq26, eq27),
                (v04x, v04y, v04z, v05x, v05y, v05z, v06x, v06y, v06z, v07x, v07y, v07z, v08x, v08y, v08z, v09x, v09y, v09z, v10x, v10y, v10z, v11x, v11y, v11z, v12x, v12y, v12z),(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))


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

    return [0,0,0,0,0,0,0,0,0,eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11,
            eq12, eq13, eq14, eq15, eq16, eq17, eq18, eq19, eq20, eq21,
            eq22, eq23, eq24, eq25, eq26, eq27]


def icosahedron_solver3():  # scipy
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

    roots = fsolve(func, init_verts)
    return roots

# x = roots
# v01 = x[0:3]
# v02 = x[3:6]
# v03 = x[6:9]
# v04 = x[9:12]
# v05 = x[12:15]
# v06 = x[15:18]
# v07 = x[18:21]
# v08 = x[21:24]
# v09 = x[24:27]
# v10 = x[27:30]
# v11 = x[30:33]
# v12 = x[33:36]
#
# # placing all verts together to form a dataframe which can be plotted
# DATA=pd.DataFrame([v01,v02,v03,v04,v05,v06,v07,v08,v09,v10,v11,v12],columns=["X","Y","Z"])

#########################
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(DATA['X'], DATA['Y'], DATA['Z'])
#
# # labeling the verts plotted
# n=[1,2,3,4,5,6,7,8,9,10,11,12]
# for i, txt in enumerate(n):
#     ax.text(DATA['X'][i], DATA['Y'][i], DATA['Z'][i],'%s' % (str(i)), size=8, zorder=1,
#     color='k')
#
# ax.set_xlabel('X Label')
# ax.set_ylabel('Y Label')
# ax.set_zlabel('Z Label')
#
# plt.show()
#########################

# tetrahedron_solver()
# tetrahedron_solver_new()
# icosahedron_solver1()
# icosahedron_solver2()
# icosahedron_solver3()
## octahedron_nsolver_new()
# octahedron_solver_new()