import numpy as np
import math
import sympy

# init 3 vertices for equilateral triangle with side length 1
v01=np.array([0,0,0])
v02=np.array([1/2,math.sqrt(3)/2,0])
v03=np.array([0,1,0])

v04x, v04y, v04z, v05x, v05y, v05z, v06x, v06y, v06z, v07x, v07y, v07z, v08x, v08y, v08z, v09x, v09y, v09z, v10x, v10y, v10z, v11x, v11y, v11z, v12x, v12y, v12z = sympy.symbols("v04x v04y v04z v05x v05y v05z v06x v06y v06z v07x v07y v07z v08x v08y v08z v09x v09y v09z v10x v10y v10z v11x v11y v11z v12x v12y v12z", real=True)



########################################################################
# Successful Tetrahedron
########################################################################
def tetrahedron_solver():
    # base
    eq1 = sympy.Eq((v02[0] - v04x)**2 + (v02[1] - v04y)**2 + (v02[2] - v04z)**2, 1.0**2)
    eq2 = sympy.Eq((v03[0] - v04x)**2 + (v03[1] - v04y)**2 + (v03[2] - v04z)**2, 1.0**2)
    # cnx to top
    eq3 = sympy.Eq((v04x - v01[0])**2 + (v04y - v01[1])**2 + (v04z - v01[2])**2, 1.0**2)

    return sympy.solve([eq1, eq2, eq3])


########################################################################
# Unsuccessful octahedron
########################################################################
def octahedron_solver():
    # hardcoded triangle side lengths 1
    # Top rim
    eq1 = sympy.Eq((v03[0] - v04x)**2 + (v03[1] - v04y)**2 + (v03[2] - v04z)**2, 1.0**2)
    eq2 = sympy.Eq((v04x - v05x)**2 + (v04y - v05y)**2 + (v04z - v05z)**2, 1.0**2)
    eq3 = sympy.Eq((v05x - v02[0])**2 + (v05y - v02[1])**2 + (v05z - v02[2])**2, 1.0**2)

    # Mid cnxn
    eq4 = sympy.Eq((v03[0] - v04x)**2 + (v03[1] - v04y)**2 + (v03[2] - v04z)**2, 1.0**2)
    eq5 = sympy.Eq((v04x - v05x)**2 + (v04y - v05y)**2 + (v04z - v05z)**2, 1.0**2)
    eq6 = sympy.Eq((v05x - v02[0])**2 + (v05y - v02[1])**2 + (v05z - v02[2])**2, 1.0**2)

    # Bottom rim
    eq7 = sympy.Eq((v02[0] - v06x)**2 + (v02[1] - v06y)**2 + (v02[2] - v06z)**2, 1.0**2)
    eq8 = sympy.Eq((v03[0] - v06x)**2 + (v03[1] - v06y)**2 + (v03[2] - v06z)**2, 1.0**2)
    eq9 = sympy.Eq((v04x - v06x)**2 + (v04y - v06y)**2 + (v04z - v06z)**2, 1.0**2)
    eq10 = sympy.Eq((v05x - v06x)**2 + (v05y - v06y)**2 + (v05z - v06z)**2, 1.0**2)

    # This one runs forever
    # return sympy.solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10])

    # This one quits after some time
    # "RecursionError: maximum recursion depth exceeded in comparison"
    sympy.nonlinsolve([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10],
              [v04x, v04y, v04z, v05x, v05y, v05z, v06x, v06y, v06z])

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

# tetrahedron_solver()
# octahedron_solver()
# icosahedron_solver1()
# icosahedron_solver2()