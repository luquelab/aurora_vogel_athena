import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import my_library as lib
import openpyxl
import solver as solv
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull

ICOS_CNXN = [[5, 0, 6], [6, 1, 7], [7, 2, 8], [8, 3, 9], [9, 4, 10],
             [12, 5, 6], [13, 6, 7], [14, 7, 8], [15, 8, 9], [16, 9, 10],
             [11, 5, 12], [12, 6, 13], [13, 7, 14], [14, 8, 15], [15, 9, 16],
             [17, 11, 12], [18, 12, 13], [19, 13, 14], [20, 14, 15], [21, 15, 16]] # NORMAL 1
# ICOS_CNXN = [[5, 0, 6], [6, 1, 7], [7, 2, 8], [8, 3, 14], [14, 3, 9],
#              [11, 6, 12], [12, 7, 13], [13, 8, 14], [14, 9, 15], [16, 9, 10], [9, 4, 10],
#              [11, 5, 6], [12, 6, 7], [13, 7, 8], [15, 9, 16],
#              [17, 11, 12], [18, 12, 13], [19, 13, 14], [20, 14, 15], [21, 15, 16]] # HIV 1

# 3j3q connection
ICOS_CNXN = [[5, 0, 6], [6, 1, 7], [7, 2, 8], [8, 3, 9], [13, 8, 9],
             [11, 6, 12], [12, 7, 13], [13, 9, 14], [14, 9, 15], [15, 4, 10], [9, 4, 15],
             [11, 5, 6], [12, 6, 7], [13, 7, 8], [15, 10, 16],
             [17, 11, 12], [18, 12, 13], [19, 13, 14], [20, 14, 15], [21, 15, 16]] # HIV 1


TETRA_CNXN = [[3, 0, 4], [4, 1, 5], [5, 2, 6], [4, 7, 5]]
OCTA_CNXN = [[2, 0, 3], [5, 1, 2], [5, 2, 6], [6, 2, 3], [6, 3, 7], [3, 7, 4], [7, 4, 8], [9, 6, 7]]
OCTA_CNXN_2 = [[4, 0, 5], [5, 1, 6], [6, 2, 7], [7, 3, 8], [4, 9, 5], [5, 10, 6], [6, 11, 7], [7, 12, 8]]

START = np.array([0, 0, 0])
END = np.array([0, 0, 0])
CORNERS = []


class PathSet:
    def __init__(self, upper_path=[], lower_path=[], distance_between=0, connection=[], triangles=[], num_subunits=[]):
        self.upper_path = upper_path
        self.lower_path = lower_path
        self.distance_between = distance_between
        self.connection = connection
        self.triangles = triangles
        self.num_subunits = num_subunits

    def assign_triangles(self, net):
        tris = []
        for triangle in self.connection:
            v1 = lib.hex_to_cart(net[triangle[0]])
            v2 = lib.hex_to_cart(net[triangle[1]])
            v3 = lib.hex_to_cart(net[triangle[2]])
            tris.append(Triangle(vertices=[v1, v2, v3]))
        self.triangles = tris

    def assign_num_subunits(self):
        total_area = 0
        for triangle in self.triangles:
            total_area += triangle.area
        unit_triangle = math.sqrt(3)/4
        self.num_subunits = round(total_area / unit_triangle * 3, 4)


class Triangle:
    def __init__(self, vertices=[]):
        self.vertices = vertices
        self.area = self.get_area()

    def get_area(self):
        v1 = self.vertices[0]
        v2 = self.vertices[1]
        v3 = self.vertices[2]
        s1 = np.linalg.norm(v1-v2)
        s2 = np.linalg.norm(v2-v3)
        s3 = np.linalg.norm(v1-v3)
        p = (s1 + s2 + s3)/2  # Half of the perimeter
        a = math.sqrt(
            abs(p*(p-s1)*(p-s2)*(p-s3))
        )  # Area calculated via Heron's formula
        return a


def get_net(move_set: PathSet, net_type):
    global START, END, CORNERS
    if net_type == "icos_TPF":
        vertex_coords = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [],
                         8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14: [], 15: [],
                         16: [], 17: [], 18: [], 19: [], 20: [], 21: []}
        top_cap_pairs = [(0, 5), (1, 6), (2, 7), (3, 8), (4, 9), 10]
        start_pt = np.array([0, 0, 0])
        CORNERS.append(start_pt - lib.rotate_sixty(move_set.upper_path[0], -1))
        for i in range(len(move_set.upper_path)):
            move = move_set.upper_path[i]
            pt1 = start_pt + move
            pt2 = pt1 + lib.rotate_sixty(lib.rotate_sixty(move, -1), -1)
            vertex_coords[top_cap_pairs[i][0]] = pt1
            vertex_coords[top_cap_pairs[i][1]] = pt2
            start_pt = pt1
        vertex_coords[top_cap_pairs[-1]] = start_pt + lib.rotate_sixty(move_set.upper_path[0], -1)
        CORNERS.append(vertex_coords[top_cap_pairs[-1]] + lib.rotate_sixty(move_set.upper_path[0], 1))

        bottom_cap_pairs = [(17, 11), (18, 12), (19, 13), (20, 14), (21, 15), 16]
        start_pt = vertex_coords[top_cap_pairs[0][1]] - move_set.distance_between - \
            lib.rotate_sixty(move_set.lower_path[0], 1)
        CORNERS.append(start_pt)
        for i in range(len(move_set.lower_path)):
            move = move_set.lower_path[i]
            pt1 = start_pt + move
            pt2 = pt1 + lib.rotate_sixty(lib.rotate_sixty(move, 1), 1)
            vertex_coords[bottom_cap_pairs[i][0]] = pt1
            vertex_coords[bottom_cap_pairs[i][1]] = pt2
            start_pt = pt1
        vertex_coords[bottom_cap_pairs[-1]] = start_pt + lib.rotate_sixty(move_set.lower_path[0], 1)
        CORNERS.append(start_pt + move_set.lower_path[0])
        print(vertex_coords)
        return vertex_coords
    elif net_type == "octa_TPF":
        vertex_coords = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [],
                         8: [], 9: []}
        top_cap_pairs = [(0, 2), (4, 3)]
        start_pt = np.array([0, 0, 0])
        CORNERS.append(start_pt + lib.rotate_sixty(move_set.upper_path[0], 1))
        # vertex 1 set b/c octahedron missing cap triangles
        vertex_coords[1] = start_pt
        for i in range(len(move_set.upper_path)):
            move = move_set.upper_path[i]
            pt1 = start_pt + move + lib.rotate_sixty(move, -1)
            pt2 = start_pt + lib.rotate_sixty(move, -1)
            vertex_coords[top_cap_pairs[i][0]] = pt1
            vertex_coords[top_cap_pairs[i][1]] = pt2
            start_pt = pt1
        CORNERS.append(vertex_coords[top_cap_pairs[-1][0]] + lib.rotate_sixty(move_set.upper_path[-1], 1))

        bottom_cap_pairs = [(9, 6), (8, 7)]
        # vertex 5 set b/c octahedron missing cap triangles
        start_pt = vertex_coords[1] - move_set.distance_between
        vertex_coords[5] = start_pt
        CORNERS.append(start_pt + lib.rotate_sixty(move_set.lower_path[0], -1))
        for i in range(len(move_set.lower_path)):
            move = move_set.lower_path[i]
            pt1 = start_pt + move + lib.rotate_sixty(move, 1)
            pt2 = start_pt + lib.rotate_sixty(move, 1)
            vertex_coords[bottom_cap_pairs[i][0]] = pt1
            vertex_coords[bottom_cap_pairs[i][1]] = pt2
            start_pt = pt1
        CORNERS.append(vertex_coords[bottom_cap_pairs[-1][0]] + lib.rotate_sixty(move, -1))
        print(vertex_coords)
        return vertex_coords
    elif net_type == "octa_TPF_2":
        vertex_coords = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [],
                         8: [], 9: [], 10: [], 11: [], 12: []}
        top_cap_pairs = [(0, 4), (1, 5), (2, 6), (3, 7), 8]
        start_pt = np.array([0, 0, 0])
        CORNERS.append(start_pt - lib.rotate_sixty(move_set.upper_path[0], -1))
        for i in range(len(move_set.upper_path)):
            move = move_set.upper_path[i]
            pt1 = start_pt + move
            pt2 = pt1 + lib.rotate_sixty(lib.rotate_sixty(move, -1), -1)
            vertex_coords[top_cap_pairs[i][0]] = pt1
            vertex_coords[top_cap_pairs[i][1]] = pt2
            start_pt = pt1
        vertex_coords[top_cap_pairs[-1]] = start_pt + lib.rotate_sixty(move_set.upper_path[0], -1)
        CORNERS.append(vertex_coords[top_cap_pairs[-1]] + lib.rotate_sixty(move_set.upper_path[0], 1) +
                       np.array([0, 4, 0]))

        bottom_cap_pairs = [(9, 4, 5), (10, 5, 6), (11, 6, 7), (12, 7, 8)]
        start_pt = vertex_coords[top_cap_pairs[0][1]] - move_set.distance_between - \
                   lib.rotate_sixty(move_set.lower_path[0], 1)
        CORNERS.append(start_pt - np.array([0, 4, 0]))
        for i in range(len(move_set.lower_path)):
            move = move_set.lower_path[i]
            pt1 = start_pt + move
            pt2 = pt1 + lib.rotate_sixty(lib.rotate_sixty(move, 1), 1)
            vertex_coords[bottom_cap_pairs[i][0]] = pt1
            vertex_coords[bottom_cap_pairs[i][1]] = pt2
            start_pt = pt1
        CORNERS.append(start_pt + lib.rotate_sixty(move_set.lower_path[-1], -1))
        print(vertex_coords)
        return vertex_coords
    elif net_type == "tetra_TPF":
        vertex_coords = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
        top_cap_pairs = [(0, 3), (1, 4), (2, 5), 6]
        start_pt = np.array([0, 0, 0])
        CORNERS.append(start_pt - lib.rotate_sixty(move_set.upper_path[0], -1))
        for i in range(len(move_set.upper_path)):
            move = move_set.upper_path[i]
            pt1 = start_pt + move
            pt2 = pt1 + lib.rotate_sixty(lib.rotate_sixty(move, -1), -1)
            vertex_coords[top_cap_pairs[i][0]] = pt1
            vertex_coords[top_cap_pairs[i][1]] = pt2
            start_pt = pt1
        vertex_coords[top_cap_pairs[-1]] = start_pt + lib.rotate_sixty(move_set.upper_path[0], -1)
        CORNERS.append(vertex_coords[top_cap_pairs[-1]] + lib.rotate_sixty(move_set.upper_path[0], 1))

        calcd_vert = calc_vert(move_set.upper_path, "tetra", "case 1")
        CORNERS.append(calcd_vert)
        CORNERS.append(calcd_vert)
        vertex_coords[7] = calcd_vert
        print(vertex_coords)
        return vertex_coords
    elif net_type == "tetra_TPF_120":
        vertex_coords = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
        top_cap_pairs = [(0, 3), (1, 4), (2, 5), 6]
        start_pt = np.array([0, 0, 0])
        CORNERS.append(start_pt - lib.rotate_sixty(move_set.upper_path[0], -1))
        for i in range(len(move_set.upper_path)):
            move = move_set.upper_path[i]
            if i == 1:
                pt1 = start_pt + move + lib.rotate_sixty(move, -1)
                pt2 = start_pt + lib.rotate_sixty(move, -1)
            else:
                pt1 = start_pt + move
                pt2 = pt1 + lib.rotate_sixty(lib.rotate_sixty(move, -1), -1)
            vertex_coords[top_cap_pairs[i][0]] = pt1
            vertex_coords[top_cap_pairs[i][1]] = pt2
            start_pt = pt1
        vertex_coords[top_cap_pairs[-1]] = start_pt + lib.rotate_sixty(move_set.upper_path[0], -1)
        CORNERS.append(vertex_coords[top_cap_pairs[-1]] + lib.rotate_sixty(move_set.upper_path[0], 1))

        calcd_vert = calc_vert(move_set.upper_path, "tetra", "case 1")
        CORNERS.append(calcd_vert)
        CORNERS.append(calcd_vert)
        vertex_coords[7] = calcd_vert
        print(vertex_coords)
        return vertex_coords


def graph_net(vertex_coords, struct_type):
    global START, END
    START, END = get_corners()
    lattice = np.array(make_lattice(start_pt=START, end_pt=END))
    shapes = make_shapes(lattice=lattice)

    # plt.figure(figsize=(8.5, 9))
    plt.figure(figsize=(9, 9)) # Not really a square..?
    plt.gca().set_aspect('equal', adjustable='box')

    x = []
    for vert in list(vertex_coords.keys()):
        x.append(lib.hex_to_cart(vertex_coords[vert]))
    x = np.array(x)

    plt.scatter(x[:, 0], x[:, 1], s=2)
    plt.scatter(lattice[:, 0], lattice[:, 1], s=1, c="black", alpha=0.1)

    # colors = ['red', 'blue', 'red', 'blue']
    colors = ['blue', 'blue', 'blue', 'blue']
    i = 0
    my_color_index = 0

    if struct_type == "icos_TPF":
        CNXN = ICOS_CNXN
        N=5
    elif struct_type == "tetra_TPF":
        CNXN = TETRA_CNXN
        N = 3
    elif struct_type == "octa_TPF":
        CNXN = OCTA_CNXN
        N = 1
    elif struct_type == "octa_TPF_2":
        CNXN = OCTA_CNXN_2
        N = 4

    for triple in CNXN:
        if i >= N:
            my_color_index += 1
            my_color_index %= 2
            i = 0
        coords = plt.Polygon(
            np.array([lib.hex_to_cart(vertex_coords[triple[0]]),
                      lib.hex_to_cart(vertex_coords[triple[1]]),
                      lib.hex_to_cart(vertex_coords[triple[2]])]),
            color=colors[my_color_index],
            alpha=0.5
        )
        i += 1

        plt.gca().add_patch(coords)

    for hex_coords in shapes:
        hexagon = plt.Polygon(hex_coords,
                              color="black",
                              alpha=0.1)

        plt.gca().add_patch(hexagon)

    plt.show()


def get_corners():
    global CORNERS
    # fix corners
    new_corners = []
    for corner in CORNERS:
        new_corner = np.array([0, 0, 0])
        new_corner[0] = corner[0] - corner[2]
        new_corner[1] = corner[1] + corner[2]
        new_corners.append(new_corner)
    CORNERS = new_corners
    # min = np.array([0, 0])
    # upper_right = CORNERS[1]
    h_start = CORNERS[1][0]
    potential_max = np.array([h_start, CORNERS[1][1], CORNERS[1][2]])
    while lib.hex_to_cart(potential_max)[0] < lib.hex_to_cart(CORNERS[3])[0]:
        h_start += 1
        potential_max = np.array([h_start, CORNERS[1][1], CORNERS[1][2]])
    max = np.array([h_start, CORNERS[1][1], CORNERS[1][2]])

    # lower_left = CORNERS[2]
    h_start = CORNERS[2][0]
    potential_min = np.array([h_start, CORNERS[1][1], CORNERS[1][2]])
    while lib.hex_to_cart(potential_min)[0] > lib.hex_to_cart(CORNERS[0])[0]:
        h_start -= 1
        potential_min = np.array([h_start, CORNERS[2][1], CORNERS[2][2]])
    min = np.array([h_start, CORNERS[2][1], CORNERS[2][2]])

    print("min "+str(min)+" max "+str(max))
    print(CORNERS)
    return min, max


def make_lattice(lat_type="hex", start_pt=np.array([0, 0]), end_pt=np.array([0, 0])):
    lattice = []
    start_cart = lib.hex_to_cart(start_pt)
    x0, y0 = start_cart
    end_cart = lib.hex_to_cart(end_pt)
    xf, yf = end_cart
    x, y = x0, y0
    while x <= xf:
        x2 = x
        while y <= yf:
            lattice.append(np.array([x2, y]))
            x2 += lib.hex_to_cart(np.array([0, 1, 1]))[0]
            y += lib.hex_to_cart(np.array([0, 1, 1]))[1]
        x += lib.hex_to_cart(np.array([1, 0, 0]))[0]
        y = y0
    x = x0 + lib.hex_to_cart(np.array([0, 1, 0]))[0]
    y = y0 + lib.hex_to_cart(np.array([0, 1, 0]))[1]
    while x <= xf:
        x2 = x
        while y <= yf:
            # print(np.array([x2, y]))
            lattice.append(np.array([x2, y]))
            x2 += lib.hex_to_cart(np.array([0, 1, 1]))[0]
            y += lib.hex_to_cart(np.array([0, 1, 1]))[1]
        x += lib.hex_to_cart(np.array([1, 0, 0]))[0]
        y = y0 + lib.hex_to_cart(np.array([0, 1, 0]))[1]

    # print(lattice)
    return lattice


def make_shapes(lat_type="hex", lattice=[]):
    shapes = []
    for pt in lattice:
        shapes.append(np.array([
            1/3 * lib.hex_to_cart(np.array([1, 1, 0])) + pt,
            1/3 * lib.hex_to_cart(np.array([0, 1, 1])) + pt,
            1/3 * lib.hex_to_cart(np.array([-1, 0, 1])) + pt,
            1/3 * lib.hex_to_cart(np.array([-1, -1, 0])) + pt,
            1/3 * lib.hex_to_cart(np.array([0, -1, -1])) + pt,
            1/3 * lib.hex_to_cart(np.array([1, 0, -1])) + pt,
        ]))
    return shapes


def make_7_5_path(p_vec, z, n=None, iteration=0, return_area=False):
    area = -1  # placeholder
    if not n:
        n_test = 0
        while 5*n_test - 4*z <= 0:
            n_test += 1
        n = n_test + iteration
    m = 5*n - 4*z
    top_path = [m*p_vec, z*p_vec, z*p_vec, z*p_vec, z*p_vec]
    bottom_path = [n*p_vec, n*p_vec, n*p_vec, n*p_vec, n*p_vec]
    distance = lib.rotate_sixty(n*p_vec, 1)
    print("We have the following 7-5 cone parameters:")
    print("z = "+str(z))
    print("n = " + str(n))
    print("m = " + str(m))
    if return_area:
        area = get_7_5_area(p_vec, z, n, m)
        print("area = " + str(3*area))
    return [top_path, bottom_path, distance, area]


def make_8_4_path(p_vec, z, n=None, iteration=0, return_area=False):
    area = -1  # placeholder
    if not n:
        n_test = 0
        a = (2 * n_test - z) * p_vec[0] + (z - n_test) * p_vec[1]
        b = (n_test - z) * p_vec[0] + (3 * n_test - 2 * z) * p_vec[1]
        while lib.hex_to_cart(lib.rotate_sixty(np.array([a, b, 0]), -1))[0] <= 0:
            n_test += 1
            a = (2 * n_test - z) * p_vec[0] + (z - n_test) * p_vec[1]
            b = (n_test - z) * p_vec[0] + (3 * n_test - 2 * z) * p_vec[1]
        n = n_test + iteration
    a = (2 * n - z) * p_vec[0] + (z - n) * p_vec[1]
    b = (n - z) * p_vec[0] + (3 * n - 2 * z) * p_vec[1]

    first_vec = np.array([a, b, 0])
    sec_vec = np.array([a + b - n*p_vec[1], n*p_vec[0] - a + n*p_vec[1], 0])
    top_path = [first_vec, sec_vec, z*p_vec, z*p_vec, z*p_vec]
    bottom_path = [n*p_vec, n*p_vec, n*p_vec, n*p_vec, n*p_vec]
    distance = lib.rotate_sixty(n*p_vec, 1)
    print("We have the following 8-4 cone parameters:")
    print("z = "+str(z))
    print("n = " + str(n))
    print("first_vec = " + str(first_vec))
    print("sec_vec = " + str(sec_vec))
    # if return_area:
    #     area = get_7_5_area(p_vec, z, n, m)
    #     print("area = " + str(3*area))
    return [top_path, bottom_path, distance, area]


def make_9_3_path(p_vec, z, n=None, iteration=0, return_area=False):
    area = -1  # placeholder
    if not n:
        n = iteration + 1
    n2 = 2*z-n
    a = z*p_vec[0] + (-z+n)*p_vec[1]
    b = (z-n)*p_vec[0] + (2*z-n)*p_vec[1]

    first_vec = np.array([a, b, 0])
    sec_vec = np.array([a + b - n*p_vec[1], -a + n*(p_vec[0] + p_vec[1]), 0])
    top_path = [z*p_vec, first_vec, sec_vec, z*p_vec, z*p_vec]
    bottom_path = [z*p_vec, z*p_vec, n2*p_vec, z*p_vec, z*p_vec]
    distance = lib.rotate_sixty(z*p_vec, 1)
    print("We have the following 9-3 cone parameters:")
    print("z = "+str(z))
    print("n = " + str(n))
    print("first_vec = " + str(first_vec))
    print("sec_vec = " + str(sec_vec))
    # if return_area:
    #     area = get_7_5_area(p_vec, z, n, m)
    #     print("area = " + str(3*area))
    return [top_path, bottom_path, distance, area]


def make_10_2_path(p_vec, z, n=None, iteration=0, return_area=False):
    area = -1  # placeholder
    if not n:
        n = iteration + 1
    a = (2*z-n)*p_vec[0] + (2*n-2*z)*p_vec[1]
    b = -2*a+2*z*p_vec[0] + n*p_vec[1]

    first_vec = np.array([a, b, 0])
    sec_vec = np.array([a + b - n * p_vec[1], a - n * p_vec[0] + n * p_vec[1], 0])
    top_path = [z * p_vec, first_vec, sec_vec, z * p_vec, z * p_vec]
    bottom_path = [z * p_vec, z * p_vec, z * p_vec, z * p_vec, z * p_vec]
    distance = lib.rotate_sixty(z * p_vec, 1)
    print("We have the following 10-2 cone parameters:")
    print("z = " + str(z))
    print("n = " + str(n))
    print("first_vec = " + str(first_vec))
    print("sec_vec = " + str(sec_vec))
    # if return_area:
    #     area = get_7_5_area(p_vec, z, n, m)
    #     print("area = " + str(3*area))
    return [top_path, bottom_path, distance, area]


def make_elongated_path(elong_type, t_vec, q_vec):
    if elong_type == "5-fold":
        top_path = [t_vec, t_vec, t_vec, t_vec, t_vec]
        bottom_path = top_path
        distance = q_vec
        return [top_path, bottom_path, distance, -1]

    if elong_type == "3-fold":
        qr_vec = lib.rotate_sixty(q_vec, -1)  # q_vec rotated -pi/3 degrees for path
        top_path = [t_vec, t_vec, qr_vec, t_vec, t_vec]
        bottom_path = [qr_vec, t_vec, t_vec, t_vec, t_vec]
        distance = lib.rotate_sixty(t_vec, 1)
        return [top_path, bottom_path, distance, -1]
    if elong_type == "2-fold":
        qr_vec = lib.rotate_sixty(q_vec, -1)  # q_vec rotated -pi/3 degrees for path
        top_path = [qr_vec, t_vec, t_vec, t_vec, t_vec]
        bottom_path = [t_vec, t_vec, t_vec, t_vec, qr_vec]
        distance = lib.rotate_sixty(t_vec, 1)
        return [top_path, bottom_path, distance, -1]


def make_tetra(top_path, option, plot="no"):
    if option == "120":
        c1 = top_path[0] + lib.rotate_sixty(top_path[0], -1)
    else:
        c1 = top_path[0]
    c2 = top_path[1]
    # Option 1
    c3_1 = np.array([-c1[1]+c2[0]+c2[1],
                  c1[0]+c1[1]-c2[0],
                     0])
    c3_2 = np.array([(-c1[0]-2*c1[1]+2*c2[0]+c2[1])/3,
                     (2*c1[0]+c1[1]-4*c2[0]-2*c2[1])/3,
                     0])
    print("Case 1:")
    print(c3_1)
    print("Case 2:")
    print(c3_2)

    top_path.append(c3_1)
    move_set = PathSet()
    move_set.upper_path = top_path
    if option == "120":
        vertex_coords = get_net(move_set, "tetra_TPF_120")
    else:
        vertex_coords = get_net(move_set, "tetra_TPF")
    if plot == "yes":
        # convert to cart
        x = []
        for vert in list(vertex_coords.keys()):
            vcart = lib.hex_to_cart(vertex_coords[vert])
            x.append(np.array([vcart[0], vcart[1], 0]))
        dat = {'v01': x[1], 'v02': x[4], 'v03': x[5],
               'r1': np.linalg.norm(x[0] - x[3]),
               'r2': np.linalg.norm(x[3] - x[4]),
               'r3': np.linalg.norm(x[5] - x[6])}
        soln_set = [solv.tetrahedron_solver_new(dat), dat['v01'],
                    dat['v02'], dat['v03']]
        face_set = np.array([list(soln_set[0].astype(float)), list(soln_set[1].astype(float)),
                            list(soln_set[2].astype(float)), list(soln_set[3].astype(float))])
        print(face_set)
        plot_my_pts(face_set, 'b')
    graph_net(vertex_coords, "tetra_TPF")


def make_octa_2(top_path, bottom_vec):
    bottom_path = [bottom_vec,
                   np.array([0, 0, 0]),
                   np.array([0, 0, 0]),
                   np.array([0, 0, 0])]
    for i in range(1, len(bottom_path)):
        bottom_path[i] = top_path[i - 1] + lib.rotate_sixty(lib.rotate_sixty(top_path[i], -1), -1) - \
                         lib.rotate_sixty(lib.rotate_sixty(bottom_path[i - 1], -1), -1)
    move_set = PathSet()
    move_set.upper_path = top_path
    move_set.distance_between = np.array([0, 0, 0])
    move_set.lower_path = bottom_path
    move_set.connection = OCTA_CNXN_2
    return move_set


def calc_vert(upper_path, struc_type, case):
    if struc_type == "tetra" and case == "case 1":
        d1 = lib.rotate_sixty(upper_path[0], 1) + \
             lib.rotate_sixty(upper_path[1], -1)
        d1_rot = lib.rotate_sixty(d1, -1)
        return lib.rotate_sixty(upper_path[0], -1) + d1 + d1_rot


def get_7_5_area(p_vec, z, n, m):  # Assumes j0 = 0
    h0 = p_vec[0]
    k0 = p_vec[1]
    m_vec = m*p_vec
    n_vec = n*p_vec
    height_vec = lib.hex_to_cart(lib.rotate_sixty(m_vec, 1) + 2*lib.rotate_sixty(n_vec, 1))
    length_vec = lib.hex_to_cart(5*n_vec)
    para_area = abs(height_vec[0]*length_vec[1]-height_vec[1]*length_vec[0]) * 4/math.sqrt(3)
    top_triangles = 4*((z*h0)**2 + (z*h0)*(n*k0) + (z*k0)**2) + (m*h0)**2 + (m*h0)*(n*k0) + (m*k0)**2
    bottom_triangles = 5*((n*h0)**2 + (n*h0)*(n*k0) + (n*k0)**2)
    return round(para_area - top_triangles - bottom_triangles, 4)


def build_path_dataset(p_vec, it_limit, cone_type):
    area_list = []
    for z in range(1, it_limit):
        for iteration in range(it_limit):
            move_set = PathSet()
            if cone_type == "7-5":
                top_path, bottom_path, distance, area = make_7_5_path(p_vec=p_vec, z=z, iteration=iteration)
            if cone_type == "8-4":
                top_path, bottom_path, distance, area = make_8_4_path(p_vec=p_vec, z=z, iteration=iteration)
            if cone_type == "9-3":
                top_path, bottom_path, distance, area = make_9_3_path(p_vec=p_vec, z=z, iteration=iteration)
            if cone_type == "10-2":
                top_path, bottom_path, distance, area = make_10_2_path(p_vec=p_vec, z=z, iteration=iteration)
            move_set.upper_path = top_path
            move_set.lower_path = bottom_path
            move_set.distance_between = distance
            move_set.connection = ICOS_CNXN

            my_net = get_net(move_set, "icos_TPF")
            move_set.assign_triangles(my_net)
            move_set.assign_num_subunits()

            area_list.append([z, iteration, move_set.num_subunits])
    df = pd.DataFrame(np.array(area_list), columns=['z', 'iteration', 'area'])
    df.to_excel("output" + str(cone_type) + ".xlsx", sheet_name='sheet_1')


def plot_my_pts(face_set, color):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    hull = ConvexHull(face_set)
    # draw the polygons of the convex hull
    color = ['b', 'r', 'r', 'r']
    c = 0
    for s in hull.simplices:
        tri = Poly3DCollection(face_set[s])
        tri.set_color(color[c])
        tri.set_edgecolor('k')
        tri.set_alpha(0.6)
        ax.add_collection3d(tri)
        c += 1
    # draw the vertices
    ax.scatter(face_set[:, 0], face_set[:, 1], face_set[:, 2], marker='o', color='purple')
    set_axes_equal(ax)
    plt.axis('off')
    plt.grid(b=None)
    # plt.show()


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

# Getting my move_set object set up with paths and a distance
move_set = PathSet()
# path_1 = [np.array([5, 5, 0]),
#           np.array([10, 10, 0]),
#           np.array([10, 10, 0]),
#           np.array([10, 10, 0]),
#           np.array([10, 10, 0])]
# path_2 = [np.array([9, 9, 0]),
#           np.array([9, 9, 0]),
#           np.array([9, 9, 0]),
#           np.array([9, 9, 0]),
#           np.array([9, 9, 0])]
# path_2 = path_1
# distance = np.array([0, 9, 9])

# top_path, bottom_path, distance, area = make_7_5_path(p_vec=np.array([1, 1, 0]), z=1, iteration=1, return_area=True)
# top_path, bottom_path, distance, area = make_8_4_path(p_vec=np.array([1, 0, 0]), z=4, n=3)
# top_path, bottom_path, distance, area = make_8_4_path(p_vec=np.array([1, 0, 0]), z=1, iteration=1)
# top_path, bottom_path, distance, area = make_9_3_path(p_vec=np.array([1, 1, 0]), z=2, n=1)
# z=2 n=3 example

# top_path, bottom_path, distance, area = make_10_2_path(p_vec=np.array([1, 0, 0]), z=3, n=1)
#
# top_path, bottom_path, distance, area = make_elongated_path(elong_type="2-fold", t_vec=np.array([1, 1, 0]),
#                                                              q_vec=np.array([0, 1, 2]))

# build_path_dataset(p_vec=np.array([2, 1, 0]), it_limit=20, cone_type="8-4")
# #
# top_path = [np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 2, 0]),
#           np.array([1, 1, 0])]
# bottom_path = [np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 2, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0])]
# distance = np.array([0, 1, 2])
# #
# #
# top_path = [np.array([1, 1, 0]),
#           np.array([1, 2, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0])]
# bottom_path = [np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 2, 0])]
# distance = np.array([-2, 2, 0])
#
#
#
# top_path = [np.array([1, 2, 0]),
#           np.array([1, 3, 0]),
#           np.array([1, 1, 0]),
#           np.array([2, 2, 0]),
#           np.array([3, 3, 0])]
# bottom_path = [np.array([1, 2, 0]),
#           np.array([3, 3, 0]),
#           np.array([2, 2, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 3, 0])]
# distance = np.array([0, 3, 3])
# # #
# top_path = [np.array([0, 3, 0]),
#           np.array([1, 3, 0]),
#           np.array([0, 6, 1]),
#           np.array([4, 5, 0]),
#           np.array([0, 3, 1])]
# bottom_path = [np.array([0, 10, 0]),
#           np.array([0, 3, 0]),
#           np.array([0, 2, 0]),
#           np.array([3, 2, 0]),
#           np.array([0, 5, 0])]
# distance = np.array([-2, 0, 1])
# #
# #
# top_path = [np.array([1, 2, 0]),
#           np.array([2, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 2, 0])]
# top_path = [np.array([1, 2, 0]), # 3-fold example
#           np.array([2, 1, 0]),
#           np.array([2, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0])]
# bottom_path = [np.array([1, 2, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([2, 1, 0]),
#           np.array([2, 1, 0])]
# distance = np.array([0, 2, 1])

# top_path = [np.array([1, 2, 0]),
#             np.array([1, 3, 0]),
#           np.array([2, 2, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 4, 0])]
# bottom_path = [np.array([1, 1, 0]),
#                np.array([3, 3, 0]),
#           np.array([1, 2, 0]),
#           np.array([1, 3, 0]),
#           np.array([0, 3, 0])]
# distance = np.array([0, 4, 4])
#
# top_path = [np.array([1, 2, 0]),
#             np.array([3, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0])]
# bottom_path = [np.array([2, 2, 0]),
#                np.array([3, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0])]
# distance = np.array([0, 0, 0])

# top_path = [np.array([2, 1, 0]),
#           np.array([2, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([3, 1, 0]),
#           np.array([3, 2, 0])]
# bottom_path = [np.array([1, 1, 0]),
#           np.array([4, 1, 0]),
#           np.array([3, 1, 0]),
#           np.array([2, 1, 0]),
#           np.array([1, 2, 0])]
# distance = np.array([0, 4, 4])
# #
# #
# # ######################################### ICOS THINGS
# top_path = [np.array([1, 2, 0]),
#           np.array([1, 3, 0]),
#           np.array([2, 3, 0]),
#           np.array([2, 2, 0]),
#           np.array([2, 1, 0])]
# bottom_path = [np.array([1, 2, 0]),
#                np.array([2, 1, 0]),
#           np.array([2, 2, 0]),
#           np.array([2, 3, 0]),
#           np.array([1, 3, 0])]
# #

# 2-fold elongated
top_path = [np.array([1, 1, 0]),
          np.array([0, 2, 0]),
          np.array([1, 2, 0]),
          np.array([2, 1, 0]),
          np.array([1, 1, 0])]
bottom_path = [np.array([1, 1, 0]),
               np.array([1, 1, 0]),
          np.array([2, 1, 0]),
          np.array([1, 2, 0]),
          np.array([0, 2, 0])]

distance = np.array([0, 1, 1])
#
# # 3-fold elongated case
# top_path = [np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 2, 0]),
#           np.array([1, 1, 0])]
# bottom_path = [np.array([1, 1, 0]),
#                np.array([1, 2, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0])]
#
# distance = np.array([0, 1, 2])

# 7-5 reflective case, no equilateral requirement
# top_path = [np.array([1, 1, 0]),
#           np.array([1, 2, 0]),
#           np.array([1, 2, 0]),
#           np.array([2, 1, 0]),
#           np.array([2, 1, 0])]
# bottom_path = [np.array([1, 1, 0]),
#                np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([3, 3, 0]),
#           np.array([1, 1, 0])]
#
# distance = np.array([0, 1, 10])

# 7-5 reflective case
top_path = [np.array([1, 1, 0]),
          np.array([1, 1, 0]),
          np.array([6, 6, 0]),
          np.array([1, 1, 0]),
          np.array([1, 1, 0])]
bottom_path = [np.array([2, 2, 0]),
               np.array([2, 2, 0]),
          np.array([2, 2, 0]),
          np.array([2, 2, 0]),
          np.array([2, 2, 0])]

distance = np.array([0, 1, 10])

# # 7-5 reflective case shifted
# top_path = [np.array([6, 6, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0])]
# bottom_path = [np.array([2, 2, 0]),
#                np.array([2, 2, 0]),
#           np.array([2, 2, 0]),
#           np.array([2, 2, 0]),
#           np.array([2, 2, 0])]
#
# distance = np.array([0, 2, 2])




# 8-4
# top_path = [np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([2, 5, 0]),
#           np.array([5, 2, 0]),
#           np.array([1, 1, 0])]
# bottom_path = [np.array([2, 2, 0]),
#                np.array([2, 2, 0]),
#           np.array([2, 2, 0]),
#           np.array([2, 2, 0]),
#           np.array([2, 2, 0])]
#
# distance = np.array([0, 2, 5])

# 9-3
# top_path = [np.array([1, 4, 0]),
#           np.array([3, 3, 0]),
#           np.array([4, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0])]
# bottom_path = [np.array([2, 2, 0]),
#                np.array([2, 2, 0]),
#           np.array([2, 2, 0]),
#           np.array([2, 2, 0]),
#           np.array([2, 2, 0])]

# top_path = [np.array([3, 3, 0]),
#           np.array([4, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 4, 0])]
# bottom_path = [np.array([2, 2, 0]),
#                np.array([2, 2, 0]),
#           np.array([2, 2, 0]),
#           np.array([2, 2, 0]),
#           np.array([2, 2, 0])]
#
# distance = np.array([0, 2, 2])

# 10-2

# top_path = [np.array([1, 3, 0]),
#           np.array([2, 3, 0]),
#           np.array([3, 2, 0]),
#           np.array([3, 1, 0]),
#           np.array([1, 1, 0])]
# bottom_path = [np.array([2, 2, 0]),
#                np.array([2, 2, 0]),
#           np.array([2, 2, 0]),
#           np.array([2, 2, 0]),
#           np.array([2, 2, 0])]
#
# distance = np.array([0, 2, 2])

# top_path = [np.array([1, 1, 0]),
#           np.array([1, 3, 0]),
#           np.array([2, 3, 0]),
#           np.array([3, 2, 0]),
#           np.array([3, 1, 0])]
# bottom_path = [np.array([2, 2, 0]),
#                np.array([2, 2, 0]),
#           np.array([2, 2, 0]),
#           np.array([2, 2, 0]),
#           np.array([2, 2, 0])]
#
# distance = np.array([0, 2, 3])

# reflective symmetry counter ex (hypoth)
# top_path = [np.array([3, 3, 0]),
#           np.array([2, 1, 0]),
#           np.array([3, 1, 0]),
#           np.array([1, 3, 0]),
#           np.array([1, 2, 0])]
# bottom_path = [np.array([3, 3, 0]),
#           np.array([2, 1, 0]),
#           np.array([1, 3, 0]),
#           np.array([3, 1, 0]),
#           np.array([1, 2, 0])]
#
# distance = np.array([0, 0, 6])

# # reflective symmetry ex
# top_path = [np.array([3, 3, 0]),
#           np.array([2, 1, 0]),
#           np.array([3, 1, 0]),
#           np.array([1, 3, 0]),
#           np.array([1, 2, 0])]
# bottom_path = [np.array([3, 1, 0]),
#           np.array([1, 3, 0]),
#           np.array([1, 2, 0]),
#           np.array([3, 3, 0]),
#           np.array([2, 1, 0])]
#
# distance = np.array([0, 1, 6])

#
# # 5-fold
# top_path = [np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0])]
# bottom_path = [np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0])]
#
# distance = np.array([0, 1, 3])

# C=((0,5,0),(0,5,0),(3,2,0),(0,2,0),(3,7,0),(0,5,0))
# d=(-4,2,0)
# D=((3,5,0),(0,9,0),(1,2,0),(1,2,0),(1,3,0),(3,5,0))

# # HIV
top_path = [np.array([0, 5, 0]),
          np.array([0, 5, 0]),
          np.array([3, 2, 0]),
          np.array([0, 2, 0]),
          np.array([3, 7, 0])]
bottom_path = [np.array([3, 5, 0]),
          np.array([0, 9, 0]),
          np.array([1, 2, 0]),
          np.array([1, 2, 0]),
          np.array([1, 3, 0])]
distance = np.array([-4, 2, 0])

# mattei 1
# top_path = [np.array([0, 9, 5]),
#           np.array([0, 10, 2]),
#           np.array([0, 2, 0]),
#           np.array([9, 4, 0]),
#           np.array([0, 0, 5])]
# bottom_path = [np.array([0, 11, 0]),
#           np.array([0, 2, 10]),
#           np.array([0, 2, 1]),
#           np.array([0, 4, 1]),
#           np.array([0, 3, 1])]
# distance = np.array([0, 0, 2])

move_set.upper_path = top_path
move_set.lower_path = bottom_path
move_set.distance_between = distance
move_set.connection = ICOS_CNXN



# # #
# # # print("top_path is " + str(top_path))
# # # print("bottom_path is " + str(bottom_path))
# # # #
# # # # # Grabbing the net coordinates
my_net = get_net(move_set, "icos_TPF")
move_set.assign_triangles(my_net)
move_set.assign_num_subunits()
print("number of subunits (via triangle sum) " + str(move_set.num_subunits))
# # # Graphing the coordinates
graph_net(my_net, "icos_TPF")
#########################################

########################################
# Tetra
#########################################
# top_path = [np.array([1, 0, 0]),
#           np.array([1, 3, 0])]
# make_tetra(top_path, option="n")

# top_path = [np.array([4, 0, 0]),
#          np.array([4, 0, 0])]
# make_tetra(top_path, option="n", plot="yes")

#########################################
# Octa
# ########################################
# top_path = [np.array([0, 2, 0]),
#           np.array([3, 0, 0])]
# bottom_path = [np.array([3, -3, 0]),
#           np.array([2, 0, 0])]
# distance = np.array([-2, 4, 0])
#
# move_set.upper_path = top_path
# move_set.lower_path = bottom_path
# move_set.distance_between = distance
# move_set.connection = OCTA_CNXN
#
# my_net = get_net(move_set, "octa_TPF")
# move_set.assign_triangles(my_net)
# move_set.assign_num_subunits()
# print("number of subunits (via triangle sum) " + str(move_set.num_subunits))
# # Graphing the coordinates
# graph_net(my_net, "octa_TPF")

#########################################
# Octa 2
#########################################
c1 = np.array([1, 1, 0])
c2 = np.array([2, 2, 0])
c3 = np.array([3, 1, 0])
c4 = np.array([2, 1, 0])
#
# # c1 = np.array([2, 0, 0])
# # c2 = np.array([3, 0, 0])
# # c3 = np.array([4, 0, 0])
# # c4 = np.array([2, 0, -2])
# # # d1 = np.array([1, 5, 0])
df = c4 + lib.rotate_sixty(c4, -1) - c1 - lib.rotate_sixty(c1, 1) + c3 + \
     lib.rotate_sixty(c3, 1) + lib.rotate_sixty(c2, 1) - lib.rotate_sixty(c2, -1)
df = np.array([df[0] + df[2], df[1] - df[2], 0])
d1 = np.array([df[1], df[0] - df[1], 0])
top_path = [c1, c2, c3, c4
            ]
# d1 = np.array([1, 1, 0])
bottom_vec = d1
#
move_set = make_octa_2(top_path, bottom_vec)
#
# f1 = lib.rotate_sixty(c1, 1) + lib.rotate_sixty(c2, -1)
# f2 = lib.rotate_sixty(c2, 1) + lib.rotate_sixty(c3, -1)
# f3 = lib.rotate_sixty(c3, 1) + lib.rotate_sixty(c4, -1)
# f4 = lib.rotate_sixty(c4, 1) + lib.rotate_sixty(c1, -1)
#
# top_path = [f1, c1, f4, d1]
# bottom_vec = f2
# move_set = make_octa_2(top_path, bottom_vec)
#
# print(move_set.lower_path)
#
my_net = get_net(move_set, "octa_TPF_2")
move_set.assign_triangles(my_net)
move_set.assign_num_subunits()
# print("number of subunits (via triangle sum) " + str(move_set.num_subunits))
# # # Graphing the coordinates
# graph_net(my_net, "octa_TPF_2")