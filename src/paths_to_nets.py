import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import my_library as lib
import openpyxl

ICOS_CNXN = [[5, 0, 6], [6, 1, 7], [7, 2, 8], [8, 3, 9], [9, 4, 10],
             [12, 5, 6], [13, 6, 7], [14, 7, 8], [15, 8, 9], [16, 9, 10],
             [11, 5, 12], [12, 6, 13], [13, 7, 14], [14, 8, 15], [15, 9, 16],
             [17, 11, 12], [18, 12, 13], [19, 13, 14], [20, 14, 15], [21, 15, 16]]
START = np.array([0, 0, 0])
END = np.array([0, 0, 0])
CORNERS = []


class PathSet:
    def __init__(self, upper_path=[], lower_path=[], distance_between=0):
        self.upper_path = upper_path
        self.lower_path = lower_path
        self.distance_between = distance_between


def get_net(move_set: PathSet):
    global START, END, CORNERS
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


def graph_net(vertex_coords):
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

    colors = ['red', 'blue']
    i = 0
    my_color_index = 0
    for triple in ICOS_CNXN:
        if i >= 5:
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
    min = np.array([0, 0])
    upper_right = CORNERS[1]
    h_start = CORNERS[1][0]
    potential_max = np.array([h_start, CORNERS[1][1], CORNERS[1][2]])
    while lib.hex_to_cart(potential_max)[0] < lib.hex_to_cart(CORNERS[3])[0]:
        h_start += 1
        potential_max = np.array([h_start, CORNERS[1][1], CORNERS[1][2]])
    max = np.array([h_start, CORNERS[1][1], CORNERS[1][2]])

    lower_left = CORNERS[2]
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
    area = -1 # placeholder
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
        print("area = " + str(area))
    return [top_path, bottom_path, distance, area]


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


def build_path_dataset(p_vec, limit):
    area_list = []
    for z in range(1, limit):
        for iteration in range(limit):
            area = make_7_5_path(p_vec=p_vec, z=z, iteration=iteration, return_area=True)[-1]
            area_list.append([z, iteration, area])
    df = pd.DataFrame(np.array(area_list), columns=['z', 'iteration', 'area'])
    df.to_excel("output.xlsx", sheet_name='sheet_1')



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

# top_path, bottom_path, distance = make_7_5_path(p_vec=np.array([1, 1, 0]), z=7, n=7)
top_path, bottom_path, distance, area = make_7_5_path(p_vec=np.array([1, 0, 0]), z=2, iteration=1, return_area=True)

build_path_dataset(p_vec=np.array([2, 1, 0]), limit=20)

top_path = [np.array([8, 14, 0]),
          np.array([14, 8, 0]),
          np.array([2, 2, 0]),
          np.array([2, 2, 0]),
          np.array([2, 2, 0])]
bottom_path = [np.array([8, 8, 0]),
          np.array([8, 8, 0]),
          np.array([8, 8, 0]),
          np.array([8, 8, 0]),
          np.array([8, 8, 0])]
distance = np.array([0, 8, 8])


move_set.upper_path = top_path
move_set.lower_path = bottom_path
move_set.distance_between = distance

# Grabbing the net coordinates
my_net = get_net(move_set)

# Graphing the coordinates
graph_net(my_net)
