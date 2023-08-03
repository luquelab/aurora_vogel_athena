import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import fsolve
from scipy.optimize import minimize
from scipy.optimize import basinhopping


ICOS_VERTS = set(range(0, 22))
VERTS_IN_USE = list()
VERTS_DICT = {}

ICOS_NET_TEMPLATE = [{"coordination": (11, 5, 12),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (5, 12, 6),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (12, 6, 13),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (6, 13, 7),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (13, 7, 14),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (7, 14, 8),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (14, 8, 15),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (8, 15, 9),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (15, 9, 16),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (9, 16, 10),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (5, 0, 6),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (6, 1, 7),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (7, 2, 8),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (8, 3, 9),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (9, 4, 10),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (11, 17, 12),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (12, 18, 13),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (13, 19, 14),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (14, 20, 15),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (15, 21, 16),
                      "distance": np.array([0, 0, 0])},
                     ]


def graph_vec(x):
    vert_list = []
    for i in range(int(len(x)/2)):
        vert_list.append(x[2*i:2*i+2])

    # placing all verts together to form a dataframe which can be plotted
    DATA = pd.DataFrame(vert_list, columns=["X", "Y"])

    fig, ax = plt.subplots()
    ax.scatter(DATA["X"], DATA["Y"])

    # labeling the verts plotted
    n = []
    for i in range(22):
        n.append(i)

    for i, txt in enumerate(n):
        ax.annotate(txt, (DATA["X"][i], DATA["Y"][i]))

    plt.show()


def dist(v1, v2):
    x_diff = (v2[0] - v1[0]) + 1 / 2 * (v2[1] - v1[1] + (-v2[2]) - (-v1[2]))
    y_diff = math.sqrt(3) / 2 * (v2[1] - v1[1] + v2[2] - v1[2])
    return math.sqrt(x_diff**2 + y_diff**2)


def hex_to_cart(v1):
    x = v1[0] + v1[1] * 1/2 + v1[2] * (-1/2)
    y = v1[1] * math.sqrt(3)/2 + v1[2] * math.sqrt(3)/2
    return np.array([x, y])


def rotate_sixty(tup, sign):
    h = tup[0]
    k = tup[1]
    j = tup[2]
    if sign == 1:
        new_tup = np.array([-j, h, k])
    else:
        new_tup = np.array([k, j, -h])
    return new_tup


def closeness_prevention(vertices):
    closeness_penalty = 0
    num_verts = int(len(vertices)/2)
    for vertex_a in range(num_verts - 1):
        for vertex_b in range(vertex_a + 1, num_verts):
            d = np.linalg.norm(vertices[2*vertex_a:2*vertex_a+2]-vertices[2*vertex_b:2*vertex_b+2])
            if d < 10**-1:
                closeness_penalty += 10**-1-d
    return closeness_penalty


def func(vertices, net, iterate, it_type):
    distance_sum = 0
    num_new_verts = int(len(vertices)/2) - 3
    for i in range(min(len(net), num_new_verts + 1)):
        shape_size = len(net[i]["coordination"])
        for j in range(len(net[i]["distance"])):
            j0 = net[i]["coordination"][j]
            jv0 = VERTS_DICT[j0]
            j1 = net[i]["coordination"][(j+1) % shape_size]
            jv1 = VERTS_DICT[j1]
            distance_sum += (net[i]["distance"][j] -
                             np.linalg.norm(vertices[jv0:jv0+2] - vertices[jv1:jv1+2]))**2
    if it_type == "basinhopping":
        distance_sum = distance_sum + closeness_prevention(vertices)
        print("closeness penalty")
        print(closeness_prevention(vertices))

        return distance_sum


def make_net(net_type, vector_list, vert_number):
    if net_type == "icosahedral":
        net = ICOS_NET_TEMPLATE
        cT = vector_list[0]
        for i in range(len(ICOS_NET_TEMPLATE)):
            d = dist(np.array([0, 0, 0]), cT)
            net[i]["distance"] = np.array([d, d, d])
    return net


def update_verts_in_use(net, verts_in_use, k):
    verts_in_use_set = set(verts_in_use)
    coordinates = net[k]["coordination"]
    for coord in coordinates:
        if coord not in verts_in_use_set:
            verts_in_use.append(coord)
    return verts_in_use


def assign_func_vertices(vertices):
    global VERTS_DICT

    assigned_vertices = np.zeros(len(VERTS_IN_USE)*2)
    for i in range(len(VERTS_IN_USE)):
        vertex = VERTS_IN_USE[i]
        assigned_vertices[2*i:2*i+2] = vertices[2*vertex:2*vertex+2]
        VERTS_DICT[vertex] = 2*i
    return assigned_vertices


def get_vertices_from_func_vertices(func_vertices):
    vertices = np.zeros(2*22)
    for i in range(len(VERTS_IN_USE)):
        vertex = VERTS_IN_USE[i]
        vertices[2*vertex:2*vertex+2] = func_vertices[2*i:2*i+2]
    return vertices


def educated_guess(func_vertices, net):
    num_new_verts = int(len(func_vertices)/2)-3
    newest_vertex = VERTS_IN_USE[-1]
    adjacent_vertex_sum = 0
    for vertex in net[num_new_verts]["coordination"]:
        if vertex != newest_vertex:
            index = VERTS_DICT[vertex]
            adjacent_vertex_sum += func_vertices[index:index+2]
    func_vertices[-2:] = adjacent_vertex_sum / 2

    return func_vertices


def iterator(vertices, net, it_type):
    global VERTS_IN_USE

    for k in range(1, int(len(vertices)/2)):
        VERTS_IN_USE = update_verts_in_use(net, VERTS_IN_USE, k)
        func_verts = assign_func_vertices(vertices)
        if it_type == "basinhopping":
            func_verts = educated_guess(func_verts, net)
            minimizer_kwargs = {"method": "BFGS"}
            new_func_vertices = basinhopping(lambda x: func(x, net, k, it_type), x0=func_verts,
                                             minimizer_kwargs=minimizer_kwargs, niter=10)
            print(new_func_vertices)
            new_func_vertices = new_func_vertices.x
            vertices = get_vertices_from_func_vertices(new_func_vertices)

        graph_vec(vertices)


cT = np.array([2, 3, 0])
my_net = make_net("icosahedral", [cT], 20)
verts = np.zeros(22*2)
verts[2*5:2*(5+1)] = hex_to_cart(rotate_sixty(cT, 1))[0], hex_to_cart(rotate_sixty(cT, 1))[1]
verts[2*11:2*(11+1)] = np.array(0), np.array(0)
verts[2*12:2*(12+1)] = hex_to_cart(cT)
VERTS_IN_USE = [11, 5, 12]
iterator(verts, my_net, "basinhopping")
