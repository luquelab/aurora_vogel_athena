import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import fsolve
from scipy.optimize import minimize

ICOS_NET_TEMPLATE = [{"coordination": (5, 0, 6),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (6, 1, 7),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (7, 2, 8),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (8, 3, 9),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (9, 4, 10),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (5, 6, 12),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (6, 7, 13),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (7, 8, 14),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (8, 9, 15),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (9, 10, 16),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (5, 11, 12),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (6, 12, 13),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (7, 13, 14),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (8, 14, 15),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (9, 15, 16),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (11, 12, 17),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (12, 13, 18),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (13, 14, 19),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (14, 15, 20),
                      "distance": np.array([0, 0, 0])},
                     {"coordination": (15, 16, 21),
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


def func(vertices, input_list, iterate, it_type):
    distance_sum = 0
    for i in range(min(len(input_list), iterate + 1)):
        for j in range(len(input_list[i]["distance"])):
            j0 = input_list[i]["coordination"][j]
            j1 = input_list[i]["coordination"][(j+1) % 3]
            distance_sum += (input_list[i]["distance"][j] - np.linalg.norm(vertices[2*j0:2*j0+2] -
                                                                              vertices[2*j1:2*j1+2]))**2
    distance_sum_vec = np.zeros(len(vertices))
    distance_sum_vec[0] = distance_sum
    if it_type == "fsolve":
        return distance_sum_vec
    if it_type == "minimize":
        return distance_sum_vec[0]


def func_wrapper(vertices, input_list, iterate):
    pass


def make_net(net_type, vector_list, vert_number):
    if net_type == "icosahedral":
        net = ICOS_NET_TEMPLATE
        cT = vector_list[0]
        for i in range(len(ICOS_NET_TEMPLATE)):
            d = dist(np.array([0, 0, 0]), cT)
            net[i]["distance"] = np.array([d, d, d])
    return net


def iterator(vertices, net, it_type):
    for k in range(1, len(vertices)):
        if it_type == "fsolve":
            new_vertices = fsolve(lambda x: func(x, net, k, it_type), x0=vertices, xtol=10**-6, full_output=True)
        if it_type == "minimize":
            new_vertices = minimize(lambda x: func(x, net, k, it_type), x0=vertices, tol=10 ** -6, method='Nelder-Mead').x
        graph_vec(new_vertices)
        vertices = new_vertices


cT = np.array([2, 3, 0])
my_net = make_net("icosahedral", [cT], 20)
verts = np.zeros(22*2)
verts[0:2] = hex_to_cart(rotate_sixty(cT, 1))[0], hex_to_cart(rotate_sixty(cT, 1))[1]
verts[2*5:2*(5+1)] = np.array(0), np.array(0)
verts[2*6:2*(6+1)] = hex_to_cart(cT)
iterator(verts, my_net, "minimize")
