import numpy as np
import math


def dist(v1, v2=np.array([0, 0, 0])):
    x_diff = (v2[0] - v1[0]) + 1 / 2 * (v2[1] - v1[1] + (-v2[2]) - (-v1[2]))
    y_diff = math.sqrt(3) / 2 * (v2[1] - v1[1] + v2[2] - v1[2])
    return math.sqrt(x_diff**2 + y_diff**2)


def hex_to_cart(v1):
    x = v1[0] + v1[1] * 1/2 + v1[2] * (-1/2)
    y = v1[1] * math.sqrt(3)/2 + v1[2] * math.sqrt(3)/2
    return np.array([x, y])


def cart_to_hex(X):
    x, y = X
    v1 = np.array([x-y/2*2/math.sqrt(3), y*2/math.sqrt(3), 0])
    return v1


def find_angle(vec1, vec2):
    cos_angle = np.dot(vec1, vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))
    return np.arccos(cos_angle)


def rotate_sixty(tup, sign):
    h = tup[0]
    k = tup[1]
    j = tup[2]
    if sign == 1:
        new_tup = np.array([-j, h, k])
    if sign == -1:
        new_tup = np.array([k, j, -h])
    return new_tup


def get_quad(tup):
    if tup[0] >= 0:
        if tup[1] >= 0:
            return 1
        else:
            return 4
    else:
        if tup[1] <= 0:
            return 3
        else:
            return 2
