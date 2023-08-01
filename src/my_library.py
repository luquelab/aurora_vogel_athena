import numpy as np
import math


class HexagonalPoint:
    def __init__(self, h=0.0, k=0.0, l=0.0):
        self.h = h
        self.k = k
        self.l = l

    def __add__(self, b):
        self.h = self.h + b.h
        self.k = self.k + b.k
        self.l = self.k + b.l

    def __sub__(self, b):
        self.h = self.h - b.h
        self.k = self.k - b.k
        self.l = self.l - b.l

    def __mul__(self, b):
        self.h = self.h * b.h
        self.k = self.k * b.k
        self.l = self.l * b.l

    def reduce(self):
        self.h = self.h - self.l
        self.k = self.k + self.l
        self.l = 0.0

    def rotate_n_pi_3(self, n):
        if n >= 0:
            h_new = self.h
            k_new = self.k
            l_new = self.l
            for i in len(range(n)):
                l_holder = l_new
                k_new = h_new
                l_new = k_new
                h_new = -l_holder
        else:
            h_new = self.h
            k_new = self.k
            l_new = self.l
            for i in len(range(n)):
                k_holder = k_new
                k_new = l_new
                l_new = -h_new
                h_new = k_holder

    def get_cartesian(self):
        return CartesianPoint(x=self.h+1/2*self.k, y=math.sqrt(3)/2*self.k)


class CartesianPoint:
    def __init__(self, x=0.0, y=0.0):
        self.x = x
        self.y = y

    def __add__(self, b):
        self.x = self.x + b.x
        self.y = self.y + b.y

    def __sub__(self, b):
        self.x = self.x - b.x
        self.y = self.y - b.y

    def __mul__(self, b):
        self.x = self.x * b.x
        self.y = self.y * b.y

    def get_hexagonal(self):
        return HexagonalPoint(h=self.x-self.y/2*2/math.sqrt(3), k=self.y*2/math.sqrt(3))

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
