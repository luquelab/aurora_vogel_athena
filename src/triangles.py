import math
import numpy as np


def convert_to_cartesian(tup):
    h = tup[0]
    k = tup[1]
    j = tup[2] # not implemented yet
    cart_tup = np.array([h+k*1/2, k*math.sqrt(3)/2])
    return cart_tup


def find_angle(vec1, vec2):
    cos_angle = np.dot(vec1, vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))
    return np.arccos(cos_angle)


def rotate_sixty(tup, sign):
    h = tup[0]
    k = tup[1]
    j = tup[2]
    if sign == 1:
        new_tup = np.array([-j, h, k])
    else:
        new_tup = np.array([k, j, -h])
    return new_tup


def check_validity(tup1, tup2):
    angle = find_angle(tup1, tup2)
    if angle > 4*math.pi/3:
        return True
    elif angle < math.pi/3:
        return True
    else:
        return False


def find_things(start_tup, start_pt=np.array([0, 0]), end_pt=np.array([1, 1])):
    rot_tup = rotate_sixty(start_tup, 1)
    h, k, j = rot_tup
    top_ratio = j/k
    tran_tup = start_tup
    diff_tup = end_pt-start_tup

    for n in range(diff_tup[0]+diff_tup[1]+1):
        for m in range(diff_tup[0]+diff_tup[1]+1):
            if m/n < rot_tup:
                while True:
                    move = 0
                    new_start_tup = start_tup+np.array([move, n, m])
                    final_tup = end_pt - (tran_tup + new_start_tup)
                    tuple_valid = check_validity(start_tup, final_tup)
                    if tuple_valid:
                        move = move + 1
                    else:
                        break


start_tup = np.array([1, 1, 0])


