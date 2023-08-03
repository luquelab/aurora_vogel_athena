import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

MOVES = dict()
UNRAVELED = None
EPSILON = 10**-8
MASTER_HEX_COORDS = []


def dist(v1, v2=np.array([0, 0, 0])):
    x_diff = (v2[0] - v1[0]) + 1 / 2 * (v2[1] - v1[1] + (-v2[2]) - (-v1[2]))
    y_diff = math.sqrt(3) / 2 * (v2[1] - v1[1] + v2[2] - v1[2])
    return math.sqrt(x_diff**2 + y_diff**2)


def hex_to_cart(v1):
    x = v1[0] + v1[1] * 1/2 + v1[2] * (-1/2)
    y = v1[1] * math.sqrt(3)/2 + v1[2] * math.sqrt(3)/2
    return np.array([x, y])


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


def check_validity(tup, line_tup):
    cart_tup = hex_to_cart(tup)
    cart_line_tup = hex_to_cart(line_tup)
    if cart_line_tup[0] != 0 and cart_line_tup[1] != 0:
        line = cart_line_tup[1] / cart_line_tup[0]
        if (cart_line_tup[0] > 0 and cart_line_tup[1] > 0) or (cart_line_tup[1] < 0 < cart_line_tup[0]):  # Quad 1/4
            return cart_tup[1] < line * cart_tup[0]
        if (cart_line_tup[0] < 0 < cart_line_tup[1]) or (cart_line_tup[0] < 0 and cart_line_tup[1] < 0):  # Quad 2/3
            return cart_tup[1] > line * cart_tup[0]
    else:
        if cart_line_tup[0] == 0:
            if cart_line_tup[1] > 0:
                return cart_tup[0] > 0
            if cart_line_tup[1] < 0:
                return cart_tup[0] < 0
        else:
            if cart_line_tup[0] > 0:
                return cart_tup[1] < 0
            if cart_line_tup[0] < 0:
                return cart_tup[1] > 0


def unravel(move_info):
    global MOVES, UNRAVELED

    # unravel
    piece_of_moves = MOVES
    pieces_unraveled = [MOVES]
    for i in range(len(move_info)):
        piece_of_moves = piece_of_moves[move_info[i]]
        pieces_unraveled.append(piece_of_moves)

    UNRAVELED = pieces_unraveled


def ravel(move_info, piece_of_moves):
    global MOVES
    pieces_unraveled = UNRAVELED

    pieces_unraveled[-1] = piece_of_moves
    for i in range(len(pieces_unraveled) - 1):
        j = (len(pieces_unraveled) - 1) - i
        pieces_unraveled[j - 1][move_info[j - 1]] = pieces_unraveled[j]
    MOVES = pieces_unraveled[0]


def graph_move_set(piece_of_moves):
    move_set = list(piece_of_moves.keys())
    move_set_cart = []
    for move in move_set:
        move_set_cart.append(hex_to_cart(np.array(move)))
    DATA = pd.DataFrame(move_set_cart, columns=["X", "Y"])

    fig, ax = plt.subplots()
    ax.scatter(DATA["X"], DATA["Y"])
    plt.show()


def get_one_path(move_info):
    global MASTER_HEX_COORDS
    piece_of_moves = MOVES
    for i in range(len(move_info)):
        piece_of_moves = piece_of_moves[move_info[i]]
    move_set = list(piece_of_moves.keys())
    if len(move_set) > 0:
        for move in move_set:
            new_move_info = move_info[:]
            new_move_info.append(move)
            get_one_path(new_move_info)
    else:
        MASTER_HEX_COORDS.append(move_info)


def graph_moves(start_tup, path_length):
    global MASTER_HEX_COORDS
    move_info = []
    i = 0
    for move in list(MOVES.keys()):
        new_move_info = move_info[:]
        new_move_info.append(move)
        get_one_path(new_move_info)
    for hex_coords in MASTER_HEX_COORDS:
        hex_coords_revised = [tuple(start_tup)] + hex_coords
        if len(hex_coords_revised) == path_length + 1:
            print(hex_coords_revised)
            prev_coord = tuple([0, 0, 0])
            for k in range(len(hex_coords_revised)):
                hex_coords_revised[k] = np.array(hex_coords_revised[k]) + np.array(prev_coord)
                prev_coord = hex_coords_revised[k]

            cart_coords = [hex_to_cart(np.array(coord)) for coord in hex_coords_revised]
            X = [coord[0] for coord in cart_coords]
            Y = [coord[1] for coord in cart_coords]
            plt.plot(X, Y, label=str(i))
            i += 1
    # plt.legend(loc='upper right')
    print("There are "+str(i)+" lines.")
    plt.show()
    print("hi")


def find_moves(start_tup, end_pt=np.array([1, 1, 0]), move_info=[], max_factor=1, final_move=False):
    rot_tup_pos = rotate_sixty(start_tup, 1)
    rot_tup_neg = rotate_sixty(rotate_sixty(start_tup, -1), -1)
    diff_tup = end_pt - start_tup
    max_move = dist(diff_tup)*max_factor
    loose_max = (abs(diff_tup[0]) + abs(diff_tup[1]) + abs(diff_tup[2]))*max_factor
    # initial stuffs
    piece_of_moves = {}

    unravel(move_info)
    if not final_move:
        for m in range(0, loose_max + 1):
            for n in range(0, loose_max + 1):
                potential_move1 = np.array([m, 0, n])
                potential_move2 = np.array([m, 0, -n])
                potential_move3 = np.array([-m, 0, n])
                potential_move4 = np.array([-m, 0, -n])
                for potential_move in [potential_move1, potential_move2]:
                    if check_validity(tup=potential_move, line_tup=rot_tup_pos) and \
                            (m != 0 or n != 0) and \
                            dist(potential_move) + EPSILON <= max_move:
                        piece_of_moves[tuple(potential_move)] = {}
                for potential_move in [potential_move3, potential_move4]:
                    if check_validity(tup=potential_move, line_tup=rot_tup_pos) and\
                            (m != 0 or n != 0) and \
                            dist(potential_move) + EPSILON <= max_move:
                        piece_of_moves[tuple(potential_move)] = {}
    else:
        if check_validity(tup=diff_tup, line_tup=rot_tup_pos) and\
                            check_validity(diff_tup, rot_tup_pos):
            piece_of_moves[tuple(diff_tup)] = {}
            print("Final Move: "+str(tuple(diff_tup)))

    # graph_move_set(piece_of_moves)

    ravel(move_info, piece_of_moves)


def get_paths(start_tup, end_pt, path_length, max_factor=1, move_info=[]):
    if path_length != 0:
        unravel(move_info)
        move_list = list(UNRAVELED[-1].keys())
        if len(move_list) > 0:
            for move in move_list:
                new_move_info = move_info[:]
                new_move_info.append(move)
                get_paths(start_tup=move, end_pt=end_pt-start_tup, path_length=path_length, max_factor=max_factor,
                          move_info=new_move_info)
        else:
            final_move = path_length == 1
            find_moves(start_tup=start_tup, end_pt=end_pt, move_info=move_info, max_factor=max_factor,
                       final_move=final_move)
            get_paths(start_tup=start_tup, end_pt=end_pt, path_length=path_length - 1, max_factor=max_factor,
                      move_info=move_info)


start_tup = np.array([1, 1, 0])
our_path_length = 2
our_max_factor = 6
end_point = np.array([5, 1, 0])
get_paths(start_tup=start_tup, end_pt=end_point, path_length=our_path_length, max_factor=our_max_factor)
print("get_paths DONE")
graph_moves(start_tup, path_length=our_path_length)

