import my_library as lib
import numpy as np
import math


class Lattice:
    def __init__(self, lattice_points=[], lattice_type="hex", bottom_left_corner=lib.CartesianPoint(),
                 top_right_corner=lib.CartesianPoint()):
        self.lattice_points = lattice_points
        self.lattice_type = lattice_type
        self.bottom_left_corner = bottom_left_corner
        self.top_right_corner = top_right_corner
        if len(self.lattice_points) == 0 and self.lattice_type and self.bottom_left_corner and self.top_right_corner:
            self.make_lattice()

    def make_lattice(self):
        if self.lattice_type == "hex":
            lattice_points = []
            bottom_left_corner_hex = self.bottom_left_corner.get_hexagonal()
            bottom_left_corner_hex.reduce()
            a3 = lib.HexagonalPoint(l=1)

            h0, k0 = bottom_left_corner_hex.h, bottom_left_corner_hex.k
            xf, yf = self.top_right_corner.x, self.top_right_corner.y

            current_lattice_point = lib.HexagonalPoint(h=math.floor(h0), k=math.floor(k0))
            # lattice_points.append(lib.HexagonalPoint(h=math.floor(h0), k=math.floor(k0)))
            while current_lattice_point.get_cartesian().x <= xf:
                parity = 0
                current_lattice_point.h += 1
                potential_lattice_point = current_lattice_point
                while current_lattice_point.get_cartesian().y <= yf:
                    h, k = potential_lattice_point.h, potential_lattice_point.k
                    if parity % 2 == 0:
                        potential_lattice_point = lib.HexagonalPoint(h=h-1, k=k+1)
                    else:
                        potential_lattice_point = lib.HexagonalPoint(h=h, k=k+1)
                    if potential_lattice_point.get_cartesian().x <= xf and \
                            potential_lattice_point.get_cartesian().y <= yf:
                        lattice_points.append(potential_lattice_point)
                        parity += 1
                    else:
                        break
            self.lattice_points = lattice_points
        else:
            print(self.lattice_type + " lattice not supported")

    def make_tiling(self, kind="numpy cartesian"):
        if kind == "numpy cartesian":
            if self.lattice_type == "hex" and len(self.lattice_points) > 0:
                tiles = []
                for point in self.lattice_points:
                    hex_pt = lib.HexagonalPoint(h=1, k=1)
                    hexagon = []
                    for i in range(6):
                        hexagon.append(1/3 * hex_pt.rotate_n_pi_3(i) + point)
                    for i in range(6):
                        hex_xy = hexagon[i].get_cartesian()
                        hexagon[i] = np.array([hex_xy.x, hex_xy.y])
                    tiles.append(np.array(hexagon))
            return tiles
        else:
            print("Kind "+str(kind)+" not supported.")
