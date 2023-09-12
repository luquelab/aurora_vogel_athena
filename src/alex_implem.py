from polyhedron import *

# PDB 3J3Y #
filt1 = [[0, 5], [1, 6], [2, 7], [3, 8], [4, 9],
         [5, 11],
         [11, 17], [12, 18], [13, 19], [14, 20], [15, 21]]
pgraph = puc.pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
                                       [0, 11], [1, 13], [7, 8], [8, 9], [9, 10], [9, 15], [10, 15],
                                       [5, 11], [6, 11], [6, 12], [6, 13], [7, 13], [7, 14], [8, 14], [8, 15], [10, 16],
                                       [11, 12], [12, 13], [13, 14], [14, 15], [10, 21],
                                       [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
                                filt=filt1)
# paths and delta
top_vector_sequence = [np.array([4, 1, 0]),
                       np.array([1, 7, 0]),
                       np.array([0, 6, 0]),
                       np.array([1, 3, 0]),
                       np.array([1, 4, 0]),
                       np.array([4, 1, 0])]
bottom_vector_sequence = [np.array([1, 2, 0]),
                          np.array([1, 2, 0]),
                          np.array([1, 3, 0]),
                          np.array([5, 4, 0]),
                          np.array([0, 9, 1]),
                          np.array([1, 2, 0])]
delta = puc.lib.HexagonalPoint(h=-9, l=-1)

upper_path_data = puc.tpf.TrianglePathData(vector_sequence=top_vector_sequence)
lower_path_data = puc.tpf.TrianglePathData(vector_sequence=bottom_vector_sequence, sigma=-1)
# # lumper
# lumper = [[0, 1, 2, 3, 4], [5, 5 + len(puc.upper_path_data.face_triangles)], [-1, -2, -3, -4, -5],
#           [-6, -6 - len(puc.upper_path_data.face_triangles)]]
# unfolding and its graph
unfolding = puc.PotentialSpacedUnfolding(upper_path_data=upper_path_data, lower_path_data=lower_path_data, delta=delta,
                                         polyhedral_graph=pgraph)
glue_instructions = [[0, 1, 2, 3, 4], [5, 10], [11, 16], [17, 18, 19, 20, 21]]
P = Polyhedron(unfolding=unfolding)
######################################################################################################
######################################################################################################
epsilon = 10**-3
delta_max = puc.math.pi/3

