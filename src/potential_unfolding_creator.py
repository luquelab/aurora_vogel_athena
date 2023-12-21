import triangle_paths as tpf
import polyhedral_graph as pg
import my_library as lib
import lattice_creator as lc
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats as st
import math
import networkx as nx

PROP_TEMP = {"coordinates": [],
             "color": 'blue',
             "associated_object": []}


class PotentialSpacedUnfolding:
    def __init__(self, vertex_properties={}, edge_properties={}, face_properties={},
                 polyhedral_graph=None, upper_path_data=None, lower_path_data=None, delta=None, struct_kind="type I"):
        self.vertex_properties = vertex_properties
        self.edge_properties = edge_properties
        self.face_properties = face_properties
        self.polyhedral_graph = polyhedral_graph
        self.upper_path_data = upper_path_data
        self.lower_path_data = lower_path_data
        self.delta = delta
        self.struct_kind = struct_kind
        if upper_path_data is None or lower_path_data is None or delta is None or polyhedral_graph is None:
            pass
        else:
            self.assign_props_with_path_data(data_kind="upper", struct_kind=self.struct_kind)
            self.assign_props_with_path_data(data_kind="lower", struct_kind=self.struct_kind)
            self.assign_face_props_with_poly_graph()
            self.assign_edge_lengths()

    def assign_props_with_path_data(self, data_kind, struct_kind):
        if data_kind == "upper" and self.struct_kind == "type I":
            vertices = self.polyhedral_graph.vertices
            face_triangles = self.upper_path_data.face_triangles
            for i in range(len(face_triangles)):
                # first vertex props assigned
                props = PROP_TEMP.copy()
                props["coordinates"] = face_triangles[i][0]
                self.vertex_properties[vertices[i]] = props
                # second vertex props assigned
                props = PROP_TEMP.copy()
                props["coordinates"] = face_triangles[i][1]
                self.vertex_properties[vertices[i + len(face_triangles)]] = props
                # third vertex props assigned
                props = PROP_TEMP.copy()
                props["coordinates"] = face_triangles[i][2]
                self.vertex_properties[vertices[i + (len(face_triangles) + 1)]] = props
        if data_kind == "lower" and self.struct_kind == "type I":
            starting_point = self.upper_path_data.face_triangles[0][1] - self.delta - \
                             self.lower_path_data.defect_triangles[0][2]
            vertices = self.polyhedral_graph.vertices
            face_triangles = self.lower_path_data.face_triangles
            for i in range(len(face_triangles)):
                # first vertex props assigned
                props = PROP_TEMP.copy()
                props["coordinates"] = face_triangles[i][0] + starting_point
                self.vertex_properties[vertices[-len(face_triangles) + i]] = props
                # second vertex props assigned
                props = PROP_TEMP.copy()
                props["coordinates"] = face_triangles[i][1] + starting_point
                self.vertex_properties[vertices[-2 * len(face_triangles) - 1 + i]] = props
                # third vertex props assigned
                props = PROP_TEMP.copy()
                props["coordinates"] = face_triangles[i][2] + starting_point
                self.vertex_properties[vertices[-2 * len(face_triangles) - 1 + (i + 1)]] = props
        else:
            "data_kind invalid or struct_kind invalid"

    def assign_face_props_with_poly_graph(self):
        for face in self.polyhedral_graph.faces:
            props = PROP_TEMP.copy()
            triangle = tuple((face[0], face[1], face[2]))
            props["coordinates"] = [self.vertex_properties[face[0]]["coordinates"],
                                    self.vertex_properties[face[1]]["coordinates"],
                                    self.vertex_properties[face[2]]["coordinates"]]
            props["associated_object"] = tpf.Triangle(vertices=props["coordinates"], labels=triangle)
            self.face_properties[triangle] = props

    def assign_edge_lengths(self):
        for edge in self.polyhedral_graph.edges:
            props = PROP_TEMP.copy()
            props["coordinates"] = [self.vertex_properties[edge.vertices[0]]["coordinates"],
                                    self.vertex_properties[edge.vertices[1]]["coordinates"]]
            props["associated_object"] = edge
            self.edge_properties[tuple((edge.vertices[0], edge.vertices[1]))] = props

    def get_edge_lengths(self):
        edge_length_list = []
        for edge in self.polyhedral_graph.edges:
            if edge.vertices not in self.polyhedral_graph.filt:
                edge_key = tuple(edge.vertices)
                v1 = self.edge_properties[edge_key]["coordinates"][0].get_numpy_cartesian()
                v2 = self.edge_properties[edge_key]["coordinates"][1].get_numpy_cartesian()
                edge_length = np.linalg.norm(v2-v1)
                edge_length_list.append(edge_length)
        return edge_length_list

    def get_bounding_box(self, slack):
        bottom_left_corner = self.vertex_properties[0]["coordinates"].get_cartesian()
        top_right_corner = self.vertex_properties[0]["coordinates"].get_cartesian()
        for vertex in self.polyhedral_graph.vertices:
            cart_coordinates = self.vertex_properties[vertex]["coordinates"].get_cartesian()
            if cart_coordinates.x < bottom_left_corner.x:
                bottom_left_corner.x = cart_coordinates.x
            if cart_coordinates.y < bottom_left_corner.y:
                bottom_left_corner.y = cart_coordinates.y
            if cart_coordinates.x > top_right_corner.x:
                top_right_corner.x = cart_coordinates.x
            if cart_coordinates.y > top_right_corner.y:
                top_right_corner.y = cart_coordinates.y
        bottom_left_corner.x = bottom_left_corner.x - slack
        bottom_left_corner.y = bottom_left_corner.y - slack
        top_right_corner.x = top_right_corner.x + slack
        top_right_corner.y = top_right_corner.y + slack
        return [bottom_left_corner, top_right_corner]

    def get_vertices(self, nump):
        if nump:
            nump_vertices = []
            vertices = self.polyhedral_graph.vertices
            for vertex in vertices:
                vertex_xy = self.vertex_properties[vertex]["coordinates"].get_cartesian()
                nump_vertices.append(np.array([vertex_xy.x, vertex_xy.y]))
            return nump_vertices
        else:
            print("Non numpy return not supported.")

    def get_total_subunits(self):
        area = 0
        for face in self.polyhedral_graph.faces:
            area += self.face_properties[tuple(face)]["associated_object"].get_area()
        unit_triangle = math.sqrt(3) / 4
        subunits = round(area / unit_triangle * 3, 4)
        return subunits


def graph_unfolding(unfolding: PotentialSpacedUnfolding, export=False, path="", name=""):
    # global START, END
    # START, END = get_corners()
    bottom_left_corner, top_right_corner = unfolding.get_bounding_box(slack=2)
    lat = lc.Lattice(bottom_left_corner=bottom_left_corner, top_right_corner=top_right_corner)
    lat_points_nump = []
    for point in lat.lattice_points:
        lat_points_nump.append(point.get_numpy_cartesian())
    lat_points_nump = np.array(lat_points_nump)
    hexagons = lat.make_tiling()

    vertices = np.array(unfolding.get_vertices(nump=True))

    plt.figure(figsize=(9, 9))  # Not really a square..?
    plt.gca().set_aspect('equal', adjustable='box')
    plt.scatter(vertices[:, 0], vertices[:, 1], s=2)
    plt.scatter(lat_points_nump[:, 0], lat_points_nump[:, 1], s=1, c="black", alpha=0.1)

    # colors = ['red', 'blue', 'red', 'blue']
    colors = ['blue', 'blue', 'blue', 'blue']
    i = 0
    my_color_index = 0

    for face in unfolding.polyhedral_graph.faces:
        face_coordinates = unfolding.face_properties[tuple(face)]["coordinates"]
        face_color = unfolding.face_properties[tuple(face)]["color"]
        face_poly = plt.Polygon(np.array([face_coordinates[0].get_numpy_cartesian(),
                                          face_coordinates[1].get_numpy_cartesian(),
                                          face_coordinates[2].get_numpy_cartesian()]),
                                color=face_color, alpha=0.3)
        plt.gca().add_patch(face_poly)

    for hexagon in hexagons:
        hexagon_poly = plt.Polygon(hexagon,
                               color="black",
                               alpha=0.1)

        plt.gca().add_patch(hexagon_poly)

    if not export:
        plt.show()
    else:
        plt.savefig(path + name)  # save the figure to file
        plt.close()


def get_stats(array, kind="edge numpy", lumper=[]):
    if kind == "degree list" and len(lumper) != 0:
        new_list = []
        for lump in lumper:
            value = 0
            for i in lump:
                value += array[i]
            new_list.append(value)
        append_value = True
        for i in range(len(array)):
            for lump in lumper:
                if (lump[0] >= 0 and i in lump) or (lump[0] < 0 and -len(array)+i in lump):
                    append_value = False
                    break
                else:
                    append_value = True
            if append_value:
                new_list.append(array[i])
        array = np.array(new_list)

    # Frequencies
    unique, counts = np.unique(array, return_counts=True)

    # Do "Global" and Individual Histogram
    # plt.hist(array)
    # plt.title("Edge Length Distribution")
    # plt.xlabel("Edge Length")
    # plt.ylabel("Number of Edges")
    # plt.savefig('path/to/save/image/to.png')  # save the figure to file
    # plt.close(fig)

    stats = {"Data": array, "Frequencies": np.asarray((unique, counts)).T, "Max": np.max(array), "Min": np.min(array),
             "Mean": np.mean(array), "Median": np.median(array), "Std": np.std(array), "Mode": st.mode(array)}
    print("DATA")
    print("Data: " + str(array))
    print("Frequencies: " + str(stats["Frequencies"]))
    print("Max: " + str(stats["Max"]))
    print("Min: " + str(stats["Min"]))
    print("Mean: " + str(stats["Mean"]))
    print("Median: " + str(stats["Median"]))
    print("Std: " + str(stats["Std"]))
    print("Mode: " + str(stats["Mode"]))
    return stats


def get_easy_stats(array):
    stats = {"Data": array, "Max": np.max(array), "Min": np.min(array), "Mean": np.mean(array),
             "Median": np.median(array), "Std": np.std(array), "Mode": st.mode(array)}
    print("DATA Easy Stats")
    print("Data: " + str(array))
    print("Max: " + str(stats["Max"]))
    print("Min: " + str(stats["Min"]))
    print("Mean: " + str(stats["Mean"]))
    print("Median: " + str(stats["Median"]))
    print("Std: " + str(stats["Std"]))
    print("Mode: " + str(stats["Mode"]))
    return stats

#######################################################################################
#######################################################################################
######################### Grabbed from                          #######################
# https://stackoverflow.com/questions/55162738/python-finding-all-possible-partitions #
# -of-a-list-of-lists-subject-to-a-size-li                                            #
#######################################################################################


# Main function
def get_all_possible_partitions(lst, c):
    yield from _get_all_possible_partitions_rec(lst, c, [False] * len(lst), [])


# Produces partitions recursively
def _get_all_possible_partitions_rec(lst, c, picked, partition):
    # If all elements have been picked it is a complete partition
    if all(picked):
        yield tuple(partition)
    else:
        # Get all possible subsets of unpicked elements
        for subset in _get_all_possible_subsets_rec(lst, c, picked, [], 0):
            # Add the subset to the partition
            partition.append(subset)
            # Generate all partitions that complete the current one
            yield from _get_all_possible_partitions_rec(lst, c, picked, partition)
            # Remove the subset from the partition
            partition.pop()


# Produces all possible subsets of unpicked elements
def _get_all_possible_subsets_rec(lst, c, picked, current, idx):
    # If we have gone over all elements finish
    if idx >= len(lst): return
    # If the current element is available and fits in the subset
    if not picked[idx] and len(lst[idx]) <= c:
        # Mark it as picked
        picked[idx] = True
        # Add it to the subset
        current.append(lst[idx])
        # Generate the subset
        yield tuple(current)
        # Generate all possible subsets extending this one
        yield from _get_all_possible_subsets_rec(lst, c - len(lst[idx]), picked, current, idx + 1)
        # Remove current element
        current.pop()
        # Unmark as picked
        picked[idx] = False
    # Only allow skip if it is not the first available element
    if len(current) > 0 or picked[idx]:
        # Get all subsets resulting from skipping current element
        yield from _get_all_possible_subsets_rec(lst, c, picked, current, idx + 1)

#######################################################################################
#######################################################################################


def brute_modularity_optim(G):
    nodes = list(G.nodes())
    nodes_mod = []
    for node in nodes:
        nodes_mod.append([node])

    # partitions = list(get_all_possible_partitions(nodes_mod, len(nodes_mod)))
    partitions = list(get_all_possible_partitions(nodes_mod, 1))
    # print("Number of partitions found: " + str(len(partitions)))
    best_mod_part = [set(nodes)]
    for part in partitions:
        # print(part)
        fixed_part = []
        for piece in part:
            # print(piece)
            fixed_piece = set()
            for num_lst in piece:
                # print(num_lst)
                fixed_piece.add(num_lst[0])
            fixed_part.append(fixed_piece)
        if nx.community.modularity(G, fixed_part, weight='weight') > nx.community.modularity(G, best_mod_part,
                                                                                             weight='weight'):
            best_mod_part = fixed_part.copy()

    print(best_mod_part)
    print(nx.community.modularity(G, best_mod_part, weight='weight'))
    return {'community': best_mod_part, 'graph': G, 'modularity': nx.community.modularity(G, best_mod_part,
                                                                                          weight='weight')}


def cluster_analysis(unfolding: PotentialSpacedUnfolding, lumper: list, data=False, calc_type='max'):
    # Cluster Analysis #1: NetworkX
    # max edge
    edge_lengths = unfolding.get_edge_lengths()
    max_edge = np.array(0)
    for edge_len in edge_lengths:
        max_edge = max(max_edge, edge_len)
    #
    vert_number = len(unfolding.polyhedral_graph.vertices)
    revised_lumper = []
    for lump in lumper:
        revised_lump = []
        for value in lump:
            if value < 0:
                revised_lump.append(vert_number+value)
            else:
                revised_lump.append(value)
        revised_lumper.append(revised_lump.copy())

    G = nx.Graph()
    initial_edges = unfolding.polyhedral_graph.edges.copy()
    edges = []
    for edge in initial_edges:
        current_edge = edge.vertices.copy()
        edge_key = tuple(current_edge)
        if edge_key in list(unfolding.edge_properties.keys()):
            v1 = unfolding.edge_properties[edge_key]["coordinates"][0].get_numpy_cartesian()
            v2 = unfolding.edge_properties[edge_key]["coordinates"][1].get_numpy_cartesian()
            edge_length = np.linalg.norm(v2 - v1)
        for lump in revised_lumper:
            if current_edge[0] in lump:
                current_edge[0] = lump[0]
            if current_edge[1] in lump:
                current_edge[1] = lump[0]
        if calc_type == "inv":
            G.add_edge(current_edge[0], current_edge[1], weight=1/edge_length, length=edge_length)
        else:
            G.add_edge(current_edge[0], current_edge[1], weight=max_edge - edge_length, length=edge_length)

    ## Brute modularity optimizer
    brute_modularity_info = brute_modularity_optim(G)
    #

    communities = nx.community.greedy_modularity_communities(G, weight='weight')
    # communities = nx.community.louvain_communities(G, weight='weight', seed=123)
    print("NetworkX greedy_modularity_communities results:")
    for community in communities:
        print("Community: " + str(sorted(community)))
    for edge in list(G.edges()):
        print("Edge: " + str(edge) + " and Length: " + str(G.edges[edge]['weight']))

    if data:
        return {"brute": brute_modularity_info, "greedy": communities, "graph": G}






## Example 1 ##
# upper_path_data = tpf.TrianglePathData(vector_sequence=[[1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0]])
# lower_path_data = tpf.TrianglePathData(vector_sequence=[[1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0]], sigma=-1)
# delta = lib.HexagonalPoint()
# pgraph = pg.PolyhedralGraph(edges=[[0, 3], [0, 4], [1, 4], [1, 5], [2, 5], [2, 6], [3, 4], [4, 5], [5, 6],
#                                    [3, 7], [4, 7], [4, 8], [5, 8], [5, 9], [6, 9]])
# unfolding = PotentialSpacedUnfolding(upper_path_data=upper_path_data, lower_path_data=lower_path_data, delta=delta,
#                                      polyhedral_graph=pgraph)
# graph_unfolding(unfolding=unfolding)


## Example 2 ##
# icosahedral graph
# filt1 = [[0, 5], [1, 6], [2, 7], [3, 8], [4, 9],
#          [5, 11],
#          [11, 17], [12, 18], [13, 19], [14, 20], [15, 21]]
# pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
#                                    [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
#                                    [5, 11], [5, 12], [6, 12], [6, 13], [7, 13], [7, 14], [8, 14], [8, 15], [9, 15], [9, 16], [10, 16],
#                                    [11, 12], [12, 13], [13, 14], [14, 15], [15, 16],
#                                    [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
#                             filt=filt1)
# # paths and delta
# upper_path_data = tpf.TrianglePathData(vector_sequence=[[1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0]])
# lower_path_data = tpf.TrianglePathData(vector_sequence=[[1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0]], sigma=-1)
# delta = lib.HexagonalPoint(k=1)
# # lumper
# lumper = [[0, 1, 2, 3, 4], [5, 5 + len(upper_path_data.face_triangles)], [-1, -2, -3, -4, -5], [-6, -6 - len(upper_path_data.face_triangles)]]
# # unfolding and its graph
# unfolding = PotentialSpacedUnfolding(upper_path_data=upper_path_data, lower_path_data=lower_path_data, delta=delta,
#                                      polyhedral_graph=pgraph)
# get_stats(unfolding.get_edge_lengths())
# print(unfolding.polyhedral_graph.get_degrees())
# get_stats(unfolding.polyhedral_graph.get_degrees(), kind="degree list", lumper=lumper)
# graph_unfolding(unfolding=unfolding)

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
## PDB 3J3Y ##
filt1 = [[0, 5], [1, 6], [2, 7], [3, 8], [4, 9],
         [5, 11],
         [11, 17], [12, 18], [13, 19], [14, 20], [15, 21]]
# pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
#                                    [0, 11], [1, 13], [7, 8], [8, 9], [9, 10], [9, 15], [10, 15],
#                                    [5, 11], [6, 11], [6, 12], [6, 13], [7, 13], [7, 14], [8, 14], [8, 15], [10, 16],
#                                    [11, 12], [12, 13], [13, 14], [14, 15], [10, 21],
#                                    [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
#                             filt=filt1)
pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
                                   [5, 6], [6, 7], [2, 13], [8, 15], [3, 15], [9, 15], [9, 10],
                                   [5, 11], [5, 12], [6, 12], [7, 12], [7, 13], [8, 13], [8, 14], [9, 16], [10, 16],
                                   [11, 12], [7, 18], [13, 14], [14, 15], [15, 16],
                                   [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
                            filt=filt1)
# # # paths and delta 1
# # top_vector_sequence = [np.array([4, 1, 0]),
# #           np.array([1, 7, 0]),
# #           np.array([0, 6, 0]),
# #           np.array([1, 3, 0]),
# #           np.array([1, 4, 0]),
# #           np.array([4, 1, 0])]
# # bottom_vector_sequence = [np.array([1, 2, 0]),
# #           np.array([1, 2, 0]),
# #           np.array([1, 3, 0]),
# #           np.array([5, 4, 0]),
# #           np.array([0, 9, 1]),
# #           np.array([1, 2, 0])]
# # delta = lib.HexagonalPoint(h=-9, l=-1)
# # paths and delta 2
# top_vector_sequence = [np.array([1, 3, 0]),
#           np.array([1, 4, 0]),
#           np.array([4, 1, 0]),
#           np.array([1, 7, 0]),
#           np.array([0, 6, 0]),
#           np.array([1, 3, 0])]
# bottom_vector_sequence = [np.array([5, 4, 0]),
#           np.array([0, 9, 1]),
#           np.array([1, 2, 0]),
#           np.array([1, 2, 0]),
#           np.array([1, 3, 0]),
#           np.array([5, 4, 0])]
# delta = lib.HexagonalPoint(h=-2, l=3)
#
# upper_path_data = tpf.TrianglePathData(vector_sequence=top_vector_sequence)
# lower_path_data = tpf.TrianglePathData(vector_sequence=bottom_vector_sequence, sigma=-1)
# # lumper
# lumper = [[0, 1, 2, 3, 4], [5, 5 + len(upper_path_data.face_triangles)], [-1, -2, -3, -4, -5], [-6, -6 - len(upper_path_data.face_triangles)]]
# # unfolding and its graph
# unfolding = PotentialSpacedUnfolding(upper_path_data=upper_path_data, lower_path_data=lower_path_data, delta=delta,
#                                      polyhedral_graph=pgraph)
# cluster_analysis(unfolding, lumper)
# get_stats(unfolding.get_edge_lengths())
# print(unfolding.polyhedral_graph.get_degrees())
# get_stats(unfolding.polyhedral_graph.get_degrees(), kind="degree list", lumper=lumper)
# print(unfolding.get_total_subunits())
# graph_unfolding(unfolding=unfolding)

## PDB 3J3Y ## ADJUSTED
filt1 = [[0, 5], [1, 6], [2, 7], [3, 8], [4, 9],
         [5, 11],
         [11, 17], [12, 18], [13, 19], [14, 20], [15, 21]]
# pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
#                                    [0, 11], [1, 13], [7, 8], [8, 9], [9, 10], [9, 15], [10, 15],
#                                    [5, 11], [6, 11], [6, 12], [6, 13], [7, 13], [7, 14], [8, 14], [8, 15], [10, 16],
#                                    [11, 12], [12, 13], [13, 14], [14, 15], [10, 21],
#                                    [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
#                             filt=filt1)
pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
                                   [5, 6], [6, 7], [2, 13], [8, 15], [3, 15], [9, 15], [9, 10],
                                   [5, 11], [5, 12], [6, 12], [7, 12], [7, 13], [8, 13], [8, 14], [9, 16], [10, 16],
                                   [11, 12], [7, 18], [13, 14], [14, 15], [15, 16],
                                   [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
                            filt=filt1)
# # paths and delta 1
# top_vector_sequence = [np.array([1, 3, 0]),
#           np.array([1, 4, 0]),
#           np.array([4, 1, 0]),
#           np.array([1, 7, 0]),
#           np.array([0, 6, 0]),
#           np.array([1, 3, 0])]
# bottom_vector_sequence = [np.array([5, 4, 0]),
#           np.array([0, 9, 1]),
#           np.array([1, 2, 0]),
#           np.array([1, 2, 0]),
#           np.array([1, 3, 0]),
#           np.array([5, 4, 0])]
# delta = lib.HexagonalPoint(h=-2, l=3)
# # paths and delta ADJ 1
top_vector_sequence = [np.array([1, 4, 0]),
          np.array([1, 4, 0]),
          np.array([4, 1, 0]),
          np.array([1, 7, 0]),
          np.array([0, 6, 0]),
          np.array([1, 4, 0])]
bottom_vector_sequence = [np.array([5, 5, 0]),
          np.array([0, 9, 1]),
          np.array([1, 2, 0]),
          np.array([1, 2, 0]),
          np.array([1, 3, 0]),
          np.array([5, 5, 0])]
delta = lib.HexagonalPoint(h=-4, k=2)
# paths and delta ADJ 2
filt1 = [[0, 5], [1, 6], [2, 7], [3, 8], [4, 9],
         [5, 11],
         [11, 17], [12, 18], [13, 19], [14, 20], [15, 21]]
pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
                                   [5, 6], [6, 7], [7, 8], [8, 9], [4, 15], [9, 15], [10, 15],
                                   [5, 11], [6, 11], [6, 12], [7, 12], [7, 13], [8, 13], [9, 13], [9, 14], [10, 16],
                                   [11, 12], [12, 13], [13, 14], [14, 15], [15, 16],
                                   [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
                            filt=filt1)
top_vector_sequence = [np.array([0, 5, 0]),
          np.array([0, 4, 0]),
          np.array([3, 2, 0]),
          np.array([0, 3, 0]),
          np.array([3, 7, 0]),
          np.array([0, 5, 0])]
bottom_vector_sequence = [np.array([3, 5, 0]),
          np.array([0, 9, 0]),
          np.array([1, 2, 0]),
          np.array([1, 2, 0]),
          np.array([1, 3, 0]),
          np.array([3, 5, 0])]
delta = lib.HexagonalPoint(h=-4, k=2)
upper_path_data = tpf.TrianglePathData(vector_sequence=top_vector_sequence)
lower_path_data = tpf.TrianglePathData(vector_sequence=bottom_vector_sequence, sigma=-1)
# lumper
lumper = [[0, 1, 2, 3, 4], [5, 5 + len(upper_path_data.face_triangles)], [-1, -2, -3, -4, -5], [-6, -6 - len(upper_path_data.face_triangles)]]
# unfolding and its graph
unfolding = PotentialSpacedUnfolding(upper_path_data=upper_path_data, lower_path_data=lower_path_data, delta=delta,
                                     polyhedral_graph=pgraph)
cluster_analysis(unfolding, lumper)
get_stats(unfolding.get_edge_lengths())
print(unfolding.polyhedral_graph.get_degrees())
get_stats(unfolding.polyhedral_graph.get_degrees(), kind="degree list", lumper=lumper)
print(unfolding.get_total_subunits())
graph_unfolding(unfolding=unfolding)


#
# # PDB 3J3Q
filt1 = [[0, 5], [1, 6], [2, 7], [3, 8], [4, 9],
         [5, 11],
         [11, 17], [12, 18], [13, 19], [14, 20], [15, 21]]
pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
                                   [5, 6], [6, 7], [7, 8], [8, 9], [4, 15], [9, 15], [10, 15],
                                   [5, 11], [5, 12], [6, 12], [7, 12], [7, 13], [8, 13], [9, 13], [9, 14], [10, 16],
                                   [11, 12], [12, 13], [13, 14], [14, 15], [15, 16],
                                   [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
                            filt=filt1)
# paths and delta 1
top_vector_sequence = [np.array([0, 5, 0]),
          np.array([0, 5, 0]),
          np.array([3, 2, 0]),
          np.array([0, 2, 0]),
          np.array([3, 7, 0]),
          np.array([0, 5, 0])]
bottom_vector_sequence = [np.array([3, 5, 0]),
          np.array([0, 9, 0]),
          np.array([1, 2, 0]),
          np.array([1, 2, 0]),
          np.array([1, 3, 0]),
          np.array([3, 5, 0])]
delta = lib.HexagonalPoint(h=-4, k=2)



upper_path_data = tpf.TrianglePathData(vector_sequence=top_vector_sequence)
lower_path_data = tpf.TrianglePathData(vector_sequence=bottom_vector_sequence, sigma=-1)
# lumper
lumper = [[0, 1, 2, 3, 4], [5, 5 + len(upper_path_data.face_triangles)], [-1, -2, -3, -4, -5], [-6, -6 - len(upper_path_data.face_triangles)]]
# unfolding and its graph
unfolding = PotentialSpacedUnfolding(upper_path_data=upper_path_data, lower_path_data=lower_path_data, delta=delta,
                                     polyhedral_graph=pgraph)
cluster_analysis(unfolding, lumper)
get_stats(unfolding.get_edge_lengths())
print(unfolding.polyhedral_graph.get_degrees())
get_stats(unfolding.polyhedral_graph.get_degrees(), kind="degree list", lumper=lumper)
print(unfolding.get_total_subunits())
graph_unfolding(unfolding=unfolding)
#
#
# ## HIV mattei et al 2016 struct 1_1 ##
# filt1 = [[0, 5], [1, 6], [2, 7], [3, 8], [4, 9],
#          [5, 11],
#          [11, 17], [12, 18], [13, 19], [14, 20], [15, 21]]
# pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
#                                    [5, 6], [6, 7], [7, 8], [3, 14], [9, 10],
#                                    [5, 11], [6, 11], [6, 12], [7, 12], [7, 13], [8, 13], [8, 14], [9, 14], [9, 15], [9, 16], [10, 16],
#                                    [11, 12], [12, 13], [13, 14], [14, 15], [15, 16],
#                                    [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
#                             filt=filt1)
# # paths and delta
# top_vector_sequence = [np.array([0, 10, 2]),
#           np.array([9, 4, 0]),
#           np.array([5, 0, -8]),
#           np.array([2, 0, 0]),
#           np.array([0, 5, 0]),
#           np.array([0, 10, 2])]
# bottom_vector_sequence = [np.array([11, 0, 0]),
#           np.array([2, 10, 0]),
#           np.array([2, 1, 0]),
#           np.array([4, 1, 0]),
#           np.array([3, 1, 0]),
#           np.array([11, 0, 0])]
# delta = lib.HexagonalPoint(l=2.0)
#
# upper_path_data = tpf.TrianglePathData(vector_sequence=top_vector_sequence)
# lower_path_data = tpf.TrianglePathData(vector_sequence=bottom_vector_sequence, sigma=-1)
# # lumper
# lumper = [[0, 1, 2, 3, 4], [5, 5 + len(upper_path_data.face_triangles)], [-1, -2, -3, -4, -5], [-6, -6 - len(upper_path_data.face_triangles)]]
# # unfolding and its graph
# unfolding = PotentialSpacedUnfolding(upper_path_data=upper_path_data, lower_path_data=lower_path_data, delta=delta,
#                                      polyhedral_graph=pgraph)
# cluster_analysis(unfolding, lumper)
# get_stats(unfolding.get_edge_lengths())
# print(unfolding.polyhedral_graph.get_degrees())
# get_stats(unfolding.polyhedral_graph.get_degrees(), kind="degree list", lumper=lumper)
# print(unfolding.get_total_subunits())
# graph_unfolding(unfolding=unfolding)
#
# ## HIV mattei et al 2016 struct 1_2 ##
# filt1 = [[0, 5], [1, 6], [2, 7], [3, 8], [4, 9],
#          [5, 11],
#          [11, 17], [12, 18], [13, 19], [14, 20], [15, 21]]
# pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
#                                    [5, 6], [6, 7], [2, 14], [8, 9], [9, 10],
#                                    [5, 11], [5, 12], [6, 12], [7, 12], [7, 13], [7, 14], [8, 14], [8, 15], [9, 15], [10, 15], [10, 16],
#                                    [11, 12], [12, 13], [13, 14], [14, 15], [15, 16],
#                                    [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
#                             filt=filt1)
# # paths and delta
# top_vector_sequence = [np.array([0, 2, 9]),
#           np.array([0, 9, 2]),
#           np.array([6, 5, 0]),
#           np.array([0, 1, 2]),
#           np.array([0, 1, 1]),
#           np.array([0, 2, 9])]
# bottom_vector_sequence = [np.array([0, 3, 1]),
#           np.array([0, 1, 3]),
#           np.array([0, 2, 0]),
#           np.array([3, 9, 0]),
#           np.array([0, 6, 7]),
#           np.array([0, 3, 1])]
# delta = lib.HexagonalPoint(h=-3.0, l=1.0)
#
# upper_path_data = tpf.TrianglePathData(vector_sequence=top_vector_sequence)
# lower_path_data = tpf.TrianglePathData(vector_sequence=bottom_vector_sequence, sigma=-1)
# # lumper
# lumper = [[0, 1, 2, 3, 4], [5, 5 + len(upper_path_data.face_triangles)], [-1, -2, -3, -4, -5], [-6, -6 - len(upper_path_data.face_triangles)]]
# # unfolding and its graph
# unfolding = PotentialSpacedUnfolding(upper_path_data=upper_path_data, lower_path_data=lower_path_data, delta=delta,
#                                      polyhedral_graph=pgraph)
# cluster_analysis(unfolding, lumper)
# get_stats(unfolding.get_edge_lengths())
# print(unfolding.polyhedral_graph.get_degrees())
# get_stats(unfolding.polyhedral_graph.get_degrees(), kind="degree list", lumper=lumper)
# print(unfolding.get_total_subunits())
# graph_unfolding(unfolding=unfolding)
#
# ## HIV mattei et al 2016 struct 1_3 ##
# filt1 = [[0, 6], [1, 7], [2, 8], [3, 9], [4, 10], [5, 11],
#          [6, 13],
#          [13, 18], [14, 19], [15, 20], [16, 21]]
# pgraph = pg.PolyhedralGraph(edges=[[0, 6], [0, 7], [1, 7], [1, 8], [2, 8], [2, 9], [3, 9], [3, 10], [4, 10], [4, 11], [5, 11], [5, 12],
#                                    [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12],
#                                    [6, 13], [7, 13], [7, 14], [8, 14], [8, 15], [9, 15], [10, 15], [10, 16], [11, 16], [11, 17], [12, 17],
#                                    [13, 14], [14, 15], [15, 16], [16, 17],
#                                    [13, 18], [14, 18], [14, 19], [15, 19], [15, 20], [16, 20], [16, 21], [17, 21]],
#                             filt=filt1)
# # paths and delta
# top_vector_sequence = [np.array([1, 2, 0]),
#           np.array([0, 4, 0]),
#           np.array([1, 4, 0]),
#           np.array([3, 2, 0]),
#           np.array([4, 0, 0]),
#           np.array([2, 1, 0]),
#           np.array([2, 0, -1])]
# bottom_vector_sequence = [np.array([0, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 0, 0])]
# delta = lib.HexagonalPoint(h=-7.0, l=6.0)
#
# upper_path_data = tpf.TrianglePathData(vector_sequence=top_vector_sequence)
# lower_path_data = tpf.TrianglePathData(vector_sequence=bottom_vector_sequence, sigma=-1)
# # lumper
# lumper = [[0, 1, 2, 3, 4, 5], [6, 6 + len(upper_path_data.face_triangles)], [-1, -2, -3, -4], [-5, -5 - len(lower_path_data.face_triangles)]]
# # unfolding and its graph
# unfolding = PotentialSpacedUnfolding(upper_path_data=upper_path_data, lower_path_data=lower_path_data, delta=delta,
#                                      polyhedral_graph=pgraph)
# cluster_analysis(unfolding, lumper)
# get_stats(unfolding.get_edge_lengths())
# print(unfolding.polyhedral_graph.get_degrees())
# get_stats(unfolding.polyhedral_graph.get_degrees(), kind="degree list", lumper=lumper)
# print(unfolding.get_total_subunits())
# graph_unfolding(unfolding=unfolding)
#
#
# ## HIV mattei et al 2016 struct 1_4 ##
# filt1 = [[0, 5], [1, 6], [2, 7], [3, 8], [4, 9],
#          [5, 11],
#          [11, 17], [12, 18], [13, 19], [14, 20], [15, 21]]
# pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
#                                    [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
#                                    [5, 11], [6, 11], [6, 12], [7, 12], [7, 13], [8, 13], [8, 14], [9, 14], [9, 15], [9, 16], [10, 16],
#                                    [6, 17], [12, 13], [13, 14], [14, 15], [15, 16],
#                                    [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
#                             filt=filt1)
# # paths and delta
# top_vector_sequence = [np.array([0, 1, 1]),
#           np.array([0, 2, 10]),
#           np.array([0, 11, 1]),
#           np.array([0, 1, 1]),
#           np.array([0, 1, 1]),
#           np.array([0, 1, 1])]
# bottom_vector_sequence = [np.array([-2, 0, 7]),
#           np.array([0, 2, 0]),
#           np.array([0, 2, 1]),
#           np.array([1, 8, 0]),
#           np.array([0, 5, 5]),
#           np.array([-2, 0, 7])]
# delta = lib.HexagonalPoint(h=-5.0, k=-2.0)
#
# upper_path_data = tpf.TrianglePathData(vector_sequence=top_vector_sequence)
# lower_path_data = tpf.TrianglePathData(vector_sequence=bottom_vector_sequence, sigma=-1)
# # lumper
# lumper = [[0, 1, 2, 3, 4], [5, 5 + len(upper_path_data.face_triangles)], [-1, -2, -3, -4, -5], [-6, -6 - len(upper_path_data.face_triangles)]]
# # unfolding and its graph
# unfolding = PotentialSpacedUnfolding(upper_path_data=upper_path_data, lower_path_data=lower_path_data, delta=delta,
#                                      polyhedral_graph=pgraph)
# cluster_analysis(unfolding, lumper)
# get_stats(unfolding.get_edge_lengths())
# print(unfolding.polyhedral_graph.get_degrees())
# get_stats(unfolding.polyhedral_graph.get_degrees(), kind="degree list", lumper=lumper)
# print(unfolding.get_total_subunits())
# graph_unfolding(unfolding=unfolding)
#
# ## HIV mattei et al 2016 struct 1_5 ##
# filt1 = [[0, 5], [1, 6], [2, 7], [3, 8], [4, 9],
#          [5, 11],
#          [11, 17], [12, 18], [13, 19], [14, 20], [15, 21]]
# pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
#                                    [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
#                                    [5, 11], [6, 11], [6, 12], [7, 12], [7, 13], [8, 13], [8, 14], [9, 14], [9, 15], [10, 15], [10, 16],
#                                    [11, 12], [12, 13], [13, 14], [14, 15], [15, 16],
#                                    [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
#                             filt=filt1)
# # paths and delta
# upper_path_data = tpf.TrianglePathData(vector_sequence=[[2, 0, 0], [1, 1, 0], [1, 2, 0], [2, 2, 0], [2, 1, 0], [2, 0, 0]])
# lower_path_data = tpf.TrianglePathData(vector_sequence=[[2, 2, 0], [1, 2, 0], [1, 1, 0], [1, 1, 0], [3, 0, 0], [2, 2, 0]], sigma=-1)
# delta = lib.HexagonalPoint(h=-1.0, l=11.0)
# # lumper
# lumper = [[0, 1, 2, 3, 4], [5, 5 + len(upper_path_data.face_triangles)], [-1, -2, -3, -4, -5], [-6, -6 - len(upper_path_data.face_triangles)]]
# # unfolding and its graph
# unfolding = PotentialSpacedUnfolding(upper_path_data=upper_path_data, lower_path_data=lower_path_data, delta=delta,
#                                      polyhedral_graph=pgraph)
# cluster_analysis(unfolding, lumper)
# get_stats(unfolding.get_edge_lengths())
# print(unfolding.polyhedral_graph.get_degrees())
# get_stats(unfolding.polyhedral_graph.get_degrees(), kind="degree list", lumper=lumper)
# print(unfolding.get_total_subunits())
# graph_unfolding(unfolding=unfolding)
#
# ## HIV mattei et al 2016 struct 1_6 ##
# filt1 = [[0, 5], [1, 6], [2, 7], [3, 8], [4, 9],
#          [5, 11],
#          [11, 17], [12, 18], [13, 19], [14, 20], [15, 21]]
# pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
#                                    [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
#                                    [5, 11], [6, 11], [6, 12], [7, 12], [7, 13], [8, 13], [8, 14], [9, 14], [9, 15], [10, 15], [10, 16],
#                                    [11, 12], [12, 13], [13, 14], [14, 15], [15, 16],
#                                    [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
#                             filt=filt1)
# # paths and delta
# upper_path_data = tpf.TrianglePathData(vector_sequence=[[2, 1, 0], [3, 2, 0], [7, 3, 0], [1, 1, 0], [2, 1, 0], [2, 1, 0]])
# lower_path_data = tpf.TrianglePathData(vector_sequence=[[3, 3, 0], [4, 2, 0], [2, 3, 0], [2, 0, 0], [4, 0, 0], [3, 3, 0]], sigma=-1)
# delta = lib.HexagonalPoint(h=-1.0, l=7.0)
# # lumper
# lumper = [[0, 1, 2, 3, 4], [5, 5 + len(upper_path_data.face_triangles)], [-1, -2, -3, -4, -5], [-6, -6 - len(upper_path_data.face_triangles)]]
# # unfolding and its graph
# unfolding = PotentialSpacedUnfolding(upper_path_data=upper_path_data, lower_path_data=lower_path_data, delta=delta,
#                                      polyhedral_graph=pgraph)
# cluster_analysis(unfolding, lumper)
# get_stats(unfolding.get_edge_lengths())
# print(unfolding.polyhedral_graph.get_degrees())
# get_stats(unfolding.polyhedral_graph.get_degrees(), kind="degree list", lumper=lumper)
# print(unfolding.get_total_subunits())
# graph_unfolding(unfolding=unfolding)


###############
##############
##############

# ## HIV mattei et al 2016 struct 1_3 ##
# filt1 = [[0, 6], [1, 7], [2, 8], [3, 9], [4, 10], [5, 11],
#          [6, 13],
#          [13, 18], [14, 19], [15, 20], [16, 21]]
# pgraph = pg.PolyhedralGraph(edges=[[0, 6], [0, 7], [1, 7], [1, 8], [2, 8], [2, 9], [3, 9], [3, 10], [4, 10], [4, 11], [5, 11], [5, 12],
#                                    [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12],
#                                    [6, 13], [7, 13], [7, 14], [8, 14], [8, 15], [9, 15], [10, 15], [10, 16], [11, 16], [11, 17], [12, 17],
#                                    [13, 14], [14, 15], [15, 16], [16, 17],
#                                    [13, 18], [14, 18], [14, 19], [15, 19], [15, 20], [16, 20], [16, 21], [17, 21]],
#                             filt=filt1)
# # paths and delta
# top_vector_sequence = [np.array([1, 2, 0]),
#           np.array([0, 4, 0]),
#           np.array([1, 4, 0]),
#           np.array([3, 2, 0]),
#           np.array([4, 0, 0]),
#           np.array([2, 1, 0]),
#           np.array([2, 0, -1])]
# bottom_vector_sequence = [np.array([0, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 1, 0]),
#           np.array([1, 0, 0])]
# delta = lib.HexagonalPoint(h=-7.0, l=6.0)
#
# upper_path_data = tpf.TrianglePathData(vector_sequence=top_vector_sequence)
# lower_path_data = tpf.TrianglePathData(vector_sequence=bottom_vector_sequence, sigma=-1)
# # lumper
# lumper = [[0, 1, 2, 3, 4, 5], [6, 6 + len(upper_path_data.face_triangles)], [-4, -3, -2, -1], [-5 - len(lower_path_data.face_triangles), -5]]
# # unfolding and its graph
# unfolding = PotentialSpacedUnfolding(upper_path_data=upper_path_data, lower_path_data=lower_path_data, delta=delta,
#                                      polyhedral_graph=pgraph)
# cluster_analysis(unfolding, lumper)
# get_stats(unfolding.get_edge_lengths())
# print(unfolding.polyhedral_graph.get_degrees())
# get_stats(unfolding.polyhedral_graph.get_degrees(), kind="degree list", lumper=lumper)
# print(unfolding.get_total_subunits())
# graph_unfolding(unfolding=unfolding)
