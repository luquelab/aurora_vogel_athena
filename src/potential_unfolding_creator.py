import triangle_paths as tpf
import polyhedral_graph as pg
import my_library as lib
import lattice_creator as lc
import numpy as np
import matplotlib.pyplot as plt

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
            starting_point = self.upper_path_data.face_triangles[0][1] - delta - \
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
            props["associated_object"] = tpf.Triangle(vertices=props["coordinates"])
            self.face_properties[triangle] = props

    def assign_edge_lengths(self):
        for edge in self.polyhedral_graph.edges:
            props = PROP_TEMP.copy()
            props["coordinates"] = [self.vertex_properties[edge.vertices[0]]["coordinates"],
                                    self.vertex_properties[edge.vertices[1]]["coordinates"]]
            props["associated_object"] = edge
            self.edge_properties[tuple((edge.vertices[0], edge.vertices[1]))] = props

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


def graph_unfolding(unfolding: PotentialSpacedUnfolding):
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
                                color=face_color, alpha=0.5)
        plt.gca().add_patch(face_poly)

    for hexagon in hexagons:
        hexagon_poly = plt.Polygon(hexagon,
                               color="black",
                               alpha=0.1)

        plt.gca().add_patch(hexagon_poly)

    plt.show()

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
pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
                                   [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
                                   [5, 11], [5, 12], [6, 12], [6, 13], [7, 13], [7, 14], [8, 14], [8, 15], [9, 15], [9, 16], [10, 16],
                                   [11, 12], [12, 13], [13, 14], [14, 15], [15, 16],
                                   [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]])
# paths and delta
upper_path_data = tpf.TrianglePathData(vector_sequence=[[1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0]])
lower_path_data = tpf.TrianglePathData(vector_sequence=[[1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0]], sigma=-1)
delta = lib.HexagonalPoint(k=1)
# unfolding and its graph
unfolding = PotentialSpacedUnfolding(upper_path_data=upper_path_data, lower_path_data=lower_path_data, delta=delta,
                                     polyhedral_graph=pgraph)
graph_unfolding(unfolding=unfolding)
