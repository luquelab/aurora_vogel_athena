import triangle_paths as tpf
import polyhedral_graph as pg
import my_library as lib

PROP_TEMP = {"coordinates": [],
             "color": [],
             "associated_object": []}


class PotentialSpacedUnfolding:
    def __init__(self, vertex_properties={}, edge_properties={}, face_properties={},
                 polyhedral_graph=None, upper_path_data=None, lower_path_data=None, delta=None, struct_kind=""):
        if upper_path_data is None or lower_path_data is None or delta is None or polyhedral_graph is None:
            self.vertex_properties = vertex_properties
            self.edge_properties = edge_properties
            self.face_properties = face_properties
            self.polyhedral_graph = polyhedral_graph
            self.upper_path_data = upper_path_data
            self.lower_path_data = lower_path_data
            self.delta = delta
            self.struct_kind = struct_kind
        else:
            self.struct_kind = struct_kind
            self.assign_props_with_path_data(data_kind="upper", struct_kind=self.struct_kind)
            self.assign_props_with_path_data(data_kind="lower", struct_kind=self.struct_kind)
            self.polyhedral_graph = polyhedral_graph
            self.assign_face_props_with_poly_graph()
            self.assign_edge_lengths()

    def assign_props_with_path_data(self, data_kind, struct_kind):
        if data_kind == "upper" and self.struct_kind == "type I":
            vertices = self.polyhedral_graph.vertices
            face_triangles = self.upper_path_data.face_triangles
            for i in range(len(face_triangles)):
                # first vertex props assigned
                props = PROP_TEMP
                props["coordinates"] = face_triangles[i][0]
                self.vertex_properties[vertices[i]] = props
                # second vertex props assigned
                props = PROP_TEMP
                props["coordinates"] = face_triangles[i][1]
                self.vertex_properties[vertices[i + len(face_triangles)]] = props
                # third vertex props assigned
                props = PROP_TEMP
                props["coordinates"] = face_triangles[i][2]
                self.vertex_properties[vertices[i + (len(face_triangles) + 1)]] = props
        if data_kind == "lower" and self.struct_kind == "type I":
            vertices = self.polyhedral_graph.vertices
            face_triangles = self.lower_path_data.face_triangles
            for i in range(len(face_triangles)):
                # first vertex props assigned
                props = PROP_TEMP
                props["coordinates"] = face_triangles[i][0]
                self.vertex_properties[vertices[-len(face_triangles) + i]] = props
                # second vertex props assigned
                props = PROP_TEMP
                props["coordinates"] = face_triangles[i][1]
                self.vertex_properties[vertices[-2*len(face_triangles)-1 + i]] = props
                # third vertex props assigned
                props = PROP_TEMP
                props["coordinates"] = face_triangles[i][2]
                self.vertex_properties[vertices[-2*len(face_triangles)-1 + (i + 1)]] = props
        else:
            "data_kind invalid or struct_kind invalid"

    def assign_face_props_with_poly_graph(self):
        for face in self.polyhedral_graph.faces:
            props = PROP_TEMP
            triangle = tuple((face[0], face[1], face[2]))
            props["coordinates"] = [self.vertex_properties[face[0]]["coordinates"],
                                    self.vertex_properties[face[1]]["coordinates"],
                                    self.vertex_properties[face[2]]["coordinates"]]
            props["associated_object"] = tpf.Triangle(vertices=props["coordinates"])
            self.face_properties[triangle] = props

    def assign_edge_lengths(self):
        for edge in self.polyhedral_graph.edges:
            props = PROP_TEMP
            props["coordinates"] = [self.vertex_properties[edge.vertices[0]]["coordinates"],
                                    self.vertex_properties[edge.vertices[1]]["coordinates"]]
            props["associated_object"] = edge
            self.edge_properties[tuple((edge.vertices[0], edge.vertices[1]))] = props

