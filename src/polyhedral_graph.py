

class PolyhedralGraph:
    def __init__(self, vertices=[], edges=[], faces=[], filt=[], kind="geodesic"):
        self.filt = filt
        if len(vertices) == 0:
            if type(edges[0]) == list:
                vertex_set = set()
                for edge in edges:
                    vertex_set.add(edge[0])
                    vertex_set.add(edge[1])
                self.vertices = list(vertex_set)
            else:
                "No implementation for populating vertices from list of non-lists"
        else:
            self.vertices = vertices
        if type(edges[0]) != Edge:
            new_edges = []
            for edge in edges:
                new_edges.append(Edge(vertices=edge))
            self.edges = new_edges
        #
        # Checking if vertex list compatible with edge list
        for edge in self.edges:
            if edge.vertices[0] in self.vertices and edge.vertices[1] in self.vertices:
                pass
            else:
                print("Invalid edge list.")
                self.edges = []
        for face in faces:
            if face[0] in self.vertices and face[1] in self.vertices and face[2] in self.vertices:
                self.faces = faces
            else:
                print("Invalid face list.")
                self.faces = []
        self.kind = kind

        if len(edges) > 0:
            self.assign_faces()

    def assign_faces(self):
        face_set = set()
        if self.kind == "geodesic":
            for i in range(len(self.edges) - 1):
                edge = self.edges[i]
                for j in range(i + 1, len(self.edges)):
                    edge_2 = self.edges[j]
                    if edge.has_common_vertex(edge_2=edge_2, retur=False):
                        v = edge.has_common_vertex(edge_2=edge_2, retur=True)
                        new_edge_verts = []
                        for ed in [edge, edge_2]:
                            for k in range(len(ed.vertices)):
                                if ed.vertices[k] != v:
                                    new_edge_verts.append(ed.vertices[k])
                        new_edge = Edge(vertices=new_edge_verts)
                        if len(new_edge.vertices) == 2 and (new_edge.is_in(self.edges) or Edge(vertices=[new_edge.vertices[-1],
                                                                              new_edge.vertices[0]]).is_in(self.edges)):
                            face_set.add(tuple(make_triangle(edge, edge_2)))
                        else:
                            pass
                            # print(str(new_edge) + " Edge non-valid or missing in edge list.")
            # converting to list
            faces = []
            for face in face_set:
                faces.append(list(face))
            self.faces = faces
        else:
            print("Assigning faces for non-geodesic polyhedral graphs is not supported.")

    def get_degrees(self):
        degrees = []
        for vertex in self.vertices:
            current_degree = 0
            for edge in self.edges:
                if vertex in edge.vertices and edge.vertices not in self.filt:
                    current_degree += 1
            degrees.append(current_degree)
        return degrees

    def __repr__(self):
        return "Vertices: " + str(self.vertices) + "\n" + "Edges: " + str(self.edges) + "\n" + "Vertices: " + str(self.faces)


class Edge:
    def __init__(self, vertices=[]):
        self.vertices = vertices

    def __repr__(self):
        return str(self.vertices)

    def is_in(self, edge_list):
        for edge in edge_list:
            if self.vertices == edge.vertices:
                return True
        return False

    def has_common_vertex(self, edge_2, vertex=None, retur=False):
        if vertex:
            for v2 in edge_2.vertices:
                if vertex - v2 == 0:
                    if retur:
                        return v2
                    else:
                        return True
                else:
                    if retur:
                        return None
                    else:
                        return False
        else:
            for v in self.vertices:
                for v2 in edge_2.vertices:
                    if v - v2 == 0:
                        if retur:
                            return v2
                        else:
                            return True
            if retur:
                return "No common vertex found."
            else:
                return False
#
#
# def find_triangles(edges, as_list=False):
#     triangles = set()
#     for i in range(len(edges)-1):
#         edge = edges[i]
#         for j in range(i+1, len(edges)):
#             edge_2 = edges[j]
#             if has_common_vertex(edge=edge, edge_2=edge_2, retur=False):
#                 v = has_common_vertex(edge=edge, edge_2=edge_2, retur=True)
#                 new_edge = []
#                 for ed in [edge, edge_2]:
#                     for i in range(len(ed)):
#                         if ed[i] != v:
#                             new_edge.append(ed[i])
#                 if len(new_edge) == 2 and (new_edge in edges or [new_edge[-1], new_edge[0]] in edges):
#                     triangles.add(tuple(make_triangle(edge, edge_2)))
#                 else:
#                     print(new_edge)
#                     return "Edge non-valid or missing in edge list."
#     if as_list:
#         new_triangles = []
#         for triangle in triangles:
#             new_triangles.append(list(triangle))
#         return new_triangles
#     else:
#         return triangles
#
#
# def has_common_vertex(edge, edge_2, vertex=None, retur=False):
#     if vertex:
#         for v2 in edge_2:
#             if vertex-v2 == 0:
#                 if retur:
#                     return v2
#                 else:
#                     return True
#     else:
#         for v in edge:
#             for v2 in edge_2:
#                 if v-v2 == 0:
#                     if retur:
#                         return v2
#                     else:
#                         return True
#         if retur:
#             return "No common vertex found."
#         else:
#             return False
#
#


def make_triangle(edge, edge_2):
    vertices = set()
    for edg in [edge, edge_2]:
        for v in edg.vertices:
            vertices.add(v)
    return vertices


# triangles = find_triangles([[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]], True)
# print(triangles)
# pg = PolyhedralGraph(edges=[[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]])
# pg = PolyhedralGraph(edges=[[0, 3], [0, 4], [1, 4], [1, 5], [2, 5], [2, 6], [3, 4], [4, 5], [5, 6],
#                             [3, 7], [4, 7], [4, 8], [5, 8], [5, 9], [6, 9]])

# filt1 = [[0, 5], [1, 6], [2, 7], [3, 8], [4, 9],
#          [5, 11],
#          [11, 17], [12, 18], [13, 19], [14, 20], [15, 21]]
#
# pgraph = PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
#                                 [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
#                                 [5, 11], [5, 12], [6, 12], [6, 13], [7, 13], [7, 14], [8, 14], [8, 15], [9, 15], [9, 16], [10, 16],
#                                 [11, 12], [12, 13], [13, 14], [14, 15], [15, 16],
#                                 [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
#                          filt=filt1)
