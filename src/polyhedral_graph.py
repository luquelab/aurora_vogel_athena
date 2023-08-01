

class PolyhedralGraph:
    def __init__(self, vertices=[], edges=[], faces=[], kind=""):
        self.vertices = vertices
        for edge in edges:
            if edge.vertices[0] in self.vertices and edge.vertices[1] in self.vertices:
                self.edges = edges
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
                            for k in range(len(ed)):
                                if ed[k] != v:
                                    new_edge_verts.append(ed[k])
                        new_edge = Edge(vertices=new_edge_verts)
                        if len(new_edge) == 2 and (new_edge in self.edges or [new_edge.vertices[-1],
                                                                              new_edge.vertices[0]] in self.edges):
                            face_set.add(tuple(make_triangle(edge, edge_2)))
                        else:
                            print(new_edge)
                            return str(new_edge) + " Edge non-valid or missing in edge list."
            # converting to list
            faces = []
            for face in face_set:
                faces.append(list(face))
            self.faces = faces
        else:
            print("Assigning faces for non-geodesic polyhedral graphs is not supported.")


class Edge:
    def __init__(self, vertices=[]):
        self.vertices = vertices

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
                for v2 in edge_2:
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
        for v in edg:
            vertices.add(v)
    return vertices


# triangles = find_triangles([[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]], True)
# print(triangles)