import potential_unfolding_creator as puc
import numpy as np


class Polyhedron:
    def __init__(self, vertices: list, edges: list, faces: list, unfolding: puc.PotentialSpacedUnfolding,
                 glue_instructs: list, vertex_properties={}, edge_properties={}, face_properties={}):
        self.vertices = vertices
        self.edges = edges
        self.faces = faces
        self.vertex_properties = vertex_properties
        self.edge_properties = edge_properties
        self.face_properties = face_properties
        self.glue_instructs = glue_instructs
        self.unfolding = unfolding
        if unfolding:
            self.populate_from_unfolding()
        if glue_instructs:
            self.glue_structure()

    def populate_from_unfolding(self):
        self.vertices = self.unfolding.polyhedral_graph.vertices
        self.edges = self.unfolding.polyhedral_graph.edges
        self.faces = self.unfolding.polyhedral_graph.faces
        self.vertex_properties = self.unfolding.vertex_properties
        self.edge_properties = self.unfolding.edge_properties
        self.face_properties = self.unfolding.face_properties

    def glue_structure(self):
        new_verts = []
        for vertex in self.vertices:
            pass

    def get_delaunay_triangulation(self, face1_key, face2_key):
        face1 = self.face_properties[face1_key]["associated_object"]
        face2 = self.face_properties[face2_key]["associated_object"]
        shared_edge = []
        not_shared = []
        for vertex in face1.vertices:
            if vertex in face2.vertices:
                shared_edge.append(vertex)
            else:
                not_shared.append(vertex)
        for vertex in face2.vertices:
            if vertex not in face1.vertices:
                not_shared.append(vertex)
        # we get gamma 1 and gamma 2
        gammas = [face1.get_angle(shared_edge[0], not_shared[0], shared_edge[1]),
                  face1.get_angle(shared_edge[0], not_shared[1], shared_edge[1])]
        if gammas[0] + gammas[1] <= puc.math.pi:
            pass
        else:
            self.edge_swap(puc.Edge(vertices=face1.get_labels(shared_edge)),
                           puc.Edge(vertices=(face1.get_labels(not_shared[0])+face2.get_labels(not_shared[1]))),
                           face1_key, face2_key, not_shared)

    def edge_swap(self, edge_old, edge_new, face1_key, face2_key, not_shared):
        edge_old_key = tuple(edge_old.vertices.copy())
        edge_new_key = tuple(edge_new.vertices.copy())
        self.edge_properties[edge_new_key] = self.edge_properties[edge_old_key]
        self.edge_properties[edge_new_key]["coordinates"] = [self.vertex_properties[edge_new.vertices[0]]["coordinates"],
                                                          self.vertex_properties[edge_new.vertices[1]]["coordinates"]]
        self.edge_properties[edge_new_key]["associated_object"] = edge_new
        del self.edge_properties[edge_old_key]
        new_edges = [edge_new]
        for edge in self.edges:
            if edge.vertices() == edge_old.vertices():
                pass
            else:
                new_edges.append(edge)
        self.edges = new_edges.copy()

        new_face1_key = tuple(not_shared[0]) + edge_new_key
        new_face2_key = tuple(not_shared[1]) + edge_new_key
        self.face_properties[new_face1_key] = self.face_properties[face1_key]
        self.face_properties[new_face1_key]["coordinates"] = [self.vertex_properties[new_face1_key[0]]["coordinates"],
                                                              self.vertex_properties[new_face1_key[1]]["coordinates"],
                                                              self.vertex_properties[new_face1_key[2]]["coordinates"]]
        self.face_properties[new_face1_key]["associated_object"] = \
            puc.tpf.Triangle(vertices=self.face_properties[new_face1_key]["coordinates"], labels=new_face1_key)
        self.face_properties[new_face2_key] = self.face_properties[face2_key]
        self.face_properties[new_face2_key]["coordinates"] = [self.vertex_properties[new_face2_key[0]]["coordinates"],
                                                              self.vertex_properties[new_face2_key[1]]["coordinates"],
                                                              self.vertex_properties[new_face2_key[2]]["coordinates"]]
        self.face_properties[new_face2_key]["associated_object"] = \
            puc.tpf.Triangle(vertices=self.face_properties[new_face2_key]["coordinates"], labels=new_face2_key)
        del self.face_properties[face1_key]
        del self.face_properties[face2_key]

