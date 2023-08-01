import my_library as lib
import math
import numpy as np


class Triangle:
    def __init__(self, vertices=[]):
        self.vertices = vertices
        self.area = self.get_area()

    def get_area(self):
        v1 = self.vertices[0]
        v2 = self.vertices[1]
        v3 = self.vertices[2]
        s1 = np.linalg.norm(v1-v2)
        s2 = np.linalg.norm(v2-v3)
        s3 = np.linalg.norm(v1-v3)
        p = (s1 + s2 + s3)/2  # Half of the perimeter
        a = math.sqrt(
            abs(p*(p-s1)*(p-s2)*(p-s3))
        )  # Area calculated via Heron's formula
        return a


class TrianglePathData:
    def __init__(self, points=[], vector_sequence=[], defect_triangles=[], face_triangles=[], sigma=1):
        self.points = points
        self.vector_sequence = vector_sequence
        self.defect_triangles = defect_triangles
        self.face_triangles = face_triangles
        self.sigma = sigma

    def assign_vector_sequence(self):
        if len(self.points) == 0:
            print("Empty set of points. Cannot calculate vector sequence.")
        else:
            vector_sequence = []
            for i in range(len(points)-1):
                vector_sequence.add(self.points[i+1]-self.points[i])
            self.vector_sequence = vector_sequence

    def assign_points(self):
        if len(self.vector_sequence) == 0:
            print("Empty vector sequence. Cannot calculate points.")
        else:
            points = []
            current_point = 0
            for i in range(len(self.vector_sequence)):
                points.append(current_point)
                current_point = self.vector_sequence[i] + current_point
            self.points = points

    def assign_defect_triangles(self):
        defect_triangles = []
        if len(self.points) == 0:
            self.assign_points()
        if len(self.vector_sequence) == 0:
            self.assign_vector_sequence()
        for i in range(len(self.vector_sequence)):
            triangle_verts = [self.points[i], self.points[i+1],
                              self.points[i+1]-self.sigma*self.vector_sequence[i].rotate_n_pi_3(1)]
            defect_tri = Triangle(vertices=triangle_verts)
            defect_triangles.append(defect_tri)

    def assign_face_triangles(self):
        face_triangles = []
        if len(self.points) == 0:
            self.assign_points()
        if len(self.vector_sequence) == 0:
            self.assign_vector_sequence()
        for i in range(len(self.vector_sequence)):
            triangle_verts = [self.points[i+1], self.points[i+1]-self.sigma*self.vector_sequence[i].rotate_n_pi_3(1),
                              self.points[i+1]+self.sigma*self.vector_sequence[i+1].rotate_n_pi_3(-1)]
            face_tri = Triangle(vertices=triangle_verts)
            face_triangles.append(face_tri)
