import my_library as lib
import math
import numpy as np


class Triangle:
    def __init__(self, vertices=[]):
        self.vertices = vertices
        self.area = self.get_area()

    def get_area(self):
        if type(self.vertices[0]) == lib.HexagonalPoint:
            v1 = self.vertices[0].get_cartesian()
            v2 = self.vertices[1].get_cartesian()
            v3 = self.vertices[2].get_cartesian()
        v1 = np.array((v1.x, v1.y))
        v2 = np.array((v2.x, v2.y))
        v3 = np.array((v3.x, v3.y))

        s1 = np.linalg.norm(v1-v2)
        s2 = np.linalg.norm(v2-v3)
        s3 = np.linalg.norm(v1-v3)
        p = (s1 + s2 + s3)/2  # Half of the perimeter
        a = math.sqrt(
            abs(p*(p-s1)*(p-s2)*(p-s3))
        )  # Area calculated via Heron's formula
        return a

    def __repr__(self):
        return "["+str(self.vertices[0])+","+str(self.vertices[1])+","+str(self.vertices[2])+"]"

    def __getitem__(self, i):
        return self.vertices[i]


class TrianglePathData:
    def __init__(self, points=[], vector_sequence=[], defect_triangles=[], face_triangles=[], sigma=1):
        self.points = points
        self.vector_sequence = vector_sequence
        self.defect_triangles = defect_triangles
        self.face_triangles = face_triangles
        self.sigma = sigma

        if len(self.vector_sequence) > 0 or len(self.points) > 0:
            if len(self.vector_sequence) > 0:
                new_vector_sequence = []
                for vector in self.vector_sequence:
                    new_vector_sequence.append(lib.HexagonalPoint(h=vector[0], k=vector[1], l=vector[2]))
                self.vector_sequence = new_vector_sequence
                self.assign_points()
            else:
                new_points = []
                for point in self.points:
                    new_points.append(lib.HexagonalPoint(h=vector[0], k=vector[1], l=vector[2]))
                self.points = new_points
                self.assign_vector_sequence()
            self.assign_defect_triangles()
            self.assign_face_triangles()

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
            current_point = lib.HexagonalPoint(h=0.0, k=0.0)
            for i in range(len(self.vector_sequence)):
                points.append(current_point)
                current_point = self.vector_sequence[i] + current_point
            points.append(current_point)
            self.points = points

    def assign_defect_triangles(self):
        defect_triangles = []
        if len(self.points) == 0:
            self.assign_points()
        if len(self.vector_sequence) == 0:
            self.assign_vector_sequence()
        for i in range(len(self.vector_sequence)):
            triangle_verts = [self.points[i], self.points[i+1],
                              self.points[i+1]-self.vector_sequence[i].rotate_n_pi_3(self.sigma)]
            defect_tri = Triangle(vertices=triangle_verts)
            defect_triangles.append(defect_tri)
        self.defect_triangles = defect_triangles

    def assign_face_triangles(self):
        face_triangles = []
        if len(self.points) == 0:
            self.assign_points()
        if len(self.vector_sequence) == 0:
            self.assign_vector_sequence()
        for i in range(len(self.vector_sequence)-1):
            triangle_verts = [self.points[i+1], self.points[i+1]-self.vector_sequence[i].rotate_n_pi_3(self.sigma),
                              self.points[i+1] + self.vector_sequence[i+1].rotate_n_pi_3(-self.sigma)]
            face_tri = Triangle(vertices=triangle_verts)
            face_triangles.append(face_tri)
        self.face_triangles = face_triangles

# dat = TrianglePathData(vector_sequence=[[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0]])