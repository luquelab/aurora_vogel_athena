import numpy as np
import math

#a = np.matrix([[1,2],[3,4]])

#def fold(p1,p2,t):
#    p = p2 - p1
#    x1 = 0
#    y1 = 0
#    z1 = 0
#    x2 = 0
#    y2 = 0
#    z2 = 0

#    u = p/np.linalg.norm(p)
#    a = 0
#    b = 0
#    c = 0


#    T = np.matrix([[1,0,0,-x1],
#                  [0,1,0,-y1],
#                  [0,0,1,-z1],
#                  [0,0,0,1]])
#    Tinv = np.matrix([1,0,0,x1],
#                  [0,1,0,y1],
#                  [0,0,1,z1],
#                  [0,0,0,1])
#    Rx = np.matrix([[1,0,0,0],
#                   [0,c/d,-b/d,0],
#                   [0,b/d,c/d,0],
#                   [0,0,0,1]])
#    Rxinv = np.matrix([[1,0,0,0],
#                       [0,c/d,b/d,0],
#                       [0,-b/d,c/d,0],
#                       [0,0,0,1]])
#    Ry = np.matrix([[d,0,-a,0],
#                    [0,1,0,0],
#                    [a,0,d,0],
#                    [0,0,0,1]])
#    Ryinv = np.matrix([[d,0,a,0],
#                       [0,1,0,0],
#                       [-a,0,d,0],
#                       [0,0,0,1]])
#    Rz = np.matrix([[math.cos(t),math.sin(t),0,0],
#                    [-math.sin(t),math.cos(t),0,0],
#                    [0,0,1,0],
#                    [0,0,0,1]])
#    reduce(np.dot,[Tinv,Rxinv,Ryinv,Rz,Ry,Rx,T,p])