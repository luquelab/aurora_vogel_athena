import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import pandas as pd
import math

########################################################################
# Distance function with hardcoded triangle side length 2 (for testing)
########################################################################
def distance_sum(verts):
    v04=vs[0:3]
    v05=vs[3:6]
    v06=vs[6:9]
    v07=vs[9:12]
    v08=vs[12:15]
    v09=vs[15:18]
    v10=vs[18:21]
    v11=vs[21:24]
    v12=vs[24:27]

    eq1 = np.array((v03[0] - v04[0])**2 + (v03[1] - v04[1])**2 + (v03[2] - v04[2])**2)-2**2
    eq2 = np.array((v04[0] - v05[0])**2 + (v04[1] - v05[1])**2 + (v04[2] - v05[2])**2)-2**2
    eq3 = np.array((v05[0] - v06[0])**2 + (v05[1] - v06[1])**2 + (v05[2] - v06[2])**2)-2**2
    eq4 = np.array((v06[0] - v02[0])**2 + (v06[1] - v02[1])**2 + (v06[2] - v02[2])**2)-2**2
    # cnxn to top rim
    eq5 = np.array((v04[0] - v01[0])**2 + (v04[1] - v01[1])**2 + (v04[2] - v01[2])**2)-2**2
    eq6 = np.array((v05[0] - v01[0])**2 + (v05[1] - v01[1])**2 + (v05[2] - v01[2])**2)-2**2
    eq7 = np.array((v06[0] - v01[0])**2 + (v06[1] - v01[1])**2 + (v06[2] - v01[2])**2)-2**2

    ## body cnxn to top
    eq8 = np.array((v07[0] - v02[0])**2 + (v07[1] - v02[1])**2 + (v07[2] - v02[2])**2)-2**2
    eq9 = np.array((v08[0] - v02[0])**2 + (v08[1] - v02[1])**2 + (v08[2] - v02[2])**2)-2**2
    eq10 = np.array((v08[0] - v03[0])**2 + (v08[1] - v03[1])**2 + (v08[2] - v03[2])**2)-2**2
    eq11 = np.array((v09[0] - v03[0])**2 + (v09[1] - v03[1])**2 + (v09[2] - v03[2])**2)-2**2
    eq12 = np.array((v09[0] - v04[0])**2 + (v09[1] - v04[1])**2 + (v09[2] - v04[2])**2)-2**2
    eq13 = np.array((v10[0] - v04[0])**2 + (v10[1] - v04[1])**2 + (v10[2] - v04[2])**2)-2**2
    eq14 = np.array((v10[0] - v05[0])**2 + (v10[1] - v05[1])**2 + (v10[2] - v05[2])**2)-2**2
    eq15 = np.array((v11[0] - v05[0])**2 + (v11[1] - v05[1])**2 + (v11[2] - v05[2])**2)-2**2
    eq16 = np.array((v11[0] - v06[0])**2 + (v11[1] - v06[1])**2 + (v11[2] - v06[2])**2)-2**2
    eq17 = np.array((v07[0] - v06[0])**2 + (v07[1] - v06[1])**2 + (v07[2] - v06[2])**2)-2**2

    # bottom rim
    eq18 = np.array((v07[0] - v08[0])**2 + (v07[1] - v08[1])**2 + (v07[2] - v08[2])**2)-2**2
    eq19 = np.array((v08[0] - v09[0])**2 + (v08[1] - v09[1])**2 + (v08[2] - v09[2])**2)-2**2
    eq20 = np.array((v09[0] - v10[0])**2 + (v09[1] - v10[1])**2 + (v09[2] - v10[2])**2)-2**2
    eq21 = np.array((v10[0] - v11[0])**2 + (v10[1] - v11[1])**2 + (v10[2] - v11[2])**2)-2**2
    eq22 = np.array((v11[0] - v07[0])**2 + (v11[1] - v07[1])**2 + (v11[2] - v07[2])**2)-2**2

    # Cnxn to bottom rim
    eq23 = np.array((v07[0] - v12[0])**2 + (v07[1] - v12[1])**2 + (v07[2] - v12[2])**2)-2**2
    eq24 = np.array((v08[0] - v12[0])**2 + (v08[1] - v12[1])**2 + (v08[2] - v12[2])**2)-2**2
    eq25 = np.array((v09[0] - v12[0])**2 + (v09[1] - v12[1])**2 + (v09[2] - v12[2])**2)-2**2
    eq26 = np.array((v10[0] - v12[0])**2 + (v10[1] - v12[1])**2 + (v10[2] - v12[2])**2)-2**2
    eq27 = np.array((v11[0] - v12[0])**2 + (v11[1] - v12[1])**2 + (v11[2] - v12[2])**2)-2**2

    # min'd around origin
    mast_eq = np.abs(eq1)+np.abs(eq2)+np.abs(eq3)+np.abs(eq4)+np.abs(eq5)+np.abs(eq6)+np.abs(eq7)+np.abs(eq8)+np.abs(eq9)+np.abs(eq10)+np.abs(eq11)+np.abs(eq12)+np.abs(eq13)+np.abs(eq14)+np.abs(eq15)+np.abs(eq16)+np.abs(eq17)+np.abs(eq18)+np.abs(eq19)+np.abs(eq20)+np.abs(eq21)+np.abs(eq22)+np.abs(eq23)+np.abs(eq24)+np.abs(eq25)+np.abs(eq26)+np.abs(eq27)
    return mast_eq

########################################################################
# Initializing stuff
########################################################################

# golden ratio
phi=(1+math.sqrt(5))/2

# regular icosahedron with triangle side length of 2
# with vertex 4 nudged from 0,-1,phi to 0.25,-1,phi
init_verts = np.array([0.25,-1,phi, # v04
               -phi,0,1, # v05
               -1,phi,0, # v06
               0,1,-phi, # v07
               phi,0,-1, # v08
               1,-phi,0, # v09
               -1,-phi,0, # v10
               -phi,0,-1, # v11
               0,-1,-phi, # v12
               ])

v01=np.array([0,1,phi])
v02=np.array([1,phi,0])
v03=np.array([phi,0,1])

# minimize to get optimal vertices, "new_verts"
new_verts_data=minimize(distance_sum, init_verts)
new_verts=new_verts.x

# assigning each vert with a variable
v04=new_verts[0:3]
v05=new_verts[3:6]
v06=new_verts[6:9]
v07=new_verts[9:12]
v08=new_verts[12:15]
v09=new_verts[15:18]
v10=new_verts[18:21]
v11=new_verts[21:24]
v12=new_verts[24:27]

# placing all verts together to form a dataframe which can be plotted
DATA=pd.DataFrame([v01,v02,v03,v04,v05,v06,v07,v08,v09,v10,v11,v12],columns=["X","Y","Z"])

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(DATA['X'], DATA['Y'], DATA['Z'])

# labeling the verts plotted
n=[1,2,3,4,5,6,7,8,9,10,11,12]
for i, txt in enumerate(n):
    ax.text(DATA['X'][i], DATA['Y'][i], DATA['Z'][i],'%s' % (str(i)), size=8, zorder=1,  
    color='k')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()