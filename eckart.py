import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as plt3d
import sys
import itertools
import DMCClusters as dmc
Wfn=dmc.wavefunction('H7O3+', 1)
xx=np.load("../coordinates/trimer/test_100H.npy")
orig = np.copy(xx)
com,evec,killL = Wfn.molecule.eckartRotate(xx)
xx-=com[:,np.newaxis,:]
print 'done'
print 'b4'
print xx[0]
xx = np.einsum('knj,kij->kni', evec.transpose(0, 2, 1), xx).transpose(0, 2, 1)
print 'af'
print xx[0]

xxp = np.copy(xx)
com, eckVecs, killList=Wfn.molecule.eckartRotate(xxp,hydro=True)
xxp-=com[:,np.newaxis,:]
# xx-=com[:,np.newaxis,:]
print 'done'
print 'b4'
print xxp[0]
xxp = np.einsum('knj,kij->kni', eckVecs.transpose(0, 2, 1), xxp).transpose(0, 2, 1)
print 'af'
print xxp[0]
ref = np.array(
            [
                 [-2.34906009e+00,  4.06869143e+00, -0.00000000e+00],
                 [ 4.69812018e+00, -0.00000000e+00,  0.00000000e+00],
                 [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                 [-2.88862583e+00,  5.00324669e+00,  1.47532198e+00],
                 [-2.88862583e+00,  5.00324669e+00, -1.47532198e+00],
                 [ 5.77725164e+00, -2.46900000e-09, -1.47532198e+00],
                 [ 5.77725164e+00, -2.46900000e-09,  1.47532198e+00],
                 [-9.12352955e-01, -1.58024167e+00, -0.00000000e+00],
                 [-9.76751990e-01,  1.69178407e+00, -0.00000000e+00],
                 [ 1.95350397e+00, -3.53000000e-09,  0.00000000e+00]])
ref=ref[[3 - 1, 9-1,8-1,10-1]]
xx=xx[:,[3 - 1, 9-1,8-1,10-1]]
xxp=xxp[:,[3 - 1, 9-1,8-1,10-1]]
xx-=xx[:,0,np.newaxis]
xxp-=xxp[:,0,np.newaxis]
fig = plt.figure()
#ref, xx[0], xx[1]
ax = fig.add_subplot(111, projection='3d')
clist = ['m','r','r','r']
walkerN = 90
for i in range(4):
    ax.scatter(xx[walkerN,i,0], xx[walkerN,i,1], xx[walkerN,i,2], c=clist[i], marker='o')
    if i > 0:
        line = plt3d.art3d.Line3D((xx[walkerN,0,0],xx[walkerN,i,0]), (xx[walkerN,0,1],xx[walkerN,i,1]), (xx[walkerN,0,2],xx[walkerN,i,2]),c='r')
        ax.add_line(line)
clist2 = ['c','b','b','b']
for i in range(4):
    ax.scatter(xxp[walkerN,i,0], xxp[walkerN,i,1], xxp[walkerN,i,2], c=clist2[i], marker='o')
    if i > 0:
        line = plt3d.art3d.Line3D((xxp[walkerN,0,0],xxp[walkerN,i,0]), (xxp[walkerN,0,1],xxp[walkerN,i,1]), (xxp[walkerN,0,2],xxp[walkerN,i,2]),c='b')
        ax.add_line(line)
clist3 = ['w','k','k','k']
for i in range(4):
    ax.scatter(ref[i,0], ref[i,1], ref[i,2], c=clist3[i], marker='o')
    if i > 0:
        line = plt3d.art3d.Line3D((ref[0,0],ref[i,0]), (ref[0,1],ref[i,1]), (ref[0,2],ref[i,2]),c='k')
        ax.add_line(line)
ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')
ax.set_zlabel('Z Axis')

plt.show()