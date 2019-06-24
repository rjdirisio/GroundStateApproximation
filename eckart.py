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
xx=np.load("../coordinates/trimer/final_allH.npy")
xx = np.array(xx[[1273435]])
xx = np.concatenate((xx,xx))
orig = np.copy(xx)

xx[:,1,2]+=1.0e-4
xx = Wfn.molecule.rotateBackToFrame(xx,3,2,1)
com,evec,killL = Wfn.molecule.eckartRotate(xx,justO=True)
xx-=com[:,np.newaxis,:]

print 'done'
print 'b4'
print xx[0]
xx = np.einsum('knj,kij->kni', evec.transpose(0, 2, 1), xx).transpose(0, 2, 1)
print 'af'
print xx[0]
xxp = np.copy(xx)
com, eckVecs, killList=Wfn.molecule.eckartRotate(xxp,hydro=True,yz=True)
xxp-=com[:,np.newaxis,:]
# xx-=com[:,np.newaxis,:]
print 'done'
print 'b4'
print xxp[0]
xxp = np.einsum('knj,kij->kni', eckVecs.transpose(0, 2, 1), xxp).transpose(0, 2, 1)
print 'af'
print xxp[0]


orig[:,1,2]-=1.0e-4
orig = Wfn.molecule.rotateBackToFrame(orig,3,2,1)
com,evec,killL = Wfn.molecule.eckartRotate(orig,justO=True)
orig-=com[:,np.newaxis,:]
print 'done'
print 'b4'
print orig[0]
orig = np.einsum('knj,kij->kni', evec.transpose(0, 2, 1), orig).transpose(0, 2, 1)
print 'af'
print orig[0]
origp = np.copy(orig)
com, eckVecsc, killList=Wfn.molecule.eckartRotate(origp,hydro=True,yz=True)
origp-=com[:,np.newaxis,:]
# orig-=com[:,np.newaxis,:]
print 'done'
print 'b4'
print origp[0]
origp = np.einsum('knj,kij->kni', eckVecsc.transpose(0, 2, 1), origp).transpose(0, 2, 1)
print 'af'
print origp[0]
alleX = np.concatenate((eckVecs.transpose((0,1,2)),eckVecsc.transpose((0,1,2))))
thp,php,xp = Wfn.molecule.extractEulers(alleX)
print thp[0],thp[2]
print php[0],php[2]
print xp[0],xp[2]
print np.degrees([thp[0],thp[2]])
print np.degrees([php[0],php[2]])
print np.degrees([xp[0],xp[2]])
print 'dun'

print (np.abs(php[0]-php[2]) > 1).sum()
print (np.abs(xp[0]-xp[2]) > 1).sum()




ref = Wfn.molecule.pullTrimerRefPos(yz=True)
ref=ref[[3 - 1, 9-1,8-1,10-1]]
mass = Wfn.molecule.get_mass()
mass = mass[[3 - 1, 9-1,8-1,10-1]]
refCOM = np.dot(mass, ref) / np.sum(mass)
ref -= refCOM
xx=xx[:,[3 - 1, 9-1,8-1,10-1]]
xxp=xxp[:,[3 - 1, 9-1,8-1,10-1]]
xx-=xx[:,0,np.newaxis]
xxp-=xxp[:,0,np.newaxis]

origp=origp[:,[3 - 1, 9-1,8-1,10-1]]
origp-=origp[:,0,np.newaxis]

fig = plt.figure()
#ref, xx[0], xx[1]
ax = fig.add_subplot(111, projection='3d')
clist = ['m','b','m','g']
walkerN = 0
for i in range(4):
    ax.scatter(xx[walkerN,i,0], xx[walkerN,i,1], xx[walkerN,i,2], c=clist[i], marker='o')
    if i > 0:
        line = plt3d.art3d.Line3D((xx[walkerN,0,0],xx[walkerN,i,0]), (xx[walkerN,0,1],xx[walkerN,i,1]), (xx[walkerN,0,2],xx[walkerN,i,2]),c='r')
        ax.add_line(line)
clist2 = ['c','b','m','g']
for i in range(4):
    ax.scatter(xxp[walkerN,i,0], xxp[walkerN,i,1], xxp[walkerN,i,2], c=clist2[i], marker='o')
    if i > 0:
        line = plt3d.art3d.Line3D((xxp[walkerN,0,0],xxp[walkerN,i,0]), (xxp[walkerN,0,1],xxp[walkerN,i,1]), (xxp[walkerN,0,2],xxp[walkerN,i,2]),c='b')
        ax.add_line(line)

clist2 = ['g','b','m','g']
for i in range(4):
    ax.scatter(origp[walkerN,i,0], origp[walkerN,i,1], origp[walkerN,i,2], c=clist2[i], marker='o')
    if i > 0:
        line = plt3d.art3d.Line3D((origp[walkerN,0,0],origp[walkerN,i,0]), (origp[walkerN,0,1],origp[walkerN,i,1]), (origp[walkerN,0,2],origp[walkerN,i,2]),c='g')
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