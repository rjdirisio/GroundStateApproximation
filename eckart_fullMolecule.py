import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as plt3d
import sys
import itertools
from matplotlib.pyplot import cm
import DMCClusters as dmc
Wfn=dmc.wavefunction('H7O3+', 1)
xx=np.load("../coordinates/trimer/final_allH.npy")
# xx = xx[[0,1,421,541]]
orig = np.copy(xx)
com,evec,killL = Wfn.molecule.eckartRotate(xx)
xx-=com[:,np.newaxis,:]
print 'done'
print 'b4'
print xx[0]
xx = np.einsum('knj,kij->kni', evec.transpose(0, 2, 1), xx).transpose(0, 2, 1)
print 'af'
print xx[0]
thp,php,xp = Wfn.molecule.extractEulers(evec.transpose(0,2,1))
print thp[0],thp[2]
print php[0],php[2]
print xp[0],xp[2]
print np.degrees([thp[0],thp[2]])
print np.degrees([php[0],php[2]])
print np.degrees([xp[0],xp[2]])
print 'dun'

print (np.abs(php[0]-php[2]) > 1).sum()
print (np.abs(xp[0]-xp[2]) > 1).sum()


ref = Wfn.molecule.pullTrimerRefPos(yz=False)
xx-=xx[:,0,np.newaxis]
ref-=ref[0,np.newaxis]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
clist = iter(cm.spectral(np.linspace(0,1,len(xx[0]))))
walkerN = 3

def pltLine(ii,wlkn,oxN,c):
    line = plt3d.art3d.Line3D((xx[walkerN,oxN-1,0],xx[walkerN,ii,0]), (xx[walkerN,oxN-1,1],xx[walkerN,ii,1]), (xx[walkerN,oxN-1,2],xx[walkerN,ii,2]),c=c)
    ax.add_line(line)
def pltRefLine(ii,oxN,c):
    line = plt3d.art3d.Line3D((ref[oxN-1,0],ref[ii,0]), (ref[oxN-1,1],ref[ii,1]), (ref[oxN-1,2],ref[ii,2]),c=c)
    ax.add_line(line)
for i in range(10):
    cc = next(clist)
    ax.scatter(xx[walkerN,i,0], xx[walkerN,i,1], xx[walkerN,i,2], c=cc, marker='o')
    ax.scatter(ref[i,0], ref[i,1], ref[i,2], c=cc, marker='o')
    if i == 4-1 or i == 5-1: #draw waters
        pltLine(i,walkerN,1,'k')
        pltRefLine(i,1,'r')
    elif i== 6-1 or i == 7-1:
        pltLine(i,walkerN,2,'k')
        pltRefLine(i,2,'r')
    elif i == 8-1 or i==9-1 or i==10-1:
        pltLine(i, walkerN, 3, 'k')
        pltRefLine(i, 3, 'r')

ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')
ax.set_zlabel('Z Axis')

plt.show()