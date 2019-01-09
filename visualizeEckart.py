import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

mpl.rcParams['legend.numpoints']=1
o=np.loadtxt("OxygensOf4")
ox = o.reshape(3,4,3)
ox=np.round(ox,8)
cls=['b','r','g']
lbl=['ReferenceGeometry','OriginalWalker','EckartedWalker']
#tetref,xx,xxnew
a=0
fig = plt.figure()
#ax = fig.gca(projection='3d')
plt.subplot(121)
for set in ox:
    for atom in set:
        plt.scatter(atom[0],atom[1],60,facecolors='none',edgecolors=cls[a],label=lbl[a])
    # com = (set[0]+set[1]+set[2]) / 3
    # plt.scatter(com[0],com[1],color=cls[a],label='com')
    #     #ax.scatter(atom[0],atom[1],atom[2],label=lbl[a],c=cls[a])
    a+=1
#ax.legend()
plt.legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0)
plt.savefig('eckart_n')
plt.show()
