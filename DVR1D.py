import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy.linalg as la
import os
import sys
import subprocess as sub
hbar = 1
mass = 1
angstr = 0.529177249
k=1
au2wn=219474.63
config=sys.argv[1]
modeN = sys.argv[2]

def initializeGrid():
    grid = np.linspace(-10,10, 1000)
    dx = grid[1] - grid[0]
    return grid,dx

def getKinetic(dim):
    ke = np.zeros((dim,dim))
    d=tcoef * (np.pi * np.pi / 3.)
    np.fill_diagonal(ke, d)
    for i in range(1,dim):
        for j in range(i):
            ke[i,j]=tcoef*\
                    (-1.)**(i-j)*\
                    (2./((i-j)**2))
    dg = np.diagonal(ke)
    s=ke+ke.T-np.diag(dg)
    return s

def diagonalizeH(H):
    energies,Wfs = la.eigh(H)
    return energies, Wfs

def scanPot(cfg):
    sub.call(['python2','toCoord_new.py',"../scan/scanFiles/scanned_"+cfg+"_"+modeN+".npy"],cwd='annes_getpot')
    potz = np.loadtxt("annes_getpot/eng_dip.dat")[:,0]
    return np.diag(potz)
#NOrmalScan
mass=1

V = scanPot(config)
print(np.diag(V))
npts = len(V)
dx=0.01*np.loadtxt("tn/tnorm_"+config+"Mode_"+modeN)
tcoef = hbar ** 2 / (2 * mass * dx * dx)
T = getKinetic(len(V))
E,Wf=diagonalizeH(T+V)
Ecm= E*au2wn
print Ecm[:5]
print 'Delta E 0-->1',Ecm[1]-Ecm[0]
np.savetxt("DVRResults/deltaE"+config+modeN,np.array([Ecm[1]-Ecm[0]]))

for i in range(3):
    plt.plot(np.diag(V))
    plt.plot(Wf[:,i]**2+(E[i]))
plt.show()
plt.savefig("DVRResults/WavePix"+config+"_""Mode_"+modeN)
np.savetxt("DVRResults/Wavefunctions_"+config+"_""Mode_"+modeN,Wf)
np.savetxt("DVRResults/EnergiesCM_"+config+"_""Mode_"+modeN,Ecm)
plt.show()
plt.close()
