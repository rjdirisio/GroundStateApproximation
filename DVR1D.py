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


def scanPotential(cfg,mdn,reff):
    sub.call(["python","getCarts.py",reff,cfg,mdn]) #'refeck or h9o4_minc
    sub.call(["python","../annes_getpot/toCoord_tet_scan.py","scan/scanned"+cfg+'Mode_'+mdn])
    sub.call(["cp", "coord.dat" ,"../annes_getpot/coord.dat"])
    sub.call(["./../annes_getpot/getpot"])
    sub.call(["cp","eng_dip.dat","DVRPotentials/eng_dip_"+cfg+'_Mode_'+mdn])

def getPotential(fl,xx=None):
    print os.getcwd()
    if xx==None:
        pot = np.loadtxt('DVRPotentials/'+fl)
        if 'eng_dip' in fl:
            pot=pot[:,0]
            if np.amin(pot) < 0:
                pot-=np.amin(pot)
            print len(pot)
            if 'rOH_scan' in fl:
                pot = np.concatenate((np.flip(pot[1001:]),np.array([pot[0]]),pot[1:1000]))
        #plt.plot(pot)
        #plt.show()
        return np.diag(pot)
    else:
        return np.diag((0.5)*k*(xx*xx))

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

#HO
"""#start
grid,dx = initializeGrid()
tcoef = hbar ** 2 / (2 * mass * dx * dx)
V = getPotential('',grid)
T = getKinetic(len(V))
plt.matshow(T)
#plt.show()
plt.close()
plt.matshow(V)
#plt.show()
plt.close()
E,Wf=diagonalizeH(T+V)
for i in range(3):
    plt.plot(grid,Wf[:,i]**2)
plt.show()
plt.close()
print E
print Wf

#rOH
#start
npts = 2001
dx=0.001
grid=np.linspace(-npts*dx,npts*dx,npts)
tcoef = hbar ** 2 / (2 * mass * dx * dx)
V = getPotential('eng_dip_scan_rOH.dat')
T = getKinetic(len(V))
plt.matshow(T)
plt.show()
plt.close()
plt.matshow(V)
plt.show()
plt.close()
E,Wf=diagonalizeH(T+V)
Ecm= E*au2wn
print Ecm
print 'Delta E 0-->1',Ecm[1]-Ecm[0]
for i in range(3):
    plt.plot(grid,Wf[:,i]**2)
plt.show()
plt.close()
#print E
#print Wf"""

#NOrmalScan
mass=1
#refScanGeom = "refEck"
refScanGeom = "h9o4_minc1"

scanPotential(config,modeN,refScanGeom)
V = getPotential("eng_dip_"+config+'_Mode_'+modeN)
npts = len(np.diag(V))
dx=0.01*np.loadtxt("tn/tnorm_"+config+"Mode_"+modeN)
tcoef = hbar ** 2 / (2 * mass * dx * dx)
T = getKinetic(len(V))
E,Wf=diagonalizeH(T+V)
Ecm= E*au2wn
print Ecm
print 'Delta E 0-->1',Ecm[1]-Ecm[0]
np.savetxt("DeltaE/DeltaE01"+refScanGeom+"_"+config+"_"+modeN,np.array([Ecm[1]-Ecm[0]]))

for i in range(3):
    plt.plot(np.diag(V))
    plt.plot(Wf[:,i]**2+(E[i]))
    #plt.ylim([0.0,0.9])
    #plt.ylim([0.0,0.2])
plt.savefig("DVRResults/WavePix"+config+refScanGeom+"_""Mode_"+modeN)
np.savetxt("DVRResults/Wavefunctions_"+config+refScanGeom+"_""Mode_"+modeN,Wf)
np.savetxt("DVRResults/EnergiesCM_"+config+refScanGeom+"_""Mode_"+modeN,Ecm)
plt.show()
plt.close()
