import numpy as np
import matplotlib as mpl
# mpl.use('Agg')
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

def scanPot(cfg,mdeN):
    sub.call(['python2', 'toCoord_new.py', "../scan/scanFiles/scanned_" + cfg + "_" + mdeN + ".npy"],
             cwd='annes_getpot')
    potz = np.loadtxt("annes_getpot/eng_dip.dat")[:,0]
    return np.diag(potz)
# NOrmalScan
mass=1
V = np.zeros(201)
V = scanPot(config,modeN)
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
# plt.show()
plt.savefig("DVRResults/WavePix"+config+"_""Mode_"+modeN)
np.savetxt("DVRResults/Wavefunctions_"+config+"_""Mode_"+modeN,Wf)
np.savetxt("DVRResults/EnergiesCM_"+config+"_""Mode_"+modeN,Ecm)
plt.close()


def expVal(psi,attr):
    """Multiplicative"""
    return psi.dot(attr*psi)

gs = np.loadtxt("DVRResults/Wavefunctions_"+config+"_Mode_"+modeN)[:,0]
fEs = np.loadtxt("DVRResults/Wavefunctions_"+config+"_Mode_"+modeN)[:,1]
# v = np.diag(scanPot(config,modeN))
# np.savetxt("DVRResults/Potential_"+config+"_Mode_"+modeN,v)

v = np.loadtxt("DVRResults/Potential_"+config+"_Mode_"+modeN)
dx=0.01*np.loadtxt("tn/tnorm_"+config+"Mode_"+modeN)
grid = np.linspace(0,dx*201,201)
avgR = expVal(gs,grid)
grid -= avgR
q = np.copy(grid)
q2 = q**2
q2p = q*q
q4p = q*q*q*q
q4 = q2**2
v0 = expVal(gs,v)
v1 = expVal(gs,q2*v)/expVal(gs,q2)
alpha =  expVal(gs,q2)/(expVal(gs,q4)-expVal(gs,q2)**2)
deltaT = 1.0**2*1.0*alpha/(2.0*1.0)
freq = deltaT + (v1-v0)
realDT = (-0.5)*np.dot(fEs,np.gradient(np.gradient(fEs,grid),grid)) - (-0.5)*np.dot(gs,np.gradient(np.gradient(gs,grid),grid))
realPE = expVal(fEs,v)-expVal(gs,v)
freqq = realDT+realPE
print('hi')
# plt.plot(grid, (Wf[:, 0] ** 2))
# plt.ylim([0,0.045])
# plt.plot(grid, (gs ** 2)*(grid**4))
# plt.plot(grid, (Wf[:, 0] ** 2)*(grid**4))
# plt.show()
print("<q2>",expVal(gs,q2))
print("<q4>",expVal(gs,q4))
print(alpha)
print(deltaT*au2wn)
print((v1-v0)*au2wn)
print('shwomp')
print(realDT*au2wn)
print(realPE*au2wn)
print((realDT+realPE)*au2wn)
print('real')
print(Ecm[1]-Ecm[0])
np.savetxt("DVRResults/GSADATA_" + config + "_""Mode_" + modeN,
           zip([expVal(gs, q2),
                expVal(gs, q2) ** 2,
                expVal(gs, q4),
                alpha,
                deltaT * au2wn,
                (v1 - v0) * au2wn,
                freq * au2wn,
                realDT * au2wn,
                realPE * au2wn,
                Ecm[1] - Ecm[0]]))

# print(v0*au2wn)
# print(v1*au2wn)
# print(au2wn*freq)