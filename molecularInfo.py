import numpy as np
import sys
from numpy import linalg as la
import os
import time
import Plot
import copy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from multiprocessing import Process, Pipe,Pool
from itertools import izip

def yaxcyc(mptshxax):
    mptshxax=mptshxax.reshape(3,3)
    yaxis = np.cross(mptshxax[-1], (0, 0, 1))
    xcomp = np.dot(mptshxax[1] - mptshxax[0], mptshxax[-1].T)
    ycomp = np.dot(mptshxax[1] - mptshxax[0], yaxis.T)
    return yaxis, xcomp, ycomp

def yaxcyc2(mptshxax):
    mptshxax=mptshxax.reshape(3,3)
    yaxis = np.cross(mptshxax[-1], (0, 0, 1))
    xcomp = np.dot(mptshxax[1] - mptshxax[0], mptshxax[-1].T) #h8-o3 . xax
    ycomp = np.dot(mptshxax[1] - mptshxax[0], yaxis.T)
    arct = np.arctan(ycomp/xcomp)
    return arct, xcomp, ycomp

def yaxcyc3(xz4):
    xz4 = xz4.reshape(3, 3)
    yaxis = np.cross(xz4[1],xz4[0]) #z cross x
    xcomp = np.dot(xz4[-1], xz4[0].T)
    ycomp = np.dot(xz4[-1], yaxis.T)
    return xcomp, ycomp

#getCompz = np.column_stack((mp, sharedH, xaxis,yaxis,zaxis))


def simpleCross(XZ):
    x=XZ[0:3]
    z = XZ[3:]
    cross = np.cross(x,z)
    return cross

def euler(zyxYX):
    z=zyxYX[0:3]
    y=zyxYX[3:6]
    x=zyxYX[6:9]
    Y=zyxYX[9:12]
    X=zyxYX[12:15]
    Z = np.array([0,0,1])
    #print 'finishThis'
    Theta = np.arccos(np.dot(z, Z) / (la.norm(z) * la.norm(Z.T)))
    if np.isnan(Theta):
        stop
    #tanPhi = np.arctan2(np.dot(y, Z.T), np.dot(x, Z.T))
    #tanChi = np.arctan2(-1 * np.dot(Y, z.T), np.dot(X, z.T))
    tanPhi = np.arctan2(np.dot(y, Z.T) /(la.norm(y)*la.norm(Z.T)), np.dot(x, Z.T)/(la.norm(x)*la.norm(Z.T)))
    tanChi = np.arctan2(-1 * np.dot(Y, z.T)/(la.norm(Y)*la.norm(z.T)), np.dot(X, z.T)/(la.norm(X)*la.norm(z.T)))
    if tanChi< 0:
        tanChi += 2 * np.pi
    if tanPhi < 0:
        tanPhi += 2 * np.pi
    return Theta,tanPhi,tanChi

def rotate1(cds):
    o2 = cds[-3:]
    z = o2[2]
    y = o2[1]
    x = o2[0]
    walkerCoords = cds[:-3].reshape(13,3)
    theta = np.arctan2(-1 * z, y)
    alpha = np.arctan2((-1 * (
            y * np.cos(theta) - np.sin(theta) * z)), x)
    r1 = np.array([[1, 0, 0],
                   [0, np.cos(theta), -1 * np.sin(theta)],
                   [0, np.sin(theta), np.cos(theta)]]
                  )
    r2 = np.array([[np.cos(alpha), -1 * np.sin(alpha), 0],
                   [np.sin(alpha), np.cos(alpha), 0],
                   [0, 0, 1]]
                  )
    rotM = np.dot(r2, r1)
    walker = rotM.dot(walkerCoords.T).T
    return walker

def rotate2(car):
    o1 = car[-3:]
    walk = car[:-3].reshape(13,3)
    z = o1[-1]
    y = o1[1]
    #x = o1[0]
    beta = np.arctan2(-1 * z, y)
    r = np.array([[1, 0, 0],
                  [0, np.cos(beta), -1 * np.sin(beta)],
                  [0, np.sin(beta), np.cos(beta)]]
                 )
    walker = r.dot(walk.T).T
    return walker

def bondLen(stuf):
    atm1 = np.array(stuf[:3])
    atm2 = np.array(stuf[3:])
    length = la.norm(atm1-atm2)
    return length

def bondAng(three):
    atm1 = np.array(three[:3])
    atm2 = np.array(three[3:6])
    atm3 = np.array(three[6:])
    left = atm1-atm2
    right = atm3-atm2
    g=la.norm(left)
    h=la.norm(right)
    i = np.dot(left,right)
    k = i/(g*h)
    j= np.arccos(i/(g*h))
    return j

def normit(OA):
    nOA = OA/la.norm(OA)
    return nOA
"""finalCoords = np.zeros(np.shape(xaxtrCoordz))
o1 = xaxtrCoordz[:, c - 1].reshape(numWalkers, 1, 3)
beta = np.zeros(numWalkers)
for walker in range(numWalkers):
    z = o1[walker, 0, 2]
    y = o1[walker, 0, 1]
    x = o1[walker, 0, 0]
    beta[walker] = np.arctan2(-1 * z, y)
    r = np.array([[1, 0, 0],
                  [0, np.cos(beta[walker]), -1 * np.sin(beta[walker])],
                  [0, np.sin(beta[walker]), np.cos(beta[walker])]]
                 )
    finalCoords[walker] = np.dot(r, xaxtrCoordz[walker].T).T
    # naa=0"""

global massH
global massO
global massD
massH=1.00782503223
massD=2.0141017778
massO=15.99491561957
massConversionFactor=1.000000000000000000/6.02213670000e23/9.10938970000e-28#1822.88839  g/mol -> a.u.
ang2bohr=1.88973
bohr2ang=1.000/ang2bohr
rad2deg=180.000/np.pi
au2wn=219474.63
global Water
Water={'H2O', 'water', 'HOH', 'h2o'}
global WaterDimer
WaterDimer = {'water dimer', 'H4O2'}
global ProtonatedWaterDimer
ProtonatedWaterDimer = {'H5O2+','HHOHOHH','H5O2+','h5o2plus','h5o2'}
global DeprotonatedWaterDimer
DeprotonatedWaterDimer = {'HOHOH','H3O2-', 'h3o2-', 'h3o2', 'H3O2'}
global ProtonatedWaterTrimer
ProtonatedWaterTrimer = {'H7O3+','O3H7+', 'H7O3plus','H7O3', 'O3H7'}
global ProtonatedWaterTetramer
ProtonatedWaterTetramer = {'H9O4+','O4H9+', 'H9O4plus','H9O4', 'O4H9'}

surfaceOptions= {'OHStretchAnti','StretchAntiIn','SharedProton','LocalOHStretch','Z-displacement','EckRotZ-displacement'}
class molecule (object):
    def __init__(self,moleculeName):
        self.name=moleculeName

        if self.name in DeprotonatedWaterDimer:
            self.nAtoms=5
            self.names=['H','O','H','O','H']
        elif self.name in ProtonatedWaterDimer:
            self.natoms=7
            self.names=['O','O','H','H','H','H','H']
        elif self.name in ProtonatedWaterTrimer:
            self.nAtoms = 10
            self.names=['O','O','O','H','H','H','H','H','H','H']

        elif  self.name in ProtonatedWaterTetramer:
            self.nAtoms = 13
            self.names = ['O', 'O', 'O', 'O','H', 'H', 'H', 'H', 'H', 'H', 'H','H','H']
        elif self.name in Water:
            self.nAtoms = 3
            self.names = ['O','H','H']

        if self.name in ProtonatedWaterDimer:
            self.potential=self.getPES()
            self.dipole=self.getDIPOLE()
        self.isotope='notDeuterated'
        self.nVibs=3*self.nAtoms-6

    def set_isotope(self,keyword):
        print 'resetting isotope to ', keyword
        if self.name in ProtonatedWaterTrimer:
            self.isotope = keyword
            if self.isotope == 'DeuteratedOnce_eigen':
                self.names = ['O','O','O','H','H','H','H','H','H','D']
            elif self.isotope == 'DeuteratedOnce_fw':
                self.names = ['O','O','O','D','H','H','H','H','H','H']
            elif self.isotope == 'notDeuterated':
                self.names = ['O','O','O','H','H','H','H','H','H','H']
            elif self.isotope == 'fullyDeuterated':
                self.names = ['O','O','O','D','D','D','D','D','D','D']
            elif self.isotope == 'notDeuteratedOnce_eigen':
                self.names = ['O','O','O','D','D','D','D','D','D','H']
            elif self.isotope == 'notDeuteratedOnce_fw':
                self.names = ['O','O','O','H','D','D','D','D','D','D']
            elif self.isotope == 'DeuteratedOnce_hydronium':
                self.names = ['O','O','O','H','H','H','H','D','H','H']
            elif self.isotope == 'notDeuteratedOnce_hydronium':
                self.names = ['O', 'O', 'O', 'D', 'D', 'D', 'D', 'H', 'D', 'H']

        elif self.name in ProtonatedWaterTetramer:
            self.isotope = keyword
            if self.isotope == 'DeuteratedOnce_eigen':
                self.names = ['O', 'O', 'O', 'O','H', 'H', 'H', 'H', 'H', 'H', 'D','H','H']
            elif self.isotope == 'DeuteratedOnce_fw':
                self.names = ['O', 'O', 'O', 'O','D', 'H', 'H', 'H', 'H', 'H', 'H','H','H']
            elif self.isotope == 'notDeuterated':
                self.names = ['O', 'O', 'O', 'O','H', 'H', 'H', 'H', 'H', 'H', 'H','H','H']
            elif self.isotope == 'fullyDeuterated':
                self.names = ['O', 'O', 'O', 'O','D', 'D', 'D', 'D', 'D', 'D', 'D','D','D']
            elif self.isotope == 'notDeuteratedOnce_eigen':
                self.names = ['O', 'O', 'O', 'O','D', 'D', 'D', 'D', 'D', 'D', 'H','D','D']
            elif self.isotope == 'notDeuteratedOnce_fw':
                self.names = ['O', 'O', 'O', 'O','H', 'D', 'D', 'D', 'D', 'D', 'D','D','D']

    def setNodalSurface(self,surfaceName,side):
        self.surfaceName=surfaceName
        if not self.surfaceName in surfaceOptions:
            print "THAT IS NOT A SURFACE! you have likely made a typo."
        self.side=side
        self.state=1

    def getPES(self):
        sys.path.insert(0,self.pathDict["potentialPath"+self.name])
        print self.pathDict["potentialPath"+self.name]
        import pes
        if self.name in DeprotonatedWaterDimer:
            pes.prepot()
            potential=pes.getpot
            print 'potential retrieved for HOHOH. Be sure to feed one walker in at a time!'
        elif self.name in ProtonatedWaterTrimer:
            print 'idk lmao'
            #potential =
            return potential


    def getDIPOLE(self):
        sys.path.insert(0,self.pathDict["dipolePath"+self.name])
        import h3o2dms2 as dms
        if self.name in DeprotonatedWaterDimer:
            dms.predip()
            dip=dms.mycalcdip
        #since all dipole moment calculations need eckart rotation...loading the reference coordinates
        self.refPos=self.loadRefCoord()
        return dip

    def loadDict(self,fileName):
        fileIn=open(fileName,'r')
        pathDict={}
        for line in fileIn:
            [key,element]=line.split()
            pathDict[key]=element
        return pathDict
    def loadRefCoord(self):
        if self.name in DeprotonatedWaterDimer:
            coordsMin=np.array([[ 0.2981678882048853 ,   -2.4557992072743176E-002,  -5.5485232545510215E-002],
                                [ -2.354423404994569 ,     0.000000000000000     ,    0.000000000000000],
                                [ -2.858918674095194 ,     1.111268022307282     ,   -1.352651141853729],
                                [  2.354423404994569 ,     0.000000000000000     ,    0.000000000000000],
                                [  2.671741580470489 ,     1.136107563104921     ,    1.382886181959795]])


            coordsSADDLE=np.array([[ 0.000000000000000 , 0.000000000000000 , 0.000000000000000],
                                   [ -2.303263755760085, 0.000000000000000 ,0.000000000000000 ],
                                   [ -2.720583162407882, 1.129745554266140 ,-1.363735721982301   ],
                                   [ 2.303263755760085, 0.000000000000000, 0.000000000000000     ],
                                   [ 2.720583162407882, 1.129745554266140, 1.363735721982301]])


        return coordsSADDLE

    def getEquilibriumEnergy(self):
        if self.name in DeprotonatedWaterDimer:
            equilibriumCoords=np.array([[ 0.2981678882048853 , -2.4557992072743176E-002,  -5.5485232545510215E-002 ],
                                        [  -2.354423404994569 ,   0.000000000000000   ,      0.000000000000000 ],
                                        [  -2.858918674095194 ,   1.111268022307282   ,     -1.352651141853729 ],
                                        [   2.354423404994569 ,   0.000000000000000  ,       0.000000000000000 ],
                                        [   2.671741580470489 ,   1.136107563104921   ,      1.382886181959795]])
            equilibriumE=self.V([equilibriumCoords])
            #print 'equilibrium position',equilibriumCoords*bohr2ang, 'and energy:',equilibriumE*au2wn,'1/cm'
        elif self.name in ProtonatedWaterTrimer:
            equilibriumE =  -9.129961286575450E-002

        elif self.name in ProtonatedWaterTetramer:
            equilibriumE = -0.12214685901476255
        return equilibriumE


    def V(self,x):
        #print 'v',self.surfaceName,self.side,
        if self.name in Water:
            print '0.0'
        v=np.array([self.potential(cart_in) for cart_in in x])
        #if self.surfaceName=='SharedProton':
        #    if self.side=='Right':
        #        r=self.calcRn(x)
        #        v[(r<0)]=1000.00
        #    elif self.side=='Left':
        #        r=self.calcRn(x)
        #        v[(r>0)]=1000.00
        #    else:
        #        donothing=1
        #elif self.surfaceName=='OHStretchAnti':
        #    #if self.side=='Both':
        #    #    r,swap=self.calcStretchAnti(x)
        #    #    if (len(swap)>0):
        #    #        print 'swap list', swap, 'v before',v[swap],
        #    #        v[np.array(swap)]=1000.00
        #    #        print 'v after', v[swap]
        #    if self.side=='Right':      
        #        #swap=self.sortProtons()
        #        r=self.calcStretchAnti(x)
        #        v[(r<0)]=1000.00
        #    elif self.side=='Left':
        #        #swap=self.sortProtons()
        #        r=self.calcStretchAnti(x)
        #        v[(r>0)]=1000.00
        #    else:
        #        donothing=1
        return v

    def calcDipole(self,x,eckartRotated=False):
        if not eckartRotated:
            eckRotx=self.eckartRotate(x)
        else:
            eckRotx=x
        dipoleVectors=self.dipole(eckRotx)

        return dipoleVectors

    """def rotateDipoleToFrame(self,coordz,dips): 
        #This is not in use
        print 'RotatingWalkers'
        print 'dip shape', np.shape(dips)
        numWalkers = coordz.shape[0]
        # print coordz[0] #a single walker
        # print coordz[0,1] #the second atom's xyz coordinates
        # translation back to Origin
        o3 = coordz[:, 3 - 1].reshape(numWalkers, 1, 3)
        #print 'o3 shape', np.shape(o3)
        trCoordz = copy.deepcopy(coordz - o3)
        o3 = o3.reshape(numWalkers,3)
        #print 'o3 shape', np.shape(o3)
        trDips = copy.deepcopy(dips-o3[:])
        #print 'trDip shape', np.shape(trDips)

        # Rotation of O2 to x axis
        # Construct rotation matrix for a single walker
        theta = np.zeros(numWalkers)
        alpha = np.zeros(numWalkers)
        dip2 = np.zeros(trDips.shape)
        xaxtrCoordz = np.zeros(np.shape(trCoordz))
        o2 = trCoordz[:, 2 - 1].reshape(numWalkers, 1, 3)

        for walker in range(numWalkers):
            z = o2[walker, 0, 2]
            y = o2[walker, 0, 1]
            x = o2[walker, 0, 0]
            theta[walker] = np.arctan2(-1 * z, y)
            alpha[walker] = np.arctan2((-1 * (
                    y * np.cos(theta[walker]) - np.sin(theta[walker]) * z)), x)
            r1 = np.array([[1, 0, 0],
                           [0, np.cos(theta[walker]), -1 * np.sin(theta[walker])],
                           [0, np.sin(theta[walker]), np.cos(theta[walker])]]
                          )
            r2 = np.array([[np.cos(alpha[walker]), -1 * np.sin(alpha[walker]), 0],
                           [np.sin(alpha[walker]), np.cos(alpha[walker]), 0],
                           [0, 0, 1]]
                          )
            rotM = np.dot(r2, r1)
            na = 0
            for coord in trCoordz[walker]:
                xaxtrCoordz[walker, na] = np.dot(rotM, coord)
                na += 1
            dip2[walker] = np.dot(rotM,trDips[walker])
        print dip2.shape
        # print 'xaxtrCoords: ',xaxtrCoordz[0]
        # Rotation of O1 to xyplane
        #finalCoords = np.zeros(np.shape(xaxtrCoordz))
        finalDips = np.zeros(np.shape(dip2))
        print finalDips.shape
        o1 = xaxtrCoordz[:, 1 - 1].reshape(numWalkers, 1, 3)
        beta = np.zeros(numWalkers)
        for walker in range(numWalkers):
            z = o1[walker, 0, 2]
            y = o1[walker, 0, 1]
            x = o1[walker, 0, 0]
            beta[walker] = np.arctan2(-1 * z, y)
            r = np.array([[1, 0, 0],
                          [0, np.cos(beta[walker]), -1 * np.sin(beta[walker])],
                          [0, np.sin(beta[walker]), np.cos(beta[walker])]]
                         )
            naa = 0
            #for coord in xaxtrCoordz[walker]: no need to rotate last walkers - just need dipole
            #    finalCoords[walker, naa] = np.dot(r, coord)
            #    naa += 1
            finalDips[walker] = np.dot(r, dip2[walker])
        # print finalCoords
        return np.round(finalDips, 12)"""

    def rotateBackToFrame(self,coordz,a,b,c):        #use the rotation matrices that I always use to reshape each coordinate back to its reference frame
        #print coordz[1]
        print 'RotatingWalkers'
        numWalkers = coordz.shape[0]
        #translation back to Origin
        o3 = coordz[:,a-1].reshape(numWalkers,1,3)
        trCoordz = copy.deepcopy(coordz-o3)
        #Rotation of O2 to x axis
        o2 = trCoordz[:,b-1,:].reshape(numWalkers,1,3)
        z = o2[:, 0, 2]
        y = o2[:, 0, 1]
        x = o2[:, 0, 0]
        theta = np.arctan2(-1 * z, y)
        alpha = np.arctan2((-1 * (
                y * np.cos(theta) - np.sin(theta) * z)), x)
        stheta = np.sin(theta)
        ctheta = np.cos(theta)
        salpha = np.sin(alpha)
        calpha = np.cos(alpha)
        r1 = np.zeros((len(trCoordz),3,3))
        r1[:,0,:]=np.tile([1,0,0],len(trCoordz)).reshape(len(trCoordz),3)
        r1[:,1,:]=np.column_stack((np.zeros(len(trCoordz)),ctheta,-1*stheta))
        r1[:,2,:]=np.column_stack((np.zeros(len(trCoordz)),stheta,ctheta))

        r2 = np.zeros((len(trCoordz),3,3))
        r2[:, 0, :] = np.column_stack((calpha,-1*salpha,np.zeros(len(trCoordz))))
        r2[:, 1, :] = np.column_stack((salpha,calpha,np.zeros(len(trCoordz))))
        r2[:, 2, :] = np.tile([0, 0, 1], len(trCoordz)).reshape(len(trCoordz),3)

        rotM = np.matmul(r2, r1)
        #print trCoordz.shape
        xaxtrCoordz = np.matmul(rotM,trCoordz.transpose(0,2,1)).transpose(0,2,1)
        #print xaxtrCoordz[0]
        #Rotation of O1 to xyplane
        o1 = xaxtrCoordz[:,c-1,:]
        z = o1[:,2]
        y = o1[:,1]
        #x= o1[:,0]
        beta = np.arctan2(-1 * z, y)
        cbeta = np.cos(beta)
        sbeta = np.sin(beta)
        r = np.zeros((len(trCoordz),3,3))
        r[:,0,:]=np.tile([1,0,0],len(trCoordz)).reshape(len(trCoordz),3)
        r[:,1,:]=np.column_stack((np.zeros(len(trCoordz)),cbeta,-1*sbeta))
        r[:,2,:]=np.column_stack((np.zeros(len(trCoordz)),sbeta,cbeta))
        mas = np.where(np.around(la.det(r),10)!=1.0)
        print mas
        finalCoords = np.matmul(r, xaxtrCoordz.transpose(0, 2, 1)).transpose(0, 2, 1)
        o1x=np.round(finalCoords[:,1-1,0],12)
        asdf = np.where(o1x<0)
        asdf2 = np.where(o1x>0)
        #print finalCoords[0]
        print np.round(finalCoords[0],12)
        return np.round(finalCoords,12)

    def bL(self,xx,atm1,atm2):
        #Rotation of O1 to xyplane
        atmO = xx[:,atm1,:]
        atmT = xx[:, atm2, :]
        lens = la.norm(atmO-atmT,axis=1)
        return lens


    #umbrell = self.ba(addedX, H2, O, 13)  # 4 11 0 now O H1 0
    def ba(self,xx,atm1,atm2,atm3): #left center right
        #Rotation of O1 to xyplane
        atmO = xx[:,atm1,:]
        atmT = xx[:, atm2, :]
        atmH = xx[:, atm3, :]
        left = atmO - atmT
        right = atmH - atmT
        return np.arccos((left*right).sum(axis=1) / (la.norm(left,axis=1)*la.norm(right,axis=1)))


    def SymInternals(self,x,rotato=True,weights=0):
        #print 'called SymInternals'
        #print 'returning values in bohr [or radians]'
        if self.name in DeprotonatedWaterDimer:
            internals= self.SymInternalsH3O2minus(x)
            #self.internalNames=internalNames
            return internals#,internalNames
        elif self.name in ProtonatedWaterTrimer:
            if (x[0,0,2]==0.0) and (x[0,1,1]==0.0) and (x[0,1,2]==0.0) and (x[0,2,0]==0.0) and (x[0,2,1]==0.0) and (x[0,2,2]==0.0):
                print 'No need to rotate.'
            else:
                #print 'Coordinates that matter: (b4): ',x[0,:3]*ang2bohr
                x = self.rotateBackToFrame(x,3,2,1)
                #print 'Coordinates that matter: (af): ', x[0, :3]*ang2bohr
                if not (x[0, 0, 2] == 0.0) and (x[0, 1, 1] == 0.0) and (x[0, 1, 2] == 0.0) and (x[0, 2, 0] == 0.0) and (
                        x[0, 2, 1] == 0.0) and (x[0, 2, 2] == 0.0):
                    rotationDidntWork
            internals = self.SymInternalsH7O3plus(x)
            return internals
        elif self.name in ProtonatedWaterTetramer:
            #print 'You just commented stuff out, fam'
            if (x[0, 0, 1] == 0.0) and (x[0, 0, 2] == 0.0) and (x[0, 1, 0] == 0.0) and (x[0, 1, 1] == 0.0) and\
                    (x[0, 1, 2] == 0.0) and (x[0, 2, 2] == 0.0):
                print 'No need to rotate.'
            else:
                # print 'Coordinates that matter: (b4): ',x[0,:3]*ang2bohr
                x = self.rotateBackToFrame(x, 2, 1, 3)
                # print 'Coordinates that matter: (af): ', x[0, :3]*ang2bohr
                if not (x[0, 0, 1] == 0.0) and (x[0, 0, 2] == 0.0) and (x[0, 1, 0] == 0.0) and (x[0, 1, 1] == 0.0) and\
                    (x[0, 1, 2] == 0.0) and (x[0, 2, 2] == 0.0):
                    rotationDidntWork
                print 'Rotating done, translate COM...'
                x-=((x[:,3-1,:]+x[:,2-1,:]+x[:,1-1,:])/3)[:,np.newaxis,:]
            internals = self.SymInternalsH9O4plus(x)

            return internals
        #elif self.name in ProtonatedWaterTetramer:
        #    return self.SymInternalsH9O4plus(x,printFlag=printFlag)
        else:
            print 'woefully unprepared to handle the calculation of the SymInternals for ', self.molecule
            #crash


    #def bondlength(self,pos,atom1=1,atom2=3):
    #    length=(pos[:,atom1,0]-pos[:,atom2,0])**2+(pos[:,atom1,1]-pos[:,atom2,1])**2+(pos[:,atom1,2]-pos[:,atom2,2])**2
    #    length=np.sqrt(length)
    #    return length

    def PltHists1D(self,cfg, thing, bound, xl, yl, overly, weits):
        if len(thing) > 90: #if len(thing) == 2: Just changed currently for plotting purposes
            theLen, xx = np.histogram(thing, bins=75, range=bound, normed=True, weights=weits)  # WEIGHTS=WEIGHTARRAY
            inin = True
            #overlay = False
        else:
            inin = False
            c = 0
            theLen = np.zeros((np.size(thing), 75))
            for i in thing:
                theLen[c], xx = np.histogram(i, bins=75, range=bound, normed=True,
                                             weights=weits[c])  # WEIGHTS=WEIGHTARRAY
                c += 1
        bnd = str(bound[0]).replace(".", "").replace("-", "") + str(bound[1]).replace(".", "").replace("-", "")
        print bnd
        mP = Plot.myPlot(cfg, '1d', bnd, bnd, xl, yl, overly, inin, theLen)
        mP.plotIt()

    def PltHists2D(self,cfg, thing1, thing2, bound1, bound2, xl, yl, overly, weits):
        if len(thing1) == 2:
            theLen, xx, yy = np.histogram2d(thing1, thing2, bins=75, range=(bound1, bound2), normed=True,
                                            weights=weits)  # WEIGHTS=WEIGHTARRAY
            inin = True
            overlay = False
        else:
            inin = False
            c = 0
            theLen = np.zeros((np.size(thing1), 75, 75))
            print np.size(theLen)
            print np.size(theLen[0])
            for i, j in zip(thing1, thing2):
                theLen[c], x2x, y2y = np.histogram2d(i, j, bins=75, range=(bound1, bound2), normed=True,
                                                     weights=weits[c])  # WEIGHTS=WEIGHTARRAY
                c += 1
        bnd1 = str(bound1[0]).replace(".", "").replace("-", "") + str(bound1[1]).replace(".", "").replace("-", "")
        bnd2 = str(bound2[0]).replace(".", "").replace("-", "") + str(bound2[1]).replace(".", "").replace("-", "")
        mP = Plot.myPlot(cfg, '2d', bnd1, bnd2, xl, yl, overly, inin, theLen)
        mP.plotIt()


    def xyzTrimerSharedHydrogens(self,atmnm,xx):
        if atmnm == 9:
            outerW = 1
        if atmnm == 10:
            outerW = 2
        mp = (xx[:,3-1]+xx[:,outerW-1])/ 2
        print mp.shape
        xaxis = np.divide((xx[:,outerW-1] - mp),la.norm(xx[:,outerW-1,] - mp,axis=1).reshape(-1,1))
        print xaxis.shape
        print xx[:,1-1].shape
        zaxis = np.cross(xx[:,1-1]-xx[:,3-1],xx[:,2-1]-xx[:,3-1],axis=1)
        yaxis = np.cross(zaxis,xaxis,axis=1)
        #xcomp11 = ((xx[:,atmnm-1]-mp)*X).sum(axis=1)
        xcomp = ((xx[:,atmnm-1] - mp)*xaxis).sum(axis=1)
        ycomp = ((xx[:,atmnm-1] - mp)*yaxis).sum(axis=1)
        zcomp = ((xx[:,atmnm-1] - mp)*zaxis).sum(axis=1)
        return xcomp, ycomp, zcomp

    def getBisectingVector(self,left, middle, right):
        bisector1 = la.norm(left - middle, axis=1).reshape(-1, 1) * (right - middle)  # |b|*a + |a|*b
        bisector2 = la.norm(right - middle, axis=1).reshape(-1, 1) * (left - middle)
        normedbisector = la.norm(bisector1 + bisector2, axis=1).reshape(-1, 1)
        bisector = (bisector1 + bisector2) / normedbisector
        #test = np.arccos(((left-middle)*bisector).sum(axis=1)/(la.norm(left)*la.norm(bisector)))
        #test2 = np.arccos(((right-middle) * bisector).sum(axis=1) / (la.norm(right) * la.norm(bisector)))
        #test3= np.arccos(((right-middle) * (left-middle)).sum(axis=1) / (la.norm(right) * la.norm(left)))
        #print 'testing'
        return bisector

    def xyzFreeHydronium(self,xx):
        xaxis = self.getBisectingVector(xx[:,9-1,:], xx[:,3 - 1,:],xx[:,10 - 1,:])
        crs = np.cross(xx[:,9-1]-xx[:,3-1],xx[:,10-1]-xx[:,3-1],axis=1)
        zaxis = crs/((la.norm(crs,axis=1))[:,np.newaxis])
        yaxis = np.cross(zaxis,xaxis,axis=1)
        xcomp = ((xx[:,8-1] - xx[:,3-1])*xaxis).sum(axis=1)
        ycomp = ((xx[:,8-1] - xx[:,3-1])*yaxis).sum(axis=1)
        zcomp = ((xx[:,8-1] - xx[:,3-1])*zaxis).sum(axis=1)
        rdistOH = la.norm(np.column_stack((xcomp,ycomp,zcomp)), axis=1)
        thetaOH = np.arccos(zcomp / rdistOH)
        phiOH = np.arctan2(ycomp,xcomp)
        phiOH[phiOH <= 0]+=(2.*np.pi)
        return rdistOH, thetaOH, phiOH

    def finalTrimerEuler(self,xx,O1, h1, h2):
        #SharedProtonCoordinateSystem
        X = np.divide((xx[:, O1 - 1, :] - xx[:,3-1]) , la.norm(xx[:, O1 - 1, :] - xx[:,3-1], axis=1).reshape(-1,1))
        crs=np.cross(xx[:, 1 - 1]-xx[:,3-1], xx[:, 2 - 1]-xx[:,3-1], axis=1)
        Z = crs / la.norm(crs,axis=1)[:,np.newaxis]
        Y = np.cross(Z, X, axis=1)

        x,y,z=self.H9GetHOHAxis(xx[:, O1 - 1], xx[:, h1 - 1], xx[:, h2 - 1])
        #I EDITED THIS ON SATURDAY AFTER I RAN THE MOST RECENT RESULTS _ CAN RUN AGAIN ON MONDAY WITH THIS ADDITION.
        exx = np.copy(x)
        x = np.copy(y)
        y = np.copy(z)
        z = np.copy(exx)
        #
        # exx = np.copy(x)
        # x = np.copy(z)
        # z = np.copy(y)
        # y = np.copy(exx)
        #
        print 'lets get weird'
        # exX = np.copy(X)
        # X = np.copy(Y)
        # Y = np.copy(Z)
        # Z = np.copy(exX)
        #
        exX = np.copy(X)
        X = np.copy(Z)
        Z = np.copy(Y)
        Y = np.copy(exX)

        Theta,tanPhi,tanChi=self.eulerMatrix(x,y,z,X,Y,Z)
        return Theta,tanPhi, tanChi


    def  H9getOOOAxis(self,xx,oW):
        # oW = xx[:,outerW-1,:]  # Coordinates of outer water Oxygen
        center = xx[:, 4 - 1, :]
        xaxis = np.divide((oW - center), la.norm(oW - center, axis=1).reshape(-1,1))  # Normalized coordinates of xaxis definition. aka vector with only x component, where x = 1
        # zaxis
        dummy = center.copy()
        dummy[:, -1] = 0.0
        OA = oW-center
        BA = dummy-center

        un = OA / la.norm(OA, axis=1)[:, None]  # unit vector pointing along OO axis
        OB = oW-dummy

        s = (un*OB).sum(axis=1)
        OC = s[:, None] * un

        AC = OA - OC
        ze = AC - BA
        zaxis = ze / la.norm(ze, axis=1)[:, None]
        zdotx= (zaxis*xaxis).sum(axis=1)
        yaxis=np.cross(zaxis,xaxis,axis=1)
        return xaxis, yaxis, zaxis

    def H9GetHOHAxis(self,o,h1,h2):
        xaxis = self.getBisectingVector(h1, o, h2)  # bisector of HOH angle
        #yaxis = np.zeros((len(h1), 3))
        #zaxis = np.zeros((len(h1), 3))
        h1on = (h1 - o) / la.norm(h1 - o, axis=1).reshape([-1, 1])
        h2on = (h2 - o) / la.norm(h2 - o, axis=1).reshape([-1, 1])
        zaxis=np.cross(h1on,h2on,axis=1)
        yaxis = np.cross(zaxis,xaxis,axis=1)
        return xaxis, yaxis, zaxis

    def calcD(self,xx,O, H1, H2, H3):
        # calc the unit vectors from O to Hi
        # unit vectors are called "a"
        #O -= 1
        #H1 -= 1
        #H2 -= 1
        #H3 -= 1

        # print 'calculating the xxition of D'
        # print xx
        O = xx[:, O]
        # Every walker's xyz coordinate for O
        H1 = xx[:, H1]
        H2 = xx[:, H2]
        H3 = xx[:, H3]
        #test = H1-O
        # print 'norm = ',la.norm(H1-O,axis=1)
        aOH1 = np.divide((H1 - O), la.norm(H1 - O, axis=1)[:, np.newaxis])  # broadcasting silliness
        aOH2 = np.divide((H2 - O), la.norm(H2 - O, axis=1)[:, np.newaxis])
        aOH3 = np.divide((H3 - O), la.norm(H3 - O, axis=1)[:, np.newaxis])
        # point in space along OH bonds that is 1 unit away from the O
        # print aOH1
        aH1 = O + aOH1
        aH2 = O + aOH2
        aH3 = O + aOH3
        # midpoint between unit vecs
        # maH1H2 = (aH1 + aH2) / 2.0
        # maH2H3 = (aH2 + aH3) / 2.0
        # maH1H3 = (aH1 + aH3) / 2.0
        # vectors between the points along the OH bonds that are 1 unit vector away from the O
        vaH1H2 = aH2 - aH1
        vaH2H3 = aH3 - aH2
        # vaH1H3 = aH3 - aH1

        # calculate vector
        # line = np.zeros((xx.shape[0], 3))
        # for i in range(xx.shape[0]):
        #     line[i] = np.cross(vaH1H2[i], vaH2H3[i].T)
        # add normalized vector to O
        line = np.cross(vaH1H2, vaH2H3, axis=1)

        #print 'max, min, ave of mag of line', np.average(la.norm(line, axis=1)), np.max(la.norm(line, axis=1)), np.min(
        #    la.norm(line, axis=1))
        #g = (line / la.norm(line, axis=1)[:, np.newaxis])
        D = O + (line / la.norm(line, axis=1)[:, np.newaxis])
        return D

    def umbrella(self,xx,O,H1,H2,H3):
        #xx*=bohr2ang
        """O,H1,H2,H3 are indices. """
        # calculate d, the trisector point of the umbrella
        D = self.calcD(xx,O, H1, H2, H3)

        # print LindseyCoords.shape, D.shape
        addedX = np.concatenate((xx, D[:, np.newaxis, :]), axis=1)  # Change D index
        getXYZ=False
        if getXYZ:
            wf = open('test_umb.xyz','w+')
            trim = ["O","O","O","O","H","H","H","H","H","H","H","H","H","F"]
            for wI, walker in enumerate(addedX):
                wf.write("14\n")
                wf.write("0.0 0.0 0.0 0.0 0.0\n")
                for aI, atm in enumerate(walker):
                    wf.write("%s %5.12f %5.12f %5.12f\n" % (trim[aI], atm[0], atm[1], atm[2]))
                wf.write("\n")
            wf.close()
        # # print 'addedLindseyCoords.shapes', addedLindseyCoords.shape
        umbrell = self.ba(addedX, H2, O, -1)  # 4 11 0 now O H1 0
        print np.degrees(umbrell)
        return umbrell

    def getfinalOOAxes(self,atmnm,xx):
        if atmnm == 11:
            outerW = 2
            at2 = 1
            at3 = 3
        elif atmnm == 12:
            outerW = 3
            at2 = 2
            at3 = 1
        elif atmnm == 13:
            outerW = 1
            at2 = 3
            at3 = 2
        oW = xx[:, outerW - 1, :]  # Coordinates of outer water Oxygen
        center = xx[:, 4 - 1, :]
        dummy = np.copy(center)
        dummy[:, -1] = 0.0
        ZBig = np.cross(xx[:,at2-1]-dummy,xx[:,at3-1]-dummy,axis=1)
        #test = la.norm(ZBig,axis=1)
        ZBig /= la.norm(ZBig,axis=1)[:,None]
        #mp = (center + oW) / 2 #no decimal
        #print mp
        #sharedH = xx[:, atmnm - 1, :]  # Coordinates of shared Hydrogen
        xaxisp = np.divide((center - oW), la.norm(center - oW, axis=1).reshape(-1,1))  # Normalized coordinates of xpaxis definition. aka vector with only x component, where x = 1
        oaHat=np.copy(xaxisp)
        OB=dummy-oW
        s=(oaHat*OB).sum(axis=1)
        OC=oaHat*s[:,np.newaxis]
        if np.all(np.around(OC,12)==np.around(OB,12)):
            ze = np.copy(ZBig)
        else:
            ze=OC-OB

        #at this point, my x axis points the 'wrong' direction.  I will flip the sign
        xaxis=np.divide((oW-center), la.norm(oW-center, axis=1).reshape(-1,1))

        zaxis = ze / la.norm(ze, axis=1)[:, None]
        sgn = np.where((ZBig*zaxis).sum(axis=1) < 0)[0]
        zaxis[sgn] = np.negative(zaxis[sgn])
        yaxis = np.cross(zaxis, xaxis, axis=1)
        return xaxis,yaxis,zaxis


    def eulerMatrix(self,x,y,z,X,Y,Z):
        #Using the
        #[X]    [. . .][x]
        #[Y] =  [. . .][y]
        #[Z]    [. . .][z]
        zdot=(z * Z).sum(axis=1) / (la.norm(z, axis=1) * la.norm(Z, axis=1))
        Yzdot=(Y * z).sum(axis=1)/(la.norm(Y,axis=1) * la.norm(z,axis=1))
        Xzdot=(X*z).sum(axis=1)/(la.norm(X,axis=1) * la.norm(z,axis=1))
        yZdot=(y*Z).sum(axis=1) / (la.norm(y,axis=1) * la.norm(Z,axis=1))
        xZdot=-(x*Z).sum(axis=1) / (la.norm(x,axis=1) * la.norm(Z,axis=1))
        Theta = np.arccos(zdot)
        tanPhi = np.arctan2(Yzdot,Xzdot)
        tanChi = np.arctan2(yZdot,xZdot) #negative baked in
        # tanChi[tanChi < 0]+=(2*np.pi)
        # tanPhi[tanPhi < 0]+=(2*np.pi)
        return Theta, tanPhi, tanChi

    def HDihedral(self,xx):
        if xx.shape[1] == 13:
            d=self.calcD(xx, 4-1, 11-1, 12-1, 13-1)
            a = 11
            b = 12
            c = 13
        elif xx.shape[1] == 10:
            d = self.calcD(xx, 3 - 1, 8 - 1, 10- 1, 9- 1)
            a = 8
            b = 10
            c = 9
        addedX = np.concatenate((xx, d[:, np.newaxis, :]), axis=1)
        di1 = self.fwiki_dihedral(addedX,b-1,c-1)
        di2 = self.fwiki_dihedral(addedX,a-1,b-1)
        di3 = self.fwiki_dihedral(addedX,c-1,a-1)
        return di1,di2,di3

    def fwiki_dihedral(self,xx,b1A,b2A):
        # https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
        b1=xx[:,b1A]-xx[:,-1]
        b2=xx[:,b2A]-xx[:,-1]
        if xx.shape[1] == 13+1:
            b3=xx[:,4-1]-xx[:,-1]
        elif xx.shape[1] == 10+1:
            b3 = xx[:, 3-1] - xx[:, -1]
        crossterm1 = np.cross(np.cross(b1,b2,axis=1),np.cross(b2,b3,axis=1),axis=1)
        term1 = (crossterm1*(b2/la.norm(b2,axis=1)[:,np.newaxis])).sum(axis=1)
        term2 = (np.cross(b1,b2,axis=1)*np.cross(b2,b3,axis=1)).sum(axis=1)
        dh = np.arctan2(term1,term2)
        #print np.degrees(dh)
        return dh

    def getHydroniumAxes(self,xx):
        X = (xx[:, 1 - 1] - xx[:, 2 - 1]) / la.norm(xx[:, 1 - 1] - xx[:, 2 - 1], axis=1)[:, np.newaxis]
        cr = np.cross(xx[:, 1 - 1] - xx[:, 2 - 1], xx[:, 3 - 1] - xx[:, 2 - 1])
        Z = cr / la.norm(cr, axis=1)[:, np.newaxis]
        Y = np.cross(Z, X)

        x = (xx[:, 13 - 1] - xx[:, 11 - 1]) / la.norm(xx[:, 13 - 1] - xx[:, 11 - 1], axis=1)[:, np.newaxis]
        cr2 = np.cross(xx[:, 13 - 1] - xx[:, 11 - 1], xx[:, 12 - 1] - xx[:, 11 - 1])
        z = cr2 / la.norm(cr2, axis=1)[:, np.newaxis]
        y = np.cross(z, x)
        return X,Y,Z,x,y,z

    def extractEulers(self,rotMs):
        # [x]    [. . .][X]
        # [y] =  [. . .][Y]
        # [z]    [. . .][Z]
        #be careful with which euler matrix you're using, body fixed space fixed stuff\
        zdot=rotMs[:,-1,-1]
        Yzdot = rotMs[:,2,1]
        Xzdot = rotMs[:,2,0]
        yZdot = rotMs[:,1,2]
        xZdot = rotMs[:,0,2]
        Theta = np.arccos(zdot)
        tanPhi = np.arctan2(Yzdot, Xzdot)
        tanChi = np.arctan2(yZdot, -xZdot)  # negative baked in
        # tanChi[tanChi < 0]+=(2*np.pi)
        # tanPhi[tanPhi < 0]+=(2*np.pi)
        return Theta,tanPhi,tanChi

    def finalPlaneShareEuler(self,xx):
        print 'get carts & eulers pre eckart'
        atmnm=11
        h1 = 8
        h2 = 7
        o  = 2
        X,Y,Z = self.getfinalOOAxes(atmnm,xx)
        x,y,z = self.H9GetHOHAxis(xx[:,o-1],xx[:,h1-1],xx[:,h2-1])
        th11,phi11,xi11 = self.eulerMatrix(x,y,z,X,Y,Z)

        atmnm=12
        h1 = 9
        h2 = 10
        o = 3
        X,Y,Z = self.getfinalOOAxes(atmnm,xx)
        x,y,z = self.H9GetHOHAxis(xx[:,o-1],xx[:,h1-1],xx[:,h2-1])
        th12,phi12,xi12 = self.eulerMatrix(x,y,z,X,Y,Z)

        atmnm=13
        h1 = 6
        h2 = 5
        o = 1
        X,Y,Z = self.getfinalOOAxes(atmnm,xx)
        x,y,z = self.H9GetHOHAxis(xx[:,o-1],xx[:,h1-1],xx[:,h2-1])
        th13,phi13,xi13 = self.eulerMatrix(x,y,z,X,Y,Z)

        umbrella=self.umbrella(xx,4-1,11-1,12-1,13-1)
        dh1,dh2,dh3=self.HDihedral(xx)
        dh1[dh1<0.0]+=2*np.pi
        dh2[dh2 < 0.0] += 2*np.pi
        dh3[dh3 < 0.0] += 2*np.pi


        #################################################################
        # XH,YH,ZH,xH,yH,zH=self.getHydroniumAxes(xx)
        # thH, phiH, xiH = self.eulerMatrix(xH, yH, zH, XH, YH, ZH)
        # mat = self.getEulerMat(thH,phiH,xiH)
        # dthH=np.degrees(thH)
        # dphiH=np.degrees(phiH)
        # dxiH=np.degrees(xiH)
        ##################################################################
        rOH11=self.bL(xx,11-1,4-1)
        rOH12=self.bL(xx,12-1,4-1)
        rOH13=self.bL(xx,13-1,4-1)
        # return dh1,dh2,dh3,xcomp11,ycomp11,zcomp11,xcomp12,ycomp12,zcomp12,xcomp13,ycomp13,zcomp13,th11,phi11,xi11,th12,phi12,xi12,th13,phi13,xi13

        print 'eckarting...'
        ocom, eVecs,kil=self.eckartRotate(xx,cart=True)
        # ocom, eVecs,kil=self.eckartRotate(xx,justO=True)

        print 'got matrix'
        xx-=ocom[:,np.newaxis,:]
        print 'done'
        print 'b4'
        print xx[0]
        xx = np.einsum('knj,kij->kni', eVecs.transpose(0, 2, 1), xx).transpose(0, 2, 1)
        print 'af'
        print xx[0]
        print 'fully rotated'
        ocomH,eVecsH,kilH=self.eckartRotate(xx,hydro=True,yz=True)
        eVecsH=eVecsH.transpose(0,2,1)
        thH,phiH,xiH=self.extractEulers(eVecsH)
        return xx[:,4-1,0],xx[:,4-1,1],xx[:,4-1,2],rOH11, rOH12, rOH13, umbrella, 2 * dh1 - dh2 - dh3, dh2 - dh3, thH, phiH, xiH, th11, phi11, xi11, th12, phi12, xi12, th13, phi13, xi13

        # umTh = self.umbrella(xx,4-1,1-1,2-1,3-1)
        # th1 = self.ba(xx,2-1,4-1,3-1)
        # th2 = self.ba(xx,3-1,4-1,1-1)
        # th3 = self.ba(xx,2-1,4-1,1-1)
        # return umTh, 2*th1-th2-th3, th2-th3, xcomp11, ycomp11, zcomp11, xcomp12, ycomp12, zcomp12, xcomp13, ycomp13, zcomp13, th11, phi11, xi11, th12, phi12, xi12, th13, phi13, xi13

        # return xx[:,4-1,0],xx[:,4-1,1],xx[:,4-1,2],xcomp11,ycomp11,zcomp11,xcomp12,ycomp12,zcomp12,xcomp13,ycomp13,zcomp13,th11,phi11,xi11,th12,phi12,xi12,th13,phi13,xi13

    def getEulerMat(self,th, ph, xi):
        a = np.array([[np.cos(ph) * np.cos(th) * np.cos(xi) - np.sin(ph) * np.sin(xi),
                       np.sin(ph) * np.cos(th) * np.cos(xi) + np.cos(ph) * np.sin(xi),
                       -np.sin(th) * np.cos(xi)]
                         ,
                      [-np.cos(ph) * np.cos(th) * np.sin(xi) - np.sin(ph) * np.cos(xi),
                       -np.sin(ph) * np.cos(th) * np.sin(xi) + np.cos(ph) * np.cos(xi),
                       np.sin(th) * np.sin(xi)]
                         ,
                      [np.cos(ph) * np.sin(th),
                       np.sin(ph) * np.sin(th),
                       np.cos(th)]
                      ])
        return np.array([[np.cos(ph) * np.cos(th) * np.cos(xi) - np.sin(ph) * np.sin(xi),
                          np.sin(ph) * np.cos(th) * np.cos(xi) + np.cos(ph) * np.sin(xi),
                          -np.sin(th) * np.cos(xi)]
                            ,
                         [-np.cos(ph) * np.cos(th) * np.sin(xi) - np.sin(ph) * np.cos(xi),
                          -np.sin(ph) * np.cos(th) * np.sin(xi) + np.cos(ph) * np.cos(xi),
                          np.sin(th) * np.sin(xi)]
                            ,
                         [np.cos(ph) * np.sin(th),
                          np.sin(ph) * np.sin(th),
                          np.cos(th)]
                         ])
    def SymInternalsH9O4plus(self,x):
        print 'Commence getting internal coordinates for tetramer'
        start = time.time()
        all = self.finalPlaneShareEuler(x)
        rOHHyd = all[3:6]
        umbDi = all[6:9]
        eulHyd = all[9:12]
        thphixi1=all[12:15]
        thphixi2=all[15:18]
        thphixi3=all[18:21]
        xyzO4 = all[0:3]
        print 'done hydronium XYZ'
        print 'time it took to get xyzs,eulers: ',str(time.time()-start)
        third = time.time()
        rOH5 = self.bL(x,5-1,1-1)
        rOH6 = self.bL(x,6-1,1-1)
        HOH516 = self.ba(x,5-1,1-1,6-1)
        rOH7 = self.bL(x,7-1,2-1)
        rOH8 = self.bL(x,8-1,2-1)
        HOH728 = self.ba(x,7-1,2-1,8-1)
        rOH9 = self.bL(x,9-1,3-1)
        rOH10 = self.bL(x,10-1,3-1)
        HOH9310 = self.ba(x,9-1,3-1,10-1)
        rO1O2 = self.bL(x,2-1,1-1)
        rO1O3 = self.bL(x,1-1,3-1)
        rO2O3 = self.bL(x,2-1,3-1)
        print 'time for rOO/rOH', str(time.time() - third)
        print 'Done with all internals'
        internal= np.array(
            (rOHHyd[0], rOHHyd[1], rOHHyd[2], umbDi[0], umbDi[1], umbDi[2], eulHyd[0], eulHyd[1], eulHyd[2],
                thphixi1[0], thphixi1[1], thphixi1[2], thphixi2[0], thphixi2[1], thphixi2[2], thphixi3[0], thphixi3[1],
                thphixi3[2], rOH5, rOH6, HOH516, rOH7, rOH8, HOH728, rOH9, rOH10, HOH9310, rO1O2, rO1O3, rO2O3,
                xyzO4[0], xyzO4[1], xyzO4[2])).T
        self.internalName = ['rOH11', 'rOH12', 'rOH13', 'umbrella', '2dihed', 'dihed-di', 'thH', 'phH', 'xiH', 'theta651',
                            'phi651', 'Chi651',
                            'theta1039', 'phi1039', 'Chi1039', 'theta728', 'phi728', 'Chi728', 'rOH5', 'rOH6',
                            'HOH516', 'rOH7', 'rOH8', 'HOH728',
                            'rOH9', 'rOH10', 'HOH9310', 'rO1O2','rO1O3','rO2O3', 'xo4', 'yo4', 'zo4']
        print 'internal shape: ',np.shape(internal)
        print 'internal[0] shape: ',np.shape(internal[0])
        return internal

    def setInternalName(self):
        if self.name in DeprotonatedWaterDimer:
            self.internalName=[]
        elif self.name in ProtonatedWaterTrimer:
            self.internalName = ['rOH9', 'rOH10', 'spHOH', 'rH8', 'thH8', 'phiH8', 'thH', 'phiH', 'xiH',
                                 'th_627', 'phi_627', 'xi_627', 'th_514', 'phi_514', 'xi_514', 'rOH_41',
                                 'rOH_51', 'aHOH_451', 'rOH_26', 'rOH_27', 'aHOH_267', 'rOO_1', 'rOO_2', 'aOOO']
        elif self.name in ProtonatedWaterTetramer:
            self.internalName = ['rOH11', 'rOH12', 'rOH13', 'umbrella', '2dihed', 'dihed-di', 'thH', 'phH', 'xiH',
                                 'theta651',
                                 'phi651', 'Chi651',
                                 'theta1039', 'phi1039', 'Chi1039', 'theta728', 'phi728', 'Chi728', 'rOH5', 'rOH6',
                                 'HOH516', 'rOH7', 'rOH8', 'HOH728',
                                 'rOH9', 'rOH10', 'HOH9310', 'rO1O2', 'rO1O3', 'rO2O3', 'xo4', 'yo4', 'zo4']

    def finalTrimerHydEuler(self,xx):
        print 'eckarting...'
        ocom, eVecs,kil=self.eckartRotate(xx,justO=True)
        print 'got Cart matrix'
        xx-=ocom[:,np.newaxis,:]
        print 'done'
        print 'b4'
        print xx[0]
        xx = np.einsum('knj,kij->kni', eVecs.transpose(0, 2, 1), xx).transpose(0, 2, 1)
        print 'af'
        print xx[0]
        print 'fully rotated'
        ocomH,eVecsH,kilH=self.eckartRotate(xx,hydro=True,yz=True)
        print 'got Hydro matrix'
        eVecsH=eVecsH.transpose(0,2,1)
        print 'ROTM FOR HYD',eVecsH[0]
        thH,phiH,xiH=self.extractEulers(eVecsH)

        # phiH[phiH<0.0]+=(2*np.pi)
        # xiH[xiH< 0.0] += (2 * np.pi)

        return thH,phiH,xiH

    def SymInternalsH7O3plus(self,x):
        print 'Commence getting internal coordinates for trimer'
        #Hydronium
        rOH9 = self.bL(x,3-1,9-1)
        rOH10 = self.bL(x,3-1,10-1)
        spHOH = self.ba(x,9-1,3-1,10-1)
        rthphi = self.xyzFreeHydronium(x)  # *ang2bohr     #Spherical free hydrogen 3xnwalkers in degrees
        print 'done with FH'
        thphixi1= self.finalTrimerEuler(x,2,7,6)
        print 'done with Euler1'
        thphixi2 = self.finalTrimerEuler(x,1,4,5)
        print 'done with Euler2'
        thH,phiH,xiH = self.finalTrimerHydEuler(x)
        rOH1 = self.bL(x,1-1,4-1)
        rOH2 = self.bL(x,1-1,5-1)
        aHOH1= self.ba(x,5-1,1-1,4-1)
        rOH3 = self.bL(x,2-1,6-1)
        rOH4 = self.bL(x,2-1,7-1)
        aHOH2= self.ba(x,6-1,2-1,7-1)
        rOO1 = self.bL(x,1-1,3-1)
        rOO2 = self.bL(x,2-1,3-1)
        aOOO = self.ba(x,1-1,3-1,2-1)
        print 'first O1H4 bond length: ', rOH1[0]*bohr2ang
        print 'first Angle: ',np.degrees(aOOO[0])
        internal = np.array((rOH9, rOH10, spHOH, rthphi[0], rthphi[1], rthphi[2], thH,phiH,xiH,
                thphixi1[0], thphixi1[1], thphixi1[2], thphixi2[0], thphixi2[1], thphixi2[2]
                , rOH1, rOH2, aHOH1, rOH3, rOH4, aHOH2, rOO1, rOO2, aOOO)).T
        self.internalName = ['rOH9', 'rOH10', 'spHOH', 'rH8', 'thH8', 'phiH8', 'thH', 'phiH', 'xiH',
                             'th_627', 'phi_627','xi_627', 'th_514', 'phi_514', 'xi_514', 'rOH_41',
                             'rOH_51', 'aHOH_451', 'rOH_26','rOH_27', 'aHOH_267','rOO_1', 'rOO_2', 'aOOO']
        return internal


    def SymInternalsH3O2minus(self,x): #get an array of all the internal coordinates associated with H3O2 minus
        #print 'calculating the internals...ver 1...'
        #internals used in jpc a paper
        #print 'called symInternalsVer1. Please only provide eckart rotated molecules. Thank you.'
        #    print 'The first walker is: \n', x[0]   

        rOH1=self.bondlength(x,atom1=1,atom2=2)
        rOH2=self.bondlength(x,atom1=3,atom2=4)
        rOO=self.bondlength(x,atom1=1,atom2=3)
        aHOO1=self.bondAngle(x,atom1=2, atom2=1, atom3=3)
        aHOO2=self.bondAngle(x,atom1=4, atom2=3, atom3=1)
        tHOOH,tRange=self.calcTorsion(x)
        HdispX,HdispY, HdispZ = self.calcCartesianSharedProtonDisplacement(x)

        #rn=self.calcRn(x)
        #NOW SYMETRIZE         

        rOH_s=np.sqrt(0.5)*(rOH1+rOH2) #symetric                                                    
        rOH_a=np.sqrt(0.5)*(rOH1-rOH2) #asym stretch                                                
        #rOH_ai=0.5*(rOH1-rOH2+rOH3-rOH4) # in phase anti sym                                       
        #rOH_ao=0.5*(rOH1-rOH2-rOH3+rOH4) #out of phase anti sym                                    

        aHOO_s=np.sqrt(0.5)*(aHOO1+aHOO2) #symetric                                                 
        aHOO_a=np.sqrt(0.5)*(aHOO1-aHOO2) #asymetric                                                

        #rearrange these       

        if rOH1.size<2:
            internal = np.array( [rOH_s, rOH_a, rOO, aHOO_s, aHOO_a,tHOOH, HdispX,HdispY,HdispZ])
            #internal = np.array( [rOH_s, rOH_a, rOO, aHOO_s, aHOO_a,tHOOH, rn,HdispY,HdispZ])      
        else:
            internal = np.array(zip(rOH_s, rOH_a, rOO, aHOO_s, aHOO_a,tHOOH, HdispX,HdispY,HdispZ))
            #internal = np.array(zip(rOH_s, rOH_a, rOO, aHOO_s, aHOO_a,tHOOH, rn,HdispY,HdispZ))    

        self.internalName=['rOH_s', 'rOH_a', 'rOO', 'rHOO_s', 'rHOO_a', 'tHOOH','HdispX','HdispY','HdispZ']
        #self.internalName=['rOH_s', 'rOH_a', 'rOO', 'rHOO_s', 'rHOO_a', 'tHOOH','rn','HdispY','HdispZ']                                                                 
        self.internalConversion=[bohr2ang,bohr2ang,bohr2ang,rad2deg,rad2deg,rad2deg,bohr2ang,bohr2ang,bohr2ang]
        return internal

    def pullTrimerRefPos(self,yz=False): #Eckart reference for the trimer is in an xyz file. Need just a 3xNatom array of reference structures. I can hard code this in
        #This one is good.
        # myBetterRef = np.array(
        #                 [
        #                     [3.15544362E-30 , 4.06869143E+00, -7.59761292E-01],
        #                     [-4.98270994E-16, -4.06869143E+00 ,-7.59761292E-01],
        #                     [-1.97215226E-30,  0.00000000E+00,  1.58929880E+00],
        #                     [ 1.47532198E+00,  5.00324669E+00, -1.29932702E+00],
        #                     [ -1.47532198E+00,  5.00324669E+00, -1.29932702E+00],
        #                     [ -1.47532198E+00, -5.00324669E+00 ,-1.29932702E+00],
        #                     [ 1.47532198E+00, -5.00324669E+00, -1.29932702E+00],
        #                     [ -3.94430453E-30, -2.22044605E-16,  3.41400471E+00],
        #                     [  7.88860905E-31, 1.69178407E+00,  6.12546816E-01],
        #                     [ -2.07183794E-16, -1.69178407E+00,  6.12546816E-01]])

        myBetterRef = np.array(
            [
                 [-2.34906009e+00,  4.06869143e+00,  0.00000000e+00],
                 [ 4.69812018e+00,  0.00000000e+00,  0.00000000e+00],
                 [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                 [-2.88862583e+00,  5.00324669e+00,  1.47532198e+00],
                 [-2.88862583e+00,  5.00324669e+00, -1.47532198e+00],
                 [ 5.77725164e+00, -2.46900000e-09, -1.47532198e+00],
                 [ 5.77725164e+00, -2.46900000e-09,  1.47532198e+00],
                 [-9.12352955e-01, -1.58024167e+00,  0.00000000e+00],
                 [-9.76751990e-01,  1.69178407e+00,  0.00000000e+00],
                 [ 1.95350397e+00, -3.53000000e-09,  0.00000000e+00]])
        if not yz:

            # sqr2 = np.sqrt(2) / 2.
            # rotM = np.array([[sqr2, -sqr2, 0.],
            #              [sqr2, sqr2, 0.],
            #              [0., 0., 1.]
            #              ])
            # myBetterRef= np.dot(rotM,myBetterRef.T).T
            # rotM = np.array([[1.,0.,0.],
            #              [0, 0., -1.],
            #              [0.,1.,0.]])
            # # rotM = np.array([[0.,0.,1.],
            # #              [0, 1, 0],
            # #              [-1.,0,0.]
            # #              ])
            # # rotM = np.array([[1., 0., 0.],
            # #                   [0, sqr2, -sqr2],
            # #                   [0., sqr2, sqr2]
            # #                   ])
            # # # rotM = np.array([[sqr2, 0., sqr2],
            # # #                   [0, 1, 0],
            # # #                   [-sqr2, 0, sqr2]
            # # #                   ])
            # return np.dot(rotM,myBetterRef.T).T
            return myBetterRef
        else:
            #rotate about y axis
            # myBetterRef[:,-1]+=1.0 #shift entire molecule above XY plane
            print 'yz'
            #Rotate y 90 deg then x 90 degkills it
            #Rotate y 90 deg then x 180

            # rotM = np.array([[0.,0.,1.],
            #              [0, 1, 0],
            #              [-1.,0,0.]
            #              ])

            rotM = np.array([[1.,0.,0.],
                         [0, 0., -1.],
                         [0.,1.,0.]])

            # # rotM2 = np.array([[1.,0.,0.],
            # #              [0, -1., 0.],
            # #              [0.,0.,-1.]])

            #
            # sqr2=np.sqrt(2) / 2.
            # rotM = np.array([[sqr2, -sqr2, 0.],
            #                  [sqr2, sqr2, 0.],
            #                  [0., 0., 1.]
            #                  ])

            # rotM = np.array([[sqr2, 0., sqr2],
            #                   [0, 1, 0],
            #                   [-sqr2, 0, sqr2]
            #                   ])
            # rotM = np.array([[1., 0., 0.],
            #                   [0, sqr2, -sqr2],
            #                   [0., sqr2, sqr2]
            #                   ])
            #z rotation.
            # rotM3=np.array([[0.70710678, -0.70710678, 0.],
            #        [0.70710678, 0.70710678, 0.],
            #        [0., 0., 1.]])
            # rotM=np.array([[-1., 0., 0.],
            #        [0, -1., 0.],
            #        [0., 0., 1.]])
            # myBetterRef = np.dot(rotM, myBetterRef.T).T
            # #rotation about x
            # rotM = np.array([[1., 0., 0.],
            #                   [0, sqr2, -sqr2],
            #                   [0., sqr2, sqr2],
            #                   ])
            myBetterRef= np.dot(rotM,myBetterRef.T).T

            return myBetterRef #myRef2

    def pullTetramerRefPos(self,yz): #Eckart reference for the trimer is in an xyz file. Need just a 3xNatom array of reference structures. I can hard code this in
        """goes O1,O2,O3,O4,..H12"""
        # myRefCOM = np.array([[0.00000000E+00,  4.81355109E+00, -4.53345972E-32],
        #                    [4.16865752E+00, -2.40677554E+00, -1.38050658E-30],
        #                    [-4.16865752E+00, -2.40677554E+00, 1.18329136E-30],
        #                    [-0.00000000E+00,  0.00000000E+00,  0.00000000E+00],
        #                    [1.79529146E-16,  5.90334467E+00, -1.46596673E+00],
        #                    [-3.94430453E-31,  5.90334467E+00,  1.46596673E+00],
        #                    [5.11244645E+00, -2.95167233E+00,  1.46596673E+00],
        #                    [5.11244645E+00, -2.95167233E+00, -1.46596673E+00],
        #                    [-5.11244645E+00, -2.95167233E+00, -1.46596673E+00],
        #                    [-5.11244645E+00, -2.95167233E+00, 1.46596673E+00],
        #                    [-1.65058312E+00, -9.52964606E-01,  3.94430453E-31],
        #                    [1.65058312E+00, -9.52964606E-01, -4.93038066E-31],
        #                    [0.00000000E+00,  1.90592921E+00, -1.72916465E-32]])
        myRefCOM = np.array([[-4.64953331e+00,  1.24583870e+00,  3.55153419e-32],
                             [ 3.40369461e+00,  3.40369461e+00, -1.29965664e-30],
                             [ 1.24583870e+00, -4.64953331e+00,  1.26414130e-30],
                             [ 3.21975274e-09, -8.62730146e-10,  8.08499391e-32],
                             [-5.70219308e+00,  1.52789803e+00, -1.46596673e+00],
                             [-5.70219308e+00,  1.52789803e+00,  1.46596673e+00],
                             [ 4.17429505e+00,  4.17429505e+00,  1.46596673e+00],
                             [ 4.17429505e+00,  4.17429505e+00, -1.46596673e+00],
                             [ 1.52789803e+00, -5.70219308e+00, -1.46596673e+00],
                             [ 1.52789803e+00, -5.70219308e+00,  1.46596673e+00],
                             [ 1.34769547e+00,  1.34769547e+00, -4.12188127e-31],
                             [4.93290781e-01, -1.84098625e+00, 4.75280392e-31],
                             [-1.84098624e+00,  4.93290777e-01,  6.35582926e-32]])
        #Rotate such that Z is along OOOOPlane
        #th = np.deg2rad(90.)
        if yz:
            # #O2 being y axis
            # th = np.deg2rad(-(75.+90))
            # rotM = np.array([[np.cos(th),-np.sin(th),0],
            #                  [np.sin(th),np.cos(th),0],
            #                  [0,0,1]])
            # myRefCOM = np.dot(rotM,myRefCOM.T).T

            #rotate 90 degrees about y axis
            rotM = np.array([[0.,0.,1.],
                         [0, 1, 0],
                         [-1.,0,0.]
                         ])
            myRefCOM= np.dot(rotM,myRefCOM.T).T
            # #nullify O2 being y axis
            # th = np.deg2rad((75. + 90))
            # rotM = np.array([[np.cos(th), -np.sin(th), 0],
            #                  [np.sin(th), np.cos(th), 0],
            #                  [0, 0, 1]])
            # myRefCOM = np.dot(rotM, myRefCOM.T).T

            print 'refStructureTurned'




        # myRefCOM, extra = self.rotateBackToFrame(np.array([myRef2, myRef2]), 2, 1, 3)  # rotate reference to OOO plane
        # print 'got Eckgeometry'
        # mass=self.get_mass()
        # com = np.dot(mass[:3], myRefCOM[:3]) / np.sum(mass[:3]) #same as overal COM
        # myRefCOM-=com
        # print 'myRefCOM',myRefCOM

        # np.savetxt('myRefCOM.txt',myRefCOM)


        # #get rotation matrix for ref geom to O2 being x axis
        # newX = (myRefCOM[2-1]-myRefCOM[4-1])/la.norm(myRefCOM[2-1]-myRefCOM[4-1])
        # oldX = (myRefCOM[1-1]-myRefCOM[2-1])/la.norm(myRefCOM[1-1]-myRefCOM[2-1])
        # th = np.deg2rad(150./2.)
        # rotM = np.array([[np.cos(th),-np.sin(th),0],
        #                  [np.sin(th),np.cos(th),0],
        #                  [0,0,1]])
        # test = np.dot(rotM,myRefCOM.T).T
        # print 'tset',test
        # stop
        # np.savetxt("rotM_EckRef",rotM)
        return myRefCOM



    def eckartRotate(self,pos,justO=False,cart=False,hydro=False,yz=False): # pos coordinates = walkerCoords numwalkersxnumAtomsx3
        """Eckart Rotate method returns the transpose of the correct matrix, meaning that when one does the dot product,
        one should transpose the matrix, or do eck.dot(___)"""
        if cart or justO or hydro:
            planar = True
        else:
            planar = False
        if self.name in ProtonatedWaterTrimer:
            self.refPos = self.pullTrimerRefPos(yz)
        else:
            if self.isotope == 'notDeuterated':
                self.refPos = self.pullTetramerRefPos(yz)
            else:
                self.refPos = self.pullTetramerRefPos(True)
        if len(pos.shape)<3:
            pos=np.array([pos])
        #Center of Mass
        print 'getting mass'
        mass=self.get_mass()
        print 'got mass, here it is', mass
        if self.name in ProtonatedWaterTetramer:
            if justO: #the OOO plane
                self.refPos=self.refPos[:3]
                refCOM = np.dot(mass[:3], self.refPos) / np.sum(mass[:3])  # same as overal COM
                com = np.dot(mass[:3],pos[:,:3])/np.sum(mass[:3])
                mass = mass[:3]
                pos = pos[:,:3,:]
            elif cart: #include central oxygen
                self.refPos = self.refPos[:4]
                com = np.dot(mass[:4], pos[:, :4]) / np.sum(mass[:4])
                refCOM =  np.dot(mass[:4], self.refPos) / np.sum(mass[:4]) #same as overal COM
                mass = mass[:4]
                pos = pos[:, :4, :]
            elif hydro:
                self.refPos = self.refPos[[4-1,13-1,11-1,12-1]]
                com = np.dot(mass[[4-1,13-1,11-1,12-1]], pos[:, [4-1,13-1,11-1,12-1]]) / np.sum(mass[[4-1,13-1,11-1,12-1]])
                refCOM = np.dot(mass[[4-1,13-1,11-1,12-1]], self.refPos) / np.sum(mass[[4-1,13-1,11-1,12-1]])  # same as overal COM
                mass = mass[[4-1,13-1,11-1,12-1]]
                pos = pos[:, [4-1,13-1,11-1,12-1],:]
            else:
                com = np.dot(mass, pos) / np.sum(mass)
        elif self.name in ProtonatedWaterTrimer:
            if justO or cart:  # the OOO plane
                self.refPos = self.refPos[:3]
                refCOM = np.dot(mass[:3], self.refPos) / np.sum(mass[:3])  # same as overal COM
                com = np.dot(mass[:3], pos[:, :3]) / np.sum(mass[:3])
                mass = mass[:3]
                pos = pos[:, :3, :]
            elif hydro:
                self.refPos = self.refPos[[3 - 1, 9-1,8-1,10-1]]
                # rotate reference so that Z axis is along OOOO Plane
                com = np.dot(mass[[3-1, 9-1,8-1,10-1]], pos[:, [3-1, 9-1,8-1,10-1]]) / np.sum(
                    mass[[3-1, 9-1,8-1,10-1]])
                refCOM = np.dot(mass[[3-1, 9-1,8-1,10-1]], self.refPos) / np.sum(
                    mass[[3-1, 9-1,8-1,10-1]])  # same as overal COM
                mass = mass[[3-1, 9-1,8-1,10-1]]
                pos = pos[:, [3-1, 9-1,8-1,10-1], :]
            else:
                com = np.dot(mass, pos) / np.sum(mass)
                refCOM = np.dot(mass,self.refPos) / np.sum(mass)
        self.refPos-=refCOM
        #First Translate:
        print 'shifting molecules'
        ShiftedMolecules=pos-com[:,np.newaxis,:]
        #Equation 3.1 in Eckart vectors, Eckart frames, and polyatomic molecules - James D. Louck and Harold W. Galbraith
        print 'starting mathy math'
        asdf = np.sum(ShiftedMolecules[:,:,:,np.newaxis]*self.refPos[np.newaxis,:,np.newaxis,:]*mass[np.newaxis,:,np.newaxis,np.newaxis],axis=1)
        myF = np.transpose(asdf,(0,2,1))
        if planar:
            print 'planar'
            # myF[:,-1] = np.cross(myF[:,0],myF[:,1])
            ind=np.where(~myF.any(axis=1))[0]
            myFp=np.copy(myF[:,:2])
            myFF = np.matmul(myFp,myFp.transpose(0,2,1))
            peval,pevecs = la.eigh(myFF)
            invRootDiagF2 = 1/np.sqrt(peval)
            axaxxa = np.where(np.isnan(invRootDiagF2))
            if len(axaxxa[0]) > 0:
                print "PLANARBADBADBAD"
                print len(axaxxa[0])
                print axaxxa[0]
                print ShiftedMolecules[axaxxa[0]]
                fileout = open("badEck_planar.xyz", "w+")
                self.printCoordsToFile(ShiftedMolecules[axaxxa[0]], fileout)
                # raise ZeroDivisionError("this is planar bad, dude")
            pevecsT = np.transpose(pevecs,(0,2,1))
            invRootF2 = np.matmul(invRootDiagF2[:, np.newaxis, :]*pevecs,pevecsT)  # -bigEvecs
            eckVecs2 = np.zeros((len(invRootDiagF2),3,3))
            eckVecs2[:,:,:2] =  np.matmul(np.transpose(myFp, (0, 2, 1)), invRootF2) #.transpose((0,2,1))
            eckVecs2[:,:,-1] = np.cross(eckVecs2[:,:,0],eckVecs2[:,:,1])
            print 'done with planar'
        else:
            print 'not planar'
            myFF = np.matmul(myF, asdf)
            bigEvals,bigEvecs=la.eigh(myFF)
            bigEvecsT=np.transpose(bigEvecs,(0,2,1))
            invRootDiagF2 = 1.0 / np.sqrt(bigEvals)
            axaxxa=np.where(np.isnan(invRootDiagF2))
            if len(axaxxa[0]) > 0:
                print "BADBADBADBADBAD"
                print len(axaxxa[0])
                print axaxxa[0]
                print ShiftedMolecules[axaxxa[0]]
                fileout = open("badEck.xyz","w+")
                self.printCoordsToFile(ShiftedMolecules[axaxxa[0]], fileout)
                # raise ZeroDivisionError("this is bad, dude")

            invRootF2=np.matmul(invRootDiagF2[:,np.newaxis,:]*-bigEvecs,-bigEvecsT,) #-bigEvecs
            eckVecs2 = np.matmul(np.transpose(myF,(0,2,1)),invRootF2)
        print 'done'
        mas = np.where(np.around(la.det(eckVecs2))==-1.0)
        if len(mas[0])!=0:
            killList2=mas
            # eckVecs2[mas,:,-1] *= -1.0
            eckVecs2[mas,-1] *= -1.0 #multiply row by -1, but actually column since this is transposed
            minus = len(mas[0])
        else:
            minus = 0
            killList2=mas[0]
        plus=len(ShiftedMolecules)-minus
        print 'Plus rotation: ',plus
        print 'Inverted Rotation: ',minus
        killList2=0
        return com, eckVecs2 , killList2

    def getInitialCoordinates(self):
        #looks up the path of the coordinates for the starting positions of the walkers and makes nWalkers copies of them and then returns that as an self.nWalkers,self.nAtoms, 3 dimensional array                      
        print 'there are ', self.nAtoms , 'atoms in ', self.name
        if self.name in DeprotonatedWaterDimer:
            coordsMin=np.array([ [   0.000000000000000 ,  0.000000000000000 ,  0.000000000000000],
                                 [  2.306185590098382 ,  0.000000000000000 ,  0.000000000000000],
                                 [  2.749724314110769 ,  -1.765018349357672 ,  0.000000000000000],
                                 [   -2.306185590098382 ,  0.000000000000000 ,  0.000000000000000],
                                 [   -2.749724314110769 ,  -1.765018349357672 ,  0.000000000000000]
                                 ])
            coordsMinRotated=np.array([[0.0000000,   0.000000000000000 ,  0.000000000000000  ],
                                       [0.0000000,  2.306185590098382 ,  0.000000000000000   ],
                                       [0.0000000,  2.749724314110769 ,  -1.765018349357672  ],
                                       [0.0000000,   -2.306185590098382 ,  0.000000000000000 ],
                                       [0.0000000,   -2.749724314110769 ,  -1.765018349357672]
                                       ])

            coordsc2v=np.array([ [   0.000000000000000 ,  0.000000000000000 ,  0.000000000000000],
                                 [  -2.306185590098382 ,  0.000000000000000 ,  0.000000000000000],
                                 [  -2.749724314110769 ,  1.765018349357672 ,  0.000000000000000],
                                 [   2.306185590098382 ,  0.000000000000000 ,  0.000000000000000],
                                 [   2.749724314110769 ,  1.765018349357672 ,  0.000000000000000]])
            coordsc2h=np.array([ [   0.000000000000000 ,   0.000000000000000       ,  0.000000000000000],
                                 [  -2.304566686034061 ,   0.000000000000000       ,  0.000000000000000],
                                 [  -2.740400260927908 ,   1.0814221449986587E-016 ,  -1.766154718409233],
                                 [   2.304566686034061 ,   0.000000000000000       ,  0.000000000000000],
                                 [  2.740400260927908  ,   1.0814221449986587E-016 ,  1.766154718409233]])
            coordsSaddle=np.array([[ 0.000000000000000   , 0.000000000000000 ,  0.000000000000000 ],
                                   [ -2.303263755760085 , 0.000000000000000 ,  0.000000000000000 ],
                                   [ -2.720583162407882 , 1.129745554266140 ,  -1.363735721982301],
                                   [ 2.303263755760085  , 0.000000000000000 ,  0.000000000000000 ],
                                   [ 2.720583162407882  , 1.129745554266140 ,  1.363735721982301 ]])

            coordsAsymStart=np.array([[ 0.2981678882048853 ,   -2.4557992072743176E-002,  -5.5485232545510215E-002],
                                      [ -2.354423404994569 ,     0.000000000000000     ,    0.000000000000000],
                                      [ -2.858918674095194 ,     1.111268022307282     ,   -1.352651141853729],
                                      [  2.354423404994569 ,     0.000000000000000     ,    0.000000000000000],
                                      [  2.771741580470489 ,     1.236107563104921     ,    1.482886181959795]])



            coordsAsymStart2=np.array([[   3.52234 ,    1.01649  ,    1.28596],
                                       [   2.35442 ,    0.00000  ,    0.00000],
                                       [   0.29817 ,   -0.02456  ,   -0.05549],
                                       [  -2.35442 ,    0.00000  ,    0.00000],
                                       [  -2.85892 ,    1.11127  ,   -1.35265]])



            coords=coordsAsymStart


        #            coords=np.array([[ 0.26591125,  0.07072797, -0.02256279],
        #WW         WW                [ 0.55610034,  2.83109547,  0.14883552],
        #  W       W                  [ 1.50122114,  1.22631416, -0.59092507],
        #   W  W  W                   [-0.11962985, -1.87021212,  0.22794889],
        #    W   W                    [ 1.20503929, -0.77837156,  0.71051114]])

        elif self.name in ProtonatedWaterTrimer:
            #Eq coords:
            #O
            #H
            #H
            #O
            #H
            #H
            #H
            #O
            #H
            #H
            eqCoords = np.array([[2.07091,-0.46191,-0.00993],
                                 [2.56267,-0.75858,0.76451],
                                 [2.70113,-0.40578,-0.73813],
                                 [0.00000,0.91527,-0.05817],
                                 [0.00000,1.67720,0.53729],
                                 [-0.87302,0.35992,0.01707],
                                 [0.87302,0.35993,0.01707],
                                 [-2.07092 ,-0.46190, -0.00993],
                                 [-2.70115, -0.40575, -0.73811],
                                 [-2.56265, -0.75862,0.76451]])
            coords = eqCoords*ang2bohr
            print 'perturb the equilibrium coordinates by a little'
        elif self.name in Water:
            eqCoords = \
                np.array([[-4.3018530574,3.6598148518, -0.081920343],
                          [-4.306042519,2.6986842642, -0.438125298],
                          [-5.2421388362,3.878384798001161,0.0488112791]]) #in angstroms already
            coords = eqCoords
        else:
            print 'This is all zeros famalama'
            eqCoords = np.array([[0, 0, 0],
                                [0, 0, 0],
                                [0, 0, 0],
                                [0, 0, 0],
                                [0, 0, 0],
                                [0, 0, 0],
                                [0, 0, 0],
                                [0, 0, 0],
                                [0, 0, 0],
                                [0, 0, 0],
                                [0, 0, 0],
                                [0, 0, 0],
                                [0, 0, 0]]
                                )
            coords = eqCoords * ang2bohr

        return coords
    def getMassSharedProton(self):
        if self.isotope=='DeuteratedTwice':
            return massH
        if self.isotope=='notDeuterated':
            return massH
        if self.isotope=='fullyDeuterated':
            return massD
    def getMassOuterProton(self):
        if self.isotope=='DeuteratedTwice':
            return massD
        if self.isotope=='notDeuterated':
            return massH
        if self.isotope=='fullyDeuterated':
            return massD


    def get_mass(self):
        mass=np.zeros((self.nAtoms))
        if self.name in ProtonatedWaterTrimer:
            mass = np.array([massO,massO,massO,massH,massH,massH,massH,massH,massH,massH])
            if self.isotope == 'DeuteratedOnce_eigen':
                mass[9] = massD
                self.names = ['O', 'O', 'O', 'H', 'H', 'H', 'H', 'H', 'H', 'D']
            elif self.isotope == 'DeuteratedOnce_fw':
                mass[3] =massD
                self.names = ['O', 'O', 'O', 'D', 'H', 'H', 'H', 'H', 'H', 'H']
            elif self.isotope == 'notDeuterated':
                print ':)'
            elif self.isotope == 'fullyDeuterated':
                mass = np.array([massO, massO, massO, massD, massD, massD, massD, massD, massD, massD])
                self.names = ['O', 'O', 'O', 'D', 'D', 'D', 'D', 'D', 'D', 'D']
            elif self.isotope == 'notDeuteratedOnce_eigen':
                self.names = ['O', 'O', 'O', 'D', 'D', 'D', 'D', 'D', 'D', 'H']
                mass = np.array([massO, massO, massO, massD, massD, massD, massD, massD, massD, massH])
            elif self.isotope == 'notDeuteratedOnce_fw':
                self.names = ['O', 'O', 'O', 'H', 'D', 'D', 'D', 'D', 'D', 'D']
                mass = np.array([massO, massO, massO, massH, massD, massD, massD, massD, massD, massD])
            elif self.isotope == 'DeuteratedOnce_hydronium':
                self.names = ['O', 'O', 'O', 'H', 'H', 'H', 'H', 'D', 'H', 'H']
                mass[7] = massD
            elif self.isotope == 'notDeuteratedOnce_hydronium':
                self.names = ['O', 'O', 'O', 'D', 'D', 'D', 'D', 'H', 'D', 'D']
                mass = np.array([massO, massO, massO, massH, massD, massD, massD, massH, massD, massD])
        elif self.name in ProtonatedWaterTetramer:
            mass = np.array([massO, massO, massO, massO, massH, massH, massH, massH, massH, massH, massH,massH,massH])
            if self.isotope == 'DeuteratedOnce_eigen':
                mass[10] = massD
                self.names = ['O', 'O', 'O', 'O','H','H','H','H','H','H','D','H','H']
            #elif self.isotope == 'DeuteratedOnce_fw':
            #    mass[3] = massD
            #    self.names = ['O', 'O', 'O', 'H', 'D', 'D', 'D', 'D', 'D', 'D']
            #elif self.isotope == 'notDeuterated':
            #    print ':)'
            elif self.isotope == 'fullyDeuterated':
                mass = np.array([massO, massO, massO,massO, massD, massD, massD, massD, massD, massD, massD,massD,massD])
            elif self.isotope == 'notDeuteratedOnce_eigen':
                self.names = ['O', 'O', 'O','O','D', 'D', 'D', 'D', 'D', 'D', 'H','D','D']
                mass = np.array([massO, massO, massO, massO,massD, massD, massD, massD, massD, massD, massH, massD, massD])
            elif self.isotope == 'notDeuteratedOnce_fw':
                self.names = ['O', 'O', 'O','O','H', 'D','D','D', 'D', 'D', 'D', 'D', 'D']
                mass = np.array([massO, massO, massO,massO, massH, massD, massD, massD, massD, massD, massD, massD, massD])
            #elif self.isotope == 'DeuteratedOnce_hydronium':
            #    self.names = ['O', 'O', 'O', 'H', 'H', 'H', 'H', 'D', 'H', 'H']
            #    mass[7] = massD
            #elif self.isotope == 'notDeuteratedOnce_hydronium':
            #    self.names = ['O', 'O', 'O', 'D', 'D', 'D', 'D', 'H', 'D', 'D']
            #    mass = np.array([massO, massO, massO, massH, massD, massD, massD, massH, massD, massD])
        elif self.name in Water:
            mass = np.array([massO,massH,massH])
        return mass*massConversionFactor

    def calcReducedmass(self,x): #Resource: Molecular Vibrations wdc - Appendix VI
        # for the shared proton stretch, page 42 in notebook #2
        massOamu=massO*massConversionFactor

        if self.name in ProtonatedWaterDimer and self.state==1 and self.surfaceName=='SharedProton':
            ##COMwat1=1.0/massWater*(massO*x[:,0,:]+massH*(x[:,3,:]+x[:,4,:]))
            ##COMwat2=1.0/massWater*(massO*x[:,1,:]+massH*(x[:,5,:]+x[:,6,:]))
            ##U=x[:,2,:]-COMwat1                                              
            ##V=x[:,2,:]-COMwat2                                       

            U=x[:,2,:]-x[:,1,:]
            V=x[:,2,:]-x[:,0,:]
            massHamu=self.getMassSharedProton()*massConversionFactor
            magU=np.sqrt(U[:,0]**2+U[:,1]**2+U[:,2]**2)
            magV=np.sqrt(V[:,0]**2+V[:,1]**2+V[:,2]**2)

            costheta= np.diag(np.dot(U,V.T))/(magU*magV)
            #mass=1.0/(2.0*((1.0/(massOamu+massHamu+massHamu))+((1-costheta)/(massHamu))))                                        
            #corresponds to calcrncom                        

            mass=1.0/(2.0*((1.0/(massOamu))+((1-costheta)/(massHamu))))
            #print 'average mass', np.average(mass)
            #Mass of water or Mass of O??                    
        elif self.name in ProtonatedWaterDimer and self.state==1 and self.surfaceName=='StretchAntiIn':

            massHamu=self.getMassOuterProton()*massConversionFactor

            U=x[:,0,:]-x[:,3,:]
            V=x[:,0,:]-x[:,4,:]

            magU=np.sqrt(U[:,0]**2+U[:,1]**2+U[:,2]**2)
            magV=np.sqrt(V[:,0]**2+V[:,1]**2+V[:,2]**2)
            costhetaWat1= np.diag(np.dot(U,V.T))/(magU*magV)

            U=x[:,1,:]-x[:,5,:]
            V=x[:,1,:]-x[:,6,:]
            magU=np.sqrt(U[:,0]**2+U[:,1]**2+U[:,2]**2)
            magV=np.sqrt(V[:,0]**2+V[:,1]**2+V[:,2]**2)
            costhetaWat2= np.diag(np.dot(U,V.T))/(magU*magV)

            g=( (1.0/massOamu)+(1.0/massHamu) )  +  1.0/2.0*((costhetaWat1/massOamu)+(costhetaWat2/massOamu))

            mass= 1.0/g

        elif self.name in DeprotonatedWaterDimer and self.surfaceName=='SharedProton' and self.state==1:
            massHamu=self.getMassSharedProton()*massConversionFactor
            U=x[:,1,:]-x[:,0,:]
            V=x[:,3,:]-x[:,0,:]
            magU=np.sqrt(U[:,0]**2+U[:,1]**2+U[:,2]**2)
            magV=np.sqrt(V[:,0]**2+V[:,1]**2+V[:,2]**2)
            costheta= np.diag(np.dot(U,V.T))/(magU*magV)
            mass=1.0/(2.0*((1.0/(massOamu))+((1.000000-costheta)/(massHamu))))
            #print 'average mass', np.average(mass)

        elif self.name in DeprotonatedWaterDimer and self.surfaceName=='Z-displacement' and self.state==1:
            massHamu=self.getMassSharedProton()*massConversionFactor
            mass=1.0/((1.0/massHamu)+(1.0/(2.0*massOamu)))

        elif  self.name in DeprotonatedWaterDimer and self.surfaceName=='LocalOHStretch' and self.state==1:
            massHamu=self.getMassOuterProton()*massConversionFactor
            mass=(massHamu*massOamu)/(massHamu+massOamu)

        elif self.name in DeprotonatedWaterDimer and self.surfaceName=='OHStretchAnti' and self.state==1:
            massHamu=self.getMassOuterProton()*massConversionFactor
            g=(1.0/massHamu)+ (1.0/massOamu)
            mass=1.0/g #(massHamu*massOamu)/(massHamu+massOamu)
            #print 'average mass', np.average(mass)
        elif self.name in ProtonatedWaterDimer and self.state==0:
            m1=self.getMassOuterProton()*massConversionFactor
            m2=2*(massO)*conversionFactor
            #m1=massH*conversionFactor
            mass=m1*m2/(m1+m2)
            print 'why are you calculating the reduced mass on the ground state?'  , end
        elif self.name in DeprotonatedWaterDimer and self.state==0:
            m2=2*(massO)*conversionFactor #Supposed to be massConversionFactor?
            m1=self.getMassOuterProton()*massConversionFactor
            mass=m1*m2/(m1+m2)
            print 'why are you calculating the reduced mass on the ground state?'  , end

        else:
            print 'not implemented for ', self.name , 'and', self.state, 'and', self.surfaceName,end

        return mass

    def sortProtons(self,x):
        #Defunct!
        #defunct
        #This is a fun function!  It was written to correct for the isomerization of H3O2-.  The isomerization was a problem
        #Because the anti symmetric stretch (and actually all of the internal coordinates) were defined by the atom positions 
        #in the coordinate array.  
        #Input is the coordinate array, x.  This function DOES NOT make a deep copy of the coordinate array passed in.  So the 
        #swapping that takes place CHANGES the original array.  For this reason, the sortProton function only needs to be called
        #once each time step, but there are many places where we'd need the protons to be sorted appropriately and automatically
        #so we'll just have to work on speed ups in here rather than frugal calling of this method.
        #print 'sorting protons'
        if self.name not in DeprotonatedWaterDimer:
            print 'failed! wrong molecule!', end
        H0O1=self.bondlength(x,atom1=0 ,atom2=1)
        H0O3=self.bondlength(x,atom1=0 ,atom2=3)
        H2O1=self.bondlength(x,atom1=2 ,atom2=1)
        H2O3=self.bondlength(x,atom1=2 ,atom2=3)
        H4O1=self.bondlength(x,atom1=4 ,atom2=1)
        H4O3=self.bondlength(x,atom1=4 ,atom2=3)
        midpoint=(x[:,1,:]+x[:,3,:])/2.0
        H0M=self.mag(x[:,0,:]-midpoint)
        H2M=self.mag(x[:,2,:]-midpoint)
        H4M=self.mag(x[:,4,:]-midpoint)
        DistanceToMidPt=np.array(zip(H0M,H2M,H4M))

        H0H2=self.bondlength(x,atom1=0 ,atom2=2)
        H0H4=self.bondlength(x,atom1=0 ,atom2=4)
        H2H4=self.bondlength(x,atom1=2 ,atom2=4)

        aveOH0=np.average(zip(H0O1,H0O3),axis=1)
        aveOH2=np.average(zip(H2O1,H2O3),axis=1)
        aveOH4=np.average(zip(H4O1,H4O3),axis=1)

        averagesAll=np.array(zip(aveOH0,aveOH2,aveOH4))

        OH1DistAll=np.array(zip(H0O1,H2O1,H4O1))
        sortByAveragesAll=np.argsort(averagesAll,axis=1)
        sortByDistToMP=np.argsort(DistanceToMidPt,axis=1)
        #checkSortNeed=np.logical_not(sortByAveragesAll[:,0]==0)





        checkSortNeed=np.logical_not(sortByDistToMP[:,0]==0)

        #print 'checkSortNeed sample', sortByAveragesAll[0:10],checkSortNeed[0:10]
        listOfSwapped=[]
        if np.any(checkSortNeed):
            #print np.where(checkSortNeed),
            #fileout=open('swappyCoordinates.xyz','a')
            for m in np.where(checkSortNeed):
                n=m[0]
                sortByAverages=sortByAveragesAll[n]
                sortByMPdist=sortByDistToMP[n]

                protonIndices=np.array([0,2,4])
                averages=averagesAll[n]
                OH1Dist=OH1DistAll[n]

                #centralProtonIndex=protonIndices[sortByAverages][0]
                centralProtonIndex=protonIndices[sortByMPdist][0]
                #outerOHIndex=protonIndices[sortByAverages][1:]
                outerOHIndex=protonIndices[sortByMPdist][1:]
                #outerOHDist=OH1Dist[sortByAverages][1:]
                outerOHDist=OH1Dist[sortByMPdist][1:]

                sortByHO1Dist=np.argsort(outerOHDist)
                H2primeIndex=outerOHIndex[sortByHO1Dist][0]
                H4primeIndex=outerOHIndex[sortByHO1Dist][1]

                positionH0prime=np.argmin(averages)*2
                positionH2prime=np.argmin(OH1Dist)*2
                positionH4prime=np.argmax(OH1Dist)*2
                if not (centralProtonIndex==0 and H2primeIndex==2 and H4primeIndex==4):
                    print 'walkN',n,
                    print 'AverageOH: ', averages,
                    print 'O1-H dis: ',OH1Dist,
                    print 'O2-H dis: ',np.array([H0O3[n],H2O3[n],H4O3[n]]),
                    print 'HH distances', H0H2[n],H0H4[n],H2H4[n],
                    print 'Distance to midpoint:', DistanceToMidPt[n],
                    print 'Order from sort by averages',sortByAverages,'order from sort by dist to MP',sortByDistToMP[n]
                    #print centralProtonIndex,H2primeIndex,HprimeIndex, '=?=',            
                    #print positionH0prime,positionH2prime, positionH4prime
                    #print x[n]
                    #self.printCoordsToFile([x[n]],fileout)
                    print 'no swap!'
                    #x[n][[0,1,2,3,4]]=x[n][[centralProtonIndex,1,H2primeIndex,3,H4primeIndex]]
                    #print 'yes swap!'
                    # self.printCoordsToFile([x[n]],fileout)
                    #print 'finally \n',x[n]
                    listOfSwapped.append(n)

        #I might not have to sort all of them!

        #fileout.close
        return listOfSwapped

    def printCoordsToFile(self,x,fileout):
        if len(x) == 1:
            x = np.array([x,x])
        au2ang=0.529177249
        if self.name in ProtonatedWaterTrimer:
            atomStr=['O','O','O','H','H','H','H','H','H','H']
        elif self.name in ProtonatedWaterTetramer:
            atomStr = ['O', 'O', 'O','O', 'H', 'H', 'H', 'H', 'H', 'H', 'H','H','H']  # needs to be updated.  Li, N, Be are for distinguishing atoms in h3o2
        for i in x:
            fileout.write("%d\nwriteout\n" % len(atomStr))
            for atmn, atm in enumerate(i):
                fileout.write("%s %5.12f %5.12f %5.12f\n" % (atomStr[atmn],atm[0],atm[1],atm[2]))
            fileout.write("\n")
        fileout.close()

    def calcAverageInternalCoordinates(self,x,dw):
        if self.name in DeprotonatedWaterDimer:
            #listOfSwapped=self.sortProtons(x)
            #print 'check swapped:', x[listOfSwapped]
            #print '\n ^^^ no really check it! ^^^ \n'
            #Local Intramolecular:
            HO1=self.bondlength(x,atom1=1, atom2=2)
            HO2=self.bondlength(x,atom1=3, atom2=4)
            magAsOH=np.absolute(HO2-HO1)
            #shared proton
            BO1=self.bondlength(x,atom1=0, atom2=1)
            BO2=self.bondlength(x,atom1=0, atom2=3)
            #Intermolecular
            OO= self.bondlength(x,atom1=1, atom2=3)
            HxH=self.bondlength(x,atom1=2, atom2=4)
            H2B=self.bondlength(x,atom1=0, atom2=2)
            H4B=self.bondlength(x,atom1=0, atom2=4)

            #PseudoRock
            Rock1=self.bondlength(x,atom1=1, atom2=4)-self.bondlength(x,atom1=1, atom2=0)
            Rock3=self.bondlength(x,atom1=3, atom2=2)-self.bondlength(x,atom1=3, atom2=0)
            #calc max, min, average, expectationValue, and  std of each
            Values=np.zeros((11,5))
            Values[:,0]=    np.min([HO1,HO2,magAsOH,BO1,Rock1,Rock3,BO2,OO,HxH,H2B,H4B],axis=1)
            Values[:,1]=    np.max([HO1,HO2,magAsOH,BO1,Rock1,Rock3,BO2,OO,HxH,H2B,H4B],axis=1)
            Values[:,2]=np.average([HO1,HO2,magAsOH,BO1,Rock1,Rock3,BO2,OO,HxH,H2B,H4B],axis=1)
            Values[:,3]=np.average([HO1,HO2,magAsOH,BO1,Rock1,Rock3,BO2,OO,HxH,H2B,H4B],axis=1,weights=dw)
            Values[:,4]=    np.std([HO1,HO2,magAsOH,BO1,Rock1,Rock3,BO2,OO,HxH,H2B,H4B],axis=1)

            internalNames=['HO1','HO2','|r2-r1|','BO1','BO2','Rock1','Rock3','OO','HxH','H2B','H4B']

        return internalNames,Values

    def calcRocks(self,x):
        BO1=self.bondlength(x,atom1=1, atom2=0)
        BO2=self.bondlength(x,atom1=3, atom2=0)
        Rock1=self.bondlength(x,atom1=1, atom2=4)-BO1
        Rock3=self.bondlength(x,atom1=3, atom2=2)-BO2
        return Rock1,Rock3
    def calcOHBondLenghts(self,x):
        if self.name in DeprotonatedWaterDimer:
            #swap=self.sortProtons(x)
            OH1,OH2=self.calcLocalOH(x)
            return OH1,OH2,self.calcStretchSym(x),self.calcStretchAnti(x)
        else:
            OH1=self.bondlength(x,atom1=1, atom2=2)
            OH2=self.bondlength(x,atom1=3, atom2=4)
            invSqrt2=1.0/np.sqrt(2.0)
            OHSym=invSqrt2*(OH1+OH2)
            OHAsym=invSqrt2*(OH1-OH2)
            return OH1, OH2, OHSym, OHAsym

    def calcSharedProtonDisplacement(self,x):
        if self.name in ProtonatedWaterDimer:
            r1=self.bondlength(x,atom1=2, atom2=1)
            r2=self.bondlength(x,atom1=2, atom2=0)
            return r2-r1
        elif self.name in DeprotonatedWaterDimer:
            r1=self.bondlength(x,atom1=0, atom2=1)
            r2=self.bondlength(x,atom1=0, atom2=3)#ha                          
            return r2-r1

    def calcStretchAntiIn(self,x):

        if self.name in ProtonatedWaterDimer:
            r1=self.bondlength(x,atom1=0, atom2=3)
            r2=self.bondlength(x,atom1=0, atom2=4)
            r3=self.bondlength(x,atom1=1, atom2=5)
            r4=self.bondlength(x,atom1=1, atom2=6)
            return 0.5*(r1+r2-r3-r4)

    def calcStretchAnti(self,x):
        #print 'calculating StretchAnti'
        if self.name in DeprotonatedWaterDimer:
            #print 'calculatingStretchAnti',
            #listOfSwapped=self.sortProtons(x)
            #listOfSwapped=[]
            r1=self.bondlength(x,atom1=1,atom2=2)
            r2=self.bondlength(x,atom1=3,atom2=4)

            #print '3 averages',np.average(r1), np.average(r2) , np.average(r1-r2), ' --- '
            #return (r1-r2), listOfSwapped
            return np.sqrt(0.5)*(r1-r2)#, listOfSwapped

    def calcStretchSym(self,x):
        if self.name in DeprotonatedWaterDimer:
            r1=self.bondlength(x,atom1=1,atom2=2)
            r2=self.bondlength(x,atom1=3,atom2=4)
            return np.sqrt(1.0/2.0)*(r1+r2)
    def calcLocalOH(self,x):
        if self.name in DeprotonatedWaterDimer:
            return self.bondlength(x,atom1=1,atom2=2),self.bondlength(x,atom1=3,atom2=4)


    def calcEckRotZdisplacement(self,x):
        newcoords=self.eckartRotate(x)
        #must be eckart rotated to a reference coordinate set that has the Os along the x axis
        return newcoords[:,0,0]


    def calcZdisplacement(self,x):
        # define midpoint
        OOMP=(x[:,1]+x[:,3])/2.0
        # define vector between O1 and MP
        OMPVec=x[:,1]-OOMP
        # define H-OOMPvector
        HMP= x[:,0]-OOMP
        # calc projection of H onto O-MP vector, and the tada!
        projection=np.zeros(x.shape[0])
        for i,(omp,hmp) in enumerate(zip(OMPVec,HMP)):
            projection[i]=np.dot(hmp,omp)
            projection[i]=projection[i]/np.linalg.norm(omp)
        return projection



    def calcCartesianSharedProtonDisplacement(self,x):

        #returns the projection of the shared proton on a molecular frame of reference.
        #Z is along the OO vector
        #Y is perpendicular to Z and X' vector (X'=the vector the bisects the OH-OH dihedral).
        #X is on the Z and X' plane and is perpendicular to the Y and Z axes.

        #  H2                      H4      X
        #   \                     /        |
        #    \          H0      /          |
        #     O1------MP-----O3            -----> Z
        #
        # Sort of a Neuman Projection:
        #
        #    Y      H2    X'
        #      \    |    /
        #        \  |  /
        #           O1 -----H4
        #
        #

        # define midpoint
        OOMP=(x[:,1]+x[:,3])/2.0
        # define vector between O1 and MP
        OMPVec=x[:,1]-OOMP
        #ZVector=OMPVec/np.linalg.norm(OMPVec,axis=1)[:,None]
        OOVec=x[:,3]-x[:,1]
        ZVector=OOVec/np.linalg.norm(OOVec,axis=1)[:,None]
        # define H-OOMPvector
        HMP= x[:,0]-OOMP
        # calc projection of H onto O-MP vector, and the tada!
        projectionZ=np.zeros(x.shape[0])
        for i,(omp,hmp) in enumerate(zip(OMPVec,HMP)):
            projectionZ[i]=np.dot(hmp,omp)
            projectionZ[i]=projectionZ[i]/np.linalg.norm(omp)
        normO1H2=(x[:,2]-x[:,1])/np.linalg.norm(x[:,2]-x[:,1],axis=1)[:,None]
        normO3H4=(x[:,4]-x[:,3])/np.linalg.norm(x[:,4]-x[:,3],axis=1)[:,None]
        XVectorPrime=(normO1H2+normO3H4)/np.linalg.norm(normO1H2+normO3H4,axis=1)[:,None]

        YVector=np.cross(XVectorPrime,ZVector)
        YVector=YVector/np.linalg.norm(YVector,axis=1)[:,None]
        #        print 'y vec normalized?',np.linalg.norm(YVector,axis=1)
        projectionY=np.zeros(x.shape[0])
        for k,(yvec,hmp) in enumerate(zip(YVector,HMP)):
            projectionY[k]=np.dot(hmp,yvec)

        XVector=np.cross(ZVector,YVector)
        XVector=XVector/np.linalg.norm(XVector,axis=1)[:,None]

        #        print 'normalized?' ,np.linalg.norm(XVector,axis=1)
        projectionX=np.zeros(x.shape[0])

        for j,(xvec,hmp) in enumerate(zip(XVector,HMP)):
            projectionX[j]=np.dot(hmp,xvec)


        ###print 'is cart projection of shared proton even working???'
        ###print 'Hdisp',HMP[0]
        ###print 'Vectors \n',XVector[0],'\n', YVector[0],'\n', ZVector[0]
        ###
        ###print 'projections:',projectionX[0],projectionY[0],projectionZ[0]
        ###
        ###print 'min',np.min(projectionX),np.min(projectionY)
        ###print 'max',np.max(projectionX),np.max(projectionY)

        return projectionX, projectionY, projectionZ

    def calcRn(self,x):
        #listOfSwapped=self.sortProtons(x)
        if self.surfaceName=='SharedProton':
            rncoord=self.calcSharedProtonDisplacement(x)
        elif self.surfaceName=='StretchAntiIn':
            rncoord=self.calcStretchAntiIn(x)
            return rncoord
        elif self.surfaceName=='OHStretchAnti':
            #print 'from calcRN',
            rncoord=self.calcStretchAnti(x)

        elif self.surfaceName=='LocalOHStretch':
            rncoord=AVERAGEGroundStateOH-self.calcLocalOH(x)

        elif self.surfaceName=='Z-displacement':
            rncoord=self.calcZdisplacement(x)
        elif self.surfaceName=='EckRotZ-displacement':
            rncoord=self.calcEckRotZdisplacement(x)
        else:
            print 'surface name is not found!', end

        return rncoord#,listOfSwapped

    def bondlength(self,pos,atom1,atom2):
        length=(pos[:,atom1,0]-pos[:,atom2,0])**2+(pos[:,atom1,1]-pos[:,atom2,1])**2+(pos[:,atom1,2]-pos[:,atom2,2])**2
        length=np.sqrt(length)
        return length

    def bondAngle(self,pos, atom1, atom2, atom3):
        #finds the angle between 3 atoms with atom2 in the center                                                                                                       
        angle=[]
        for molecule in pos:
            a=molecule[atom1,:]-molecule[atom2,:]
            b=molecule[atom3,:]-molecule[atom2,:]
            angle.append(np.arccos(np.dot(b,a)/(self.mag(np.array([a]))*self.mag(np.array([b])))))
        angle=np.array(angle)
        #print 'ANGLE SHAPE, ', np.shape(angle.reshape(angle.shape[0]))
        return angle.reshape(angle.shape[0])

    """def bondAng(three):
        atm1 = np.array(three[:3])
        atm2 = np.array(three[3:6])
        atm3 = np.array(three[6:])
        left = atm1 - atm2
        right = atm3 - atm2
        return np.arccos(np.dot(left, right) / (la.norm(left)) * (la.norm(right)))"""

    def calcTorsion(self,posList,pointList=None):
        if pointList==None:
            if self.name in DeprotonatedWaterDimer:
                atom0, atom1, atom2,atom3=2,1,3,4
            elif self.name in ProtonatedWaterDimer:
                atom0, atom1, atom2,atom3=3,0,1,5
            else:
                atom0, atom1, atom2,atom3=0,1,2,3
        else:
            atom0=pointList[0]
            atom1=pointList[1]
            atom2=pointList[2]
            atom3=pointList[3]
            print 'calculating the torsion angle between ', atom0,atom1,atom2, atom3
        torsion=np.zeros(posList.shape[0])
        for i,pos in enumerate(posList):
            a=pos[atom2]-pos[atom1]
            b=pos[atom1]-pos[atom0]
            c=pos[atom3]-pos[atom2]
            a=a/np.sqrt(np.sum(a**2))
            b=b/np.sqrt(np.sum(b**2))
            c=c/np.sqrt(np.sum(c**2))
            n1=np.cross(b,a)
            n2=np.cross(a,c)
            m1=np.cross(n1,a)
            x1=np.dot(n1,n2)
            y1=np.dot(m1,n2)
            torsion[i]=np.arctan2(y1,x1)
            if torsion[i]<0.0:
                torsion[i]=2*np.pi+torsion[i]
            if torsion[i]>np.pi:
                torsion[i]=2*np.pi-torsion[i]
        #print 'torsion max, min, ave', np.max(torsion*rad2deg), np.min(torsion*rad2deg), np.average(torsion*rad2deg)
        return torsion,[0.0, np.pi]

    def symmetrizeCoordinates(self,x,dw):
        if self.name in DeprotonatedWaterDimer:
            xSym=np.concatenate((x, #1     
                                 self.exchange(x,[(1,3),(2,4)]))) #2
            dwSym=np.concatenate((dw,dw))
        return xSym,dwSym

    def exchange(self, x, listOfExchanges):
        #This function takes in x and returns xprime where xprime has the atoms exchanged according to the listOfExchanges which is a list of tuples that indicates which atoms should be exchanged.                                                                                                                                                                        
        xprime=1.0*x #deep copy of x                                                                                                                                    
        for (exAtom1,exAtom2) in listOfExchanges:
            temp=1.0*xprime[:,exAtom1]
            xprime[:,exAtom1]=xprime[:,exAtom2]
            xprime[:,exAtom2]=temp
        return xprime


    def mag(self,xList):
        magnitude=np.zeros(xList.shape[0])
        for i,x in enumerate(xList):
            magnitude[i]=np.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
        return magnitude

