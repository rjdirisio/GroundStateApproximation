import numpy as np
import sys
from numpy import linalg as la
import os
import time
import Plot
import copy
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

        finalCoords = np.matmul(r, xaxtrCoordz.transpose(0, 2, 1)).transpose(0, 2, 1)
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


    def SymInternals(self,x,rotato=False,weights=0):
        #print 'called SymInternals'
        #print 'returning values in bohr [or radians]'
        if self.name in DeprotonatedWaterDimer:
            internals= self.SymInternalsH3O2minus(x)
            #self.internalNames=internalNames

            return internals#,internalNames

        #elif self.name in ProtonatedWaterDimer:
        #    return self.SymInternalsH5O2plus(x)
        elif self.name in ProtonatedWaterTrimer:
            if rotato:
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
            if rotato:
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


    def xyzSharedHydrogens(self,atmnm,xx):
        #X = Outer water - OO Midpoint to outer water. Basically, along O3O1 'bond'
        #Y =
        #Z =
        #Lindsey's index: [walker#,atom#,(x,y,z)]
        if atmnm == 9:
            outerW = 1
        if atmnm == 10:
            outerW = 2
        #print 'Shape Coordinates: ',np.shape(xx)

        oW = xx[:,outerW-1,:]  # Coordinates of outer water Oxygen
        center = xx[:,3-1,:]
        mp = (center + oW )/ 2
        sharedH = xx[:,atmnm - 1,:]  # Coordinates of shared Hydrogen
        #print 'mp=',mp
        a = la.norm(oW-mp,axis=1)
        xaxis = np.divide((oW - mp), (la.norm((oW - mp),axis=1).reshape(-1,1)))  # Normalized coordinates of xaxis definition. aka vector with only x component, where x = 1
        #print "xaxis", xaxis
        # print np.shape(xaxis)
        # print xaxis
        yaxis = np.zeros((np.shape(xx)[0], 3))
        xcomp = np.zeros(np.shape(xx)[0])
        ycomp = np.zeros(np.shape(xx)[0])
        # print sharedH[0]
        # print xaxis[0]
        # print mp[0]
        mpsharedHxaxis=np.column_stack((mp,sharedH,xaxis))
        pool = Pool(12)
        res=pool.map(yaxcyc,(mpsharedHxaxis),chunksize=1000)
        pool.close()
        pool.join()
        nres=np.array(res)
        xcomp = nres[:,1]
        ycomp = nres[:,2]
        #for i in range(len(mp)):
        #    yaxis[i] = np.cross(xaxis[i], (0, 0, 1))
        #    xcomp[i] = np.dot(sharedH[i] - mp[i], xaxis[i].T)
        #    ycomp[i] = np.dot(sharedH[i] - mp[i], yaxis[i].T)
        # print "yaxis\n",yaxis
        # print "xcomp\n",xcomp
        # print "ycomp\n",ycomp
        #print yaxis[0],xcomp[0],ycomp[0]
        #stop
        #print  'np.average(xcomp)',np.average(xcomp)
            #xcomp*=-1

        zComp = xx[:,atmnm-1,-1].copy()
        return xcomp, ycomp, zComp

    def getBisectingVector(self,left, middle, right):
        # print np.shape(left)
        # print np.shape(right)
        # print np.shape(middle)
        # print np.shape(la.norm(left-middle,axis=1).reshape(-1,1))
        bisector1 = la.norm(left - middle, axis=1).reshape(-1, 1) * (right - middle)  # |b|*a + |a|*b
        bisector2 = la.norm(right - middle, axis=1).reshape(-1, 1) * (left - middle)
        # print bisector1.shape()
        normedbisector = la.norm(bisector1 + bisector2, axis=1).reshape(-1, 1)
        # print normedbisector.shape()
        bisector = (bisector1 + bisector2) / normedbisector
        return bisector

    def xyzFreeHydronium(self,xx):

        xaxis = self.getBisectingVector(xx[:,1-1,:], xx[:,3 - 1,:],xx[:,2 - 1,:])
        zcomp = xx[:,8-1,-1].copy()
        h8 = xx[:,8-1,:].copy()
        #print 'Shape of H8 coordinates', np.shape(h8)
        #print 'Len of h8',len(h8)
        o3 = xx[:,3-1,:].copy()
        yaxis = np.zeros((len(h8),3))
        xcomp = np.zeros(len(h8))
        ycomp = np.zeros(len(h8))
        #print np.shape(xaxis)
        #print np.shape(h8-o3)

        h8o3xaxis = np.column_stack((o3,h8,xaxis))
        pool = Pool(12)
        res = pool.map(yaxcyc2, (h8o3xaxis),chunksize=1000) #h8 = sharedH
        pool.close()
        pool.join()

        nres = copy.deepcopy(np.array(res))
        #print type(nres)
        phiOH = nres[:,0]
        #print phiOH[0:10]
        xcomp = nres[:, 1]
        ycomp = nres[:, 2]
        #print 'testing phiOH vs phiOH'
        #print ycomp[0:10]

        #xcomp = xcomp.copy()
        #xcomp2 = np.empty_like(xcomp)
        #np.copyto(xcomp2, xcomp)
        #xcomp = xcomp2
        #ycomp = np.array(ycomp.tolist())
        #for i in range(len(zcomp)):
        #    yaxis[i] = np.cross(xaxis[i], (0, 0, 1))
        #    xcomp[i] = np.dot(h8[i] - o3[i], xaxis[i].T)
        #    ycomp[i] = np.dot(h8[i] - o3[i], yaxis[i].T)
        rOH = np.zeros((len(h8), 3))
        rOH[:, 0] = xcomp
        rOH[:, 1] = ycomp
        rOH[:, 2] = zcomp
        rdistOH = la.norm(rOH, axis=1)
        #phiOH2 = np.arctan(ycomp/xcomp)  # To have it peak at 90
        thetaOH = np.arccos(zcomp / rdistOH)
        # thetaOH = np.degrees(np.arccos(zcomp / rdistOH)) #/(la.norm(zcomp)*la.norm(rdistOH)))
        # phiOH = np.degrees(np.arctan(ycomp/xcomp))  # To have it peak at 90

        return rdistOH, thetaOH, phiOH

    def getOOOEulerAxes(self,Ox, centr):
        xaxis = np.divide((Ox - centr), la.norm(Ox - centr, axis=1).reshape(-1,1))  # Normalized coordinates of xaxis definition. aka vector with only x component, where x = 1 #OO vector
        #print 'xaxis for OOO: ', xaxis[0]
        yaxis = np.zeros((len(Ox), 3))
        zaxis = np.array([0, 0, 1])  # in our definition
        zax = np.tile(zaxis,(len(xaxis),1))
        #print 'xax',xaxis
        XZ = np.column_stack((xaxis,zax))
        #np.cross(xaxis,zax)
        #print 'XZ!', XZ
        #print 'xz[0]',XZ[0]
        #print 'xz[1]',XZ[1]
        pool = Pool(12)
        res = pool.map(simpleCross, XZ,chunksize=1000)  #
        pool.close()
        pool.join()
        yaxis = np.array(res)
        #print yaxis[1]
        #for i in range(len(Ox)):
        #    yaxis[i] = np.cross(xaxis[i], (0, 0, 1))
        #print yaxis[1]
        #stop
        #zaxis = np.array([0, 0, 1])  # in our definition
        return xaxis, yaxis, zaxis

    def getHOHEulerAxes(self,o, h1, h2):
        xax = self.getBisectingVector(h1, o, h2)  # bisector of HOH angle
        yax = np.zeros((len(h1), 3))
        zax = np.zeros((len(h1), 3))
        #print (h1-o).shape
        #print la.norm(h1-o,axis=1).shape
        h1on = (h1-o)/ la.norm(h1-o,axis=1).reshape([-1,1])
        #print h1on
        h2on = (h2-o)/ la.norm(h2-o,axis=1).reshape([-1,1])
        h1onh2on = np.column_stack((h1on,h2on))
        #print h1onh2on
        pool = Pool(12)
        res = pool.map(simpleCross,h1onh2on,chunksize=1000)
        pool.close()
        pool.join()
        zax = np.array(res)
        #print zax
        zx = np.column_stack((xax,zax))
        pool = Pool(12)
        res2 = pool.map(simpleCross, zx,chunksize=1000)
        pool.close()
        pool.join()
        yax = np.array(res2)
        #print 'compare'
        #print yax
        #print yax,zax
        """WARNING - THIS GIVES A DIFFERENT ANSWER THAN MULTIPROCESSING"""
        #for i in range(len(h1)):
        #    zax[i] = np.cross(h1[i] - o[i] / (la.norm(h1[i] - o[i])), h2[i] - o[i] / (
        #        la.norm(h2[i] - o[i])))  # o and h1 for testing, first instance. normalize this vector, since it's not normed like normal in this treatment.
        #    yax[i] = np.cross(xax[i], zax[i])
        #print yax
        """OLD DEFINITION - END"""

        return xax, yax, zax

    def eulerH2O__V2(self,xx,h1, O1, h2):
        O = xx[:,O1 - 1,:].copy()
        center = xx[:,3 - 1,:].copy()
        hL = xx[:,h1 - 1,:].copy()
        hR = xx[:,h2 - 1,:].copy()
        X, Y, Z = self.getOOOEulerAxes(O, center)  # gets X,Y,Z axes based on water
        x, y, z = self.getHOHEulerAxes(O, hL, hR)  # THESE ARE ALL WRONG
        # Phi = angle btw x axes
        # Theta = angle btw Z axes
        # Chi = angle btw line of node and y axis
        # or (the same definitions)
        # Phi = rotation about OOO's Z axis
        # Theta = rotation about line of nodes
        # Chi = roataion about H2O's z axis
        # [x] = [c(ph)c(th)c(ch)-s(ph)s(ch)           . . .                 . . .] [X]
        # [y] = [        . . .             -s(ph)c(th)s(ch)+c(ph)c(ch)      . . .] [Y]
        # [z] = [        . . .                        . . .                 c(th)] [Z]
        # NEW takeaway:
        # alpha = sin(th)sin(ch)
        # beta = -sin(th)cos(ch)
        # tan x = -alpha/beta

        #print np.shape(z)
        #print np.shape(Z)
        # N = np.zeros((numWalkers,3)) #Line of node axis
        tanPhi = np.zeros(len(O))
        Theta = np.zeros(len(O))
        tanChi = np.zeros(len(O))

        #print z
        zyxYX = np.column_stack((z,y,x,Y,X))
        #print zyxYX
        #print zyxYX.shape
        pool = Pool(12)
        res = pool.map(euler, zyxYX,chunksize=1000)
        pool.close()
        pool.join()

        nres = np.array(res)
        Theta = nres[:, 0]
        tanPhi = nres[:,1]
        tanChi = nres[:,2]

        #testing
        #for i in range(len(O)):
        #    Theta[i] = np.arccos(np.dot(z[i], Z.T)/(la.norm(z[i])*la.norm(Z.T)))
        #    if np.isnan(Theta[i]):
        #        stop
        #    tanPhi[i] = np.arctan2(np.dot(y[i], Z.T) , np.dot(x[i], Z.T))
        #    tanChi[i] = np.arctan2(-1 * np.dot(Y[i], z[i].T) , np.dot(X[i], z[i].T))
        #    if tanChi[i] < 0:
        #        tanChi[i]+=2*np.pi
        #    if tanPhi[i] < 0:
        #        tanPhi[i] += 2*np.pi
        return Theta,tanPhi, tanChi

        #return np.degrees(Theta), np.degrees(tanPhi), np.degrees(tanChi)
#!!!!!!!!!!!!!!!!!!!!!!!!!h9o4
    def H9xyzSharedHydrogens(self,atmnm,xx):
        # X = Outer water - OO Midpoint to outer water. Basically, along O3O1 'bond'
        # Y =
        # Z =
        # Lindsey's index: [walker#,atom#,(x,y,z)]
        if atmnm == 11:
            outerW = 2
        elif atmnm == 12:
            outerW = 3
        elif atmnm == 13:
            outerW = 1

        # print 'Shape Coordinates: ',np.shape(xx)

        oW = xx[:, outerW - 1, :]  # Coordinates of outer water Oxygen
        center = xx[:, 4 - 1, :]
        mp = (center + oW) / 2

        sharedH = xx[:, atmnm - 1, :]  # Coordinates of shared Hydrogen
        #xaxis = np.divide((oW - mp), la.norm(oW - mp, axis=1).reshape(-1,1))  # Normalized coordinates of xaxis definition. aka vector with only x component, where x = 1
        xaxis = np.divide((mp - oW), la.norm(oW - mp, axis=1).reshape(-1,1))  # Normalized coordinates of xaxis definition. aka vector with only x component, where x = 1

        dummy = center.copy()
        dummy[:, -1] = 0.0
        oaHat=np.copy(xaxis)
        OB=dummy-oW
        s=(oaHat*OB).sum(axis=1)
        OC=oaHat*s[:,np.newaxis]
        ze=OC-OB
        #at this point, my x axis points the 'wrong' direction.  I will flip the sign
        xaxis=np.negative(xaxis)
        zaxis = ze / la.norm(ze, axis=1)[:, None]
        zdotx= (zaxis*xaxis).sum(axis=1)
        yaxis=np.cross(zaxis,xaxis,axis=1)
        # # zaxis
        # dummy = center.copy()
        # dummy[:, -1] = 0.0
        # OA = center - oW
        # BA = center - dummy
        # st = time.time()
        # un = OA / la.norm(OA, axis=1)[:, None]  # unit vector pointing along OO axis
        #
        # #print 'UN',un
        # # print la.norm(un,axis=1)
        # OB = dummy - oW
        # OC = np.zeros((len(OB), 3))
        # for i in range(len(oW)):
        #     s = np.dot(un[i], OB[i])
        #     OC[i] = s * un[i]
        # AC = OA - OC
        # ze = AC - BA
        # zaxis = ze / la.norm(ze, axis=1)[:, None]
        # yaxis = np.cross(zaxis,xaxis,axis=1)
        """ I think I made a mistake earlier on this... Paralellization was not account for xax"""
        xcomp = ((sharedH-mp)*xaxis).sum(axis=1)
        ycomp = ((sharedH-mp)*yaxis).sum(axis=1)
        zcomp = ((sharedH-mp)*zaxis).sum(axis=1)
        return xcomp, ycomp, zcomp

    def H9getOOOAxis_v1(self,xx,oW):
        # oW = xx[:,outerW-1,:]  # Coordinates of outer water Oxygen
        center = xx[:, 4 - 1, :]
        xaxis = np.divide((center -oW), la.norm(center- oW, axis=1).reshape(-1,1))  # Normalized coordinates of xaxis definition. aka vector with only x component, where x = 1
        # zaxis
        dummy = center.copy()
        dummy[:, -1] = 0.0
        oaHat=np.copy(xaxis)
        OB=dummy-oW
        s=(oaHat*OB).sum(axis=1)
        OC=oaHat*s[:,np.newaxis]
        ze=OC-OB
        #at this point, my x axis points the 'wrong' direction.  I will flip the sign
        xaxis=np.negative(xaxis)
        zaxis = ze / la.norm(ze, axis=1)[:, None]
        zdotx= (zaxis*xaxis).sum(axis=1)
        yaxis=np.cross(zaxis,xaxis,axis=1)
        return xaxis, yaxis, zaxis

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
        yaxis = np.zeros((len(h1), 3))
        zaxis = np.zeros((len(h1), 3))
        h1on = (h1 - o) / la.norm(h1 - o, axis=1).reshape([-1, 1])
        h2on = (h2 - o) / la.norm(h2 - o, axis=1).reshape([-1, 1])
        zaxis=np.cross(h1on,h2on,axis=1)
        yaxis = np.cross(zaxis,xaxis,axis=1)
        return xaxis, yaxis, zaxis

    def H9eulerH2O__V2(self,xx,h1,O1,h2):
        O = xx[:, O1 - 1, :].copy()
        center = xx[:, 4 - 1, :].copy()
        hL = xx[:, h1 - 1, :].copy()
        hR = xx[:, h2 - 1, :].copy()
        X, Y, Z = self.H9getOOOAxis(xx, O)  # gets X,Y,Z axes based on water
        x, y, z = self.H9GetHOHAxis(O, hL, hR)  # THESE ARE ALL WRONG
        # testing
        print z.shape
        Theta = np.arccos((z*Z).sum(axis=1)/(la.norm(z,axis=1) * la.norm(Z,axis=1)))
        #Theta = np.arccos((z*Z).sum(axis=1))
        #tanPhi = np.arctan2((y * Z).sum(axis=1)/(la.norm(y,axis=1) * la.norm(Z,axis=1)),
        #                    (x*Z).sum(axis=1)/(la.norm(x,axis=1) * la.norm(Z,axis=1)))
        #tanChi = np.arctan2(-1 * (Y * z).sum(axis=1) / (la.norm(Y, axis=1) * la.norm(z, axis=1)),
        #                    (X * z).sum(axis=1) / (la.norm(X, axis=1) * la.norm(z, axis=1)))
        print 'hello'
        tanPhi = np.arctan2((Y * z).sum(axis=1)/(la.norm(Y,axis=1) * la.norm(z,axis=1)),
                            (X*z).sum(axis=1)/(la.norm(X,axis=1) * la.norm(z,axis=1)))
        tanChi = np.arctan2(-(y*Z).sum(axis=1) / (la.norm(y,axis=1) * la.norm(Z,axis=1)),
                            (x*Z).sum(axis=1) / (la.norm(x,axis=1) * la.norm(Z,axis=1)))

        tanChi[tanChi < 0]+=(2*np.pi)
        tanPhi[tanPhi < 0]+=(2*np.pi)
        return Theta, tanPhi, tanChi

    def planarXYZ(self,xx):
        """Deviation from COM"""
        #com,allEckVecs,killList
        start = time.time()
        print 'eckarting justOOO...'
        oCom,oEckVectors,killList= self.eckartRotate(xx,justO=True)
        print 'done'
        print 'that took ', time.time()-start,' seconds'

        #xxNew=xx-oCom[:,np.newaxis,:]
        xxNew = np.zeros(xx.shape)
        for walker in range(len(xx)):
            xxNew[walker] = np.dot(xx[walker], oEckVectors[walker])

        """tetRef = self.pullTetramerRefPos()
        #rotTetRef,rotTetRef2 = self.rotateBackToFrame(np.array([tetRef,tetRef]),2,1,3) #the reference geometry
        rotTetRef = self.pullTetramerRefPos()# the reference geometry

        file = open("OxygensOf3","w+")
        walkern=2
        for i in range(3):
            file.write("%5.12f %5.12f %5.12f\n" % (rotTetRef[i,0],rotTetRef[i,1],rotTetRef[i,2]))
        for j in range(3):
            file.write("%5.12f %5.12f %5.12f\n" % (xx[walkern,j,0],xx[walkern,j,1],xx[walkern,j,2]))
        for k in range(3):
            file.write("%5.12f %5.12f %5.12f\n" % (xxNew[walkern,k,0],xxNew[walkern,k,1],xxNew[walkern,k,2]))
        file.close()
        stop"""
        o4 = xx[:,4-1,:] #central atom
        COM = (xx[:,1-1,:]+xx[:,2-1,:]+xx[:,3-1,:])/3

        xaxis = (xxNew[:,2-1,:]-COM)/la.norm(COM-xxNew[:,2-1,:],axis=1)[:,np.newaxis] #define x axis relative to eckarted walker.
        zaxis = np.tile([0,0,1],(len(xaxis),1))
        zaxis2 = o4[:,-1]/np.abs(o4[:,-1])
        yaxis = np.cross(zaxis, xaxis,axis=1)  # z cross x
        disp=o4-COM
        xcomp = (disp*xaxis).sum(axis=1)
        ycomp = (disp*yaxis).sum(axis=1)
        zcomp = o4[:,-1]
        #return disp[:,0],disp[:,1],disp[:,2]
        #print xcomp[0],ycomp[0],zcomp[0]
        #print xcomp[6],ycomp[6],zcomp[6]
        #stop
        return xcomp,ycomp,zcomp


    def oldPlanarXYZ(self,xx):
        O2 = xx[:,2-1,:] #(0,0,0)
        O1 = xx[:,1-1,:] #(x,0,0)
        O4 = xx[:,4-1,:]
        #print 'norming xaxis...'
        xaxis = (O1-O4)/la.norm(O1-O4,axis=1)[:,None] #should just be 1,0,0
        #print 'tiling zaxis...'
        zaxis = np.tile([0,0,1],(len(xaxis),1))
        yaxis = np.zeros((len(zaxis),3))
        xcomp = np.zeros(len(zaxis))
        ycomp = np.zeros(len(zaxis))
        zcomp = np.zeros(len(zaxis))
        #print 'prepping and running xcomp and ycomp'
        yaxis = np.cross(zaxis, xaxis,axis=1)  # z cross x
        xcomp = (O4*xaxis).sum(axis=1)
        ycomp = (O4*yaxis).sum(axis=1)
        zcomp = xx[:, 4 - 1, -1]
        #print 'done with planarXYZ'
        #for at in range(len(xaxis)):
        #    yaxis[at] = np.cross(zaxis[at],xaxis[at])
        #    xcomp[at] = np.dot(O4[at],xaxis[at].T)
        #    ycomp[at] = np.dot(O4[at],yaxis[at].T)
        return xcomp,ycomp,zcomp

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
        maH1H2 = (aH1 + aH2) / 2.0
        maH2H3 = (aH2 + aH3) / 2.0
        maH1H3 = (aH1 + aH3) / 2.0
        # vectors between the points along the OH bonds that are 1 unit vector away from the O
        vaH1H2 = aH2 - aH1
        vaH2H3 = aH3 - aH2
        vaH1H3 = aH3 - aH1

        # calculate vector
        line = np.zeros((xx.shape[0], 3))
        for i in range(xx.shape[0]):
            line[i] = np.cross(vaH1H2[i], vaH2H3[i].T)
        # add normalized vector to O

        #print 'max, min, ave of mag of line', np.average(la.norm(line, axis=1)), np.max(la.norm(line, axis=1)), np.min(
        #    la.norm(line, axis=1))
        g = (line / la.norm(line, axis=1)[:, np.newaxis])
        D = O + (line / la.norm(line, axis=1)[:, np.newaxis])
        return D

    def umbrella(self,xx,O,H1,H2,H3):
        xx*=bohr2ang
        """O,H1,H2,H3 are indices. """
        # calculate d, the trisector point of the umbrella
        D = self.calcD(xx,O, H1, H2, H3)

        # print LindseyCoords.shape, D.shape
        addedX = np.concatenate((xx, D[:, np.newaxis, :]), axis=1)  # Change D index
        # print 'addedLindseyCoords.shapes', addedLindseyCoords.shape
        umbrell = self.ba(addedX, H2, O, 13)  # 4 11 0 now O H1 0

        return umbrell


    def finalPlanarXyz(self,xx):
        #For any geometry, use reference to determine xyz compz
        xx=xx[:,:4]
        #mass=self.get_mass()[:3]
        #ocom = (mass[:3]*xxp).sum(axis=1)/ np.sum(mass[:3])
        #xxp-= ocom[:,np.newaxis]
        #xx-=ocom[:,np.newaxis]

        ocom, eVecs, kill=self.eckartRotate(xx[:,:3],True)
        for q in eVecs:
            if q[-1,-1]==0.:
                q[-1,-1]=1.
        xx-=ocom[:,np.newaxis,:]
        xxNew = np.zeros(xx.shape)
        asdf = np.copy(xxNew)
        for walker in range(len(xx)):
            xxNew[walker] = np.dot(xx[walker], eVecs[walker])
            #asdf[walker] =  np.dot(eVecs[walker], xx[walker])

        """rotTetRef = self.pullTetramerRefPos()# the reference geometry
        file = open("OxygensOf4","w+")
        walkern=1
        for i in range(4):
            file.write("%5.12f %5.12f %5.12f\n" % (rotTetRef[i,0],rotTetRef[i,1],rotTetRef[i,2]))
        for j in range(4):
            file.write("%5.12f %5.12f %5.12f\n" % (xx[walkern,j,0],xx[walkern,j,1],xx[walkern,j,2]))
        for k in range(4):
            file.write("%5.12f %5.12f %5.12f\n" % (xxNew[walkern,k,0],xxNew[walkern,k,1],xxNew[walkern,k,2]))
        file.close()"""
        print xxNew[0]
        print xxNew[1]
        print 'idk'
        return xxNew[:,4-1,0],xxNew[:,4-1,1],xxNew[:,4-1,2],

    def getfinalOOAxes(self,atmnm,xx):
        print 'fssdf'
        if atmnm == 11:
            outerW = 2
        elif atmnm == 12:
            outerW = 3
        elif atmnm == 13:
            outerW = 1
        oW = xx[:, outerW - 1, :]  # Coordinates of outer water Oxygen
        center = xx[:, 4 - 1, :]
        mp = (center + oW) / 2 #no decimal
        print mp
        sharedH = xx[:, atmnm - 1, :]  # Coordinates of shared Hydrogen
        xaxis = np.divide((mp - oW), la.norm(oW - mp, axis=1).reshape(-1,1))  # Normalized coordinates of xaxis definition. aka vector with only x component, where x = 1
        dummy = center.copy()
        dummy[:, -1] = 0.0
        oaHat=np.copy(xaxis)
        OB=dummy-oW
        s=(oaHat*OB).sum(axis=1)
        OC=oaHat*s[:,np.newaxis]
        ze=OC-OB
        #at this point, my x axis points the 'wrong' direction.  I will flip the sign
        xaxis=np.negative(xaxis)
        zaxis = ze / la.norm(ze, axis=1)[:, None]
        #I don't think this is correct, I think I should be taking X x Y to get Z
        #On second thought, I think this is okay
        negZ=np.where(xx[:,4-1,-1]<0)
        zaxis[negZ,-1]=np.negative(zaxis[negZ,-1])
        yaxis = np.cross(zaxis, xaxis, axis=1)
        return xaxis,yaxis,zaxis


    def eulerMatrix(self,x,y,z,X,Y,Z):
        Theta = np.arccos((z*Z).sum(axis=1)/(la.norm(z,axis=1) * la.norm(Z,axis=1)))
        tanPhi = np.arctan2((Y * z).sum(axis=1)/(la.norm(Y,axis=1) * la.norm(z,axis=1)),
                            (X*z).sum(axis=1)/(la.norm(X,axis=1) * la.norm(z,axis=1)))
        tanChi = np.arctan2((y*Z).sum(axis=1) / (la.norm(y,axis=1) * la.norm(Z,axis=1)),
                            -(x*Z).sum(axis=1) / (la.norm(x,axis=1) * la.norm(Z,axis=1)))

        tanChi[tanChi < 0]+=(2*np.pi)
        tanPhi[tanPhi < 0]+=(2*np.pi)
        return Theta, tanPhi, tanChi

    def finalPlaneShareEuler(self,xx):
        #For any geometry, use reference to determine xyz compz
        xxp=np.copy(xx[:,:4])
        ocom, eVecs, kill=self.eckartRotate(xxp[:,:3],True)
        xx-=ocom[:,np.newaxis,:]
        #xxNew = np.zeros(xx.shape)
        #asdf = np.copy(xxNew)
        yy=np.copy(xx)
        for walker in range(len(yy)):
            for atm in range(len(yy[walker])):
                yy[walker,atm]=np.dot(yy[walker,atm],eVecs[walker])


        for walker in range(len(xx)):
            xx[walker] = np.dot(xx[walker], eVecs[walker])
        planarXYZ=np.copy(np.array([xx[:,4-1,0],xx[:,4-1,1],xx[:,4-1,2]]))
        ##########Got O4 - Begin SharedProtons#############
        #zOfEck4=planarXYZ[-1]
        #Do normal stuff, but flip z based on sign of eckarted shit
        ######615930827
        print 'central O computed.  Moving to eulers and shared XYZs'
        atmnm=11
        h1 = 8
        h2 = 7
        o  = 2
        mp = 0.5*(xx[:,o-1]+xx[:,4-1])
        X,Y,Z = self.getfinalOOAxes(atmnm,xx)
        xcomp11 = ((xx[:,atmnm-1]-mp)*X).sum(axis=1)
        ycomp11 = ((xx[:,atmnm-1]-mp)*Y).sum(axis=1)
        zcomp11 = ((xx[:,atmnm-1]-mp)*Z).sum(axis=1)
        x,y,z = self.H9GetHOHAxis(xx[:,o-1],xx[:,h1-1],xx[:,h2-1])
        th11,phi11,xi11 = self.eulerMatrix(x,y,z,X,Y,Z)

        atmnm=12
        h1 = 9
        h2 = 10
        o = 3
        mp = 0.5 * (xx[:, o - 1] + xx[:, 4 - 1])
        X,Y,Z = self.getfinalOOAxes(atmnm,xx)
        xcomp12 = ((xx[:,atmnm-1]-mp)*X).sum(axis=1)
        ycomp12 = ((xx[:,atmnm-1]-mp)*Y).sum(axis=1)
        zcomp12 = ((xx[:,atmnm-1]-mp)*Z).sum(axis=1)
        x,y,z = self.H9GetHOHAxis(xx[:,o-1],xx[:,h1-1],xx[:,h2-1])
        th12,phi12,xi12 = self.eulerMatrix(x,y,z,X,Y,Z)

        atmnm=13
        h1 = 6
        h2 = 5
        o = 1
        mp = 0.5 * (xx[:, o - 1] + xx[:, 4 - 1])
        X,Y,Z = self.getfinalOOAxes(atmnm,xx)
        xcomp13 = ((xx[:,atmnm-1]-mp)*X).sum(axis=1)
        ycomp13 = ((xx[:,atmnm-1]-mp)*Y).sum(axis=1)
        zcomp13 = ((xx[:,atmnm-1]-mp)*Z).sum(axis=1)
        x,y,z = self.H9GetHOHAxis(xx[:,o-1],xx[:,h1-1],xx[:,h2-1])
        th13,phi13,xi13 = self.eulerMatrix(x,y,z,X,Y,Z)


        return planarXYZ[0],planarXYZ[1],planarXYZ[2],xcomp11,ycomp11,zcomp11,xcomp12,ycomp12,zcomp12,xcomp13,ycomp13,zcomp13,th11,phi11,xi11,th12,phi12,xi12,th13,phi13,xi13


    def SymInternalsH9O4plus(self,x):
        print 'Commence getting internal coordinates for tetramer'
        start = time.time()
        all = self.finalPlaneShareEuler(x)
        xyz11 = all[3:6]
        xyz12 = all[6:9]
        xyz13 = all[9:12]
        thphixi1=all[12:15]
        thphixi2=all[15:18]
        thphixi3=all[18:21]
        xyzO4 = all[0:3]


        #xyz11 = self.H9xyzSharedHydrogens(11,x)
        #xyz12 = self.H9xyzSharedHydrogens(12, x)
        #xyz13 = self.H9xyzSharedHydrogens(13, x)
        print 'done hydronium XYZ'
        print 'time it took to get xyzs: ',str(time.time()-start)
        second = time.time()
        #thphixi1 = self.H9eulerH2O__V2(x,6,5,1)
        #thphixi1 = self.H9eulerH2O__V2(x,6,1,5)
        #thphixi2 = self.H9eulerH2O__V2(x, 9, 3,10)
        #thphixi3 = self.H9eulerH2O__V2(x, 8, 2, 7)
        print 'done euler angles'
        print 'time it took to eulers: ', str(time.time()-second)
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
        fourth = time.time()
        #xyzO4 = self.finalPlanarXyz(x)
        print 'time for planar stuff', str(time.time() - fourth)
        print 'Done with all internals'

        internal = np.array(
            zip(xyz11[0], xyz11[1], xyz11[2], xyz12[0], xyz12[1], xyz12[2], xyz13[0], xyz13[1], xyz13[2],
                thphixi1[0], thphixi1[1], thphixi1[2], thphixi2[0], thphixi2[1], thphixi2[2], thphixi3[0], thphixi3[1],
                thphixi3[2],rOH5, rOH6, HOH516, rOH7, rOH8, HOH728, rOH9, rOH10, HOH9310, rO1O2,rO1O3,rO2O3,xyzO4[0],xyzO4[1],xyzO4[2]))

        print 'internal shape: ',np.shape(internal)
        print 'internal[0] shape: ',np.shape(internal[0])
        #RYAN COMMENTED THIS OUT self.internalConversion=[bohr2ang,bohr2ang,bohr2ang,rad2deg,rad2deg,rad2deg,bohr2ang,bohr2ang,bohr2ang]
        """self.internalName = ['xH11','yH11','zH11','xH12','yH12','zH12','xH13','yH13','zH13','theta651','phi651','Chi651',
                             'theta1039', 'phi1039', 'Chi1039','theta728','phi728','Chi728','rOH5','rOH6','HOH516','rOH7','rOH8','HOH728',
                             'rOH9', 'rOH10', 'HOH9310','rO4O1','rO4O3','aO1O2O3','xO4','yO4','zO4']"""

        self.internalName = ['xH11', 'yH11', 'zH11', 'xH12', 'yH12', 'zH12', 'xH13', 'yH13', 'zH13', 'theta651',
                             'phi651', 'Chi651',
                             'theta1039', 'phi1039', 'Chi1039', 'theta728', 'phi728', 'Chi728', 'rOH5', 'rOH6',
                             'HOH516', 'rOH7', 'rOH8', 'HOH728',
                             'rOH9', 'rOH10', 'HOH9310', 'rO1O2','rO1O3','rO2O3', 'xO4', 'yO4', 'zO4']
                            #'rOH9', 'rOH10', 'HOH9310', 'rO4O1', 'rO4O2', 'rO4O3', 'Oumbrella', 'thetaOx', 'thetaOx2']

        return internal

    def SymInternalsH7O3plus(self,x):
        #3N-6=24: cartesian shared proton (6) ; H8 r,th,phi (3) ; Euler angles for flanking waters: th,ph,xi (6) ; OOdists+OOOtheta (3) ; OH dists + HOH theta flanking (6)
        #print 'Changing to angstrom coordinates'
        #rOH1 = self.bondlength(x, atom1=1 - 1, atom2=4 - 1)
        #print 'before conversion', rOH1[0]
        #x *= bohr2ang

        #rOH1 = self.bondlength(x, atom1=1 - 1, atom2=4 - 1)
        #print 'after conversion', rOH1[0]
        #xyzO9=parmap(self.xyzSharedHydrogens,(self,9,x))-
        
        print 'Commence getting internal coordinates for trimer'
        xyzO9 = self.xyzSharedHydrogens(9, x)  # *ang2bohr  #Cartesian shared hydrogen 3xnwalkers
        print 'done with xyzO9'
        xyz10 = self.xyzSharedHydrogens(10, x)  # *ang2bohr #Cartesian shared hydrogen 3xnwalkers
        print 'done with xyz10'
        rthphi = self.xyzFreeHydronium(x)  # *ang2bohr     #Spherical free hydrogen 3xnwalkers in degrees
        print 'done with FH'
        thphixi1 = self.eulerH2O__V2(x, 6, 2,7)  # Euler angles - 3xnwalkers - going to need to do some adjustments with +/-360
        print 'done with Euler1'
        thphixi2 = self.eulerH2O__V2(x, 5, 1,4)  # Euler angles - 3xnwalkers - going to need to do some adjustments with +/-360
        print 'done with Euler2'
        rOH1 = self.bondlength(x, atom1=1 - 1, atom2=4 - 1)  # *ang2bohr
        rOH2 = self.bondlength(x, atom1=1 - 1, atom2=5 - 1)  # *ang2bohr
        aHOH1 = self.bondAngle(x, atom1=5 - 1, atom2=1 - 1, atom3=4 - 1)  # RIGHT SIDE OF TRIMER r1,r2,theta
        rOH3 = self.bondlength(x, atom1=2 - 1, atom2=6 - 1)  ##*ang2bohr
        rOH4 = self.bondlength(x, atom1=2 - 1, atom2=7 - 1)  # *ang2bohr
        aHOH2 = self.bondAngle(x, atom1=6 - 1, atom2=2 - 1, atom3=7 - 1)  # LEFT SIDE OF TRIMER r1,r2,theta
        rOO1 = self.bondlength(x, atom1=1 - 1, atom2=3 - 1)  # *ang2bohr
        rOO2 = self.bondlength(x, atom1=2 - 1, atom2=3 - 1)  # *ang2bohr
        aOOO = self.bondAngle(x, atom1=1 - 1, atom2=3 - 1, atom3=2 - 1)  # OOO ANGLE r1,r2,theta

        #Symmetric & asymmetric stretch
        #h9/h10 sym-asymm
        symmHB = (1/np.sqrt(2)) * (xyzO9[0]+xyz10[0])
        asymmHB = (1/np.sqrt(2)) * (xyzO9[0]-xyz10[0])
        print 'first symm', symmHB[0]
        print 'first asymm', asymmHB[0]
        print 'first O1H4 bond length: ', rOH1[0]
        print 'first Angle: ',aOOO[0]
        #self.PltHists1D('allH',aOOO.flatten(),(70,180),'OOOAngle','allHTesting/Probability Density',False,weits)
        #self.PltHists1D('allH', rthphi[1].flatten(), (0, 180), 'Theta (H8)', 'allHTesting/Probability Density', False, weits)
        #self.PltHists1D('allH', rthphi[2].flatten(), (-180, 180), 'Phi (H8)', 'allHTesting/Probability Density', False,weits)

        #SYMMETRIZE?
        #print "Shouldn't these be symmetric and antisymmetric stretches? Bends and rocks? Not just cartesian coordinates."

        #print 'xyzO9: ',np.shape(xyzO9)
        #print 'rthphi: ',np.shape(rthphi[0])
        #internal = np.array(zip(xyzO9[0],xyzO9[1],xyzO9[2],xyz10[0],xyz10[1],xyz10[2],rthphi[0],rthphi[1],rthphi[2],
        #                        thphixi1[0],thphixi1[1],thphixi1[2],thphixi2[0],thphixi2[1],thphixi2[2]
        #                        ,rOH1,rOH2,aHOH1,rOH3,rOH4,aHOH2,rOO1,rOO2,aOOO))
        internal = np.array(
            zip(symmHB, xyzO9[1], xyzO9[2], asymmHB, xyz10[1], xyz10[2], rthphi[0], rthphi[1], rthphi[2],
                thphixi1[0], thphixi1[1], thphixi1[2], thphixi2[0], thphixi2[1], thphixi2[2]
                , rOH1, rOH2, aHOH1, rOH3, rOH4, aHOH2, rOO1, rOO2, aOOO))
        #print 'internal shape: ',np.shape(internal)
        #print 'internal[0] shape: ',np.shape(internal[0])
        #RYAN COMMENTED THIS OUT self.internalConversion=[bohr2ang,bohr2ang,bohr2ang,rad2deg,rad2deg,rad2deg,bohr2ang,bohr2ang,bohr2ang]
        self.internalName = ['symmStretchHB', 'yH9', 'zH9', 'asymmStretchHB', 'yh10', 'zh10', 'rH8', 'thH8', 'phiH8','th_627','phi_627',
                             'xi_627','th_514','phi_514','xi_514','rOH_41','rOH_51','aHOH_451','rOH_26','rOH_27','aHOH_267',
                             'rOO_1','rOO_2','aOOO']
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

    def pullTrimerRefPos(self): #Eckart reference for the trimer is in an xyz file. Need just a 3xNatom array of reference structures. I can hard code this in
        myRef = np.array([ [-2.677869210066,  3.854059894133, 0.000000000000],
                           [4.692722823509, - 0.000000000000,  0.000000000000],
                           [0.000000000000,  0.000000000000,  0.000000000000],
                           [-3.300280868789,  4.749851143725, 1.466631117729],
                           [ -3.285993479891,  4.729288353663,  -1.484977070838],
                           [5.758560142898,  -0.000000005149,  -1.484919577747],
                           [5.783382924818,  -0.000000001709,  1.466733420402],
                           [-0.847110300023,  -1.617745856796,  -0.000000000000],
                           [-1.109445938627,  1.609532820363,  0.010222069526],
                           [1.954938728445,  0.006207766673,  0.010188860037] ]) #goes O1,O2,O3,H4,..H10

        myRef2 = np.array([[4.15875766E+00, -7.07188287E-01,  9.71623738E-04 ],
                           [-4.15865808E+00, -7.07217393E-01,  9.60925166E-04],
                           [-1.44995671E-04,  1.46727440E+00, -4.05794501E-04],
                           [5.12540100E+00, -1.21166877E+00,  1.46792260E+00],
                           [5.10321557E+00, -1.20192966E+00, -1.48369235E+00],
                           [-5.10316080E+00, -1.20203703E+00, -1.48364794E+00],
                           [-5.12516161E+00, -1.21167684E+00,  1.46801170E+00],
                           [9.08440280E-04,  3.29338915E+00, -1.55804935E-03],
                           [1.72882181E+00,  5.55063992E-01,  1.03940985E-02],
                           [-1.72966108E+00,  5.55909445E-01,  1.03558993E-02]]) #Not rotated to xy plane

        myBetterRef = np.array(
                        [
                            [3.15544362E-30 , 4.06869143E+00, -7.59761292E-01],
                            [-4.98270994E-16, -4.06869143E+00 ,-7.59761292E-01],
                            [-1.97215226E-30,  0.00000000E+00,  1.58929880E+00],
                            [ 1.47532198E+00,  5.00324669E+00, -1.29932702E+00],
                            [ -1.47532198E+00,  5.00324669E+00, -1.29932702E+00],
                            [ -1.47532198E+00, -5.00324669E+00 ,-1.29932702E+00],
                            [ 1.47532198E+00, -5.00324669E+00, -1.29932702E+00],
                            [ -3.94430453E-30, -2.22044605E-16,  3.41400471E+00],
                            [  7.88860905E-31, 1.69178407E+00,  6.12546816E-01],
                            [ -2.07183794E-16, -1.69178407E+00,  6.12546816E-01]])


        return myBetterRef #myRef2

    def pullTetramerRefPos(self): #Eckart reference for the trimer is in an xyz file. Need just a 3xNatom array of reference structures. I can hard code this in
        """goes O1,O2,O3,O4,..H12"""
        myRef2 = np.array([[0.00000000E+00,  4.81355109E+00, -4.53345972E-32],
                           [-4.16865752E+00, -2.40677554E+00,  1.18329136E-30],
                           [4.16865752E+00, -2.40677554E+00, -1.38050658E-30],
                           [-0.00000000E+00,  0.00000000E+00,  0.00000000E+00],
                           [1.79529146E-16,  5.90334467E+00, -1.46596673E+00],
                           [-3.94430453E-31,  5.90334467E+00,  1.46596673E+00],
                           [-5.11244645E+00, -2.95167233E+00, -1.46596673E+00],
                           [-5.11244645E+00, -2.95167233E+00,  1.46596673E+00],
                           [5.11244645E+00, -2.95167233E+00,  1.46596673E+00],
                           [5.11244645E+00, -2.95167233E+00, -1.46596673E+00],
                           [-1.65058312E+00, -9.52964606E-01,  3.94430453E-31],
                           [1.65058312E+00, -9.52964606E-01, -4.93038066E-31],
                           [0.00000000E+00,  1.90592921E+00, -1.72916465E-32]])

        myRefCOM,extra = self.rotateBackToFrame(np.array([myRef2,myRef2]),2,1,3) #rotate reference to OOO plane
        mass = self.get_mass()
        com = np.dot(mass[:3], myRefCOM[:3]) / np.sum(mass[:3])

        myRefCOM-=com
        #rotateAboutZ
        x,y,z=myRefCOM[2-1]-com
        th=np.arctan2(y,x)
        rz=np.array([[np.cos(th),np.sin(th),0],
                     [-np.sin(th),np.cos(th),0],
                     [0,0,1]])
        b4=np.copy(myRefCOM)
        for atm in range(len(myRefCOM)):
            myRefCOM[atm]=rz.dot(myRefCOM[atm])
        return myRefCOM



    def eckartRotate(self,pos,justO=False,specialCond=False): # pos coordinates = walkerCoords numwalkersxnumAtomsx3
        """Eckart Rotate method returns the transpose of the correct matrix, meaning that when one does the dot product,
        one should """
        nMolecules=pos.shape[0]
        allEckVecs = np.zeros((nMolecules,3, 3))
        if self.name in ProtonatedWaterTrimer:
            self.refPos = self.pullTrimerRefPos()
        else:
            self.refPos = self.pullTetramerRefPos()
        if len(pos.shape)<3:
            pos=np.array([pos])
        newCoord=np.zeros(pos.shape)
        killList = []
        #Center of Mass
        mass=self.get_mass()
        #com=np.dot(mass,pos)/np.sum(mass)
        if justO: #the OOO plane
            self.refPos=self.refPos [:3]
            com = np.dot(mass[:3],pos[:,:3])/np.sum(mass[:3])
            mass = mass[:3]
            pos = pos[:,:3,:]
        else:
            com = np.dot(mass, pos) / np.sum(mass)

        #First Translate:
        ShiftedMolecules=pos-com[:,np.newaxis,:]

        #Equation 3.1 in Eckart vectors, Eckart frames, and polyatomic molecules - James D. Louck and Harold W. Galbraith
        start = time.time()
        myFF = np.zeros((len(ShiftedMolecules),3,3))
        myF = np.zeros((len(ShiftedMolecules),3,3))
        st=time.time()
        asdf = np.sum(ShiftedMolecules[:,:,:,np.newaxis]*self.refPos[np.newaxis,:,np.newaxis,:]*mass[np.newaxis,:,np.newaxis,np.newaxis],axis=1)

        myF = np.transpose(asdf,(0,2,1))
        myFF = np.matmul(myF,asdf)
        bigEvals,bigEvecs=la.eigh(myFF)
        #bigEvals=np.sort(bigEvals,axis=1)
        bvec = copy.deepcopy(bigEvecs)
        #bigEvecs=bigEvecs[:,:,(2,1,0)]
        bigEvecsT=np.transpose(bigEvecs,(0,2,1))

        msk=np.where(bigEvals[:,0]==0.0)
        if len(msk[0])==0:
            invRootDiagF2 = 1.0 / np.sqrt(bigEvals)
        else:
            invRootDiagF2 = np.zeros(bigEvals.shape)
            invRootDiagF2[:,1:]=1.0 / np.sqrt(bigEvals[:,1:])

        invRootF2=np.matmul(invRootDiagF2[:,np.newaxis,:]*-bigEvecs,-bigEvecsT,) #-bigEvecs
        print myF
        eckVecs2 = np.matmul(np.transpose(myF,(0,2,1)),invRootF2)
        eckVecs2[:,-1]=np.cross(eckVecs2[:,0],eckVecs2[:,1])
        #print 'did it work'
        plus=0
        minus=0
        mas = np.where(np.around(la.det(eckVecs2),10)==-1.0)
        print 'wlks neg for mine'
        print mas
        if len(mas[0])!=0:
            killList2=mas
            #eckVecs2[mas] = np.negative(eckVecs2[mas])
            minus = len(mas[0])

        else:
            killList2=mas[0]

        plus=len(ShiftedMolecules)-minus
        print 'Plus rotation: ',plus
        print 'Inverted Rotation: ',minus
        return com, eckVecs2, killList2

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """
        start2=time.time()
        plus = 0
        minus = 0
        for moli,molecule in enumerate(ShiftedMolecules):
            Fvec=np.zeros((3,3))
            for atom,massa,eckatom in zip(molecule,mass,self.refPos):
                Fvec=Fvec+massa*np.outer(eckatom,atom)
            #F from eqn 3.4b             - vectorsthat connect the dotz
            FF=np.dot(Fvec,Fvec.transpose())
            #Diagonalize FF, R: which is actually F in the paper              
            sortEigValsF,sortEigVecF=np.linalg.eigh(FF)
            sortEigVecFT=-sortEigVecF.transpose()
            #sortEigVecFT=sortEigVecF.transpose() #RYAN: SIGN DOESN"T CHANGE ANYTHING
            if len(np.where(sortEigValsF<=0)[0])!=0:
                #sortEigVecFT=np.abs(sortEigVecFT)                                                            
                sortEigValsF=np.abs(sortEigValsF)
                invRootDiagF=copy.deepcopy(sortEigValsF)
                for e,element in enumerate(sortEigValsF):
                    if element>0:
                        invRootDiagF[e]=1.0/np.sqrt(element)
            #Get the inverse sqrt of diagonalized(FF)                                                         
            else:
                invRootDiagF=1.0/np.sqrt(sortEigValsF)
            # F^{-1/2}
            invRootF=np.dot(invRootDiagF[np.newaxis,:]*-sortEigVecF,sortEigVecFT)

            #invRootF = np.dot(invRootDiagF[np.newaxis, :] * sortEigVecF, sortEigVecFT)
        
            #3.4C? ryan
            eckVecs=np.dot(Fvec.transpose(),invRootF)

            allEckVecs[moli] = eckVecs

            if moli==6:
                print 'asdf'
            if not justO:
                #newCoordRD[moli] = np.dot(eckVecs, molecule.T).T wrong!
                newCoord[moli] = np.dot(molecule, eckVecs)
                detEck=la.det(eckVecs)
                if np.around(detEck)==-1.:
                    killList.append(moli)
                    minus+=1
                elif np.around(detEck)==1.:
                    plus+=1

            else:
                detEck = la.det(eckVecs)
                if np.around(detEck) == -1.:
                    killList.append(moli)
                    print moli
                    minus+=1
                else:
                    plus+=1
        
            if len(np.where(np.isnan(newCoord[moli]))[0])!=0:
                print 'whaaaaaT?! nan',np.where(np.isnan(newCoord[moli])),'\ncoords:\n',newCoord[moli]
                print '   molecule number:',moli,'\n   sortEigValsF: \n', sortEigValsF,'\n   molecule: \n', molecule,
                print '\n   eckVecs \n', eckVecs
                octopus
        #print 'Lindsey: ',time.time()-start
        
        print "Whew ! Done with eckart."
        print 'plus',plus
        print 'minus',minus
        #print 'first allEckVecs', allEckVecs[0]
        #print allEckVecs
        #ff = open('allHTesting/eckartRotatedMolecule',"w+")
        #elf.printCoordsToFile(newCoord,ff)
        #ff.close()
        #print 'recorded new eckart coordinates'
        print np.array_equal(np.around(allEckVecs,5),np.around(eckVecs2,5))
        if self.name in ProtonatedWaterTrimer or self.name in ProtonatedWaterTetramer:
            return com,allEckVecs,killList
        else:
            return newCoord """


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
        au2ang=0.529177249
        fakeMoleculeNames=['O','O','O','H','H','H','H','H','H','H','H'] #needs to be updated.  Li, N, Be are for distinguishing atoms in h3o2
        #for particle in x:
        fileout.write(str(self.nAtoms)+' \n'+' \n')
        for atomName, atom in zip(fakeMoleculeNames,x):
            print atom
            fileout.write(str(atomName)+'   '+str(au2ang*atom[0])+"   "+str(au2ang*atom[1])+"   "+str(au2ang*atom[2])+"\n")
        fileout.write("\n")
        return

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

