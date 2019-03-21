import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import sys
import DMCClusters as dmc
import subprocess as sub
angstr = 0.529177
au2wn=219474.63

coordinateSet = sys.argv[1]
coordModel = sys.argv[2]
modeNum = sys.argv[3]
Wfn=dmc.wavefunction('H9O4+', 1)
if 'allD' in coordinateSet:
    Wfn.setIsotope('fullyDeuterated')
if '1He' in coordinateSet:
    Wfn.setIsotope('notDeuteratedOnce_eigen')
if '1Hw' in coordinateSet:
    Wfn.setIsotope('notDeuteratedOnce_fw')
path='../spectra/'
gPath = '../Gmats/tetramer/'
GfileName = gPath+coordinateSet+'.gmat'
numWalkers = 201
def writeNewWalkers(walkz,title):
    coordAr = ['O','O','O','O','H','H','H','H','H','H','H','H','H']
    walkz2 = walkz*angstr
    writeFile = open(title,"w")
    for walker in walkz2:
        writeFile.write("13\n \n")
        for k,atm in enumerate(walker):
            writeFile.write("%s %5.12f  %5.12f  %5.12f\n" % (coordAr[k],atm[0],atm[1],atm[2]))
        writeFile.write("\n")
    writeFile.close()

def getEqGeom():
    RefGeom =np.array([[0.00000000E+00,  4.81355109E+00, -4.53345972E-32],
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
           [0.00000000E+00,  1.90592921E+00, -1.72916465E-32]]) #must be COM of Os at origin and O2=(0,0,0)
    eqGeom = np.array([[-1.55099190 ,  1.94067311 ,  0.14704161],
             [-0.91203521,  -2.30896272,   0.14764850],
             [2.46079102,   0.36718848,   0.14815394],
             [0.00253217,   0.00164013,  -0.47227522],
             [-1.96589559,   2.46292466,  -0.54312627],
             [-2.13630186,   1.99106023,   0.90604777],
             [-1.18003749,  -2.92157144,  -0.54090532],
             [-0.65410291,  -2.84939169,   0.89772271],
             [ 2.79828182,   0.87002791,   0.89281564],
             [ 3.12620054,   0.43432898,  -0.54032031],
             [-0.31106354,  -0.91215572,  -0.20184621],
             [ 0.95094197,   0.18695800,  -0.20259538],
             [-0.63272209,   0.72926470,  -0.20069859]])
    return eqGeom


def getEulerMat(th,ph,xi):
    a=np.array([[np.cos(ph)*np.cos(th)*np.cos(xi)-np.sin(ph)*np.sin(xi),
                      np.sin(ph) * np.cos(th) * np.cos(xi) + np.cos(ph) * np.sin(xi),
                      -np.sin(th)*np.cos(xi)]
                        ,
                     [-np.cos(ph)*np.cos(th)*np.sin(xi)-np.sin(ph)*np.cos(xi),
                      -np.sin(ph) * np.cos(th) * np.sin(xi) + np.cos(ph) * np.cos(xi),
                      np.sin(th)*np.sin(xi)]
                        ,
                     [np.cos(ph)*np.sin(th),
                      np.sin(ph)*np.sin(th),
                      np.cos(th)]
                     ])
    return np.array([[np.cos(ph)*np.cos(th)*np.cos(xi)-np.sin(ph)*np.sin(xi),
                      np.sin(ph) * np.cos(th) * np.cos(xi) + np.cos(ph) * np.sin(xi),
                      -np.sin(th)*np.cos(xi)]
                        ,
                     [-np.cos(ph)*np.cos(th)*np.sin(xi)-np.sin(ph)*np.cos(xi),
                      -np.sin(ph) * np.cos(th) * np.sin(xi) + np.cos(ph) * np.cos(xi),
                      np.sin(th)*np.sin(xi)]
                        ,
                     [np.cos(ph)*np.sin(th),
                      np.sin(ph)*np.sin(th),
                      np.cos(th)]
                     ])


def RotateToOOSystem(geol,outerO,shft,geomn=np.zeros(1)):

    k=len(geomn)
    if k!=1:
        shit=True
    else:
        shit=False
    geol-=shft #o4: shift OOOCOM to origin
    shftp=[0,0,0]
    if k!=1:
        xprime=geomn-shftp
    else:
        xprime = geol[outerO-1]-geol[4-1]#my x axis defined in my general coordinates
    x,y,z=xprime[0],xprime[1],xprime[2]
    theta=np.arctan2(z,x) #o4, theta=0, don't rotate so z=0.
    r1=np.array([[np.cos(theta),0,np.sin(theta)],
                [0,1,0],
                [-np.sin(theta),0,np.cos(theta)]
                 ])
    geop = np.zeros((geol.shape))
    for atm in range(len(geol)):
        geop[atm] = r1.dot(geol[atm])
    x,y=geop[outerO-1,0],geop[outerO-1,1]
    phi = np.arctan2(-y, x)
    r2=np.array([[np.cos(phi),-np.sin([phi]),0],
                [np.sin(phi),np.cos(phi),0],
                 [0, 0, 1]
                 ])
    for atm in range(len(geol)):
        geop[atm] = r2.dot(geop[atm])
    # if outerO == 2 and k==1:
    #     geop[:,0]=-1*geop[:,0]
    #     geop[:,1] =-1 *geop[:,1]
    #
    # if outerO == 3:
    #     geop[:, 0] = -1 * geop[:, 0]
    #     geop[:, 1] = -1 * geop[:, 1]
    # if shit:
    #     geop[:,-1]=np.negative(geop[:,-1])
    #     th=np.pi
    #     r25=np.array([[1,0,0],
    #                  [0,np.cos(th),-np.sin(th)],
    #                  [0,np.sin(th),np.cos(th)]]
    #                 )
    #     for atm in range(len(geop)):
    #         geop[atm]=r25.dot(geop[atm])
    geopp=geop-geop[4-1]
    return geop

def eRotateToOOSystem(geoq,outerO,shft):
    geoq-=shft
    xprime = geoq[outerO-1]-geoq[4-1]#my x axis defined in my general coordinates
    x,y,z=xprime[0],xprime[1],xprime[2]
    theta=np.arctan2(z,x)
    #if theta > np.pi or theta < 0:
    #    theta+=2.*np.pi
    #print 'theta', np.degrees(theta)
    r1 = np.array([[np.cos(theta), 0, np.sin(theta)],
                   [0, 1, 0],
                   [-np.sin(theta), 0, np.cos(theta)]
                   ])

    geop = np.zeros((geoq.shape))
    for atm in range(len(geoq)):
        geop[atm] = r1.dot(geoq[atm])
    x,y=geop[outerO-1,0],geop[outerO-1,1]
    phi = np.arctan2(-y, x)
    #print 'phi',np.degrees(phi)
    r2=np.array([[np.cos(phi),-np.sin([phi]),0],
                [np.sin(phi),np.cos(phi),0],
                 [0, 0, 1]
                 ])
    for atm in range(len(geoq)):
        geop[atm] = r2.dot(geop[atm])

    if np.round(geop[1-1,-1]) > 0 or np.round(geop[2-1,-1]) > 0 or np.round(geop[3-1,-1]) > 0:
        th=np.pi
        r25=np.array([[1,0,0],
                     [0,np.cos(th),-np.sin(th)],
                     [0,np.sin(th),np.cos(th)]]
                    )
        for atm in range(len(geop)):
            geop[atm]=r25.dot(geop[atm])
         #geop[:,-1]=np.negative(geop[:,-1])
    if not (np.round(geop[1-1,-1])<=0. and np.round(geop[2-1,-1])<=0. and np.round(geop[3-1,-1])<=0.):
        #rotate about x axis until whichever one is
        #print 'rotatata'
        th=np.pi+np.pi/2.
        r3=np.array([[1,0,0],
                     [0,np.cos(th),-np.sin(th)],
                     [0,np.sin(th),np.cos(th)]]
                    )
        for atm in range(len(geop)):
            geop[atm]=r3.dot(geop[atm])
        #print 'rotatata'

    # if outerO == 2 and k==1:
    #     geop[:,0]=-1*geop[:,0]
    #     geop[:,1] =-1 *geop[:,1]
    #
    # if outerO == 3:
    #     geop[:, 0] = -1 * geop[:, 0]
    #     geop[:, 1] = -1 * geop[:, 1]

    geopp=geop-geop[4-1]
    return geop


def transformHOH(geom,r1,r2,ang,outerO,h1,h2,eulermat):
    HOH=np.zeros((3,3))

    #geom-=geom[outerO-1]
    ####added -1 to x comp
    ##6 X 5 , 9 X 10, 8 X 7 (Right X Left
    if h1 == 6 or h1==8:  #'Right' if looking down OO
        HOH[1,:]=np.array([r1 * np.cos(ang / 2), r1 * np.sin(ang / 2), 0.])
        HOH[2,:]=np.array([r2*np.cos(ang/2),-r2*np.sin(ang/2),0.])
    else:
        HOH[1, :] = np.array([r1 * np.cos(ang / 2), r1 * np.sin(ang / 2), 0.])
        HOH[2, :] = np.array([r2 * np.cos(ang / 2), -r2 * np.sin(ang / 2), 0.])
    # HOH[1, :] = np.array([r1 * np.cos(ang / 2.), r1 * np.sin(ang / 2.), 0.])
    # HOH[2, :] = np.array([r2 * np.cos(ang / 2.), -r2 * np.sin(ang / 2.), 0.])
    # try
    #try placing HOH with bisector on OO vector
    #HOH-=geom[4-1]
    #sanityCheckR=RotateToOOSystem(sanityCheck,outerO,sanityCheck[4-1])
    #mpo=(outerO+geom[4-1])/2
    if outerO==2:
        geom = Wfn.molecule.rotateBackToFrame(np.array([geom, geom]), 3, 2, 1)[0]
    if outerO==3:
        geom=Wfn.molecule.rotateBackToFrame(np.array([geom,geom]),1,3,2)[0]
    geomp=eRotateToOOSystem(geom,outerO,geom[outerO-1])
    #geomp = RotateToOOSystem(geom, outerO, mpo)
    geomp-=geomp[outerO-1]
    geomp[h1 - 1] = eulermat.T.dot(HOH[1])
    geomp[h2 - 1] = eulermat.T.dot(HOH[2])

    fuck=Wfn.molecule.rotateBackToFrame(np.array([geomp,geomp]),2,1,3)[0]
    return fuck

def transformSharedHs(geom,atmnm,mp,carts,outerO,geomn=np.zeros(1)):
    # if outerO==2:
    #     geom = Wfn.molecule.rotateBackToFrame(np.array([geom, geom]), 3, 2, 1)[0]
    # if outerO==3:
    #     geom=Wfn.molecule.rotateBackToFrame(np.array([geom,geom]),1,3,2)[0]
    geomp=RotateToOOSystem(geom,outerO,mp,geomn) #the oxygen that corresponds with the sharedH
    #geomp-=mp
    geomp[atmnm-1]=carts
    #geomp-=geomp[2-1]
    return Wfn.molecule.rotateBackToFrame(np.array([geomp,geomp]),2,1,3)[0]


def getOuters(age,outerMomentsR,outerMomentsE):
    """
    9 ('theta', 1.5707963267947298), 11
    10('phi', 3.917022130300791), 11
    11 ('Chi', 3.141764102329608), 11
    12 ('theta', 1.570796326795287), 12
    13 ('phi', 3.1415926535902767), 12
    14('Chi', 3.141592653591482), 12
    15 ('theta', 1.570796326795087), 13
    16('phi', 3.1415926535899237), 13
    17 ('Chi', 3.141592653590946), 13

    18('rOH5', 1.8486050034911365),
    19 ('rOH6', 1.8486050034914412),
    20('HOH516', 1.8505061097571658),
    21('rOH7', 1.848605003491622),
    22('rOH8', 1.848605003491682),
    23('HOH728', 1.8505061097565314),
    24('rOH9', 1.848605003491019),
    25 ('rOH10', 1.8486050034910018),
    26 ('HOH9310', 1.8505061097582232),"""
    # Get outers might as well do it right here
    #H6O1H5 - H13
    theta = outerMomentsE[6]
    phi = outerMomentsE[7]  # 3.14159265358988
    chi = outerMomentsE[8]
    eulermat = getEulerMat(theta, phi, chi)
    r1=outerMomentsR[0]
    r2=outerMomentsR[1]
    ang=outerMomentsR[2]
    age = rotateToOOAndPlaceOuters(age, r1, r2, ang, 1, 6, 5, eulermat)

    # H2O7H8
    theta = outerMomentsE[0]
    phi = outerMomentsE[1]
    chi = outerMomentsE[2]
    #print 'theta,phi,chi'
    #print np.degrees(theta)
    #print np.degrees(phi)
    #print np.degrees(chi)
    eulermat = getEulerMat(theta, phi, chi)
    r1=outerMomentsR[3]
    r2=outerMomentsR[4]
    ang=outerMomentsR[5]
    age = rotateToOOAndPlaceOuters(age, r1, r2, ang, 2, 8, 7, eulermat)

    theta = outerMomentsE[3]
    phi = outerMomentsE[4]
    chi = outerMomentsE[5]
    #print 'theta,phi,chi'
    #print np.degrees(theta)
    #print np.degrees(phi)
    #print np.degrees(chi)
    eulermat = getEulerMat(theta, phi, chi)
    r1=outerMomentsR[6]
    r2=outerMomentsR[7]
    ang=outerMomentsR[8]
    age = rotateToOOAndPlaceOuters(age, r1, r2, ang, 3, 9, 10, eulermat)
    #print age
    return age


# def rotate2D(x1,x2,xc,yc):
#     theta = np.arccos(np.dot(x1,x2)/(la.norm(x1)*la.norm(x2)))
#     #theta2= np.arctan2(-XYZ[1],XYZ[0])
#     xcomp = xc*np.cos(theta) - yc*np.sin(theta)
#     ycomp = yc*np.cos(theta) + xc*np.sin(theta)
#     return xcomp,ycomp

def getTriangle(ge,triMoments,sharedMoments):
    print 'begin getTriangle'
    """Gets you the coordinates of All oxygens and Hydrogens in core"""
    """
    27 ('rO1O2', 8.262888698477584), 
	28('rO1O3', 8.262888698473365), 
	29('rO2O3', 8.262888698476258), 
	30('xO4', 0.005158589118831137), 
	31('yO4', -6.291090302810263e-16),
	32 ('zO4', 0.8436517294463016)]"""
    #O2 is at origin, so that's 0,0,0
    ge[1-1,0]= triMoments[0] #rO1O2
    #Pg 33 of lab notebook.  Bx (x comp of O3) = rO2O3^2-rO1O3^2+Cx^2 / 2Cx, Cx is x comp of O1
    ge[3-1,0]= (triMoments[2]**2 - triMoments[1]**2 + ge[1-1,0]**2) / (2*ge[1-1,0])
    ge[3-1,1]=np.sqrt(triMoments[2]**2-ge[3-1,0]**2)
    Wfn.molecule.rotateBackToFrame(np.array([ge,ge]),2,1,3)
    oCom, oEckVectors, killList = Wfn.molecule.eckartRotate(np.array([ge[:4],ge[:4]]), justO=True)
    ge-=oCom[0,np.newaxis,:]
    gomo=np.dot(ge,oEckVectors[0])
    #PLACE O4!############
    gomo[4-1]=triMoments[3:6]
    ######################
    #asdf=Wfn.molecule.rotateBackToFrame(np.array([gomo,gomo]),2,1,3)[0]
    #export(Wfn.molecule.rotateBackToFrame(np.array([gomo,gomo]),2,1,3)[0],'newnewnew')
    age=np.copy(gomo)

    """0[('xH11', -0.5321874975974219), 
    1('yH11', 2.7067953923930632e-17),
    2 ('zH11', 0.05484504812449089),
    3 ('xH12', -0.5321874975979741),
    4 ('yH12', 8.743866449023069e-16),
    5 ('zH12', 0.05484504812451695), 
    6('xH13', -0.5321874975976805),
    7 ('yH13', -2.977884463189319e-16),
    8 ('zH13', 0.05484504812446696),"""
    """EQ
    -5.267726156514476177e-01
    -8.452510647260627619e-02
    -5.439273391335344709e-02
    -5.269572170506972020e-01
    -8.277303830664685391e-02
    -5.261157613133353450e-02
    -5.257875868187682489e-01
    -7.711274431814060804e-02
    -5.690159542498002265e-02"""

    age = rotateToOOAndPlaceHs(age, 11, sharedMoments[:3], 2)
    age = Wfn.molecule.rotateBackToFrame(np.array([age, age]), 2, 1, 3)[0]
    oCom, oEckVectors, killList = Wfn.molecule.eckartRotate(np.array([age[:4], age[:4]]), justO=True)
    age -= oCom[0, np.newaxis, :]
    age = np.dot(age, oEckVectors[0])

    age = rotateToOOAndPlaceHs(age,13,sharedMoments[6:],1)
    age=Wfn.molecule.rotateBackToFrame(np.array([age, age]), 2, 1, 3)[0]
    oCom, oEckVectors, killList = Wfn.molecule.eckartRotate(np.array([age[:4], age[:4]]), justO=True)
    age -= oCom[0, np.newaxis, :]
    age= np.dot(age, oEckVectors[0])

    age = rotateToOOAndPlaceHs(age, 12, sharedMoments[3:6], 3)
    age=Wfn.molecule.rotateBackToFrame(np.array([age, age]), 2, 1, 3)[0]
    oCom, oEckVectors, killList = Wfn.molecule.eckartRotate(np.array([age[:4], age[:4]]), justO=True)
    age -= oCom[0, np.newaxis, :]
    age= np.dot(age, oEckVectors[0])


    return age

def constructEuler(X,Y,Z,x,y,z):
    theta,phi,xi=Wfn.molecule.eulerMatrix(np.array([x,x]),np.array([y,y]),np.array([z,z]),
                                          np.array([X, X]),np.array([Y,Y]),np.array([Z,Z]))
    em=getEulerMat(theta[0],phi[0],xi[0])
    return em

def rotateToOOAndPlaceHs(x,sharedHn,xyzCompz,outerOn):
    orX =np.array([1.,0.,0.])
    orZ =np.cross(x[2-1],x[1-1])/la.norm(np.cross(x[2-1],x[1-1]))
    orY = np.cross(orZ,orX)/la.norm(np.cross(orZ,orX))
    xax,yax,zax=Wfn.molecule.getfinalOOAxes(sharedHn,np.array([x,x]))
    emat=constructEuler(orX,orY,orZ,xax[0]/la.norm(xax[0]),yax[0]/la.norm(yax[0]),zax[0]/la.norm(zax[0]))

    #idx=np.where(emat==-0.0)
    #emat[idx] = 0.0

    #print la.det(emat)
    #rotate
    y=np.copy(x)
    newInts = emat.T.dot(xyzCompz)
    mpy = (y[outerOn-1]+y[4-1])*0.500000000
    xyz=mpy+newInts
    for atm in range(len(x)):
        x[atm]=emat.dot(x[atm])
    asdf = (x[outerOn-1]-x[4-1])/la.norm(x[outerOn-1]-x[4-1])
    mp = (x[outerOn-1]+x[4-1])*0.500000000
    x-=mp
    x[sharedHn-1]=xyzCompz
    x+=mp
    for atm in range(len(x)):
        x[atm]=emat.T.dot(x[atm])
    xp = np.copy(x)*angstr

    y[sharedHn-1]=xyz
    y*=angstr
    return x

def rotateToOOAndPlaceOuters(x, r1, r2, ang, o, hl, hr, ema):
    if o == 1:
        sh = 13
    elif o == 2:
        sh = 11
    else:
        sh = 12

    orX =np.array([1.,0.,0.])
    orZ =np.cross(x[2-1],x[1-1])/la.norm(np.cross(x[2-1],x[1-1]))
    orY = np.cross(orZ,orX)/la.norm(np.cross(orZ,orX))
    xax,yax,zax=Wfn.molecule.getfinalOOAxes(sh,np.array([x,x]))

    ematOO=constructEuler(orX,orY,orZ,xax[0]/la.norm(xax[0]),yax[0]/la.norm(yax[0]),zax[0]/la.norm(zax[0]))
    for atm in range(len(x)):
        x[atm]=ematOO.dot(x[atm])

    asdf = (x[o-1]-x[4-1])/la.norm(x[o-1]-x[4-1])

    HOH = np.zeros((3,3))
    HOH[2,:]=np.array([r1 * np.cos(ang / 2), r1 * np.sin(ang / 2), 0.])
    HOH[1,:]=np.array([r2*np.cos(ang/2),-r2*np.sin(ang/2),0.])
    # ##########################
    # xxp = np.copy(x)
    # for i in range(len(HOH)):
    #     xxp[i]=np.dot(ema,xxp[i])
    # shh=np.copy(xxp[o-1])
    # xxp-=shh
    # xxp[hl-1]=np.copy(HOH[1])
    # xxp[hr - 1] = np.copy(HOH[2])
    # xxp+=shh
    # for i in range(len(xxp)):
    #     xxp[i]=np.dot(ema.T,xxp[i])
    # xxp*=angstr
    # ################################
    for i in range(len(HOH)):
        HOH[i]=(ema.T).dot(HOH[i])
    shft=np.copy(x[o-1])
    x-=shft
    x[hl-1]=HOH[1]
    x[hr-1]=HOH[2]
    x+=shft
    for atm in range(len(x)):
        x[atm]=ematOO.T.dot(x[atm])
    y = np.copy(x)*angstr
    return x
    #now x is in OO coordinate system


def export(gem,extr):
    gem2=gem*angstr
    fl = open('testExtract'+extr+'.xyz','w')
    j=0
    k=0
    for i in gem2:
        if j != 1:
            if not np.array_equal(i,np.array([0.,0.,0.])):
                if j<4:
                    k+=1
                else:
                    k+=1
        else:
            k+=1
        j+=1
    fl.write(str(k)+'\nwut\n')
    j=0
    for i in gem2:
        if j != 1:
            if not np.array_equal(i,np.array([0.,0.,0.])):
                if j<4:
                    fl.write("O %5.12f, %5.12f, %5.12f\n" % (i[0],i[1],i[2]))
                else:
                    fl.write("H %5.12f, %5.12f, %5.12f\n" % (i[0],i[1],i[2]))

        else:
            fl.write("O %5.12f, %5.12f, %5.12f\n" % (i[0],i[1],i[2]))
        j+=1
    fl.close()
                
#def outerWaters(lastMoments,coordinates):
def scanAlongCoordinate(modeNumber,config):
    tmatr = np.loadtxt(Tmatname)
    sharedH1vec = tmatr[modeNumber]
    print 'sharedHvec', la.norm(sharedH1vec)
    np.savetxt("tn/tnorm_"+config+"Mode_"+str(modeNumber), [la.norm(sharedH1vec)])
    normedSharedH1 = sharedH1vec / la.norm(sharedH1vec)
    dx = 0.01
    scannedGeoms = np.zeros((numWalkers, 13, 3))
    #export(finalConstructedGeometry, 'testWell')
    print 'numWalkers/2', numWalkers / 2
    scannedGeoms[numWalkers / 2] = np.load("../coordinates/tetramer/"+coordinateSet+".npy")  # assign central point.

    for i in range(numWalkers / 2):
        avGeom = np.zeros(np.shape(eqG))
        g = (dx * (i + 1)) * normedSharedH1
        averageMomentChange = averageMoments + g
        shit=np.copy(averageMomentChange)
        partiallyConstructedGeometry = getTriangle(avGeom, averageMomentChange[27:], averageMomentChange[:9])
        scannedGeoms[(numWalkers / 2) + 1 + i] = getOuters(partiallyConstructedGeometry, averageMomentChange[18:27],
                                                           averageMomentChange[9:18])
    f = -dx * normedSharedH1
    for j in range(numWalkers / 2):
        avGeom = np.zeros(np.shape(eqG))
        averageMomentChange = averageMoments + (j + 1) * f
        partiallyConstructedGeometry = getTriangle(avGeom, averageMomentChange[27:], averageMomentChange[:9])
        scannedGeoms[(numWalkers / 2) - 1 - j] = getOuters(partiallyConstructedGeometry, averageMomentChange[18:27],
                                                           averageMomentChange[9:18])
    np.save("scan/scanned"+config+"Mode_"+str(modeNumber)+".npy", scannedGeoms)
    writeNewWalkers(scannedGeoms, "scan/scanned"+config+"Mode_"+str(modeNumber)+".xyz")
    #print 'start'
    #sub.call(["sh", "move.sh"])
    print 'done'
    #potz = np.loadtxt("../annes_getpot/eng_dip.dat")
    #np.savetxt("../../RyanDVR/1DDVR/Potentials/eng_dip_"+config+"_Mode"+str(modeNumber),potz)
    #plt.plot(au2wn * potz[:, 0])
    #plt.savefig("ScannedPotential"+config+"_Mode"+str(modeNumber)+".png")


eqG = getEqGeom()
eqG,extra = Wfn.molecule.rotateBackToFrame(np.array([eqG,eqG]),2,1,3) #rotate reference to OOO plane
mass = Wfn.molecule.get_mass()
#Now, start getting r-<r> by using transformation matrix.
avGeom = np.zeros(np.shape(eqG)) #Where we will fill in our coordinates
Tmatname ='TransformationMatrix'+coordModel+'.data'
averageMoments = np.loadtxt('averageInternalsWithNewEckart_'+coordinateSet)
averageMoments = np.around(averageMoments,12)
partiallyConstructedGeometry=getTriangle(avGeom,averageMoments[27:],averageMoments[:9])
#averagedGeometry=getOuters(partiallyConstructedGeometry,averageMoments[18:27],averageMoments[9:18])
# finalConstructedGeometry=getOuters(partiallyConstructedGeometry,averageMoments[18:27],averageMoments[9:18])
# xx=Wfn.molecule.rotateBackToFrame(np.array([finalConstructedGeometry,finalConstructedGeometry]),2,1,3)[0]
# export(xx,'desting')
scanAlongCoordinate(int(modeNum),coordModel)