import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import sys
import DMCClusters as dmc
import subprocess as sub
angstr = 0.529177
au2wn=219474.63

coordinateSet = sys.argv[1]
testName = sys.argv[2]
kill = sys.argv[3]
Wfn=dmc.wavefunction('H7O3+', 1)
path='../spectra/'
gPath = '../Gmats/trimer/'
GfileName = gPath+coordinateSet+'.gmat'
numWalkers = 201
internalName = ['rOH8', 'thH8', 'phiH8', 'rOH9', 'thH9', 'phiH9', 'rOH10', 'thH10', 'phiH10',
                         'th_627', 'phi_627', 'xi_627', 'th_514', 'phi_514', 'xi_514', 'rOH_41',
                         'rOH_51', 'aHOH_451', 'rOH_26', 'rOH_27', 'aHOH_267', 'rOO_1', 'rOO_2', 'aOOO']
def writeNewWalkers(walkz,title):
    coordAr = ['O','O','O','H','H','H','H','H','H','H']
    walkz2 = walkz*angstr
    writeFile = open(title,"w")
    for walker in walkz2:
        writeFile.write("%d\n \n" % walkz.shape[1])
        for k,atm in enumerate(walker):
            writeFile.write("%s %5.12f  %5.12f  %5.12f\n" % (coordAr[k],atm[0],atm[1],atm[2]))
        writeFile.write("\n")
    writeFile.close()

def getEqGeom():
    eqGeom = np.array(
        [
            [-1.55099190 ,  1.94067311 ,  0.14704161],
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
             [-0.63272209,   0.72926470,  -0.20069859]
        ])
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
    return a

def oxygenFrame(r1,r2,theta):
    o1 = np.array([r1 * np.cos(theta / 2) , r1 * np.sin(theta / 2), 0.0])
    o2 = np.array([r2*np.cos(theta/2) , -r2*np.sin(theta/2),0.0])
    o3 = [0.0 , 0.0 , 0.0]
    # o1 = [r1 , 0.0 , 0.0]
    # rz = np.array([
    #     [np.cos(theta),-1.0*np.sin(theta),0.0],
    #     [np.sin(theta), np.cos(theta), 0.0],
    #     [0,0,1]
    # ])
    # o2 = rz.dot([r2,0.0,0.0]) #place r2 on x axis and then rotation to xy plane
    return np.array([o1,o2,o3])

def sphToCart(r,theta,phi):
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return np.array([x,y,z])
def sharedHs(internz):
    """
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    """
    h8 = internz[:3]
    xyzH8 = sphToCart(h8[0],h8[1],h8[2])
    h9 = internz[3:6]
    xyzH9 = sphToCart(h9[0],h9[1],h9[2])
    h10 = internz[6:9]
    xyzH10 = sphToCart(h10[0],h10[1],h10[2])
    return np.array([xyzH8,xyzH9,xyzH10])

def placeWater(xx,thPhXi_r1r2th,oind,Hind1,Hind2):
    if oind == 1:
        other = 2
    else:
        other = 1
    thE,phE,xiE = thPhXi_r1r2th[0:3]
    r1,r2,theta = thPhXi_r1r2th[3:]
    #rotate molecule such that OO is on x axis,
    X = (xx[oind - 1] - xx[3 - 1])/ la.norm(xx[oind - 1] - xx[3 - 1])
    crs = np.cross(xx[1 - 1] - xx[3 - 1], xx[2 - 1] - xx[3 - 1])
    Z = crs / la.norm(crs)
    Y = np.cross(Z, X)
    xx = Wfn.molecule.rotateBackToFrame(np.array([xx,xx]),3,oind,other)[0]
    Xn = (xx[oind - 1] - xx[3 - 1]) / la.norm(xx[oind - 1] - xx[3 - 1])
    crsn = np.cross(xx[1 - 1] - xx[3 - 1], xx[2 - 1] - xx[3 - 1])
    Zn = crsn / la.norm(crsn)
    Yn = np.cross(Zn, Xn)


    # if Zn[-1] == -1: #rotate 180 to get the axes aligned the correct way
    #     rx = np.array([
    #         [1,0,0],
    #         [0,np.cos(np.pi),-np.sin(np.pi)],
    #         [0,np.sin(np.pi),np.cos(np.pi)]
    #     ])
    #     for i in range(len(xx)):
    #         xx[i] = rx.dot(xx[i])


    #here on, it's rotate to "oo axis"
    #well, let's get weird
    """
    exx = np.copy(x)
        x = np.copy(y)
        y = np.copy(z)
        z = np.copy(exx)
        print 'lets get weird' #big X Y Z are the OO plane
        exX = np.copy(X)
        X = np.copy(Z)
        Z = np.copy(Y)
        Y = np.copy(exX)
    """
    flipAx = np.array([
        [0,1,0],
        [0,0,1],
        [1,0,0.]
    ])
    for i in range(len(xx)):
            xx[i] = flipAx.dot(xx[i])
    #/done with OO axis
    xx -= xx[oind-1]
    Xp = (xx[oind - 1] - xx[3 - 1]) / la.norm(xx[oind - 1] - xx[3 - 1])
    crsp = np.cross(xx[1 - 1] - xx[3 - 1], xx[2 - 1] - xx[3 - 1])
    Zp = crsp / la.norm(crsp)
    Yp = np.cross(Zp, Xp)
    #place HOH when bisector on x axis
    HOH = np.zeros((3,3))
    # HOH[1,:]=np.array([r1 * np.cos(theta / 2), r1 * np.sin(theta / 2), 0.])
    # HOH[2,:]=np.array([r2*np.cos(theta/2),-r2*np.sin(theta/2),0.])
    HOH[1,:]=np.array([r1 * np.cos(theta / 2) , 0.0 , r1 * np.sin(theta / 2)])
    HOH[2,:]=np.array([r2*np.cos(theta/2) , 0.0 , -r2*np.sin(theta/2)])
    # HOH[1,:]=np.array([0.0, r1 * np.cos(theta / 2) , r1 * np.sin(theta / 2)])
    # HOH[2,:]=np.array([0.0, r2*np.cos(theta/2) , -r2*np.sin(theta/2)])

    # prelim = np.copy(xx)
    # prelim[Hind1-1] = HOH[1]
    # prelim[Hind2 - 1] = HOH[2]
    # xx[Hind1-1] = HOH[1]
    # xx[Hind2-1] = HOH[2]

    # writeNewWalkers(np.array([prelim,prelim]),'testing_prelim.xyz')
    #for now, we are not going to worry about eulers
    # thE = np.pi
    # phE = np.pi/2
    # xiE *= -1.0
    weird = np.array([
        [0,0,1],
        [1,0,0],
        [0,1,0.]
    ])
    for i in range(len(HOH)):
        HOH[i] = np.dot(weird, HOH[i])
    hRotM = getEulerMat(thE,phE,xiE)
    # hRotM = weird.dot(hRotM)
    # test = Wfn.molecule.extractEulers(np.array([hRotM,hRotM])) #looks ok
    for i in range(len(HOH)):
        HOH[i] = np.dot(hRotM, HOH[i])
    xx[Hind1-1] = HOH[1]
    xx[Hind2-1] = HOH[2]
    return xx

def genCarts(internals):
    """
    :param internals: 3N-6 coordinates
    :return: Cartesian coordinates
    """
    fullCart = np.zeros((10,3))
    roo1,roo2,aooo = internals[-3:]
    beginX = oxygenFrame(roo1,roo2,aooo)

    beginX = np.array([beginX,beginX])
    com, eckVecs, killList  = Wfn.molecule.eckartRotate(beginX,planar=True,lst=[1-1,2-1,3-1],dip=True)
    beginX -= com[:, np.newaxis, :]
    beginX = np.einsum('knj,kij->kni', eckVecs.transpose(0, 2, 1), beginX).transpose(0, 2, 1)[0]

    Hs = sharedHs(internals[:9])
    fullCart[0:3] = beginX
    fullCart -= fullCart[3-1]
    fullCart[8-1] = Hs[0]
    fullCart[9 - 1] = Hs[1]
    fullCart[10 - 1] = Hs[2]
    fullCart += fullCart[3-1]

    fullCart = placeWater(fullCart,internals[[9,10,11,18,19,20]],2,6,7)
    fullCart = np.array([fullCart,fullCart])
    com, eckVecs, killList = Wfn.molecule.eckartRotate(fullCart, planar=True, lst=[1 - 1, 2 - 1, 3 - 1], dip=True)
    fullCart -= com[:, np.newaxis, :]
    fullCart = np.einsum('knj,kij->kni', eckVecs.transpose(0, 2, 1), fullCart).transpose(0, 2, 1)[0]
    # writeNewWalkers(np.array([fullCart,fullCart]),'test_firstWater.xyz')

    fullCart = placeWater(fullCart,internals[[12,13,14,15,16,17]], 1,4,5)
    fullCart = np.array([fullCart, fullCart])
    com, eckVecs, killList = Wfn.molecule.eckartRotate(fullCart, planar=True, lst=[1 - 1, 2 - 1, 3 - 1], dip=True)
    fullCart -= com[:, np.newaxis, :]
    fullCart = np.einsum('knj,kij->kni', eckVecs.transpose(0, 2, 1), fullCart).transpose(0, 2, 1)[0]
    writeNewWalkers(np.array([fullCart,fullCart]),'test.xyz')
    # bnd = np.zeros((10,10))
    # for i in range(10):
    #     for j in range(10):
    #         bnd[i,j] = la.norm(fullCart[i]-fullCart[j])
    #
    # for i in bnd:
    #     print(i)
    # print('hi')
    return fullCart

def scanNormalModes(Tmatname,modeNumber,numWalkers,fll):
    averageMoments = np.loadtxt("averageInternalsWithNewEckart_"+fll)
    #for eq cds
    averageMoments = np.loadtxt("averageInternalsWithNewEckart_" + "eq_allH_rn_spc_eqH")
    fll = "eq_allH_rn_spc_eqH"
    #for eq cds
    tmatr = np.loadtxt(Tmatname)
    normedVec = tmatr[modeNumber]
    np.savetxt("tn/tnorm_"+fll+"Mode_"+str(modeNumber), [la.norm(normedVec)])
    normedVec /= la.norm(normedVec)
    dx = 0.01
    scannedGeoms = np.zeros((numWalkers, 10, 3))
    print 'numWalkers/2', numWalkers / 2

    rotRef = genCarts(averageMoments)
    rotRef = np.array([rotRef, rotRef])
    com, eckVecs, killList = Wfn.molecule.eckartRotate(rotRef, planar=True, lst=[1 - 1, 2 - 1, 3 - 1], dip=True)
    rotRef -= com[:, np.newaxis, :]
    rotRef = np.einsum('knj,kij->kni', eckVecs.transpose(0, 2, 1), rotRef).transpose(0, 2, 1)[0]
    scannedGeoms[numWalkers / 2] = rotRef  # assign central point, saddle pt geom.

    for i in range(numWalkers / 2):
        g = (dx * (i + 1)) * normedVec
        internalChange = averageMoments + g
        scannedGeoms[(numWalkers / 2) + 1 + i] = genCarts(internalChange)
    f = -dx * normedVec
    for j in range(numWalkers / 2):
        internalChange = averageMoments + (j + 1) * f
        scannedGeoms[(numWalkers / 2) - 1 - j] = genCarts(internalChange)
    np.save("scan/scanFiles/scanned_"+fll+"_"+str(modeNumber)+".npy", scannedGeoms)
    writeNewWalkers(scannedGeoms,"scan/scanFiles/scanned_"+fll+"_"+str(modeNumber)+".xyz")

    print 'done'

def finDiffFCs(Tmatname,modeNumber,fll):
    stence = 5
    averageMoments = np.loadtxt("averageInternalsWithNewEckart_" + fll)
    # averageMoments = np.loadtxt("averageInternalsWithNewEckart_" + "eq_allH_rn_spc_eqH")
    tmatr = np.loadtxt(Tmatname)
    normedVec = tmatr[modeNumber]
    np.savetxt("tn/tnorm_" + fll + "Mode_" + str(modeNumber), [la.norm(normedVec)])
    normedVec /= la.norm(normedVec)
    dx = 0.01
    scannedGeoms = np.zeros((stence, 10, 3))
    print 'stence/2', stence / 2

    rotRef = genCarts(averageMoments)
    rotRef = np.array([rotRef, rotRef])
    com, eckVecs, killList = Wfn.molecule.eckartRotate(rotRef, planar=True, lst=[1 - 1, 2 - 1, 3 - 1], dip=True)
    rotRef -= com[:, np.newaxis, :]
    rotRef = np.einsum('knj,kij->kni', eckVecs.transpose(0, 2, 1), rotRef).transpose(0, 2, 1)[0]
    scannedGeoms[stence / 2] = rotRef  # assign central point, saddle pt geom.

    for i in range(stence / 2):
        g = (dx * (i + 1)) * normedVec
        internalChange = averageMoments + g
        scannedGeoms[(stence / 2) + 1 + i] = genCarts(internalChange)
    f = -dx * normedVec
    for j in range(stence / 2):
        internalChange = averageMoments + (j + 1) * f
        scannedGeoms[(stence / 2) - 1 - j] = genCarts(internalChange)
    np.save("scan/finDiff/scanned_" + fll + "_" + str(modeNumber) + ".npy", scannedGeoms)
    writeNewWalkers(scannedGeoms, "scan/finDiff/scanned_" + fll + "_" + str(modeNumber) + ".xyz")
    
    
# eqG = getEqGeom()
# refG = Wfn.molecule.pullTrimerRefPos(dip=True)
Tmatname ='TransformationMatrix'+coordinateSet+'_'+testName+'_'+kill+'.datatest'
fll = coordinateSet+'_'+testName+'_'+kill
for i in range(24):
    scanNormalModes(Tmatname,i,201,fll)
# averageMoments = np.loadtxt('averageInternalsWithNewEckart_'+fll)
# cartCds = genCarts(averageMoments)
