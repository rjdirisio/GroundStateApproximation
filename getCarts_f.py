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
numWalkers = 201
if 'tet' not in coordinateSet and 'h9o4' not in coordinateSet:
    Wfn=dmc.wavefunction('H7O3+', 1)
    internalName = ['rOH8', 'thH8', 'phiH8', 'rOH9', 'thH9', 'phiH9', 'rOH10', 'thH10', 'phiH10',
                             'th_627', 'phi_627', 'xi_627', 'th_514', 'phi_514', 'xi_514', 'rOH_41',
                             'rOH_51', 'aHOH_451', 'rOH_26', 'rOH_27', 'aHOH_267', 'rOO_1', 'rOO_2', 'aOOO']
else:
    Wfn=dmc.wavefunction('H9O4+', 1)
    internalName = ['rOH11', 'rOH12', 'rOH13', 'thH11', 'thH12', 'thH13', 'phH11', 'phH12', 'phH13',
                    'theta651', 'phi651', 'Chi651', 'theta1039', 'phi1039', 'Chi1039', 'theta728', 'phi728',
                    'Chi728', 'rOH5', 'rOH6', 'HOH516', 'rOH7', 'rOH8', 'HOH728', 'rOH9', 'rOH10', 'HOH9310', 'rO1O2',
                    'rO1O3', 'rO2O3','xo4', 'yo4', 'zo4']
def writeNewWalkers(walkz,title,syst):
    if syst == 'trimer':
        coordAr = ['O','O','O','H','H','H','H','H','H','H']
    else:
        coordAr = ['O','O', 'O', 'O', 'H', 'H','H', 'H', 'H', 'H', 'H', 'H', 'H']
    walkz2 = walkz*angstr
    writeFile = open(title,"w")
    for walker in walkz2:
        writeFile.write("%d\n \n" % walkz.shape[1])
        for k,atm in enumerate(walker):
            writeFile.write("%s %5.12f  %5.12f  %5.12f\n" % (coordAr[k],atm[0],atm[1],atm[2]))
        writeFile.write("\n")
    writeFile.close()

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

def placeWater_trimer(xx,thPhXi_r1r2th,oind,Hind1,Hind2):
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
    #well, let's get weird
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
    HOH[1,:]=np.array([r1 * np.cos(theta / 2) , 0.0 , r1 * np.sin(theta / 2)])
    HOH[2,:]=np.array([r2*np.cos(theta/2) , 0.0 , -r2*np.sin(theta/2)])
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

def genCarts_trimer(internals):
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
    fullCart = placeWater_trimer(fullCart,internals[[9,10,11,18,19,20]],2,6,7)
    fullCart = np.array([fullCart,fullCart])
    com, eckVecs, killList = Wfn.molecule.eckartRotate(fullCart, planar=True, lst=[1 - 1, 2 - 1, 3 - 1], dip=True)
    fullCart -= com[:, np.newaxis, :]
    fullCart = np.einsum('knj,kij->kni', eckVecs.transpose(0, 2, 1), fullCart).transpose(0, 2, 1)[0]
    fullCart = placeWater_trimer(fullCart,internals[[12,13,14,15,16,17]], 1,4,5)
    fullCart = np.array([fullCart, fullCart])
    com, eckVecs, killList = Wfn.molecule.eckartRotate(fullCart, planar=True, lst=[1 - 1, 2 - 1, 3 - 1], dip=True)
    fullCart -= com[:, np.newaxis, :]
    fullCart = np.einsum('knj,kij->kni', eckVecs.transpose(0, 2, 1), fullCart).transpose(0, 2, 1)[0]
    # writeNewWalkers(np.array([fullCart,fullCart]),'test.xyz',system='trimer')
    return fullCart

def placeWater_tetramer(xx,thPhXi_r1r2th,oind,Hind1,Hind2):
    thE,phE,xiE = thPhXi_r1r2th[0:3]
    r1,r2,theta = thPhXi_r1r2th[3:]
    if oind == 1:
        other = 3
    elif oind == 2:
        other = 1
    elif oind == 3:
        other = 2
    #rotate molecule such that OO is on x axis,
    # X = (xx[oind - 1] - xx[4 - 1])/ la.norm(xx[oind - 1] - xx[4 - 1])!
    # crs = np.cross(xx[1 - 1] - xx[4 - 1], xx[2 - 1] - xx[4 - 1])!
    # Z = crs / la.norm(crs)
    # Y = np.cross(Z, X)
    xx = Wfn.molecule.rotateBackToFrame(np.array([xx,xx]),4,oind,other)[0]
    # Xn = (xx[oind - 1] - xx[4 - 1]) / la.norm(xx[oind - 1] - xx[4 - 1])
    # crsn = np.cross(xx[1 - 1] - xx[4 - 1], xx[2 - 1] - xx[4 - 1])!
    # Zn = crsn / la.norm(crsn)
    # Yn = np.cross(Zn, Xn)
    #well, let's get weird
    flipAx = np.array([
        [0,1,0],
        [0,0,1],
        [1,0,0.]
    ])
    for i in range(len(xx)):
            xx[i] = flipAx.dot(xx[i])
    #/done with OO axis
    xx -= xx[oind-1]
    # Xp = (xx[oind - 1] - xx[4 - 1]) / la.norm(xx[oind - 1] - xx[4 - 1])
    # crsp = np.cross(xx[1 - 1] - xx[4 - 1], xx[2 - 1] - xx[4 - 1])
    # Zp = crsp / la.norm(crsp)
    # Yp = np.cross(Zp, Xp)
    #place HOH when bisector on x axis
    HOH = np.zeros((3,3))
    HOH[1,:]=np.array([r1 * np.cos(theta / 2) , 0.0 , r1 * np.sin(theta / 2)])
    HOH[2,:]=np.array([r2*np.cos(theta/2) , 0.0 , -r2*np.sin(theta/2)])
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


def outerOxygenFrame(oo12,oo13,oo23):
    """For tetramer"""
    # A = o2,B = o1, C = o3
    # a = rBC = oo13
    # b = rAC = oo23
    # c = rAB = oo12
    # o2 = [0.0,0.0,0.0]
    # o1 = [oo12,0.0,0.0]
    # cBAC = (oo23 ** 2 + oo12 ** 2 - oo13 ** 2) / (2 * oo23 * oo12)
    # yc = np.sqrt(oo23**2 - cBAC**2)
    # o3 = [oo23*cBAC, yc,0.0]
    # switch o2 and o3
    #B = O3
    #A = O1
    #C = 02
    a = oo23
    b = oo12
    c = oo13
    A = [0.0,0.0,0.0]
    B = [c,0.0,0.0]
    xC = b*((b**2+c**2-a**2)/(2*b*c))
    yC = np.sqrt(b**2 - xC**2)
    C = np.array([xC,yC,0.0])
    return np.array([A,C,B])

def placeHyd(xx,oxSkele,hyd,dispXYZ):
    """For tetramer"""
    xx[[0,1,2]] = oxSkele
    xx[[4-1,11-1,12-1,13-1]] = hyd
    xx[[4-1,11-1,12-1,13-1]] += dispXYZ
    return xx
def genCarts_tetramer(internals):
    """
        internalName = ['rOH11', 'rOH12', 'rOH13', 'thH11', 'thH12', 'thH13', 'phH11', 'phH12', 'phH13',
                    'theta651', 'phi651', 'Chi651', 'theta1039', 'phi1039', 'Chi1039', 'theta728', 'phi728',
                    'Chi728', 'rOH5', 'rOH6', 'HOH516', 'rOH7', 'rOH8', 'HOH728', 'rOH9', 'rOH10', 'HOH9310', 'rO1O2',
                    'rO1O3', 'rO2O3','xo4', 'yo4', 'zo4']

        internalName = ['rOH11', 'rOH12', 'rOH13', 'thH11','thH12','thH13','phH11','phH12','phH13',
                     'theta728','phi728', 'Chi728','theta1039', 'phi1039', 'Chi1039', 'theta651', 'phi651',
                     'Chi651', 'rOH5', 'rOH6','HOH516', 'rOH7', 'rOH8', 'HOH728','rOH9', 'rOH10', 'HOH9310', 'rO1O2', 'rO1O3', 'rO2O3',
                     'xo4', 'yo4', 'zo4']
    """
    fullCart = np.zeros((13,3))
    ro1o2, ro1o3, ro2o3 = internals[-6:-3]
    beginX = outerOxygenFrame(ro1o2, ro1o3, ro2o3)
    beginX = np.array([beginX, beginX])
    com, eckVecs, killList = Wfn.molecule.eckartRotate(beginX, planar=True, lst=[1 - 1, 2 - 1, 3 - 1], dip=True)
    beginX -= com[:, np.newaxis, :]
    beginX = np.einsum('knj,kij->kni', eckVecs.transpose(0, 2, 1), beginX).transpose(0, 2, 1)[0]

    hydronium = np.zeros((4,3))
    hydronium[1:] = sharedHs(internals[[0,3,6,1,4,7,2,5,8]])
    lst = [4 - 1, 11 - 1, 12 - 1, 13 - 1]
    mass = Wfn.molecule.get_mass()
    mass = mass[lst]
    ocomH = np.dot(mass, hydronium) / np.sum(mass)
    hydronium -= ocomH[np.newaxis,:]
    fullCart = placeHyd(fullCart, beginX,hydronium,internals[-3:])
    writeNewWalkers(np.array([fullCart,fullCart]),'test.xyz',syst='tetramer')
    #place waters
    #water 516
    water1Ints = internals[[15,16,17,18,19,20]]
    fullCart = placeWater_tetramer(fullCart,water1Ints,1,5,6)
    fullCart = np.array([fullCart,fullCart])
    com, eckVecs, killList = Wfn.molecule.eckartRotate(fullCart, planar=True, lst=[1 - 1, 2 - 1, 3 - 1,4-1], dip=True)
    fullCart -= com[:, np.newaxis, :]
    fullCart = np.einsum('knj,kij->kni', eckVecs.transpose(0, 2, 1), fullCart).transpose(0, 2, 1)[0]
    #water 728
    water1Ints = internals[[9,10,11,21,22,23]]
    fullCart = placeWater_tetramer(fullCart,water1Ints,2,7,8)
    fullCart = np.array([fullCart,fullCart])
    com, eckVecs, killList = Wfn.molecule.eckartRotate(fullCart, planar=True, lst=[1 - 1, 2 - 1, 3 - 1,4-1], dip=True)
    fullCart -= com[:, np.newaxis, :]
    fullCart = np.einsum('knj,kij->kni', eckVecs.transpose(0, 2, 1), fullCart).transpose(0, 2, 1)[0]
    #water 9310
    water1Ints = internals[[12,13,14,24,25,26]]
    fullCart = placeWater_tetramer(fullCart,water1Ints,3,9,10)
    fullCart = np.array([fullCart,fullCart])
    com, eckVecs, killList = Wfn.molecule.eckartRotate(fullCart, planar=True, lst=[1 - 1, 2 - 1, 3 - 1,4-1], dip=True)
    fullCart -= com[:, np.newaxis, :]
    fullCart = np.einsum('knj,kij->kni', eckVecs.transpose(0, 2, 1), fullCart).transpose(0, 2, 1)[0]
    writeNewWalkers(np.array([fullCart,fullCart]),'test2.xyz','tetramer')
    return fullCart


def scanNormalModes(Tmatname,modeNumber,numWalkers,fll):
    testing=False
    if testing:
        print('testing')
    else:
        tmatr = np.loadtxt(Tmatname)
        normedVec = tmatr[modeNumber]
        np.savetxt("tn/tnorm_"+fll+"Mode_"+str(modeNumber), [la.norm(normedVec)])
        normedVec /= la.norm(normedVec)
    dx = 0.01
    print 'numWalkers/2', numWalkers / 2
    if 'tet' not in fll:
        scannedGeoms = np.zeros((numWalkers, 10, 3))
        print('trimer scan')
        averageMoments = np.loadtxt("averageInternalsWithNewEckart_" + fll)
        rotRef = genCarts_trimer(averageMoments)
        rotRef = np.array([rotRef, rotRef])
        com, eckVecs, killList = Wfn.molecule.eckartRotate(rotRef, planar=True, lst=[1 - 1, 2 - 1, 3 - 1], dip=True)
        rotRef -= com[:, np.newaxis, :]
        rotRef = np.einsum('knj,kij->kni', eckVecs.transpose(0, 2, 1), rotRef).transpose(0, 2, 1)[0]
        sysStr = 'trimer'
    else:
        scannedGeoms = np.zeros((numWalkers, 13, 3))
        print('tetramer scan')
        if testing:
            averageMoments = np.loadtxt("averageInternalsWithNewEckart_h9o4_minc")
        else:
            averageMoments = np.loadtxt("averageInternals_" + fll)
        rotRef = genCarts_tetramer(averageMoments)
        sysStr = 'tetramer'
    scannedGeoms[numWalkers / 2] = rotRef  # assign central point, saddle pt geom.

    for i in range(numWalkers / 2):
        g = (dx * (i + 1)) * normedVec
        internalChange = averageMoments + g
        if 'tet' in fll:
            scannedGeoms[(numWalkers / 2) + 1 + i] = genCarts_tetramer(internalChange)
        else:
            scannedGeoms[(numWalkers / 2) + 1 + i] = genCarts_trimer(internalChange)
    f = -dx * normedVec
    for j in range(numWalkers / 2):
        internalChange = averageMoments + (j + 1) * f
        if 'tet' in fll:
            scannedGeoms[(numWalkers / 2) - 1 - j] = genCarts_tetramer(internalChange)
        else:
            scannedGeoms[(numWalkers / 2) - 1 - j] = genCarts_trimer(internalChange)
    np.save("scan/scanFiles/scanned_"+fll+"_"+str(modeNumber)+".npy", scannedGeoms)
    writeNewWalkers(scannedGeoms,"scan/scanFiles/scanned_"+fll+"_"+str(modeNumber)+".xyz",sysStr)
    print 'done'


Tmatname ='TransformationMatrix'+coordinateSet+'_'+testName+'_'+kill+'.datatest'
fll = coordinateSet+'_'+testName+'_'+kill
if 'tet' in fll:
    nvibs = 33
else:
    nvibs = 24
for i in range(nvibs):
    scanNormalModes(Tmatname,i,201,fll)
