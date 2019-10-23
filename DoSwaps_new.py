import os, numpy as np
import glob
import sys
import DMCClusters as dmc
Wfn=dmc.wavefunction('H7O3+', 1)
# xx=np.load("../coordinates/tetramer/rachelDips.npy")*ang2bohr
# xx=np.load("partrig_structures.npy")*ang2bohr
# dips = np.load("partrig_dips_noshift.npy")


print "../coordinates/trimer/"+sys.argv[1]
fileList = glob.glob("../coordinates/trimer/"+sys.argv[1])
numAtoms = 10
numWalkers=0
angstr=0.529177
#initalizing arrays to store atom coordinates + walker information

def SwapTwoAtoms(walkers,atm1,atm2):
    atm1-=1
    atm2-=1
    origAtom1 = np.copy(walkers[:,atm1])
    walkers[:,atm1] = walkers[:,atm2]
    walkers[:,atm2] = origAtom1
    return walkers

def swapChunk(walkers,atmG1,atmG2):
    atmG1-=1
    atmG2-=1
    oldChunk = np.copy(walkers[:,atmG1])
    walkers[:,atmG1] = walkers[:,atmG2]
    walkers[:, atmG2] = oldChunk
    return walkers

def swapChunkComplex(walkers,pairs):
    atmG1 = np.column_stack((np.tile(9,len(walkers)),np.tile(1,len(walkers)),np.tile(4,len(walkers)),np.tile(5,len(walkers))))
    pairs2 = np.copy(pairs)
    pairs2[pairs==6] = 999
    pairs2[pairs==7] = 6
    pairs2[pairs2==999] = 7
    print pairs.shape,pairs2.shape
    atmG2 = np.column_stack((np.tile(10,len(walkers)),np.tile(2,len(walkers)),(pairs),(pairs2)))
    atmG1-=1
    atmG2-=1
    print atmG1,atmG1.shape
    atmG1 = np.int_(atmG1)
    atmG2 = np.int_(atmG2)
    #oldChunk = np.copy(walkers[:,atmG1.astype(int)])
    for i in range(len(walkers)):
        #print(i)
        oldChunk = np.copy(walkers[i,atmG1[i]])
        walkers[i,atmG1[i]] = walkers[[i],atmG2[i]]
        walkers[i, atmG2[i]] = oldChunk
    return walkers

def swapStuff(myWalkers,dw):
    origWalkers = np.copy(myWalkers)
    if 'llH' in config or 'llD' in config or '1Hh' in config or '1Dh' in config or 'test' or 'ref' in config:
        #H4H5
        newWalkers45 = SwapTwoAtoms(np.copy(myWalkers),4,5)
        newWalkers67 = SwapTwoAtoms(np.copy(myWalkers), 6, 7)
        newWalkers4567 = SwapTwoAtoms(np.copy(newWalkers45), 6, 7)
        big4 = np.concatenate((origWalkers,newWalkers45,newWalkers67,newWalkers4567),axis=0)
        # pairS = np.sign(big4[:,[4],-1].flatten())*np.sign(big4[:,[6],-1].flatten() )
        # pairS[pairS==-1.0] = 7
        # pairS[pairS==1.0] = 6
        # swap4 = swapChunkComplex(np.copy(big4),pairS)
        swap4 = swapChunk(np.copy(big4),np.array([9,1,4,5]),np.array([10,2,6,7]))
        print('caw')
        return np.concatenate((big4,swap4),axis=0), np.tile(dw,8)
    #
    # elif '1He' in config or '1De' in config:
    #    print "Eigen Symmetrized"
    #    # H4H5Only
    #    H4H5 = swapTwoCoordinates(4,5,myWalkers)
    #    # H6H7Only
    #    H6H7 = swapTwoCoordinates(6,7,myWalkers)
    #    #SwapTwoAtOnce
    #    H4H5H6H7 = swapTwoCoordinates(6,7,H4H5)
    #    big4 = [myWalkers,H4H5,H6H7,H4H5H6H7]
    #    return big4+negZ(copy.deepcopy(big4))
    #
    # elif '1Hw' in config or '1Dw' in config:
    #     # H4H5Only
    #     print "Free Water Symmetrized"
    #     H4H5 = swapTwoCoordinates(4, 5, myWalkers)
    #     big2 = [myWalkers, H4H5]
    #     return big2+negZ(copy.deepcopy(big2))


for config in fileList:
    print config

    # numWalkers,xyzCoords, weightArray = putWalkersInArray2(config)
    cds=np.load(config)
    ocom, evecs, kil = Wfn.molecule.eckartRotate(cds, planar=True, lst=[0, 1, 2], dip=True)
    cds-=ocom[:,np.newaxis]
    cds = np.einsum('knj,kij->kni', evecs.transpose(0, 2, 1), cds).transpose(0, 2, 1)
    # cds[[4-1,5-1]]=cds[[5-1,4-1]]
    dw = np.load(config[:-4]+"_dw.npy")

    newCds,newDw = swapStuff(cds,dw)
    # newCds = Wfn.molecule.rotateBackToFrame(newCds,3,2,1)
    # flout = open("refCds.xyz", "w+")
    # Wfn.molecule.printCoordsToFile(newCds, flout)
    newCdsp = np.copy(newCds)
    newCdsp[:,:,-1]*=-1.0
    newCds = np.concatenate((newCds,newCdsp),axis=0)
    ocom, evecs, kil = Wfn.molecule.eckartRotate(newCds, planar=True, lst=[0, 1, 2], dip=True)
    newCds-=ocom[:,np.newaxis]
    newCds = np.einsum('knj,kij->kni', evecs.transpose(0, 2, 1), newCds).transpose(0, 2, 1)

    np.save("../coordinates/trimer/ffinal_"+config[-8:],newCds)
    np.save("../coordinates/trimer/ffinal_"+config[-8:-4]+"_dw",np.concatenate((newDw,newDw)))

    print "Done"
