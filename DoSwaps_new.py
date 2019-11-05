import os, numpy as np
import glob
import sys
import DMCClusters as dmc
trimer = False
tetramer = True
if trimer:
    Wfn=dmc.wavefunction('H7O3+', 1)
    print "../coordinates/trimer/" + sys.argv[1]
    fileList = glob.glob("../coordinates/trimer/" + sys.argv[1])
    numAtoms = 10
else:
    Wfn = dmc.wavefunction('H9O4+', 1)
    print "../coordinates/tetramer/" + sys.argv[1]
    fileList = glob.glob("../coordinates/tetramer/" + sys.argv[1])
    numAtoms = 13
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

def swapStuff(myWalkers,dw):
    origWalkers = np.copy(myWalkers)
    if tetramer:
        print 'tetramer'
        if 'llH' in config or 'llD' in config or 'test' in config or 'ref' in config:
            #H5H6Only
            H6H5 = SwapTwoAtoms(np.copy(origWalkers),5, 6)
            print "allH/allD"
            # H7H8Only
            H8H7 = SwapTwoAtoms(np.copy(origWalkers),7, 8)
            # H9H10Only
            H10H9 = SwapTwoAtoms(np.copy(origWalkers),9, 10)
            # SwapTwoAtOnce
            # H6H5,H7H8
            H6H5H8H7 = SwapTwoAtoms(np.copy(H6H5),7, 8)
            # H6H5,H9H10
            H6H5H10H9 = SwapTwoAtoms(np.copy(H6H5),9, 10)
            # H7H8,H9H10
            H8H7H10H9 = SwapTwoAtoms(np.copy(H8H7),9, 10)
            # AllSwapped
            H6H5H8H7H10H9 = SwapTwoAtoms(np.copy(H6H5H8H7),9, 10)
            big8 = np.concatenate((origWalkers,H6H5,H8H7,H10H9,H6H5H8H7,H6H5H10H9,H8H7H10H9,H6H5H8H7H10H9),axis=0)

            swap8_1 = swapChunk(np.copy(big8),np.array([5, 6, 1, 13]), np.array([10, 9, 3, 12]))
            swap8_2 = swapChunk(np.copy(big8),np.array([5, 6, 1, 13]), np.array([7, 8, 2, 11]))
            swap8_3 = swapChunk(np.copy(big8),np.array([7, 8, 2, 11]), np.array([10, 9, 3, 12]))
            swap8_4 = swapChunk(np.copy(swap8_1),np.array([5, 6, 1, 13]), np.array([7, 8, 2, 11]))
            swap8_5 = swapChunk(np.copy(swap8_2),np.array([5, 6, 1, 13]), np.array([10, 9, 3, 12]))

            # new8_56_910[swappedBoy] = swapLobe([5, 6, 1, 13], [10, 9, 3, 12], big8[swappedBoy])
            # new8_56_78[swappedBoy] = swapLobe([5, 6, 1, 13], [7, 8, 2, 11], big8[swappedBoy])
            # new8_78_910[swappedBoy] = swapLobe([7, 8, 2, 11], [10, 9, 3, 12], big8[swappedBoy])
            # new8_56_910_65with78 = copy.deepcopy(new8_56_910)
            # new8_56_78_65with910 = copy.deepcopy(new8_56_78)
            # new8_56_910_65with78[swappedBoy] = swapLobe([5, 6, 1, 13], [7, 8, 2, 11], new8_56_910[swappedBoy])
            # new8_56_78_65with910[swappedBoy] = swapLobe([5, 6, 1, 13], [10, 9, 3, 12], new8_56_78[swappedBoy])

            oswCds = np.concatenate((big8, swap8_1,swap8_2,swap8_3,swap8_4,swap8_5), axis=0)
            oswCdsZ = np.copy(oswCds)
            oswCdsZ[:,:,-1] *= -1.0
            print('caw')
            return np.concatenate((oswCds,oswCdsZ)), np.tile(dw,96)

        elif '1He' in config or '1De' in config:
            print('1He')
            newWalkers45 = SwapTwoAtoms(np.copy(myWalkers),4,5)
            newWalkers67 = SwapTwoAtoms(np.copy(myWalkers), 6, 7)
            newWalkers4567 = SwapTwoAtoms(np.copy(newWalkers45), 6, 7)
            big4 = np.concatenate((origWalkers,newWalkers45,newWalkers67,newWalkers4567),axis=0)
            oswCdsZ = np.copy(big4)
            oswCdsZ[:,:,-1] *= -1.0
            print('caw')
            return np.concatenate((big4,oswCdsZ)), np.tile(dw,8)

        elif '1Hw' in config or '1Dw' in config:
            print('1Hw')
            newWalkers67 = SwapTwoAtoms(np.copy(myWalkers), 6, 7)
            big2 = np.concatenate((origWalkers,newWalkers67),axis=0)
            oswCdsZ = np.copy(big2)
            oswCdsZ[:,:,-1] *= -1.0
            print('caw')
            return np.concatenate((big2,oswCdsZ)), np.tile(dw,4)
    else:
        if 'llH' in config or 'llD' in config:
            fix

            print('caw')
            return np.concatenate((big4, swap4), axis=0), np.tile(dw, 16)
        elif '1He' in config or '1De' in config:
            print('1He')
            newWalkers45 = SwapTwoAtoms(np.copy(myWalkers), 4, 5)
            newWalkers67 = SwapTwoAtoms(np.copy(myWalkers), 6, 7)
            newWalkers4567 = SwapTwoAtoms(np.copy(newWalkers45), 6, 7)
            big4 = np.concatenate((origWalkers, newWalkers45, newWalkers67, newWalkers4567), axis=0)
            oswCdsZ = np.copy(big4)
            oswCdsZ[:, :, -1] *= -1.0
            print('caw')
            return np.concatenate((big4, oswCdsZ)), np.tile(dw, 8)
        elif '1Hw' in config or '1Dw' in config:
            print('1Hw')
            newWalkers67 = SwapTwoAtoms(np.copy(myWalkers), 6, 7)
            big2 = np.concatenate((origWalkers, newWalkers67), axis=0)
            oswCdsZ = np.copy(big2)
            oswCdsZ[:, :, -1] *= -1.0
            print('caw')
            return np.concatenate((big2, oswCdsZ)), np.tile(dw, 4)

for config in fileList:
    print config
    cds=np.load(config)
    if trimer:
        ocom, evecs, kil = Wfn.molecule.eckartRotate(cds, planar=True, lst=[0, 1, 2], dip=True)
    else:
        ocom, evecs, kil = Wfn.molecule.eckartRotate(cds, planar=True, lst=[0,1,2,3], dip=True)
    cds-=ocom[:,np.newaxis]
    cds = np.einsum('knj,kij->kni', evecs.transpose(0, 2, 1), cds).transpose(0, 2, 1)
    dw = np.load(config[:-4]+"_dw.npy")
    newCds,newDw = swapStuff(cds,dw)
    if trimer:
        ocom2, evecs2, kil = Wfn.molecule.eckartRotate(newCds, planar=True, lst=[0, 1, 2], dip=True)
    else:
        ocom2, evecs2, kil = Wfn.molecule.eckartRotate(newCds, planar=True, lst=[0, 1, 2, 3], dip=True)
    newCds-=ocom2[:,np.newaxis]
    newCds = np.einsum('knj,kij->kni', evecs2.transpose(0, 2, 1), newCds).transpose(0, 2, 1)
    if trimer:
        np.save("../coordinates/trimer/ffinal_"+config[-8:],newCds)
        np.save("../coordinates/trimer/ffinal_"+config[-8:-4]+"_dw",newDw)
    else:
        np.save("../coordinates/tetramer/tetffinal_" + config[-8:], newCds)
        np.save("../coordinates/tetramer/tetffinal_" + config[-8:-4] + "_dw", newDw)
    print "Done"
