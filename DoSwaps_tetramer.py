import os, numpy as np
import sys
from numpy import linalg as la
import glob
import copy
print "../coordinates/tetramer/"+sys.argv[1]
fileList = glob.glob("../coordinates/tetramer/"+sys.argv[1])

numAtoms = 13
numWalkers=0
angstr=0.529177
#initalizing arrays to store atom coordinates + walker information
class Walker:
    def __init__(self):
        self.coords = np.zeros((numAtoms,3))
        self.weight = 0.0
        self.WalkerNo = 0

def flipCoords(walkerCrds):
    newCs = np.negative(walkerCrds)
    return newCs

def giveWalkersTheirCoordinates(cartCoords,myWalkers):
    u = range(numWalkers*3)
    j = 0
    u = u[::3]
    for v in u:
        a = np.copy(cartCoords[:, v:v + 3])
        myWalkers[j].coords=a
        j+=1

def giveWalkersTheirWeight(wtArray):
    for b in range(numWalkers):
        myWalkers[b].weight = wtArray[b]
        myWalkers[b].WalkerNo = b

def writeTheWalkers(someCfww,wfn):
    wfn = wfn+config[-4:]
    wfi = open(wfn,"w")
    for each in someCfww:
        for d in range(len(each)):
            wfi.write('13\n')
            wfi.write('%5.12f 0 0 0 0\n' % each[d].weight)
            wfi.write('O %5.12f %5.12f %5.12f\n' % (each[d].coords[0][0],each[d].coords[0][1],each[d].coords[0][2]))
            wfi.write('O %5.12f %5.12f %5.12f\n' % (each[d].coords[1][0],each[d].coords[1][1],each[d].coords[1][2]))
            wfi.write('O %5.12f %5.12f %5.12f\n' % (each[d].coords[2][0],each[d].coords[2][1],each[d].coords[2][2]))
            wfi.write('O %5.12f %5.12f %5.12f\n' % (each[d].coords[3][0],each[d].coords[3][1],each[d].coords[3][2]))
            wfi.write('H %5.12f %5.12f %5.12f\n' % (each[d].coords[4][0],each[d].coords[4][1],each[d].coords[4][2]))
            wfi.write('H %5.12f %5.12f %5.12f\n' % (each[d].coords[5][0],each[d].coords[5][1],each[d].coords[5][2]))
            wfi.write('H %5.12f %5.12f %5.12f\n' % (each[d].coords[6][0],each[d].coords[6][1],each[d].coords[6][2]))
            wfi.write('H %5.12f %5.12f %5.12f\n' % (each[d].coords[7][0],each[d].coords[7][1],each[d].coords[7][2]))
            wfi.write('H %5.12f %5.12f %5.12f\n' % (each[d].coords[8][0],each[d].coords[8][1],each[d].coords[8][2]))
            wfi.write('H %5.12f %5.12f %5.12f\n' % (each[d].coords[9][0],each[d].coords[9][1],each[d].coords[9][2]))
            wfi.write('H %5.12f %5.12f %5.12f\n' % (each[d].coords[10][0],each[d].coords[10][1],each[d].coords[10][2]))
            wfi.write('H %5.12f %5.12f %5.12f\n' % (each[d].coords[11][0],each[d].coords[11][1],each[d].coords[11][2]))
            wfi.write('H %5.12f %5.12f %5.12f\n' % (each[d].coords[12][0],each[d].coords[12][1],each[d].coords[12][2]))
            wfi.write('\n')
    wfi.close()


def putWalkersInArray2(cfg):
    f = open(cfg , "r")
    linez = f.readlines()
    nWalkers = len(linez)/16
    wtAr = np.zeros(nWalkers)
    filledCoords = np.zeros((numAtoms,nWalkers*3))
    ind=0
    globe=0
    globexyz = 0
    for line in linez:
        splitz = line.split()
        #If weight line.
        if len(splitz) == 5: #len(splitz) == 1
            wtAr[globe] = float(splitz[0])
            globe+=1
        elif len(splitz) == 4:
            splitz2 = [float(x) for x in splitz if x != 'D' and x !='H' and x!='O']
            #print ind, globexyz, filledCoords.size, splitz, filledCoords[ind,globexyz:globexyz+3]
            filledCoords[ind,globexyz:globexyz+3] = splitz2
            ind+=1
        elif line=='\n':
            ind=0
            globexyz += 3
    f.close()
    #filledCoords = filledCoords * angstr
    return nWalkers,filledCoords,wtAr


def swapTwoCoordinates(ind1,ind2,walkersToSwap):
    ind1-=1
    ind2-=1
    newWalkers = copy.deepcopy(walkersToSwap)
    for walk in newWalkers:
        hold = np.copy(walk.coords[ind1][:])
        walk.coords[ind1][:] = np.copy(walk.coords[ind2][:])
        walk.coords[ind2][:] = hold
    return newWalkers

def swapLobe(group1,group2,walkersToSwap):
    #Group1 = [6,5,1,13] , Group2 = [9,10,3,12] so a list of indices
    group1[:] = [o-1 for o in group1]
    group2[:] = [q-1 for q in group2]
    newWalkers = copy.deepcopy(walkersToSwap)
    for ind in range(len(group1)):
        for walk in newWalkers:
            hold = np.copy(walk.coords[group1[ind]][:])
            walk.coords[group1[ind]][:] = np.copy(walk.coords[group2[ind]][:])
            walk.coords[group2[ind]][:] = hold
    group1[:] = [o+1 for o in group1]
    group2[:] = [q+1 for q in group2]
    return newWalkers

def negZ(walkers):
    for set in walkers:
        for walk in set:
            #walk.coords[:,:] = np.negative(walk.coords[:,:])
            walk.coords[:,-1] = np.negative(walk.coords[:,-1])
    return walkers


def swapStuff():
    if 'llH' in config or 'llD' in config:
        # H5H6Only
        H6H5 = swapTwoCoordinates(5, 6,myWalkers)
        print "allH/allD",H6H5[0].coords
        # H7H8Only
        H8H7 = swapTwoCoordinates(7, 8,myWalkers)
        # H9H10Only
        H10H9 = swapTwoCoordinates(9, 10,myWalkers)

        #SwapTwoAtOnce
        #H6H5,H7H8
        H6H5H8H7 = swapTwoCoordinates(7,8,H6H5)
        #H6H5,H9H10
        H6H5H10H9 = swapTwoCoordinates(9,10,H6H5)
        #H7H8,H9H10
        H8H7H10H9 = swapTwoCoordinates(9,10,H8H7)
        #AllSwapped
        H6H5H8H7H10H9 = swapTwoCoordinates(9,10,H6H5H8H7)
        #original = myWalkers

        big8 = [myWalkers,H6H5,H8H7,H10H9,H6H5H8H7,H6H5H10H9,H8H7H10H9,H6H5H8H7H10H9]
        new8_56_910 = copy.deepcopy(big8)
        new8_56_78 = copy.deepcopy(big8)
        new8_78_910 = copy.deepcopy(big8)

        #SwapLobes
        #Just 56<-->910 Chunk
        for swappedBoy in range(len(big8)):
            new8_56_910[swappedBoy] = swapLobe([5,6,1,13],[10,9,3,12],big8[swappedBoy])
        # Just 56<-->78 Chunk
        for swappedBoy in range(len(big8)):
            new8_56_78[swappedBoy] = swapLobe([5, 6, 1, 13], [7, 8, 2, 11], big8[swappedBoy])
        # Just 78<-->910 Chunk
        for swappedBoy in range(len(big8)):
            new8_78_910[swappedBoy] = swapLobe([7, 8, 2, 11], [10,9, 3, 12], big8[swappedBoy])

        # SwapLobes
        new8_56_910_65with78 = copy.deepcopy(new8_56_910)
        new8_56_78_65with910 = copy.deepcopy(new8_56_78)
        #new8_78_910_78with65 = copy.deepcopy(new8_78_910) redundant.
        #After first flip, taked those flipped and flip them once more to get a unique combination that was not accessible through original or one flip.
        for swappedBoy in range(len(big8)):
            new8_56_910_65with78[swappedBoy] = swapLobe([5,6,1,13], [7, 8, 2, 11], new8_56_910[swappedBoy])
        # Just 56<-->78 Chunk
        for swappedBoy in range(len(big8)):
            new8_56_78_65with910[swappedBoy] = swapLobe([5,6,1,13], [10,9,3,12], new8_56_78[swappedBoy])
        # Just 78<-->910 Chunk
        #for swappedBoy in range(len(big8)):
        #    new8_78_910_78with65[swappedBoy] = swapLobe([7, 8, 2, 11], [5,6,1,13], new8_78_910[swappedBoy])


        return big8+new8_56_910+new8_56_78+new8_78_910+new8_56_910_65with78+new8_56_78_65with910+\
               negZ(big8) + negZ(new8_56_910) + negZ(new8_56_78) + negZ(new8_78_910) + negZ(new8_56_910_65with78) + negZ(new8_56_78_65with910)

    elif '1He' in config or '1De' in config:
        # H5H6Only
        H6H5 = swapTwoCoordinates(5, 6,myWalkers)
        print "1he/1de",H6H5[0].coords
        # H7H8Only
        H8H7 = swapTwoCoordinates(7, 8,myWalkers)
        # H9H10Only
        H10H9 = swapTwoCoordinates(9, 10,myWalkers)

        #SwapTwoAtOnce
        #H6H5,H7H8
        H6H5H8H7 = swapTwoCoordinates(7,8,H6H5)
        #H6H5,H9H10
        H6H5H10H9 = swapTwoCoordinates(9,10,H6H5)
        #H7H8,H9H10
        H8H7H10H9 = swapTwoCoordinates(9,10,H8H7)
        #AllSwapped
        H6H5H8H7H10H9 = swapTwoCoordinates(9,10,H6H5H8H7)
        #original = myWalkers

        big8 = [myWalkers,H6H5,H8H7,H10H9,H6H5H8H7,H6H5H10H9,H8H7H10H9,H6H5H8H7H10H9]
        new8_56_910 = copy.deepcopy(big8)
        #new8_56_78 = copy.deepcopy(big8)
        #new8_78_910 = copy.deepcopy(big8)
        #SwapLobes
        #Just 56<-->910 Chunk
        for swappedBoy in range(len(big8)):
            new8_56_910[swappedBoy] = swapLobe([5,6,1,13],[10,9,3,12],big8[swappedBoy])
        return big8+new8_56_910+negZ(big8)+negZ(new8_56_910)

    elif '1Hw' in config or '1Dw' in config:
        # H7H8Only
        H8H7 = swapTwoCoordinates(7, 8,myWalkers)
        # H9H10Only
        H10H9 = swapTwoCoordinates(9, 10,myWalkers)
        #H7H8,H9H10
        H8H7H10H9 = swapTwoCoordinates(9,10,H8H7)
        big4=[myWalkers,H8H7,H10H9,H8H7H10H9]
        #new4_56_910 = copy.deepcopy(big4)
        #new4_56_78 = copy.deepcopy(big4)
        new4_78_910 = copy.deepcopy(big4)
        # Just 78<-->910 Chunk
        for swappedBoy in range(len(big4)):
            new4_78_910[swappedBoy] = swapLobe([7, 8, 2, 11], [10,9, 3, 12], big4[swappedBoy])
        return big4+new4_78_910


for config in fileList:
    print config
    numWalkers, xyzCoords, weightArray = putWalkersInArray2(config)
    #np.save("weights"+sys.argv[1],weightArray)
    #tpo
    #np.save("wtAr_100.npy",weightArray)
    #stop
    myWalkers = [Walker() for i in range(numWalkers)]
    giveWalkersTheirCoordinates(xyzCoords, myWalkers)
    giveWalkersTheirWeight(weightArray)
    print myWalkers[0].coords
    SymWalkers = swapStuff()
    # name = 'Sym'+config[-10:]
    name = '../coordinates/tetramer/RSwapped_'
    writeTheWalkers(SymWalkers, name)
    print "Done"