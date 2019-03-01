import os, numpy as np
from numpy import linalg as la
import glob
import copy
import sys
print "../coordinates/trimer/"+sys.argv[1]
fileList = glob.glob("../coordinates/trimer/"+sys.argv[1])
numAtoms = 10
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
    for walkInstance in someCfww:
        for each in walkInstance:
            wfi.write('10\n')
            wfi.write('%5.12f 0 0 0 0\n' % each.weight)
            wfi.write('O %5.12f %5.12f %5.12f\n' % (each.coords[0][0],each.coords[0][1],each.coords[0][2]))
            wfi.write('O %5.12f %5.12f %5.12f\n' % (each.coords[1][0],each.coords[1][1],each.coords[1][2]))
            wfi.write('O %5.12f %5.12f %5.12f\n' % (each.coords[2][0],each.coords[2][1],each.coords[2][2]))
            wfi.write('H %5.12f %5.12f %5.12f\n' % (each.coords[3][0],each.coords[3][1],each.coords[3][2]))
            wfi.write('H %5.12f %5.12f %5.12f\n' % (each.coords[4][0],each.coords[4][1],each.coords[4][2]))
            wfi.write('H %5.12f %5.12f %5.12f\n' % (each.coords[5][0],each.coords[5][1],each.coords[5][2]))
            wfi.write('H %5.12f %5.12f %5.12f\n' % (each.coords[6][0],each.coords[6][1],each.coords[6][2]))
            wfi.write('H %5.12f %5.12f %5.12f\n' % (each.coords[7][0],each.coords[7][1],each.coords[7][2]))
            wfi.write('H %5.12f %5.12f %5.12f\n' % (each.coords[8][0],each.coords[8][1],each.coords[8][2]))
            wfi.write('H %5.12f %5.12f %5.12f\n' % (each.coords[9][0],each.coords[9][1],each.coords[9][2]))
            wfi.write('\n')
    wfi.close()


def putWalkersInArray2(cfg):
    f = open(cfg , "r")
    linez = f.readlines()
    nWalkers = len(linez)/13
    wtAr = np.zeros(nWalkers)
    filledCoords = np.zeros((numAtoms,nWalkers*3))
    ind=0
    globe=0
    globexyz = 0
    for line in linez:
        splitz = line.split()
        #If weight line.
        if len(splitz) == 1:
            wtAr[globe] = float(splitz[0])
            globe+=1
        elif len(splitz) == 4:
            splitz2 = [float(x) for x in splitz if x != 'D' and x !='H' and x!='O']
            #print ind, globexyz, filledCoords.size, splitz, filledCoords[ind,globexyz:globexyz+3]
            #print globe
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

def negZ(walkers):
    for set in walkers:
        for walk in set:
            walk.coords[:,-1] = np.negative(walk.coords[:,-1])
    return walkers


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
    return newWalkers

def swapStuff():
    if 'llH' in config or 'llD' in config or '1Hh' in config or '1Dh' in config:
        print 'Full symmetrize'
        # H4H5Only
        H4H5 = swapTwoCoordinates(4,5,myWalkers)
        # H6H7Only
        H6H7 = swapTwoCoordinates(6,7,myWalkers)
        #SwapTwoAtOnce
        H4H5H6H7 = swapTwoCoordinates(6,7,H4H5)

        big4 = [myWalkers,H4H5,H6H7,H4H5H6H7]

        new4_45_67 = copy.deepcopy(big4)

        # SwapLobes
        for swappedBoy in range(len(big4)):
            new4_45_67[swappedBoy] = swapLobe([5,4,1,9],[6,7,2,10],big4[swappedBoy])
        return big4+new4_45_67+negZ(copy.deepcopy(big4))+negZ(copy.deepcopy(new4_45_67))

    elif '1He' in config or '1De' in config:
        print "Eigen Symmetrized"
        # H4H5Only
        H4H5 = swapTwoCoordinates(4,5,myWalkers)
        # H6H7Only
        H6H7 = swapTwoCoordinates(6,7,myWalkers)
        #SwapTwoAtOnce
        H4H5H6H7 = swapTwoCoordinates(6,7,H4H5)

        big4 = [myWalkers,H4H5,H6H7,H4H5H6H7]
        return big4+negZ(copy.deepcopy(big4))

    elif '1Hw' in config or '1Dw' in config:
        # H4H5Only
        print "Free Water Symmetrized"
        H4H5 = swapTwoCoordinates(4, 5, myWalkers)
        big2 = [myWalkers, H4H5]
        return big2+negZ(copy.deepcopy(big2))


for config in fileList:
    print config
    numWalkers,xyzCoords, weightArray = putWalkersInArray2(config)
    myWalkers = [Walker() for i in range(numWalkers)]
    giveWalkersTheirCoordinates(xyzCoords,myWalkers)
    giveWalkersTheirWeight(weightArray)
    print myWalkers[0].coords
    SymWalkers = swapStuff()
    #name = 'Sym'+config[-10:]
    name = '../coordinates/trimer/Rotated_Symmetrized_'
    writeTheWalkers(SymWalkers,name)
    print "Done"
