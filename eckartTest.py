def eckartRotate(self, pos, justO=False, specialCond=False):  # pos coordinates = walkerCoords numwalkersxnumAtomsx3
    nMolecules = pos.shape[0]
    allEckVecs = np.zeros((nMolecules, 3, 3))
    if self.name in ProtonatedWaterTrimer:
        self.refPos = self.pullTrimerRefPos()
    else:
        self.refPos = self.pullTetramerRefPos()
    if len(pos.shape) < 3:
        pos = np.array([pos])
    Fvec = np.zeros((3, 3))
    Fvec2 = np.zeros((3, 3))
    newCoord = np.zeros(pos.shape)
    newCoordRD = np.zeros(pos.shape)
    rHC = np.zeros((nMolecules, 3))
    rHCprime = np.zeros((nMolecules, 3))
    killList = []
    # Center of Mass
    mass = self.get_mass()
    com = np.dot(mass, pos) / np.sum(mass)

    if justO:
        self.refPos = self.refPos[:3]
        mass = mass[:3]
        pos = pos[:, :3, :]
    # com_ref = np.dot(mass,self.refPos)/np.sum(mass)
    ###########TEST
    # self.refPos = self.refPos-com_ref This did nothing.
    ###########/TEST
    rs = 0.000
    rsprime = 0.000

    # First Translate:
    ShiftedMolecules = pos - com[:, np.newaxis, :]
    plus = 0
    elze = 0
    # Equation 3.1 in Eckart vectors, Eckart frames, and polyatomic molecules - James D. Louck and Harold W. Galbraith
    for moli, molecule in enumerate(ShiftedMolecules):
        Fvec = np.zeros((3, 3))
        for atom, massa, eckatom in zip(molecule, mass, self.refPos):
            Fvec = Fvec + massa * np.outer(eckatom, atom)
        # F from eqn 3.4b             - vectorsthat connect the dotz
        FF = np.dot(Fvec, Fvec.transpose())
        # Diagonalize FF
        sortEigValsF, sortEigVecF = np.linalg.eigh(FF)
        sortEigVecFT = -sortEigVecF.transpose()
        if specialCond:
            print 'special condition activated!'
            print 'eigenvals \n', sortEigValsF
            print 'vect \n', sortEigVecFT
        if len(np.where(sortEigValsF <= 0)[0]) != 0:
            # sortEigVecFT=np.abs(sortEigVecFT)
            sortEigValsF = np.abs(sortEigValsF)
            invRootDiagF = sortEigValsF
            for e, element in enumerate(sortEigValsF):
                if element > 0:
                    invRootDiagF[e] = 1.0 / np.sqrt(element)
        # Get the inverse sqrt of diagonalized(FF)
        else:
            invRootDiagF = 1.0 / np.sqrt(sortEigValsF)
        # F^{-1/2}
        invRootF = np.dot(invRootDiagF[np.newaxis, :] * -sortEigVecF, sortEigVecFT)
        eckVecs = np.dot(Fvec.transpose(), invRootF)
        # print np.trace(eckVecs)
        # if moli == 999:
        #    stopit
        allEckVecs[moli] = eckVecs
        # print 'eckvecs & molecule'
        # print eckVecs.shape
        # print molecule.shape
        # print 'trace', np.trace(eckVecs)
        if not justO:
            # newCoordRD[moli] = np.dot(eckVecs, molecule.T).T wrong!
            newCoord[moli] = np.dot(molecule, eckVecs)
            detEck = la.det(eckVecs)
            if np.around(detEck) == -1.:
                killList.append(moli)
            elif np.around(detEck) == 1.:
                plus += 1
            else:
                print detEck
                elze += 1
            """o = open('orig.xyz',"w")
            nc = open('nc.xyz',"w")
            ncRD =  open('ncRD.xyz',"w")
            print 'ref'
            self.printCoordsToFile(self.refPos,o)
            o.close()
            print 'new'
            self.printCoordsToFile(newCoord[0], nc)
            nc.close()
            print 'newRD'
            self.printCoordsToFile(newCoordRD[0], ncRD)
            ncRD.close()
            stop"""
            # print 'O1 x y z: ',pos[0, 0, :]
            # print 'O2 x y z: ',pos[0, 1, :]
            # x = pos[0, 1, :]-pos[0, 2, :]
            # y = pos[0, 0, :]-pos[0, 2, :]
            # print 'x',x
            # print 'y',y
            # zs = np.dot(np.cross(x,y),(0,0,1)) #x of O2 and y of O1
            # print 'cross and dot: ', zs #x of O2 and y of O1

            # print 'new O1 x y z: ', newCoord[0, 0, :]
            # print 'new O2 x y z: ', newCoord[0, 1, :]
            # x = newCoord[0, 1, :]-newCoord[0, 2, :]
            # y = newCoord[0, 0, :]-newCoord[0, 2, :]
            # print 'new x',x
            # print 'new y',y
            # newZs = np.dot(np.cross(x,y),(0,0,1)) #x of O2 and y of O1
            # print 'new cross and dot: ', newZs #x of O2 and y of O1
            # print 'trace',np.trace(eckVecs)

            # if np.trace(eckVecs) < 0 :
            #    trace
            # if np.sign(zs) > 0 and np.sign(newZs) > 0:
            #    nm=0
            #    #print 'good'
            # else:
            #    if np.sign(zs) != np.sign(newZs):
            #        print 'bad'
            #        stoop
            #    else:
            #        print 'weird'
            #        stoop

            # newCoord[moli]= np.dot(eckVecs,molecule) I believe this to be incorrect.
        # else:
        # if np.around(np.trace(eckVecs)) != 2.:
        #    print np.trace(eckVecs)
        #    stop
        # for n,oneCoord in enumerate(molecule):

        #    print 'oneCoord',oneCoord,oneCoord.shape
        #    print 'eckVecs',eckVecs,eckVecs.shape
        #    print 'dotpdt', np.dot(oneCoord,eckVecs)
        #    print 'dotpdt', np.dot(eckVecs,oneCoord)

        #    newCoordRD[moli,n]=np.dot(oneCoord,eckVecs)

        # newCoordRD[moli] = np.dot(eckVecs,molecule.T).T
        # print newCoord[moli]
        # print np.dot(eckVecs,molecule.T).T

        if len(np.where(np.isnan(newCoord[moli]))[0]) != 0:
            print 'whaaaaaT?! nan', np.where(np.isnan(newCoord[moli])), '\ncoords:\n', newCoord[moli]
            print '   molecule number:', moli, '\n   sortEigValsF: \n', sortEigValsF, '\n   molecule: \n', molecule,
            print '\n   eckVecs \n', eckVecs
            octopus

    # print 'REference Eckarted Coords Shifted by COM', self.refPos-com_ref
    # print 'Old Coord - Shifted', pos[0]-com[0]
    # print 'Lindsey New Coord', newCoord[0]
    # print 'Ryan New Coord', newCoordRD[0]
    # print 'L Difference: ', (pos[0]-com[0])-newCoord[0]
    # print 'R Difference: ', (pos[0]-com[0]) - newCoordRD[0]
    # stop

    print "Whew ! Done with eckart."
    print 'plus', plus
    print 'else', elze
    print 'first allEckVecs', allEckVecs[0]
    # ff = open('allHTesting/eckartRotatedMolecule',"w+")
    # elf.printCoordsToFile(newCoord,ff)
    # ff.close()
    # print 'recorded new eckart coordinates'
    if self.name in ProtonatedWaterTrimer or self.name in ProtonatedWaterTetramer:
        return com, allEckVecs, killList
    else:
        return newCoord