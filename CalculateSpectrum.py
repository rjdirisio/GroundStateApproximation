import numpy as np
import os
import usefulFunctions as use
import time
import scipy.linalg as sla
import numpy.linalg as la
import itertools
import gc
import sys
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

ProtonatedWaterTrimer = {'H7O3+','O3H7+', 'H7O3plus','H7O3', 'O3H7'}
ProtonatedWaterTetramer = {'H9O4+','O4H9+', 'H9O4plus','H9O4', 'O4H9'}
global ProtonatedWaterTrimer
au2wn=219474.63
au2ang=0.529177249
massConversionFactor=1.000000000000000000/6.02213670000e23/9.10938970000e-28
ang2bohr=1.88973
bohr2ang=1.000/ang2bohr
global massH
global massO
global massD
massH=1.00782503223
massD=2.0141017778
massO=15.99491561957
massO*=massConversionFactor


verbose=False



class HarmonicApproxSpectrum(object):
    def __init__(self,wfn,coords,dw,path,testName): #CHANGE THIS
        self.wfn=wfn
        self.coords=coords
        self.dw=dw
        self.nVibs=self.wfn.molecule.nVibs
        self.path=path
        self.testname = testName

    def LoadG(self,GfileName,calcAllGs=False):
        print 'does ', GfileName, 'exist?'
        if calcAllGs:
            print 'Doesnt matter! Calculate ALL GS!!!!!! >>:)'
            allGs, gnm = self.calculateG_all(self.coords, self.dw)  # tested and gives same results
            gspl = GfileName.split("/")
            walkSet, _ = gspl[-1].split(".")
            np.save("allGs/allGM" + walkSet + ".npy", allGs)
            return gnm
        else:
            if not os.path.isfile(GfileName):
                print 'no!'
                gnm = self.calculateG(self.coords,self.dw)
                if 'test' not in GfileName:
                    np.savetxt(GfileName,gnm)
            else:
                print 'yes! ^-^'
                gnm=np.loadtxt(GfileName)
            return gnm

    def calculateG(self,eckartRotatedCoords,descendantWeights):
        # dx = 1e-6
        dx = 1e-4
        # dx = 1e-1
        start = time.time()
        mass = self.wfn.molecule.get_mass()
        print 'Start calculating G'
        # threwOut=0
        print 'summing up the descendants', np.sum(descendantWeights)
        sumDescendants = np.sum(descendantWeights)
        gnm = np.zeros((self.nVibs,self.nVibs))
        # mwpartialderv_all = np.zeros((len(descendantWeights), self.nVibs, self.nVibs))
        for atom in range(self.wfn.molecule.nAtoms):
            for coordinate in range(3):
                #WHERE INTERNAL COORDINATES ARE USED
                print 'dx number',atom*3+(coordinate+1), 'atom:',atom, 'coordinate',coordinate
                deltax=np.zeros((eckartRotatedCoords.shape))
                # deltax[:,1,2]+=dx #perturbs the x,y,z coordinate of the atom of interest - testing
                # deltax[:,5,1]+=dx #perturbs the x,y,z coordinate of the atom of interest - testing
                deltax[:,atom,coordinate]=deltax[:,atom,coordinate]+dx #perturbs the x,y,z coordinate of the atom of interest

                coordPlus=self.wfn.molecule.SymInternals(eckartRotatedCoords+deltax,False)
                coordMinus=self.wfn.molecule.SymInternals(eckartRotatedCoords-deltax,False)
                # coordPlus=np.where(coordPlus == np.isnan(coordPlus),coordPlus,0.)
                # coordMinus = np.where(coordMinus== np.isnan(coordMinus), coordMinus, 0.)
                pnan=np.where(np.isnan(coordPlus))
                mnan=np.where(np.isnan(coordMinus))
                # print len(pnan[0])/3., 'nans in coordPlus'
                # print len(mnan[0])/3., 'nans in coordMinus'
                # if len(pnan[0]) > 0 or len(mnan[0]) > 0:
                #     print 'setting nans to zero'
                #     descendantWeights[pnan[0]]=0.0
                #     descendantWeights[mnan[0]] = 0.0
                #     coordPlus[pnan]=0.0
                #     coordMinus[mnan]=0.0
                bigIdx = np.where(np.abs(coordPlus-coordMinus) > 1.)
                #Add 2pi
                if self.wfn.molecule.name in ProtonatedWaterTrimer:
                    badIdx = np.where(bigIdx[1] == 6)
                    # bigIdx = np.where(bigIdx[1]!=6)
                else:
                    badIdx = np.where(bigIdx[1] == 9)
                    # fixAdjustedEckartStuffInTetramer

                #For 2pi stuff going from 0 to 2pi
                # qual1 = np.abs(coordPlus - coordMinus) > 1.0
                # print np.sum(qual1),'walkers with weird derivatives'
                qual1 = np.where(np.abs(coordPlus - coordMinus) > 1.0)
                # if np.any(qual1):
                #     qual2 = np.logical_or(qual1[1]==10 , qual1[1] == 13)
                # else:
                #     qual2 = []
                #if difference is negative, add pi/2.  Otherwise,
                # coordPlus[qual1[0][qual2],qual1[1][qual2]] += (np.pi/2.0) * np.sign(coordMinus[qual1[0][qual2],qual1[1][qual2]] - coordPlus[qual1[0][qual2],qual1[1][qual2]])
                coordPlus[np.abs(coordPlus - coordMinus) > 1.0] += (-1.0 * 2. * np.pi) * np.sign(coordPlus[np.abs(coordPlus - coordMinus) > 1.0])
                # descendantWeights[(np.abs(coordPlus-coordMinus) > 1.)[:,0]]=0.0

                partialderv=(coordPlus-coordMinus)/(2.0*dx) #Discretizing stuff - derivative with respect to our perturbation
                excessCount = len(bigIdx[0])
                # np.savetxt("pderv_"+str(atom*3+(coordinate+1)),partialderv[:,[3,6]])
                ######Commented Out
                print excessCount, 'walkers have bad stuff (other than theta)'
                if len(bigIdx[0])!=0:
                    for badWalk in np.column_stack(bigIdx):
                        print badWalk[0],badWalk[1], partialderv[badWalk[0],badWalk[1]],coordPlus[badWalk[0],badWalk[1]],coordMinus[badWalk[0],badWalk[1]]
                if len(badIdx)!=0:
                    print 'Theta H was bad in ',len(badIdx[0]), 'walkers'

                if len(coordPlus[np.abs(coordPlus-coordMinus) > 1.])>0 :
                    badfl = open("badFile.xyz","w+")
                    self.wfn.molecule.printCoordsToFile(eckartRotatedCoords[np.where(np.abs(coordPlus-coordMinus) > 1.0)[0]],badfl)
                    add2pifuckyouuuu
                ###################


                for i,pd in enumerate(partialderv):
                    mwpd2 = (partialderv[i,:,np.newaxis]*partialderv[i,np.newaxis,:])/mass[atom]
                    gnm+=mwpd2*descendantWeights[i]
        gnm/=sumDescendants
        return gnm

    def calculateG_all(self,eckartRotatedCoords,descendantWeights):
        #RYAN - THE COORDINATES ARE NOT ACTUALLY ECKART ROTATED
        #Input is x which is a NAtoms x 3(coordinates) sized array
        #input is also dx, the perturbation size, usually .001                                                
        #output is the G matrix, which is a self.nVibs*self.nVibs sized array (there are self.nVibs internals)
        #(3N-6 x 3N-6)
        print 'CALCULATE ALL G MATRICES'
        dx=1e-4
        start=time.time()
        mass=self.wfn.molecule.get_mass()
        print 'Start calculating G'
        # threwOut=0
        print 'summing up the descendants', np.sum(descendantWeights)
        sumDescendants=np.sum(descendantWeights)
        mwpartialderv_all = np.zeros((len(descendantWeights),self.nVibs,self.nVibs))
        for atom in range(self.wfn.molecule.nAtoms):
            for coordinate in range(3):
                #WHERE INTERNAL COORDINATES ARE USED
                cycleTime = time.time()
                print 'dx number',atom*3+(coordinate+1), 'atom:',atom, 'coordinate',coordinate
                deltax=np.zeros((eckartRotatedCoords.shape))
                # deltax[:,1,2]+=dx #perturbs the x,y,z coordinate of the atom of interest
                deltax[:,atom,coordinate]=deltax[:,atom,coordinate]+dx #perturbs the x,y,z coordinate of the atom of interest
                coordPlus=self.wfn.molecule.SymInternals(eckartRotatedCoords+deltax,False)
                coordMinus=self.wfn.molecule.SymInternals(eckartRotatedCoords-deltax,False)
                # coordPlus=np.where(coordPlus == np.isnan(coordPlus),coordPlus,0.)
                # coordMinus = np.where(coordMinus== np.isnan(coordMinus), coordMinus, 0.)
                pnan=np.where(np.isnan(coordPlus))
                mnan=np.where(np.isnan(coordMinus))
                print len(pnan[0])/3., 'nans in coordPlus'
                print len(mnan[0])/3., 'nans in coordMinus'
                if len(pnan[0]) > 0 or len(mnan[0]) > 0:
                    print 'setting nans to zero'
                descendantWeights[pnan[0]]=0.0
                descendantWeights[mnan[0]] = 0.0
                coordPlus[pnan]=0.0
                coordMinus[mnan]=0.0
                bigIdx = np.where(np.abs(coordPlus-coordMinus) > 1.)
                #Add 2pi
                if self.wfn.molecule.name in ProtonatedWaterTrimer:
                    badIdx = np.where(bigIdx[1] == 6)
                    # bigIdx = np.where(bigIdx[1]!=6)
                else:
                    badIdx = np.where(bigIdx[1] == 9)
                    # fixAdjustedEckartStuffInTetramer

                # qual = (np.abs(coordPlus - coordMinus) > 1.0) * (np.abs(coordPlus-coordMinus) < (np.pi+0.5))
                # print len(coordPlus[qual])
                # coordPlus[qual] += (-1.0 * np.pi) * np.sign(coordPlus[qual])
                # coordPlus[np.abs(coordPlus-coordMinus) > (np.pi+0.5)]+=(-1.0*2.*np.pi)*np.sign(coordPlus[np.abs(coordPlus-coordMinus) > (np.pi+0.5)])
                # print len(coordPlus[np.abs(coordPlus-coordMinus) > 1.0]), 'potential theta'

                #For 2pi stuff going from 0 to 2pi
                # qual1 = np.abs(coordPlus - coordMinus) > 1.0
                # print np.sum(qual1),'walkers with weird derivatives'
                coordPlus[np.abs(coordPlus - coordMinus) > 1.0] += (-1.0 * 2. * np.pi) * np.sign(coordPlus[np.abs(coordPlus - coordMinus) > 1.0])
                # descendantWeights[(np.abs(coordPlus-coordMinus) > 1.)[:,0]]=0.0

                partialderv=(coordPlus-coordMinus)/(2.0*dx) #Discretizing stuff - derivative with respect to our perturbation
                excessCount = len(bigIdx[0])
                print excessCount, 'walkers have bad stuff (other than theta)'
                if len(bigIdx[0])!=0:
                    for badWalk in np.column_stack(bigIdx):
                        print badWalk[0],badWalk[1], partialderv[badWalk[0],badWalk[1]],coordPlus[badWalk[0],badWalk[1]],coordMinus[badWalk[0],badWalk[1]]
                if len(badIdx)!=0:
                    print 'Theta H was bad in ',len(badIdx[0]), 'walkers'
                if len(coordPlus[np.abs(coordPlus-coordMinus) > 1.])>0 :
                    badfl = open("badFile.xyz","w+")
                    self.wfn.molecule.printCoordsToFile(eckartRotatedCoords[np.where(np.abs(coordPlus-coordMinus) > 1.0)[0]],badfl)
                    add2pifuckyouuuu
                mwpd2 = (partialderv[:,:,np.newaxis]*partialderv[:,np.newaxis,:])/mass[atom]
                mwpartialderv_all += mwpd2
                print 'dx timing: ', str(time.time()-cycleTime), 'secs'
        print 'timing for G matrix', time.time()-start
        print 'dividing by ',sumDescendants

        gmat = np.average(mwpartialderv_all, axis=0, weights=descendantWeights)
        return mwpartialderv_all,gmat


    def diagonalizeRootG(self,G):
        w,v=np.linalg.eigh(G)
        #print w,v
        #octomom
        #vinv=np.linalg.inv(v)
        invRootDiagG=np.diag(1.0/np.sqrt(w))  #1.0/rootDiagG  #MAYBE THIS IS OK??
        # for i,ValinvRootDiagG in enumerate(invRootDiagG):#range(self.nVibs):
        #     if verbose: print 'alp[',i,']',ValinvRootDiagG[i]
        invRootG=np.dot(v,np.dot(invRootDiagG,v.transpose()))
        #invG=np.dot(invRootG,invRootG.transpose())
        #checkG=np.linalg.inv(invG)
        #self.GHalfInv=invRootG
        return  invRootG

    def calculateSecondMoments(self,x,dw,setOfWalkers,kil,ek):
        #calculate average internals
        print 'Calculating moments . . .'
        internals=self.wfn.molecule.SymInternals(x)
        averageInternals=np.average(internals,weights=dw,axis=0)
        np.savetxt('averageInternals_'+setOfWalkers+'_'+ek+'_'+kil,averageInternals)
        #calculate moments
        print 'averageInternals: ',zip(self.wfn.molecule.internalName,averageInternals)
        moments=internals-averageInternals
        np.save("moments_"+setOfWalkers+'_'+ek+'_'+kil+".npy",moments)
        print 'Moments: ', moments
        #calculate second moments
        walkerSize =len(dw)
        print 'calculate second moments'
        # dmap = np.memmap('dmap'+setOfWalkers+'_'+ek+'_'+kil, dtype='float64', mode='w+', shape=(int(walkerSize)))
        # dmap[:] = dw[:]
        # dmap.flush()
        # del dw
        # del dmap
        # gc.collect()
        sh = np.memmap('smap'+setOfWalkers+'_'+ek+'_'+kil, dtype='float64', mode='w+', shape=(int(walkerSize), int(self.nVibs),int(self.nVibs)))
        print 'walkerSize',walkerSize
        if walkerSize > 1000000:
            chunkSize=500000
        else:
            chunkSize=100
        cycle = np.arange(0,walkerSize,chunkSize)
        last=False
        a=time.time()
        for stt in cycle:
            print time.time()-a
            print stt
            k=stt+chunkSize
            if k > walkerSize:
                last=True
            if not last:
                sh[stt:k,:,:]=moments[stt:k,:,np.newaxis]*moments[stt:k,np.newaxis,:]
            else:
                sh[stt:,:,:]=moments[stt:,:,np.newaxis]*moments[stt:,np.newaxis,:]
        sh.flush()
        del sh
        #del mh
        print 'returning moments...'
        return moments

    def overlapMatrix(self,q,dw,poE,walkerSet,kil,ek):
        nwalkers = len(dw)
        # lg = open("overlapLog.txt","w+")
        potE = np.copy(poE)*au2wn
        #Csontruct constants and on diagonal elements
        a = np.average(q * q * q / np.average(q * q, axis=0, weights=dw), axis=0, weights=dw) * (1 / np.average(q * q, axis=0, weights=dw))
        b = -1 / np.average(q * q, axis=0, weights=dw)
        bq2aq1 = 1 + a * q + b * q * q
        #Construct Overlap Matrix
        #Construct diagonal elements
        nvibs2 = self.nVibs * 2
        engageKineticCoupling=False
        if (not os.path.isfile('ezOvMat_'+walkerSet+"noKEInFunds"+'_'+ek+'_'+kil+".npy")) or engageKineticCoupling:
            print 'ez part of overlap matrix doesnt exist'
            dgnl2=[]
            dgnl2.append(1)
            #dgnl.append(1) #ground state with itself
            for y in range(self.nVibs*2):
                if y < self.nVibs: #Fundamentals
                    #dgnl.append(np.average(q2[:,y],weights=dw))
                    dgnl2.append(np.average((q*q)[:,y],weights=dw))
                else: #Overtones
                    a1sqrt=1+a[y-self.nVibs]*(q[:,y-self.nVibs])+b[y-self.nVibs]*(q[:,y-self.nVibs])*(q[:,y-self.nVibs])
                    a1 =  np.average(a1sqrt*a1sqrt,weights=dw)
                    #dgnl.append(a1)
                    dgnl2.append(a1)
            del a1sqrt
            ########combos for overlap 2
            for (x,y) in itertools.combinations(q.T,2):
                a2p=np.average(x*x*y*y,weights=dw)
                dgnl2.append(a2p)
            #np.savetxt("dgnl2",dgnl2)

            #############Overlap matrix constructed######################
            overlap2=np.diag(dgnl2)
            ham2 = np.diag(np.zeros(len(dgnl2)))
            #fig = plt.figure()
            #plt.matshow(overlap2)
            #plt.colorbar()
            #plt.savefig("overlapBegin.png")
            #plt.close()
            ############Off diagonal Elements#############################
            #q
            # bq^2+aq+1 ASDF
            del a
            del b
            gc.collect()
            overlap2[0,1:self.nVibs+1]=np.average(q,axis=0,weights=dw)
            ham2[0, 1:self.nVibs + 1] = np.average(q*potE[:,np.newaxis], axis=0, weights=dw)
            overlap2[0,self.nVibs+1:self.nVibs*2+1] = np.average(bq2aq1, axis=0, weights=dw)
            ham2[0, self.nVibs + 1:self.nVibs * 2 + 1] = np.average(bq2aq1*potE[:,np.newaxis], axis=0, weights=dw)
            #combinations
            nvibs2 = self.nVibs*2
            vbo=np.flip(np.arange(self.nVibs+1))
            curI = nvibs2+1
            nxtI = curI+self.nVibs-1
            for combo in range(self.nVibs):
                overlap2[0,curI:nxtI]=np.average(q[:,combo,np.newaxis]*q[:,(combo+1):],axis=0,weights=dw)
                ham2[0, curI:nxtI] = np.average(q[:, combo, np.newaxis]*q[:, (combo + 1):]*potE[:,np.newaxis], axis=0, weights=dw)
                curI = nxtI
                nxtI += vbo[combo+1]-1
            ##Fundamentals with other funds - already calculated
            pst=np.triu_indices_from(overlap2[1:self.nVibs+1,1:self.nVibs+1],k=1)
            af=pst[0]+1
            bf=pst[1]+1
            overlap2[tuple((af,bf))] = np.copy(overlap2[0,nvibs2+1:])
            # ham2[tuple((af,bf))] = np.copy(ham2[0,nvibs2+1:])
            if engageKineticCoupling:

                ########KINETIC COUPLING#######
                gmats = np.load("allGs/allGM" + walkerSet + '_' + ek + '_' + kil+".npy")
                ghinv = np.load("allGs/ghinv_" + walkerSet+ '_' + ek + '_' + kil+".npy")
                vecc = np.load("allGs/vecs_T_L" + walkerSet+ '_' + ek + '_' + kil+".npy")
                # tmat = np.loadtxt("TransformationMatrix" + walkerSet +'_'+ek+'_'+kil+ ".datatest")
                gmatz = np.matmul(vecc.T, np.matmul(ghinv, np.matmul(gmats, np.matmul(ghinv, vecc))))
                # gmatz2 = np.matmul(la.inv(vecc),np.matmul(ghinv,np.matmul(gmats,np.matmul(ghinv,vecc))))
                plt.matshow(gmats[0])
                plt.colorbar()
                plt.savefig("gmats_orig")
                plt.close()
                plt.matshow(gmatz[0])
                plt.colorbar()
                plt.savefig("gmatz")
                plt.close()
                # plt.matshow(gmatzTest[0])
                # plt.colorbar()
                # plt.savefig("gmatzTest")
                # plt.close()
                for combo in range(self.nVibs-1):
                    print combo
                    # <2,0|gab|0,2>
                    # <0,0|gab|0,0>
                    # <1,1|gab|1,1>
                    # <2,0|gab|0,0>
                    # <0,0|gab|0,2>
                    # <1,0| V |0,1>
                    if combo == 1:
                        aa1=np.average(bq2aq1[:,combo,np.newaxis]*gmatz[:,combo,(combo+1):]*bq2aq1[:, (combo + 1):],axis=0,weights=dw)
                        aa2=np.average(gmatz[:,combo,(combo+1):],axis=0,weights=dw)
                        aa3=np.average((q[:, combo, np.newaxis])**2 * (q[:, (combo + 1):])**2 *gmatz[:,combo,(combo+1):],axis=0,weights=dw)
                        aa4=np.average(bq2aq1[:, combo, np.newaxis] * gmatz[:, combo, (combo + 1):], axis=0, weights=dw)
                        aa5=np.average(gmatz[:, combo, (combo + 1):] * bq2aq1[:, (combo + 1):], axis=0, weights=dw)
                        aa6=np.average(q[:, combo, np.newaxis] * q[:, (combo + 1):] * potE[:, np.newaxis]/au2wn, axis=0, weights=dw)
                        print '<2,0|gab|0,2>',aa1*au2wn
                        print '<0,0|gab|0,0>',aa2*au2wn
                        print '<1,1|gab|1,1>',aa3*au2wn
                        print '<2,0|gab|0,0> (subtract this #)',aa4*au2wn
                        print '<0,0|gab|0,2> (subtract this #)',aa5*au2wn
                        print '1/4(2,2 + 0,0 + 1,1 - 2,0 - 0,2)',0.25*(aa1+aa2+aa3-aa4-aa5)*au2wn
                        print 'potential energy avg',aa6*au2wn
                        print 'total',(0.25*(aa1+aa2+aa3-aa4-aa5)+aa6)*au2wn
                        del aa1,aa2,aa3,aa4,aa5,aa6
                    ham2[combo + 1, combo + 2:self.nVibs + 1] = \
                        0.25*(np.average(bq2aq1[:,combo,np.newaxis]*gmatz[:,combo,(combo+1):]*bq2aq1[:, (combo + 1):],axis=0,weights=dw)\
                        +np.average(gmatz[:,combo,(combo+1):],axis=0,weights=dw)\
                        +np.average((q[:, combo, np.newaxis])**2 * (q[:, (combo + 1):])**2 *gmatz[:,combo,(combo+1):],axis=0,weights=dw)\
                        -np.average(bq2aq1[:,combo,np.newaxis]*gmatz[:,combo,(combo+1):],axis=0, weights=dw)\
                        -np.average(gmatz[:,combo,(combo+1):]*bq2aq1[:,(combo+1):], axis=0,weights=dw))\
                        +np.average(q[:, combo, np.newaxis] * q[:, (combo + 1):] * potE[:, np.newaxis] /au2wn, axis=0, weights=dw) #pot E is in wavenumbers, KE is in a.u.
                    ham2[combo + 1, combo + 2:self.nVibs + 1]*=au2wn #change all to wavenumbers
                walkerSet = walkerSet + 'KEInFunds'
                ###############
                if not os.path.isfile('ezOvMat_' + walkerSet + "noKEInFunds" + '_' + ek + '_' + kil + ".npy"):
                    print 'Have spectrum calc continue through'
                    walkerSet = walkerSet + "_KE"
                else:

                    dov = np.diagonal(overlap2)
                    overlap2 = overlap2 / np.sqrt(dov[:, None])
                    overlap2 = overlap2 / np.sqrt(dov)
                    ham2 = ham2 / np.sqrt(dov[:, None])
                    ham2 = ham2 / np.sqrt(dov)
                    overlap2 = overlap2 + overlap2.T - np.eye(len(overlap2))
                    print 'Hammie and Overlap Man are Constructed. Nearing the end...'
                    print 'abs overlap Max Min', np.amax(np.abs(overlap2)), np.amin(np.abs(overlap2))
                    ham2 = ham2 + ham2.T  # nothing on diagonal
                    np.savetxt(
                        'JUSTFUNDSWITHKEINCLUDEDoffDiagonalCouplingsInPotential2_' + walkerSet + ek + kil + '.dat',
                        ham2)
                    print 'finished.'
                    youarefinished
            else:
                walkerSet = walkerSet + 'noKEInFunds'
                ham2[tuple((af,bf))] = np.copy(ham2[0,nvibs2+1:])

            ##Funds with overtones
            ##MEMMAPS
            # overtones with other overtones
            print "Ovs with Ovs"
            for combo in range(self.nVibs):
                print combo
                #asdf=np.average(bq2aq1[:, combo, np.newaxis] * bq2aq1[:, (combo + 1):], axis=0, weights=dw)
                overlap2[self.nVibs + combo + 1, self.nVibs + combo + 2:nvibs2 + 1] = np.average(bq2aq1[:, combo, np.newaxis] * bq2aq1[:, (combo + 1):], axis=0, weights=dw)
                ham2[self.nVibs + combo + 1, self.nVibs + combo + 2:nvibs2 + 1] = np.average(bq2aq1[:, combo, np.newaxis]*bq2aq1[:, (combo + 1):]*potE[:,np.newaxis], axis=0, weights=dw)
            print 'done. overs with themselves'

            #Funds with Overtones
            hav = np.zeros((self.nVibs,self.nVibs))
            hamhav = np.zeros((self.nVibs,self.nVibs))
            for ind in range(self.nVibs):
                hav[ind,:]=np.average(q[:,ind,np.newaxis]*bq2aq1,weights=dw,axis=0)
                hamhav[ind,:]=np.average(q[:,ind,np.newaxis]*bq2aq1*potE[:,np.newaxis],weights=dw,axis=0)
            overlap2[1:self.nVibs + 1, self.nVibs + 1:nvibs2 + 1] = hav
            ham2[1:self.nVibs + 1, self.nVibs + 1:nvibs2 + 1] = hamhav
            if 'test' not in walkerSet:
                np.save("ezOvMat_"+walkerSet+'_'+ek+'_'+kil+".npy",overlap2)
                np.save("ezHamMat_"+walkerSet+'_'+ek+'_'+kil+".npy",ham2)
        else:
            print 'ez overlap and ham already exist! loading...'
            overlap2 = np.load('ezOvMat_'+walkerSet+'noKEInFunds'+'_'+ek+'_'+kil+".npy")
            ham2 = np.load('ezHamMat_'+walkerSet+'noKEInFunds'+'_'+ek+'_'+kil+".npy")
        # funds with combos and overtones with combos
        ####PUT BACK HERE

        #HEFTY BOI
        bigMem =False
        if bigMem:
            print 'bigMemActivated'
            lnsize = int((self.nVibs * self.nVibs - self.nVibs) / 2.)
            if not os.path.isfile('lm_' + walkerSet+'_'+ek+'_'+kil + ".npy"):
                print 'no lm , calculating...'
                lm = np.zeros((len(q), lnsize))
                for combo in range(self.nVibs):
                    if combo == 0:
                        prev = 0
                    print prev
                    print prev + self.nVibs - combo - 1
                    lm[:, prev:(prev + self.nVibs - combo - 1)] = q[:, combo, np.newaxis] * q[:, (combo + 1):]
                    prev += self.nVibs - 1 - combo
                np.save('lm_' + walkerSet+'_'+ek+'_'+kil + ".npy",lm)
            else:
                print 'lm exists, laoding...'
                lm=np.load('lm_' + walkerSet+'_'+ek+'_'+kil + ".npy")

                #########################################
            sumDw = np.sum(dw)
            if self.wfn.molecule.name in ProtonatedWaterTrimer:
                if 'allH' in walkerSet:
                    splitArs = 96
                else:
                    splitArs = len(q)  # smallest, least memory, fastest?
            else:
                splitArs = len(q)
            qsize = q.shape[1]
            print 'splitting arrays'

            if len(dw) % splitArs != 0.0:
                print 'not divisible~!!!'
                octopus

            # Reshaping for chopping instead of array splits
            q = np.reshape(q, (splitArs, nwalkers / splitArs, -1))
            bq2aq1 = np.reshape(bq2aq1, (splitArs, nwalkers / splitArs, -1))
            potE = np.reshape(potE, (splitArs, nwalkers / splitArs))
            dw = np.reshape(dw, (splitArs, nwalkers / splitArs))
            lm = np.reshape(lm, (splitArs, nwalkers / splitArs, -1))

            # q = np.array_split(q, splitArs)
            # bq2aq1 = np.array_split(bq2aq1, splitArs)
            # potE = np.array_split(potE, splitArs)
            # dw = np.array_split(dw, splitArs)
            # lm = np.array_split(lm, splitArs)

            cyc = 0
            fuco = np.zeros((qsize, lnsize))
            ovco = np.zeros((qsize, lnsize))
            hfuco = np.zeros((qsize, lnsize))
            hovco = np.zeros((qsize, lnsize))
            coco = np.zeros((lnsize, lnsize))
            hcoco = np.zeros((lnsize, lnsize))

            for qq, bb, pp, dd, ll in itertools.izip(q, bq2aq1, potE, dw, lm):

                print 'cycle', cyc, 'out of ', splitArs
                if cyc == 0:
                    st = time.time()
                fuco += np.sum(qq[:, :, np.newaxis] * ll[:, np.newaxis, :] * dd[:, np.newaxis, np.newaxis], axis=0)
                ovco += np.sum(bb[:, :, np.newaxis] * ll[:, np.newaxis, :] * dd[:, np.newaxis, np.newaxis], axis=0)
                hfuco += np.sum(
                    qq[:, :, np.newaxis] * ll[:, np.newaxis, :] * pp[:, np.newaxis, np.newaxis] * dd[:, np.newaxis,np.newaxis],axis=0)
                hovco += np.sum(
                    bb[:, :, np.newaxis] * ll[:, np.newaxis, :] * pp[:, np.newaxis, np.newaxis] * dd[:, np.newaxis,np.newaxis],axis=0)
                coco += np.sum(ll[:, :, np.newaxis] * ll[:, np.newaxis, :] * dd[:, np.newaxis, np.newaxis], axis=0)
                hcoco += np.sum(ll[:, :, np.newaxis] * ll[:, np.newaxis, :] * pp[:, np.newaxis, np.newaxis] * dd[:, np.newaxis,np.newaxis],axis=0)
                cyc += 1
                print time.time()-st , 'secs'
            ovlas = np.triu_indices_from(overlap2[nvibs2 + 1:, nvibs2 + 1:], k=1)
            g = ovlas[0] + nvibs2 + 1
            h = ovlas[1] + nvibs2 + 1
            asco = np.triu_indices_from(np.zeros((lnsize, lnsize)), k=1)
            overlap2[tuple((g, h))] = coco[asco] / sumDw
            overlap2[1:self.nVibs + 1, self.nVibs * 2 + 1:] = fuco / sumDw  # FC
            overlap2[self.nVibs + 1:2 * self.nVibs + 1, self.nVibs * 2 + 1:] = ovco / sumDw
            ham2[tuple((g, h))] = hcoco[asco] / sumDw
            ham2[1:self.nVibs + 1, self.nVibs * 2 + 1:] = hfuco / sumDw
            ham2[self.nVibs + 1:2 * self.nVibs + 1, self.nVibs * 2 + 1:] = hovco / sumDw
            #########################################
        else:
            print 'smol Mem Activated'
            lst = []
            for o, t in itertools.combinations(np.arange(self.nVibs), 2):
                lst.append((o, t))
            ct = 0
            print len(lst)
            print 'funds and overs with combos'
            for (i, j) in lst:
                print ct
                qqiqj = q * q[:, i, np.newaxis] * q[:, j, np.newaxis]
                bqq = bq2aq1 * q[:, i, np.newaxis] * q[:, j, np.newaxis]
                overlap2[1:self.nVibs + 1, self.nVibs * 2 + 1 + ct] = np.average(qqiqj, weights=dw, axis=0)
                overlap2[self.nVibs + 1:2 * self.nVibs + 1, self.nVibs * 2 + 1 + ct] = np.average(bqq, weights=dw, axis=0)
                ham2[1:self.nVibs + 1, self.nVibs * 2 + 1 + ct] = np.average(qqiqj* potE[:, np.newaxis], weights=dw, axis=0)
                ham2[self.nVibs + 1:2 * self.nVibs + 1, self.nVibs * 2 + 1 + ct] = np.average(bqq* potE[:, np.newaxis], weights=dw, axis=0)
                ct += 1
            print 'hello'
            print 'combos with combos'
            # combos with combos
            ovlas = np.triu_indices_from(overlap2[nvibs2 + 1:, nvibs2 + 1:], k=1)
            g = ovlas[0]
            h = ovlas[1]
            g += nvibs2 + 1
            h += nvibs2 + 1
            ct = 0
            for [(w, x), (y, z)] in itertools.combinations(lst, 2):
                if ct % 1000 == 0:
                    print(ct)
                qwxyz = q[:, w] * q[:, x] * q[:, y] * q[:, z]
                overlap2[tuple((g[ct], h[ct]))] = np.average(qwxyz, weights=dw)
                ham2[tuple((g[ct], h[ct]))] = np.average(qwxyz * potE, weights=dw)
                ct += 1
        dov = np.diagonal(overlap2)
        np.save("overlap2NOTTRANSPOSED_notNormed" + walkerSet + ek + kil + '.dat', overlap2)
        np.save("ham2NOTTRANSPOSED_notNormed" + walkerSet + ek + kil + '.dat', ham2)
        overlap2 = overlap2 / np.sqrt(dov[:, None])
        overlap2 = overlap2 / np.sqrt(dov)
        ham2 = ham2 / np.sqrt(dov[:, None])
        ham2 = ham2 / np.sqrt(dov)
        overlap2 = overlap2 + overlap2.T - np.eye(len(overlap2))
        print 'Hammie and Overlap Man are Constructed. Nearing the end...'
        print 'abs overlap Max Min', np.amax(np.abs(overlap2)), np.amin(np.abs(overlap2))
        ham2 = ham2 + ham2.T  # nothing on diagonal
        return ham2,overlap2

    def calculateSpectrum(self, coords,dw,GfileName,pe,dips,setOfWalkers,testName,kill,ecked,diPath):
        sumDw = np.sum(dw)
        #Equations are referenced from McCoy, Diken, and Johnson. JPC A 2009,113,7346-7352
        #print 'I better be calculating the spectrum from ECKART ROTATED coordinates!'
        """Commence getting our modified descendent weights"""
        print 'Commence Calculate Spectrum Function'
        #dipoleMoments = np.zeros(np.shape(dips))

        # com, eckVecs = self.wfn.molecule.eckartRotate(coords)
        #print dips
        # if testName == 'oxEck':
        #     justO = True
        # else:
        #     justO = False
        if not ecked:
            if self.wfn.molecule.name in ProtonatedWaterTrimer:
                com, eckVecs, killList = self.wfn.molecule.eckartRotate(coords, planar=True,lst=[0,1,2],dip=True)
            else:
                com, eckVecs, killList = self.wfn.molecule.eckartRotate(coords, planar=True, lst=[0,1,2,3],dip=True)
            #
            # print 'getting eckarted dipole moments. . .'
            # if setOfWalkers == 'fSymtet_allD':
            #     com, eckVecs, killList = self.wfn.molecule.eckartRotate(coords, justO=False,yz=True)
            # else:
            #     com, eckVecs, killList = self.wfn.molecule.eckartRotate(coords, justO=True)
            dips = dips - com  # added this to shift dipole to center of mass before eckart roatation - translation of dipole should NOT matter
            #print 'killList = '+str(len(killList))+' Walkers out of ' + str(len(dw))
            dipoleMoments = np.zeros(np.shape(dips))
            print 'dips shifted to COM: ', dips[0]
            for den in range(dips.shape[0]):
                dipoleMoments[den] = np.dot(dips[den],eckVecs[den]) #because eckVecs is .T?
                #dipoleMoments2[b] = np.dot(eckVecs[b],dips[b]) #This definitely isn't right
            # dipz=np.dot(dips,eckVecs)
            print 'dipole moment - after eckart rotation: ', dipoleMoments[0]
            np.save(diPath+'eng_dip_'+setOfWalkers+'_eckart.npy',np.column_stack((pe,dipoleMoments)))
            del dips
        else:
            dipoleMoments = dips
            print 'dips loaded from eckart',dipoleMoments[0]
            del dips
        # stoptopus
        #dips = dips*ang2bohr
        #First, What is the G Matrix for this set of walkers based on the SymInternals coordinates
        calcAllGs = False
        if not os.path.isfile('q_'+setOfWalkers+'_'+testName+'_'+kill+'.npy'):
            self.G=self.LoadG(GfileName,calcAllGs=calcAllGs)
            self.wfn.molecule.setInternalName()
            if not os.path.isfile("smap"+setOfWalkers+'_'+testName+'_'+kill):
                moments=self.calculateSecondMoments(coords, dw,setOfWalkers,kill,testName)
            else:
                print "moments_"+setOfWalkers
                if not os.path.isfile("moments_"+setOfWalkers+'_'+testName+'_'+kill+".npy"):
                    if not os.path.isfile("mu2ave_"+setOfWalkers+'_'+testName+'_'+kill+".npy"):
                        print 'no moments...calculating'
                        internals = self.wfn.molecule.SymInternals(coords)
                        averageInternals = np.average(internals, weights=dw, axis=0)
                        momentsmoments = internals - averageInternals
                        if 'test' not in setOfWalkers:
                            np.save("moments_"+setOfWalkers+'_'+testName+'_'+kill+".npy",moments)
                    else:
                        moments = np.zeros((self.nVibs,len(coords)))
                else:
                    print 'moments already calculated'
                    moments = np.load("moments_"+setOfWalkers+'_'+testName+'_'+kill+".npy")
            q,q2=self.calculateQCoordinates(moments,dw,GfileName,setOfWalkers,kill,testName)
            print 'done with normal modes'
            #q4ave=np.average(q4,axis=0,weights=dw)
            q2ave=np.average(q2,axis=0,weights=dw)
            qave =np.average(q,axis=0,weights=dw)
            #print q.shape
            print '/\/\/\/\/\/\/\/\/\/\ '
            print 'some averages',
            print 'q\n',qave
            print 'q^2 \n',q2ave #average of q-squared
            #print 'q^4 \n',q4ave
            print '/\/\/\/\/\/\/\/\/\/\ '
            if 'test' not in setOfWalkers and 'top' not in setOfWalkers:
                np.save("q_"+setOfWalkers+'_'+testName+'_'+kill+".npy",q)
                np.save("q2_"+setOfWalkers+'_'+testName+'_'+kill+".npy",q2)
        else:
            if calcAllGs:
                self.G = self.LoadG(GfileName,calcAllGs=calcAllGs)
            print 'loading qs from file'
            q=np.load('q_'+setOfWalkers+'_'+testName+'_'+kill+'.npy')
            q2=np.load('q2_'+setOfWalkers+'_'+testName+'_'+kill+'.npy')
            q2ave = np.average(q2, axis=0, weights=dw)
        #Now calculate the Potential energy
        print 'calculating PE'
        potentialEnergy=self.calculatePotentialEnergy(coords,pe)
        print 'Potential Energy', potentialEnergy
        overlapTime=True
        if overlapTime:
            ham2,overlap2=self.overlapMatrix(q,dw,potentialEnergy,setOfWalkers,kill,testName)
            overlapMs = self.path + 'redH/'
            np.savetxt(overlapMs + 'overlapMatrix2_' + setOfWalkers + testName + kill + '.dat', overlap2)
            np.savetxt(overlapMs + 'offDiagonalCouplingsInPotential2_' + setOfWalkers + testName + kill + '.dat', ham2)
        else:
            overlapMs = self.path + 'redH/'
        #V_0=<0|V|0>
        print 'overlap matrix done or skipped'
        print 'Potential Energy',potentialEnergy
        V_0=np.average(potentialEnergy[:,None],axis=0,weights=dw)

        # Vq= <1_q|V|1_q> bra:state with one quanta in mode q, ket: state with one quanta in mode q
        Vq=np.average(potentialEnergy[:,None]*q2,axis=0,weights=dw)

        #Vq2d=<q^2Vq^2> (<> is an average of the descendant weights)
        Vq2d=np.zeros((self.nVibs,self.nVibs)) # one quanta in one mode and 1 quanta in another mode

        #q2ave2d=<q^2 q^2> (<> is an average of the descendant weights)  
        q2ave2d=np.zeros((self.nVibs,self.nVibs))

        for i in range(coords.shape[0]):
            Vq2d=Vq2d+np.outer(q2[i],potentialEnergy[i]*q2[i])*dw[i]
            q2ave2d=q2ave2d+np.outer(q2[i],q2[i])*dw[i]
        Vq2d=Vq2d/np.sum(dw)
        q2ave2d=q2ave2d/np.sum(dw)        
        
        #Now calculate the kinetic energy
        print 'ZPE: average v_0',V_0*au2wn
        np.save(overlapMs+"v_0"+setOfWalkers + testName + kill,V_0*au2wn)
        print 'Vq', Vq*au2wn

        q4 = q*q*q*q
        q8 = q4*q4
        alpha=q2ave/(np.average(q4,weights=dw,axis=0)-q2ave**2) # Equation #11
        #alphaPrime=0.5/q2ave   #Equation in text after #8

        alphaPrime = np.sqrt(8. * np.average(q4,weights=dw,axis=0) / (np.average(q8,weights=dw,axis=0) - np.average(q4,weights=dw,axis=0) ** 2))

        #        print 'how similar are these?', zip(alpha,alphaPrime) Still a mystery to me why there were 2 dfns of alpha
        #kineticEnergy= hbar**2 nquanta alpha/(2 mass)
        Tq=1.0**2*1.0*alpha/(2.0*1.0) #Equation #10

        TqP = 1.0**2*2.0*alphaPrime/(2.0*1.0) #N=2 !!!!!!!!!!

        Tq2d=np.zeros((self.nVibs,self.nVibs)) #Tij=Ti+Tj
        #Finish calculate the potential and kinetic energy for the combination and overtone bands 
        #put Equation 8 into equation 4, algebra...
        for ivibmode in range(self.nVibs): #combination band calculations + overtones (i=j)
            for jvibmode in range(ivibmode):
                Vq2d[ivibmode,jvibmode]=Vq2d[ivibmode,jvibmode]/q2ave2d[ivibmode,jvibmode]
                Vq2d[jvibmode,ivibmode]=Vq2d[jvibmode,ivibmode]/q2ave2d[jvibmode,ivibmode]
                Tq2d[ivibmode,jvibmode]=Tq[ivibmode]+Tq[jvibmode]
                Tq2d[jvibmode,ivibmode]=Tq2d[ivibmode,jvibmode]
            #Vq2d[ivibmode,ivibmode]=Vq2d[ivibmode,ivibmode]*4.0*alpha[ivibmode]**2-Vq[ivibmode]*4.0*alpha[ivibmode]+V_0
            #Vq2d[ivibmode,ivibmode]=Vq2d[ivibmode,ivibmode]/(4.0*alpha[ivibmode]**2*q2ave2d[ivibmode,ivibmode]-4.0*alpha[ivibmode]*q2ave[ivibmode]+1.0)
            #Vq2d[ivibmode, ivibmode] = Vq2d[ivibmode, ivibmode] * 4.0 * alpha[ivibmode] ** 2 - Vq[ivibmode] * 4.0 *alpha[ivibmode] + V_0
            Tq2d[ivibmode,ivibmode]=Tq[ivibmode]*2.0

        T2P = np.copy(Tq2d)
        np.fill_diagonal(T2P,TqP)

        #Let's get overtones set up
        a = np.average(q*q*q / np.average(q*q,axis=0,weights=dw),axis=0,weights=dw) * (1/np.average(q*q,axis=0,weights=dw))
        b = -1/np.average(q*q,axis=0,weights=dw)
        vovers = np.average((1+a*q+b*q*q)*potentialEnergy[:,None]*(1+a*q+b*q*q),axis=0,weights=dw)
        vovers /= (np.average(np.square(1+a*q+b*q*q),axis=0,weights=dw))
        np.fill_diagonal(Vq2d,vovers)
        # print 'overtone combo potential',Vq2d
        #Energy is V+T, don't forget to subtract off the ZPE
            #Eq2d=np.zeros((self.nVibs,self.nVibs))
        Eq2d=(Vq2d+Tq2d)-V_0 #related to frequencies later on
        Eq2dP = (Vq2d + T2P) - V_0  # related to frequencies later on

        Eq=(Vq/q2ave+Tq)-V_0
        #noKE_funds = Vq/q2ave - V_0
        #noKE_combos = Vq2d-V_0
        #Lastly, the dipole moments to calculate the intensities! Ryan - Do you even have to eckart the molecule, why not just eckart rotate the dipole moment?

        #dipoleMoments=self.wfn.molecule.dipole(coords)
        #notRotateddipoleMoments = dips
        #dips =self.wfn.molecule.rotateDipoleToFrame(coords, dips) #NEED TO HAVE THIS SO THAT MY ECKART ROTATION OF THE COORDINATES (WHICH ARE IN XY PLANE) MATCH THE ECKART ROTATION OF THE DIPOLE

        #dips =self.wfn.molecule.loadDips(dips) #get rotated dipoles



        #stop
        #print 'getting eckarted dipole moments. . .'
        #dipoleMoments = self.wfn.molecule.rotateDipoleToFrame(coords,dips)
        #print dipoleMoments

        #print 'Dipole Moments!: ',dipoleMoments
        averageDipoleMoment_0=np.average(dipoleMoments,axis=0,weights=dw)
        print 'Average Dipole Moment!: ',averageDipoleMoment_0
        #aMu_i=<1_i|dipoleMoment|0>=Sum_n^walkers(Dipole_n*q_i*dw_n)/Sum(dw) #Equation 4 for \psi_n=\psi_1
        aMu=np.zeros((3,self.nVibs))  
        #aMu2d_{i,j}=<1_i,1_j|dipoleMoment|0>=Sum_n^walkers(Dipole_n*q_i*q_j*dw_n)/Sum(dw) #Equation 4 for \psi_n=\psi_{1,1}
        aMu2d=np.zeros((3,self.nVibs,self.nVibs))
        print 'dip:'
        print dipoleMoments.shape
        #q = nWalkers x 24
        #dw = nWalkers
        #dipoleMoments = nWalkers x 3
        #aMu = dipole (x) -Q coord - 24x3
        for P,R,d in zip(dipoleMoments,q,dw): #P = dipole moment (nx3), R = q coordinate (3n-6x1), d = descendent wight (n)
            # wtf = np.outer(P,-R)
            # wtf2 = wtf*d
            aMu=aMu+np.outer(P,-R)*d  #off by a sign change for some reason so I changed it to be -R instead of R...I guess P could be backwards?
            #aMu=aMu+np.outer(P,R)*d  #Ryan - this doesn't change the results, as of right now (7/19/18)
            temp=np.outer(R,R)        #aMu += Dipole.coordinate * DQ
            temp=temp.reshape((self.nVibs*self.nVibs,1))
            aMu2d=aMu2d+(np.outer(P,temp).reshape((3,self.nVibs,self.nVibs)))*d
        #aMu = 24x3 dipole contribution in that coordinate?

        # test = (-q[:, np.newaxis, :] * dipoleMoments[:, :, np.newaxis] * dw[:, np.newaxis, np.newaxis]).sum(axis=0)
        # qout = q[:,:,np.newaxis]*q[:,np.newaxis,:]
        # test2 = dipoleMoments[:,np.newaxis,:,np.newaxis]*qout[:,:,np.newaxis,:]*dw[:,np.newaxis,np.newaxis,np.newaxis]
        # test2 = np.sum(test2,axis=0).transpose(1,0,2)
        # print test.shape
        # test/=np.sum(dw)
        # test2/=np.sum(dw)

        aMu=aMu/np.sum(dw)
        aMu2d=aMu2d/np.sum(dw)






        # testI = np.sum(np.square(test)/q2ave,axis=0)
        # testI2 = np.square(la.norm(test.T,axis=1))/q2ave
        # testI22 = la.norm(np.square(test/q2ave),axis=0)
        #np.savetxt("WtvsMu-qs-a.data",zip(dw,np.sum(dipoleMoments**2,axis=1),q[:,5],q[:,6]))

        aMu=aMu.transpose() #3x24 instead of 24x3
        aMu2d=aMu2d.transpose()  #switches around the axes 24x24x3

        magAvgMu=np.zeros(self.nVibs)
        magAvgMu2=np.zeros(self.nVibs)

        #Finding |<1_i|dipoleMoment|0>|^2
        """Getting intensity!"""
        #print 'aMu: ', aMu.shape,aMu
        #np.savetxt('aMuDividedBySqrt_q2',aMu/np.sqrt(q2ave))
        for m,mode in enumerate(aMu): #eq 4
            magAvgMu[m]=magAvgMu[m]+(mode[0]**2+mode[1]**2+mode[2]**2)/(q2ave[m])  #     += mu_x^2+mu_y^2+mu_z^2 /<q^2>
        print magAvgMu

        magMu2d=np.zeros((self.nVibs,self.nVibs))
        #Finding |<1_i,1_j|dipoleMoment|0>|^2
        for ivibmode in range(self.nVibs):
            #aMu2d[ivibmode,ivibmode]=2*alpha[ivibmode]*aMu2d[ivibmode,ivibmode]-averageDipoleMoment_0
            #magMu2d[ivibmode,ivibmode]=np.sum(aMu2d[ivibmode,ivibmode]**2)
            #magMu2d[ivibmode,ivibmode]=magMu2d[ivibmode,ivibmode]/(q2ave2d[ivibmode,ivibmode]*4.0*alpha[ivibmode]**2-4.0*alpha[ivibmode]*q2ave[ivibmode]+1.0)
            #magMu2d[ivibmode, ivibmode] = magMu2d[ivibmode, ivibmode]/(q2ave2d[ivibmode, ivibmode]*4.0*alpha[ivibmode]**2-4.0*alpha[ivibmode]*q2ave[ivibmode]+1.0)
            for jvibmode in range(ivibmode):
                #print ivibmode,jvibmode
                magMu2d[ivibmode,jvibmode]=np.sum(aMu2d[ivibmode,jvibmode]**2)
                magMu2d[ivibmode,jvibmode]=magMu2d[ivibmode,jvibmode]/q2ave2d[ivibmode,jvibmode]
                magMu2d[jvibmode,ivibmode]=magMu2d[ivibmode,jvibmode]

        #First: get values I want, which are not the magnitude squared, just the:
        #Fundamentals: <i=1 | u | 0>
        mu10 = aMu/ np.sqrt(q2ave)[:,None]
        overlapMs = self.path+'redH/'
        np.savetxt(overlapMs+'1mu0' + setOfWalkers + testName + kill,mu10)

        #Overtone intensities are currently incorrect, let's do it from scratch
        overpoly = (b*q2+a*q+1)
        rMu20 = np.zeros((3,self.nVibs))
        for walker in range(len(dipoleMoments)):
            rMu20 += np.outer(dipoleMoments[walker],overpoly[walker])*dw[walker]
        rMu20 /= np.sum(dw)
        rMu20 /= np.sqrt(np.average(np.square(b*q2+a*q+1),weights=dw,axis=0))
        rMu20=rMu20.T
        magRMu20 = np.sum(np.square(rMu20),axis=1)
        np.savetxt(overlapMs+'2mu0' + setOfWalkers + testName + kill,rMu20)
        np.fill_diagonal(magMu2d,magRMu20)

        #Combinations: <i=1,j=1|u|0>
        tupz = []
        for i in range(self.nVibs):
            for j in range(i):
                tupz.append([i,j])
        fMu = np.zeros((len(tupz),3))
        fMup= np.zeros((len(tupz),3))
        ct=0
        avg2pairs = np.zeros(len(tupz))
        avg2pairsp = np.zeros(len(tupz))

        for glomp in tupz:
            fMu[ct,:] = aMu2d[glomp[0],glomp[1]]
            avg2pairs[ct] = np.sqrt(q2ave2d[glomp[0],glomp[1]])
            ct+=1
        ct=0
        for (x,y) in itertools.combinations(np.arange(self.nVibs),2):
            fMup[ct,:] = aMu2d[x,y]
            avg2pairsp[ct] = np.sqrt(q2ave2d[x,y])
            ct+=1
        mu110 = fMu/avg2pairs[:,None]
        mu110p = fMup/avg2pairsp[:,None]

        np.savetxt(overlapMs+'11mu0'+ setOfWalkers + testName + kill,mu110)
        np.savetxt(overlapMs+'11mu0p'+ setOfWalkers + testName + kill,mu110p)
        # comoz=np.square(la.norm(mu110p,axis=1))
        # intFromAnalysis = []
        # for (x,y) in itertools.combinations(np.arange(self.nVibs),2):
        #         #print Eq2d[x,y]*au2wn , magMu2d[x,y],'combination bands' , x,y
        #         intFromAnalysis.append(magMu2d[x,y])
        # #asdf=magMu2d[np.triu_indices(magAvgMu.shape[0],k=1)]
        # ar = np.array(intFromAnalysis)
        # print np.around(ar,5) == np.around(comoz,5)
        print 'done'
        
        #######################################################################################################

        print "m alpha   alpha'=0.5/<q2>        Vq            Vq/q2ave      Tq            Eq          |<1|dipole|0>|"
        for i in range(9):
            print i, alpha[i], alphaPrime[i],au2wn*Vq[i], au2wn*Vq[i], Vq[i]/q2ave[i]*au2wn, Eq[i]*au2wn, magAvgMu[i]

        print 'energies!', zip(Vq2d[0]*au2wn,Tq2d[0]*au2wn,Eq2d[0]*au2wn)

        print 'Spectrum Info'
        #fundamentalFileName=self.path+'fundamentalsJUSTPOTENTIAL'+testName+'.data'
        comboFileName=self.path+'combinationOverrtone_'+setOfWalkers+testName+kill+'.data'
        #fundamentalFile=open(fundamentalFileName,'w')
        comboFile=open(comboFileName,'w')
        comboFile2=open(comboFileName+'2','w')
        #justVComboFile = open(comboFileName+"_justV",'w')
        """magAvgMu = Intensity for fundamentals
           magMu2d  = Intensity for combination and overtone bands"""
        for (x,y) in itertools.combinations(np.arange(self.nVibs),2):
                #print Eq2d[x,y]*au2wn , magMu2d[x,y],'combination bands' , x,y
                comboFile2.write(str( Eq2d[x,y]*au2wn)+"   "+str(magMu2d[x,y])+"       "+str(x)+" "+str(y)+"\n")

        for i in range(self.nVibs):
            print Eq[i]*au2wn, magAvgMu[i], 'mode ',i
            #fundamentalFile.write(str(noKE_funds[i]*au2wn)+"   "+str(magAvgMu[i])+"       "+str(i)+"\n")
        for i in range(self.nVibs):
            for j in range(i):
                print Eq2d[i,j]*au2wn , magMu2d[i,j],'combination bands' , i,j
                comboFile.write(str( Eq2d[i,j]*au2wn)+"   "+str(magMu2d[i,j])+"       "+str(i)+" "+str(j)+"\n")
                #justVComboFile.write(str(noKE_combos[i,j]*au2wn)+"   "+str(magMu2d[i,j])+"       "+str(i)+" "+str(j)+"\n")
        for i in range(self.nVibs):
            print Eq2d[i,i]*au2wn , magMu2d[i,i],'overtone bands' , i,i
            comboFile.write(str( Eq2d[i,i]*au2wn)+"   "+str(magMu2d[i,i])+"       "+str(i)+" "+str(i)+"\n")
            #comboFile2.write(str( Eq2d[i,i]*au2wn)+"   "+str(magMu2d[i,i])+"       "+str(i)+" "+str(i)+"\n")
            comboFile2.write(str(Eq2dP[i, i] * au2wn) + "   " + str(magMu2d[i, i]) + "       " + str(i) + " " + str(i)+"\n")
        comboFile.close()
        comboFile2.close()
        return Eq*au2wn,magAvgMu, Eq2dP*au2wn,magMu2d


    def calculatePotentialEnergy(self,coords,pe):
        equilibriumEnergy=self.wfn.molecule.getEquilibriumEnergy()
        electronicEnergy = pe
        relativePotentialEnergy=electronicEnergy-equilibriumEnergy
        #print 'average relative potential energy', np.average(relativePotentialEnergy)
        #print 'descendant weighted average relative potential energy', np.average(relativePotentialEnergy,weights=dw)
        return relativePotentialEnergy

    def calculateQCoordinates(self,moments, dw,gmf,setOfWalkers,kil,eckt):
        print 'calculating Normal coordinates'
        walkerSize=len(dw)
        if not os.path.isfile("mu2ave_" + setOfWalkers + '_' + eckt + '_' + kil + ".npy"):
            print 'no mu2ave, calculating...'
            sh = np.memmap('smap' + setOfWalkers + '_' + eckt + '_' + kil, dtype='float64', mode='r',
                           shape=(int(walkerSize), int(self.nVibs), int(self.nVibs)))
            if 'tet' in setOfWalkers:
                # dh = np.memmap('dmap' + setOfWalkers + '_' + eckt + '_' + kil, dtype='float64', mode='r',
                #                shape=int(walkerSize))
                if walkerSize > 1000000:
                    chunkSize = 1000000
                else:
                    chunkSize = 10
                cycle = np.arange(0, walkerSize, chunkSize)
                last = False
                mu2AveR = np.zeros((self.nVibs, self.nVibs))
                for j in cycle:
                    k = j + chunkSize
                    if k > walkerSize:
                        last = True
                    if not last:
                        mu2AveR += np.sum(sh[j:k] * dw[j:k, np.newaxis, np.newaxis], axis=0)
                    else:
                        mu2AveR += np.sum(sh[j:] * dw[j:, np.newaxis, np.newaxis], axis=0)
                mu2Ave = np.divide(mu2AveR, np.sum(dw))
            else:
                mu2Ave = np.average(sh, axis=0, weights=dw)
            np.save('mu2ave_' + setOfWalkers + '_' + eckt + '_' + kil + '.npy', mu2Ave)
        else:
            print 'loading second moments matrix'
            mu2Ave=np.load("mu2ave_"+setOfWalkers+'_'+eckt+'_'+kil+'.npy')
        # mu2Ave_test = np.load("mu2ave_Sfinal_allH_rn_spc_finalVersion_addPiToChi.npy")

        # testing=False
        #
        # if testing:
        #     # self.wfn.molecule.internalName=['xH11', 'yH11', 'zH11', 'xH12', 'yH12', 'zH12', 'xH13', 'yH13', 'zH13', 'theta651',
        #     #                  'phi651', 'Chi651',
        #     #                  'theta1039', 'phi1039', 'Chi1039', 'theta728', 'phi728', 'Chi728', 'rOH5', 'rOH6',
        #     #                  'HOH516', 'rOH7', 'rOH8', 'HOH728',
        #     #                  'rOH9', 'rOH10', 'HOH9310', 'rO1O2','rO1O3','rO2O3', 'xO4', 'yO4', 'zO4']
        #     #####testing####
        #     # print 'eulers'
        #     # a = np.copy(mu2Ave[:9,:9])
        #     # b = np.copy(mu2Ave[:9, 18:])
        #     # c = np.copy(mu2Ave[18:,:9])
        #     # d = np.copy(mu2Ave[18:,18:])
        #     # ab=np.hstack((a,b))
        #     # cd=np.hstack((c,d))
        #     # mu2Ave =np.vstack((ab,cd))
        #     #
        #     #
        #     # a = np.copy(self.G[:9, :9])
        #     # b = np.copy(self.G[:9, 18:])
        #     # c = np.copy(self.G[18:, :9])
        #     # d = np.copy(self.G[18:, 18:])
        #     # ab = np.hstack((a, b))
        #     # cd = np.hstack((c, d))
        #     # self.G = np.vstack((ab, cd))
        #     # self.wfn.molecule.internalName = np.concatenate(
        #     #     (self.wfn.molecule.internalName[:9], self.wfn.molecule.internalName[18:]))
        #     print 'hi'
        #     print 'xyzO4 omitted'
        #     # mu2Ave=mu2Ave[:-3,:-3]
        #     # self.G=self.G[:-3,:-3]
        #     # self.wfn.molecule.internalName =self.wfn.molecule.internalName[:-3]
        #
        #     #just X & Y omitted
        #     # print 'x&y'
        #     # a = np.copy(mu2Ave[:29,:29])
        #     # b = np.copy(mu2Ave[:29, 31:])
        #     # c = np.copy(mu2Ave[31:,:29])
        #     # d = np.copy(mu2Ave[31:,31:])
        #     # ab=np.hstack((a,b))
        #     # cd=np.hstack((c,d))
        #     # mu2Ave =np.vstack((ab,cd))
        #     #
        #     # a = np.copy(self.G[:29, :29])
        #     # b = np.copy(self.G[:29, 31:])
        #     # c = np.copy(self.G[31:, :29])
        #     # d = np.copy(self.G[31:, 31:])
        #     # ab = np.hstack((a, b))
        #     # cd = np.hstack((c, d))
        #     # self.G = np.vstack((ab, cd))
        #     # self.wfn.molecule.internalName = np.concatenate(
        #     #     (self.wfn.molecule.internalName[:29], self.wfn.molecule.internalName[31:]))
        #     # Just Z omitted
        #     # print 'just z'
        #     # mu2Ave=mu2Ave[:-1,:-1]
        #     # self.G=self.G[:-1,:-1]
        #     # self.wfn.molecule.internalName =self.wfn.molecule.internalName[:-1]
        #     print 'off diagonal couplingggs'
        #     self.G[:29,29:]=0.0
        #     self.G[29:,:29]=0.0
        ########Testing########################
        # Gomit = [1,2,4,5,7,8,9,10,11,12,13,14,17,20,23]
        # Gomit = [0,1,2,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
        # Gomit = [15 16 18 19]
        # Gomit = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,17,20,21,22,23]
        # PinternalName = ['rOH_41','rOH_51', 'rOH_26', 'rOH_27']

        # # """self.internalName = ['rOH8', 'thH8', 'phiH8', 'rOH9', 'thH9', 'phiH9', 'rOH10', 'thH10', 'phiH10',
        # #                          'th_627', 'phi_627', 'xi_627', 'th_514', 'phi_514', 'xi_514', 'rOH_41',
        # #                          'rOH_51', 'aHOH_451', 'rOH_26', 'rOH_27', 'aHOH_267', 'rOO_1', 'rOO_2', 'aOOO']"""
        # PinternalName = ['rOH9', 'rOH10']
        # # PinternalName = ['rOH8', 'rOH9',  'rOH10', 'rOH_41','rOH_51', 'rOH_26', 'rOH_27', 'rOO_1', 'rOO_2']
        #
        # PinternalName = ['rOH8', 'thH8', 'phiH8', 'rOH9', 'thH9', 'phiH9', 'rOH10', 'thH10', 'phiH10', 'rOH_41',
        #                          'rOH_51', 'aHOH_451', 'rOH_26', 'rOH_27', 'aHOH_267', 'rOO_1', 'rOO_2', 'aOOO']
        # PinternalName = ['rOH8', 'thH8', 'phiH8', 'rOH9', 'thH9', 'phiH9', 'rOH10', 'thH10', 'phiH10', 'rOH_41',
        #                 'rOH_51', 'aHOH_451', 'rOH_26', 'rOH_27', 'aHOH_267', 'rOO_1', 'rOO_2', 'aOOO']
        # Gomit = [9,10,11,12,13, 14]
        # Gomit = [10, 11, 13, 14]
        # Gomit = [9,12]
        # # # # # # Gomit = [1,  2,  4,  5,  7,  8]+[9, 10, 11, 12, 13, 14]
        # # # # # # Gomit2 = [10, 11, 13, 14]
        # # # # # # Ginclude = [0,1,2,3,4,5,6,7,8,15,16,17,18,19,20,21,22,23]
        # # # # # # Ginclude2 = [0,1,2,3,4,5,6,7,8,9,12,15,16,17,18,19,20,21,22,23]
        # self.G[Gomit,:]=0.0
        # self.G[:,Gomit] = 0.0
        # mu2Ave[Gomit, :] = 0.0
        # mu2Ave[:,Gomit] = 0.0
        # self.G = self.G[self.G!=0.0].reshape(self.nVibs-len(Gomit),self.nVibs-len(Gomit))
        # mu2Ave = mu2Ave[mu2Ave!=0.0].reshape(self.nVibs-len(Gomit),self.nVibs-len(Gomit))
        # badTest = len(Gomit)
        # PinternalName = self.wfn.molecule.internalName[:9]+self.wfn.molecule.internalName[15:]
        # PinternalName = self.wfn.molecule.internalName[:10]+[self.wfn.molecule.internalName[12]]+self.wfn.molecule.internalName[14:]
        # PinternalName = self.wfn.molecule.internalName[:9] + self.wfn.molecule.internalName[10:12] + self.wfn.molecule.internalName[13:]

        ########Testing############################
        # mu2Ave = np.round(mu2Ave,11)
        # plt.matshow(mu2Ave)
        # plt.colorbar()
        # plt.show()
        # self.G = np.round(self.G,11)
        ########Testing############################
        GHalfInv=self.diagonalizeRootG(self.G)
        mu2AvePrime=np.dot(GHalfInv,np.dot(mu2Ave,GHalfInv)) # mass weights g^-1/2 . sm . g^-1/2
        eigval,vects=np.linalg.eigh(mu2AvePrime)
        seigval,svects = sla.eigh(mu2Ave,self.G)
        # print self.G == self.G.T
        #eigvalt, vectt = sla.eigh(mu2Ave,self.G)
        # self.wfn.molecule.internalName = ['rOH11', 'rOH12', 'rOH13', 'umbrella', '2dihed', 'dihed-di', 'thH', 'phH', 'xiH', 'theta651',
        #                     'phi651', 'Chi651',
        #                     'theta1039', 'phi1039', 'Chi1039', 'theta728', 'phi728', 'Chi728', 'rOH5', 'rOH6',
        #                     'HOH516', 'rOH7', 'rOH8', 'HOH728',
        #                     'rOH9', 'rOH10', 'HOH9310', 'rO1O2','rO1O3','rO2O3', 'xo4', 'yo4', 'zo4']
        #print 'test'
        """testVect = np.dot(self.G,vects)
        testest = np.dot(np.dot(vects.conj(),self.G),vects)
        testest2 = np.dot(np.dot(vects, mu2Ave), vects)

        s1=np.sum(np.square(vects1),axis=1)
        s2=np.sum(np.square(vects),axis=1)
        vects3=vects/np.linalg.norm(vects,axis=1)[:,np.newaxis]
        s3 = np.sum(np.square(vects3),axis=1)"""

        print 'diagonalized <mu^2>'
        #for i in range(self.wfn.molecule.nVibs):
        #    print 'v[',i,',]:',eigval[i]
        #    print vects[:,i]

        #GAAH                                                              
        TransformationMatrix=np.dot(vects.transpose(),GHalfInv)
        # TransformationMatrix = svects.T
        TmatRy = np.dot(GHalfInv,vects)
        np.save("allGs/vecs_T_L"+setOfWalkers+'_'+eckt+'_'+kil,vects.T)
        np.save("allGs/ghinv_"+setOfWalkers+'_'+eckt+'_'+kil,GHalfInv)

        # oppTmat = np.dot(GHalfInv,vects)
        # whatHappens = np.dot()
        #print 'RESULTS'
        Tnorm=1.0*TransformationMatrix # NEED TO DEEP COPY :-0
        alpha=1.0*eigval #AS in alpha_j in equation 5 of the JPC A 2011 h5o2 dier paper                                          
        #save the transformation matrix for future reference                                                                      
        TMatFileName='TransformationMatrix'+setOfWalkers+'_'+eckt+'_'+kil+'.data'
        #stop
        #TMatFileName = 'allHTesting/spectra/TransformationMatrix.data'
        np.savetxt(TMatFileName+"test",TransformationMatrix)
        # gmf = gmf[30:-5]
        # print 'assignment file name', setOfWalkers+gmf
        # if 'Eck' in gmf:
        #     gmf=setOfWalkers
        # if not testing:
        #     assignF = open(self.path+'assignments_'+setOfWalkers+"_"+eckt+kil,'w+')
        # else:
        assignF = open(self.path+'assignments_'+setOfWalkers+"_"+eckt+kil,'w+')
        for i,(vec,q2) in enumerate(zip(TransformationMatrix,eigval)):
            assignF.write("%d\n" % i)
            if verbose: print '\n',i,':<q^2 >=',q2
            #        alpha[i]=alpha[i]/(firstterm-(alpha[i]*alpha[i]))                                                            
            if verbose: print 'vec    ',vec
            Tnorm[i]=np.abs(vec)/np.max(np.abs(vec))
            sortvecidx=np.argsort(abs(vec))[::-1]
            #print 'sortvecidx: '
            sortedvec=vec[sortvecidx]

            sortedNames=use.sortList(abs(vec),self.wfn.molecule.internalName)[::-1]
            # sortedNames=use.sortList(abs(vec),PinternalName)[::-1]

            print 'assignments\n',zip(sortedNames,sortvecidx,sortedvec)[0:9]
            # for hu in range(self.nVibs - badTest):
            for hu in range(self.nVibs):
                assignF.write("%s %d %5.12f " % (sortedNames[hu],sortvecidx[hu],sortedvec[hu]))
            assignF.write('\n \n')
        assignF.close()

        #stop

        #calculate q from the moments and the transformation matrix
        # q=[]
        # #print moments.shape
        # for s in moments:
        #     q.append(np.dot(TransformationMatrix,s))
        #
        # q=np.array(q)
        #
        # pairz = [(9,12),(10,13),(11,14)]
        # u=0
        # for internA in mu2Ave:
        #     print(u)
        #     for pair in pairz:
        #         print('pair')
        #         print(internA[pair[0]], internA[pair[1]])
        #         print(np.abs(internA[pair[0]]),np.abs(internA[pair[1]]))
        #         print(np.abs(internA[pair[0]])-np.abs(internA[pair[1]]))
        #
            # u+=1
        q = np.matmul(TransformationMatrix, moments.T).T
        q2=q**2
        # q4=q**4
        return q, q2  # ,q4
