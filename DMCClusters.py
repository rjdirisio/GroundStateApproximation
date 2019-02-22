import sys
import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import molecularInfo
import usefulFunctions as use
import subprocess

au2wn=219474.63
au2ang=0.529177249
massConversionFactor=1.000000000000000000/6.02213670000e23/9.10938970000e-28#1822.88839 


class wavefunction (object):

    def __init__(self,moleculeName, nWalkers, potentialName='Ground State',dtau=10.0):
        self.molecule=molecularInfo.molecule(moleculeName)

        print 'initialized', nWalkers,' coordinates for ', self.molecule.name
        self.nAtoms=self.molecule.nAtoms
        self.nDim=3

        #based on the molecule find the coordinates
        self.x=np.zeros((nWalkers,self.nAtoms,self.nDim))
        self.x[:]=self.molecule.getInitialCoordinates()
        self.D=.5
        self.set_dtau(dtau=dtau)
        self.mass=self.molecule.get_mass()
        self.sigma_dx=(2.0000*self.D*self.dtau/self.mass)**0.5
        self.alpha=.50/self.dtau
        self.recrossing=False

    def setNodalSurface(self,surfaceName,side='Both'):
        self.surfaceName=surfaceName
        self.side=side # Right: x>node
        self.molecule.setNodalSurface(self.surfaceName,self.side)
        self.recrossing=True 
    def setIsotope(self,keyword):
        self.molecule.set_isotope(keyword)
        self.mass=self.molecule.get_mass()
        self.sigma_dx=(2.0000*self.D*self.dtau/self.mass)**0.5
        print 'isotope made! new masses: ', self.mass
        
    def set_dtau(self,dtau=10.0):
        self.dtau=dtau        
        print 'dtau is ', self.dtau, 'imaginary atomic time units'

    def exchange(self,x, listOfExchanges): #Symmetry ops. exch(wv)
        xprime=1.0*x
        for (exAtom1,exAtom2) in listOfExchanges:
            print 'exchanging atom #',exAtom1,'with atom #',exAtom2
            temp=1.0*xprime[:,exAtom1]
            xprime[:,exAtom1]=xprime[:,exAtom2]*1.0
            xprime[:,exAtom2]=temp
        return xprime

    def diffuse(self): #random displacement
        dx=np.zeros((self.currentPop,self.nAtoms,self.nDim))
        for atom in range(self.nAtoms):
            centerOfGaussian=0.000000
            dx[:,atom:atom+1,:]=np.random.normal(centerOfGaussian,self.sigma_dx[atom],(self.currentPop,1,self.nDim))

        return dx

    def exportCoords(self,x,fileName,nDesc): #xyz
        fileout=open(fileName,'w')
        for particle, note in zip(x,nDesc):
            fileout.write(str(self.nAtoms)+' \n'+str(note)+' \n')
            for atomName, atom in zip(self.molecule.names,particle):
                fileout.write(str(atomName)+'   '+str(au2ang*atom[0])+"   "+str(au2ang*atom[1])+"   "+str(au2ang*atom[2])+"\n")
            fileout.write("\n")
        fileout.close()
        return

    def getPot(self,x):
        self.exportCoords(x,'walkers_intTimeStep'+self.molecule.name)
        process = subprocess.Popen(["sub.batch"],stdin=subprocess.PIPE) #call batch file
        process.wait()
        potential = np.loadtxt('vvalue') #import vvalue file with new potential energy
        return potential


    def loadCoords(self,fileName):  #xyz coords -> coords and dWs
        print 'Loading: ', fileName
        print 'Commence loadCoords'
        print 'fileName = ',fileName
        print 'numberOfAtoms, ',self.nAtoms

        filein=open(fileName,'r')
        filedata=filein.readlines()
        nMolecules=len(filedata)/(self.nAtoms+3)
        repeatUnit=self.nAtoms+3
        coord=np.zeros((nMolecules,self.nAtoms,3))
        descCoord=np.zeros((nMolecules))
        #Ryan Edits
        #pe = np.zeros(nMolecules)
        #dip = np.zeros((nMolecules,3))
        for ln in range(nMolecules):
            #print ln
            coordstemp=[]
            if 'eckart' in fileName:
                print 'doing nothing . . .'
            #else:
                #pe[ln] = float(filedata[ln*repeatUnit+1].split()[1].replace(',',''))
                #dip[ln,:] = float(filedata[ln*repeatUnit+1].split()[2].replace(',','')),float(filedata[ln*repeatUnit+1].split()[3].replace(',','')),float(filedata[ln*repeatUnit+1].split()[4].replace(',',''))
            for an in range(self.nAtoms):
                #print an
                #if an == 9:
                #    print 'hi'
                coordstemp.append(float(filedata[ln*repeatUnit+2+an].split()[1]))
                coordstemp.append(float(filedata[ln*repeatUnit+2+an].split()[2]))
                coordstemp.append(float(filedata[ln*repeatUnit+2+an].split()[3]))
            coord[ln,:,:]=np.reshape(coordstemp,(self.nAtoms,3))
            try:
                descCoord[ln]=float(filedata[ln*repeatUnit+1].split()[0])
            except:
                descCoord[ln]=1.0
        #coord=coord/au2ang #FIGURE OUT IF THIS IS NECESSARY - I don't think it is, fam
        return coord,descCoord #,pe,dip

    def propagate(self,x,nSteps,setV_ref=False,ConstantV_ref=0,printCensus=True,initialPop=0,testing=False):
        #print 'ready to propagate for',nSteps, 'steps on x (shaped:',x.shape,') '
        if testing:
            np.random.seed(0)

        if initialPop==0:#if it is the default value of zero...it needs to be set.
            initialPop=x.shape[0]
        self.currentPop=x.shape[0]
        if setV_ref:
            v_ref=ConstantV_ref
        else:
            #v_ref=np.average(self.molecule.V(x))+(self.alpha*(1.0-float(self.currentPop)/float(initialPop)))!!
            v_ref = np.average(self.getPot(x)) + (self.alpha * (1.0 - float(self.currentPop) / float(initialPop)))
        population=[]
        population.append(self.currentPop)
        
        vRefList=[]
        vRefList.append(v_ref)

        descendants=np.zeros(self.currentPop)
        whoYaFrom=np.arange(self.currentPop)
        
        printRate=100

        for step in range(nSteps):
            if step%printRate==0 and printCensus: print 'step #',step, 
            dx=self.diffuse()
            x=x+dx

            #v=self.molecule.V(x)!!!!!!!!!!!!!!!!!!!!
            v = self.getPot(x,nSteps)
            #Elimination of a random selection of walkers 
            #in the classically forbidden region
            N_r=np.random.random(self.currentPop)
            P_d=1-np.exp(-(v-v_ref)*self.dtau)
            Diff=N_r-P_d
            mask_survive=(Diff>0)
            nDeaths=np.sum(np.array(Diff<0).astype(int))

            if step%printRate==0 and printCensus: print 'Census: Current Population:',self.currentPop,'Deaths:',nDeaths,
            
            #LASSIE MODE
            #DEATH BY PES HOLES
            mask_not_holes=(v>-0.0000005) #mask
            inHole=(v<=-0.0000005)
            mask_survive=np.logical_and(mask_survive,mask_not_holes)
            if step%printRate==0 and printCensus: print 'Death by holes',np.sum(np.array(v<0.0).astype(int)),
            # removal of all walkers that cross and a random selection of walkers that are too close to the node

            if self.recrossing:
                #oldDist,swapped=self.molecule.calcRn(x-dx)
                #newDist,newswapped=self.molecule.calcRn(x)
                oldDist=self.molecule.calcRn(x-dx)
                newDist=self.molecule.calcRn(x)
                crossed=(oldDist*newDist<0)

                newReducedMass=self.molecule.calcReducedmass(x)
                oldReducedMass=self.molecule.calcReducedmass(x-dx)
                P_recrossDeath=np.exp(-2.0*(oldDist)*newDist*np.sqrt(oldReducedMass*newReducedMass)/self.dtau) ##mass is reduced mass!
                #P_recrossDeath=np.exp(-2.00000000*(oldDist)*newDist*newReducedMass/(self.dtau)) #
                #P_recrossDeath[newswapped]=10.0

                #caution, this is clever, you are combining the probability of crossing and the probability of recrossing
                N_r_recross=np.random.random(self.currentPop)
                Diff=N_r_recross-P_recrossDeath
                mask_survive_recross=(Diff>0.000000)
                mask_died_by_recrossing=(Diff<0.00000)
                tempRecrossCensus=np.sum(np.array(Diff<0.0000000).astype(int))
                mask_survive=np.logical_and(mask_survive,mask_survive_recross)
                #if len(newswapped)>0: print step,'     new swapped walkers that survived:',newswapped,'crossed:',crossed[newswapped],'died by recrossing:',mask_died_by_recrossing[newswapped],'survived:',mask_survive[newswapped]
                #if len(swapped)>0: print step,'     swapped walkers that survived',swapped, 'crossed:',crossed[swapped],'died by recrossing:',mask_died_by_recrossing[swapped],'survived:',mask_survive[swapped]
                #if len(swapped)>0 or len(newswapped)>0:print 'P_rec death: ',P_recrossDeath[swapped],P_recrossDeath[newswapped],
                #if len(swapped)>0 or len(newswapped)>0 : print 'did they survive?',mask_survive[swapped], mask_survive[newswapped],v[swapped],v[newswapped]
                if step%printRate==0 and  printCensus: print 'Deaths by recrossing: ',tempRecrossCensus, 'crossed',np.sum(crossed.astype(int))



            survivors=x[mask_survive]            

            # Creation of a random selection of walkers in the classically allowed region and who did not die by recrossing, so the survivors
        
            P_exp_b=np.exp(-(v-v_ref)*self.dtau)-1.0
            P_exp_b[np.logical_not(mask_survive)]=0.000000 #BECAUSE THE DEAD CANNOT REPRODUCE!! NO ZOMBIE MOTHERS!!
            
            #no really, this is important! it can cause at least a 3% difference in the energy!
            #if self.recrossing:
            #    #P_exp_b[mask_died_by_recrossing]=0.00000  #if it crossed, the probability that it gives birth should be zero
            #    P_exp_b[crossed]=0.00000  #if it crossed, the probability that it gives birth should be zero   
            #P_exp_b[inHole]=0.00000  #if it fell in a hole

            weight_P_b=P_exp_b.astype(int)
            P_b=P_exp_b-weight_P_b            #classically allowed region                            
            
            Diff=N_r-P_b
            mask_b = (Diff<0)
            #if self.recrossing:
                #if len(swapped)>0 or len(newswapped)>0 : print 'were any swapped walkers parents?', np.any(mask_b[swapped])
            next_gen=x[mask_b]
            new_pop_whoYaFrom=whoYaFrom[mask_b]
            nBirths=np.sum(np.array(Diff<0).astype(int))
            addBirthtot=0
            new_pop=next_gen
            
           #for the additional births                                                   
            for n,(particle,weight) in enumerate(zip(x,weight_P_b)):
                if weight>0: #i.e. the dead can't reproduce                              
                    if weight>10:
                        #this really shouldn't happen                                   
                        print 'weight of ',n,' is too big, resetting to 0'
                        print v[n]*au2wn,'<',v_ref*au2wn, weight, '\n',x[n]
                        print 'is it in a hole?',inHole[n], P_exp_b[n]
                        print 'this really should not be happening is a sign of a larger issue!'
                        weight=0

                    addBirthtot=addBirthtot+weight
                    temp=np.array([particle])
                    temp_whoYaFrom=np.array([whoYaFrom[n]])
                    for i in range(weight-1):
                        temp=np.concatenate((temp,np.array([particle])))
                        temp_whoYaFrom=np.concatenate((temp_whoYaFrom,np.array([whoYaFrom[n]])))

                    new_pop=np.concatenate((new_pop,temp))
                    new_pop_whoYaFrom=np.concatenate((new_pop_whoYaFrom,temp_whoYaFrom))

            next_gen=new_pop
            next_gen_whoYaFrom=new_pop_whoYaFrom

            #collect survivors and next generation                                       
            new_population=np.concatenate((survivors,next_gen))
            whoYaFrom=np.concatenate((whoYaFrom[mask_survive],next_gen_whoYaFrom))

            x=new_population
            
            #let's make sure none of the walkers crossed to the other side...
            #position_on_coord=self.molecule.calcRn(x)
            #print position_on_coord[0],position_on_coord.shape
            #if step%printRate==0 and printCensus:print 'wrong side!', position_on_coord[(position_on_coord<0)],np.average(position_on_coord)
            #if step%printRate==0 and printCensus:print 'from these array locations!',np.where(position_on_coord<0)
#            if np.where(position_on_coord>0)[0].shape[0]>1:
#               end
            
            self.currentPop=x.shape[0]
            if step%printRate==0 and printCensus: print 'Births:',nBirths, "Extra Births: ", addBirthtot,

            #readjust V_ref
            #v_average=np.average(self.molecule.V(new_population))!!!!!!!!!!!!!!!!!!!!!!
            v_average = np.average(self.getPot(new_population))
            if not setV_ref:
                v_ref=v_average+(self.alpha*(1.00-float(self.currentPop)/float(initialPop)))
            if step%printRate==0 and printCensus: print 'v_ref',v_ref, '=', v_average,'+', (self.alpha*(1-float(self.currentPop)/float(initialPop)))
            vRefList.append(v_ref)
            population.append(self.currentPop)
            if self.currentPop<5:
                print 'massive walker die off!', end

        for anc in whoYaFrom:
            descendants[anc]=descendants[anc]+1
        sys.stdout.flush()
        return  vRefList, population, x, descendants





#wavefunction
#initialize

#propagate

#calculate V


