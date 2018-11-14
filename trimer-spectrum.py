#!/usr/bin/python
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import DMCClusters as dmc
import time
import sys
import os
import glob
import usefulFunctions as use
import Plot
import CalculateSpectrum

angstr = 0.529177

def PltHists1D(cfg, thing, bound, xl, yl, overly, weits):
    theLen, xx = np.histogram(thing, bins=75, range=bound, normed=True, weights=weits)  # WEIGHTS=WEIGHTARRAY
    inin = True
    overlay = False
    bnd = str(bound[0]).replace(".", "").replace("-", "") + str(bound[1]).replace(".", "").replace("-", "")
    print bnd
    mP = Plot.myPlot(cfg, '1d', bnd, bnd, xl, yl, overly, inin, theLen)
    mP.plotIt()




# H E R M I T E  P O L Y N O M I A L  A P P R O X I M A T I O N
au2wn=219474.63
nBins=51

if len(sys.argv)<4:
    print 'Usage: ./HOHOH-groundState.py N_size nReps descendantSteps nRepsDW'
    end


starttime=time.time()

stateGround='stateGround'
state='stateGround'
DWstate='DWGround'
molecule='H7O3'
dTau=10

N_size=int(sys.argv[1])
nReps=int(sys.argv[2])
descendantSteps=int(sys.argv[3])
nRepsDW=int(sys.argv[4])
nStart=int(sys.argv[5])
coordinateSet = sys.argv[6]
eckG = sys.argv[7]
testName = sys.argv[8]
kill = sys.argv[9]
#figure out which nrep number we're at in the directory of interest
#fileParameterName=molecule+'-'+state+'-'+DWstate+'-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)
path='allHTesting/'
#print 'path: ', path
#preExistingFiles=glob.glob(path+'*'+fileParameterName+'*')
#print 'the files that already exist are:', preExistingFiles

#outputFile=open('resultsH7O3/testing.data','w')
#outputFile.write('the files that already exist are: '+str(preExistingFiles)+'\n')

equilibrationSteps=500
propagationSteps=500
averaged_vref=[]
list_of_pop_list=[] #not used

Wfn=dmc.wavefunction('H7O3+', N_size) #define our wavefunction
if 'allD' in coordinateSet:
    Wfn.setIsotope('fullyDeuterated')
if '1He' in coordinateSet:
    Wfn.setIsotope('notDeuteratedOnce_eigen')
if '1H_w' in coordinateSet:
    Wfn.setIsotope('notDeuteratedOnce_fw') #good
if '1H_h' in coordinateSet:
    Wfn.setIsotope('notDeuteratedOnce_hydronium')

Wfn.setNodalSurface('OHStretchAnti','Both') #define one nodal surface, here is Antisymm stretch
gatheredSymEckRotCoords=[] #we are going to get a series of eckart rotated coordinates, and
gatheredSymDW=[]           #descendent weights (symmetrized?)
nWalkersTotal=0
#for each iwfn in nReps, 


head = '/home/netid.washington.edu/rjdiri/'
#Lindsey's code: Be careful.  Here, we calculate the HOA spectrum on each individual wavefunction, but then also on the overall xyz file.  
#Because of this, in my first step, I really only need to do it for one 'wavefunction' aka my (Anne's) entire simulation
for iwfn in range(nStart,nReps):
    #print '   REPETITION NUMBER: ', iwfn
    #groundPath='data'+molecule+'/hermite/'

    if 'full' in coordinateSet:
        cds = head+'Documents/h7o3/Hermite_others/SymmetrizedWalkers_r_s_r/finalFull_'+coordinateSet[-4:]  # big boi that I symmetrized before rotating.
        GfileName = 'allHTesting/spectra/TheGMatrixFinalFull_'+coordinateSet[-4:]+'.gmat'
        dipF = 'eng_dip_full_'+coordinateSet[-4:]+'.dat'

    #old
    #elif 'symOdd' in coordinateSet:
    #    cds = head + 'Documents/h7o3/Hermite_others/SymmetrizedWalkers_r_s_r/symOdd_allD'# big boi that I symmetrized before rotating.
    #    GfileName = 'allHTesting/spectra/TheGMatrixFinalFull_symOdd_allD.gmat'
    #    dipF = 'eng_dip_symOdd_allD.dat'

    #elif 'symEven' in coordinateSet:
    #    cds = head + 'Documents/h7o3/Hermite_others/SymmetrizedWalkers_r_s_r/symEven_allD'# big boi that I symmetrized before rotating.
    #    GfileName = 'allHTesting/spectra/TheGMatrixFinalFull_symEven_allD.gmat'
    #    dipF = 'eng_dip_symEven_allD.dat'

    elif 'odd_allH' in coordinateSet:
        cds = head + 'Documents/h7o3/Hermite_others/SymmetrizedWalkers_r_s_r/symOdd_allH'# big boi that I symmetrized before rotating.
        GfileName = 'allHTesting/spectra/TheGMatrixFinalFull_symOdd_allH.gmat'
        dipF = 'eng_dip_symOdd_allH.dat'

    elif 'even_allH' in coordinateSet:
        cds = head + 'Documents/h7o3/Hermite_others/SymmetrizedWalkers_r_s_r/symEven_allH'# big boi that I symmetrized before rotating.
        GfileName = 'allHTesting/spectra/TheGMatrixFinalFull_symEven_allH.gmat'
        dipF = 'eng_dip_symEven_allH.dat'

    elif 'odd_allD2' in coordinateSet:
        cds = head + 'Documents/h7o3/Hermite_others/SymmetrizedWalkers_r_s_r/symOdd_allD2'# big boi that I symmetrized before rotating.
        GfileName = 'allHTesting/spectra/TheGMatrixFinalFull_symOdd_allD2.gmat'
        dipF = 'eng_dip_symOdd_allD2.dat'

    elif 'even_allD2' in coordinateSet:
        cds = head + 'Documents/h7o3/Hermite_others/SymmetrizedWalkers_r_s_r/symEven_allD2'# big boi that I symmetrized before rotating.
        GfileName = 'allHTesting/spectra/TheGMatrixFinalFull_symEven_allD2.gmat'
        dipF = 'eng_dip_symEven_allD2.dat'


    elif coordinateSet == 'test_1000':
        GfileName = 'allHTesting/spectra/TheGMatrix_1000new.gmat'
        cds = head+'Documents/h7o3/Hermite_others/test_1000_allH'  # 1000 test file
        dipF = 'test_1000_allHDips'

    elif coordinateSet == 'test_10k':
        GfileName = 'allHTesting/spectra/TheGMatrix_10knew.gmat'
        cds = head+'Documents/h7o3/Hermite_others/test_10k'  # 1000 test file
        dipF = 'eng_dip_test_10k.dat'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Testing with my new spectrum

    elif coordinateSet == 'herm':
        cds = head+'Documents/h7o3/flavorsOfCoordinates/rotatedOrig_allH' #not symmetrized
        GfileName = 'allHTesting/spectra/TheGMatrixnew_herm_allH.gmat'
        dipF = 'eng_dip_rotatedOrig_allH.dat'

    elif coordinateSet == 'full_noZ_allH':
        cds = head+'Documents/h7o3/flavorsOfCoordinates/full_noZ_allH' #not symmetrized
        GfileName = 'allHTesting/spectra/TheGMatrixnew_full_noZ_allH.gmat'
        dipF = 'eng_dip_full_noZ_allH.dat'

    elif coordinateSet == 'justFWSwap_allH':
        cds = head+'Documents/h7o3/flavorsOfCoordinates/justFWSwap_allH' #not symmetrized
        GfileName = 'allHTesting/spectra/TheGMatrixnew_justFWSwap_allH.gmat'
        dipF = 'eng_dip_justFWSwap_allH.dat'

    elif coordinateSet == 'justFWSwap_pmZ_allH':
        cds = head+'Documents/h7o3/flavorsOfCoordinates/justFWSwap_pmZ_allH' #not symmetrized
        GfileName = 'allHTesting/spectra/TheGMatrixnew_justFWSwap_pmZ_allH.gmat'
        dipF = 'eng_dip_justFWSwap_pmZ_allH.dat'

    elif coordinateSet == 'onlyPMZ_allH':
        cds = head+'Documents/h7o3/flavorsOfCoordinates/onlyPMZ_allH' #not symmetrized
        GfileName = 'allHTesting/spectra/TheGMatrixnew_onlyPMZ_allH.gmat'
        dipF = 'eng_dip_onlyPMZ_allH.dat'




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Testing
    """elif coordinateSet == 'test':
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/in100AllH/xyz_trC_cutH'  # 100 walker test file
        GfileName = 'allHTesting/spectra/TheGMatrix_100.gmat'
    elif coordinateSet == 'test_1000':
        GfileName = 'allHTesting/spectra/TheGMatrix_1000.gmat'
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite/xyz_trC_thouH'  # 1000 test file
        #dipF = 'test_1000Dipoles'

    elif coordinateSet == 'allD_test_1000':
        GfileName = 'allHTesting/spectra/TheGMatrixallD_1000.gmat'
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite_others/allD_test_1000'  # 1000 test file
        #dipF = 'allD_test_1000Dipoles'

    elif coordinateSet == 'symtest_10000':
        GfileName = 'allHTesting/spectra/TheGMatrixsym_10kAllH_Z.gmat'
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite/sym10kAllH_Z'  # 1000 test file
        #dipF = 'symtest_10000'

    elif coordinateSet == 'one':
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite/oneWalker'  # one walker for very detailed testing
        GfileName = 'allHTesting/spectra/TheGMatrix_One.gmat'

    elif coordinateSet == 'rtest':
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/testFileForRotation/notRotated.xyz'  # rotation test file
    elif coordinateSet == 'test_10000':
        GfileName = 'allHTesting/spectra/TheGMatrix10000.gmat'
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite/xyz_trC_13thou' #10,000 walker test file
        #dipF = 'test_10000Dipoles'
    elif coordinateSet == 'refEck':
        GfileName = 'allHTesting/spectra/TheGMatrix_EckRef.gmat'
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite/final_eckartReference.xyz'
    elif coordinateSet == 'eckNotRotated':
        GfileName = 'allHTesting/spectra/TheGMatrix_notRoteckart_Dontuse.gmat'
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/getEckartFrame/eckart.xyz'
        #dipF = 'test_10000Dipoles'

    elif coordinateSet == 'symNOZ_10k': #test without +/- Z
        GfileName = 'allHTesting/spectra/TheGMatrixSymNOZ10k.gmat'
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite/symNOZ_10k' #10k * 8 permutations
        #dipF = 'symNOZ_10kDipoles'

    #    elif coordinateSet == 'symtest_10000NOZ': #test without +/- Z
#        GfileName = 'allHTesting/spectra/TheGMatrixSymNOZ10k.gmat'
#        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite/10kSymmed_allH_noZ' #10k * 8 permutations
#        dipF = 'symNOZTest10kDips' #currently crap

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Deuterations

    elif coordinateSet == 'herm_allD':
        GfileName = 'allHTesting/spectra/TheGMatrix_allDHerm.gmat'
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite_others/Hermite_allD'
        #dipF = 'hermDipolesallD'

    elif coordinateSet == 'herm_1H_w':
        GfileName = 'allHTesting/spectra/TheGMatrix_1H_wHerm.gmat'
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite_others/Hermite_1H_w'
        #dipF = 'hermDipoles1H_w'


    elif coordinateSet == 'herm_1H_h':
        GfileName = 'allHTesting/spectra/TheGMatrix_1H_hHerm.gmat'
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite_others/Hermite_1H_h'
        #dipF = 'hermDipoles1H_h'


    elif coordinateSet == 'herm_1He':
        GfileName = 'allHTesting/spectra/TheGMatrix_1HeHerm.gmat'
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite_others/Hermite_1He'
        #dipF = 'hermDipoles_1He'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Full deuteration sets - Symmetrized, no +/- Z in symmetrization
    elif coordinateSet == 'fullNOZ':
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite/symNOZFull'
        GfileName = 'allHTesting/spectra/TheGMatrixFullNOZ.gmat'
        #dipF = 'fullNOZDipoles'
    elif coordinateSet == 'fullNOZ_allD':
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite_others/symHermiteNOZ_allD'
        GfileName = 'allHTesting/spectra/TheGMatrixFullNOZ_allD.gmat'
        #dipF = 'fullNOZDipoles_allD'
    elif coordinateSet == 'fullNOZ_1He':
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite_others/symHermiteNOZ_1He'
        GfileName = 'allHTesting/spectra/TheGMatrixFullNOZ_1He.gmat'
        #dipF = 'fullNOZDipoles_1He'
    elif coordinateSet == 'fullNOZ_1H_w':
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite_others/symHermiteNOZ_1H_w'
        GfileName = 'allHTesting/spectra/TheGMatrixFullNOZ_1H_w.gmat'
        #dipF = 'fullNOZDipoles_1H_w'
    elif coordinateSet == 'fullNOZ_1H_h':
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite_others/symHermiteNOZ_1H_h'
        GfileName = 'allHTesting/spectra/TheGMatrixFullNOZ_1H_h.gmat'
        #dipF = 'fullNOZDipoles_1H_h'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~C2 rotation , nor +/-Z included
    elif coordinateSet == 'symNOC2_allD':
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite_others/symHermiteNOC2_allD'
        GfileName = 'allHTesting/spectra/TheGMatrixsymNOC2_allD.gmat'
        #dipF = 'Dipoles_symNOC2_allD'

    elif coordinateSet == 'symNOC2_allH':
        cds = '/home/netid.washington.edu/rjdiri/Documents/h7o3/Hermite_others/symNOC2_allH'
        GfileName = 'allHTesting/spectra/TheGMatrixsymNOC2_allH.gmat'
        #dipF = 'Dipoles_symNOC2_allH'"""

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Fully symmetrized, all deuterations

    print GfileName
    if eckG == 't':
        GfileName = 'allHTesting/spectra/TheGMatrix_EckRef.gmat'
    if kill == 'k':
        GfileName = GfileName[:-4]+'killd.gmat'
    if 'nz' in kill:
        print ':3 +/- Z , same gmatrix as regular'
        #Same as regular

    #dipF = "dipoles/"+dipF
    dipF = "newDipolesAndPE/" + dipF

    #alreadyEckarted = '/home/netid.washington.edu/rjdiri/udrive/h7o3/LindseyDMC/allHTesting/eckartRotatedMolecule'
    #load in the ground state
    #groundStateWfnName='Wfn-'+str(iwfn)+'-'+molecule+'-stateGround-DWGround-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.xyz'

    #symCoords, symDW, pe, dip = Wfn.loadCoords(testing100)  # MODIFIED - PE AND DIPOLE MOMENT RETURNED
    #symCoords,symDW,pe,dip=Wfn.loadCoords(testingUnSymAllHName) #MODIFIED - PE AND DIPOLE MOMENT RETURNED

    #symCoords, symDW, pe, dip = Wfn.loadCoords(cds) getting PE and DW from before code running
    symCoords, symDW = Wfn.loadCoords(cds)

    if coordinateSet=='refEck':
        symCoords /= angstr
    print 'Got symCoords!'
    pdip = np.loadtxt(dipF)
    pe = pdip[:,0]
    dip = pdip[:,1:]
    print 'Shape of dipole: ', np.shape(dip)
    print 'NUMBER OF WALKERS IN allH: ',symCoords.shape[0]

    #symCoords2, symDW, pe, dip = Wfn.loadCoords(alreadyEckarted)  # MODIFIED - PE AND DIPOLE MOMENT RETURNED
    #print 'NUMBER OF WALKERS IN allH: ', symCoords2.shape[0]

    #This stuff was for testing eckart -----------------
    #anan=open('allHTesting/angstromedShiz',"w")
    #Wfn.molecule.printCoordsToFile(symCoords,anan)
    #anan.close()
    #print 'ANAN DONE'
    #print 'Internals without eckart: '
    #internals = Wfn.molecule.SymInternalsH7O3plus(symCoords, symDW) #TO CHECK MY OLD CODE AND TO MAKE SURE ECKART IS DONE PROPERLY
    #print 'DONE with internals without eckart'
    #---------------------------------------------------

    #print 'symCoords: ', symCoords
    #eckRotcoords=Wfn.molecule.eckartRotate(GroundCoords) #Eckart rotate ground state
    #nwalkers=GroundCoords.shape[0]

    #symCoords,symDW=Wfn.molecule.symmetrizeCoordinates(GroundCoords,groundDW) #symmetrize ground state (Not eckarted)
    #print 'Commence eckart rotation of symmetrized coordinates'
    #~~~~~~~~~~~~~Eckart NOT NECESSARY AT THE MOMENT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~symEckRotCoords=Wfn.molecule.eckartRotate(symCoords) #symmetrize ground state (eckarted)

    symEckRotCoords = symCoords #*angstr
    nWalkersTotal+=symCoords.shape[0]

    if coordinateSet == 'rtest':
        internals = Wfn.molecule.SymInternals(symEckRotCoords,True)
        doneWithTest


    if iwfn==nStart:
        gatheredSymEckRotCoords=symEckRotCoords
        gatheredSymDW=symDW
    else:
        gatheredSymEckRotCoords=np.append(gatheredSymEckRotCoords,symEckRotCoords,axis=0) #gatheredSymEckRotCoords and dw collect the coordinates and descendent weights after we have propagated
        gatheredSymDW=np.append(gatheredSymDW,symDW,axis=0)
    print "all Walkers: ",nWalkersTotal
    #Wfn.molecule.getZsOfEckart()
    #fileOut=open('symEckRotCoords.xyz', 'w')
    #Wfn.molecule.printCoordsToFile(symEckRotCoords,fileOut)
    #fileOut.close()
    #print 'two symmeterized walkers!'
    #i=26
    #print symEckRotCoords[i],'\n',symEckRotCoords[i+nwalkers]
    #ryantest=True
    #if ryantest:
    #    HOASpectrum = CalculateSpectrum.HarmonicApproxSpectrum(Wfn, symEckRotCoords, symDW, path,testName)  # Initialize class
    #    symEckRotCoords = HOASpectrum.calculateSpectrum(symEckRotCoords, symDW, GfileName, pe, dip, coordinateSet, testName, kill)




    iwantToPlotStuff=False

    if iwantToPlotStuff:
        print 'INTERNAL COORDINATES :-O'
        internals = Wfn.molecule.SymInternalsH7O3plus(symEckRotCoords)
        print np.shape(internals[:, 9])
        print np.shape(symDW)
        print 'Theta6', np.amin(internals[:, 9]), np.amax(internals[:, 9])
        print 'Phi6', np.amin(internals[:, 10]), np.amax(internals[:, 10])
        print 'Xi6', np.amin(internals[:, 11]), np.amax(internals[:, 11])
        print 'Theta4', np.amin(internals[:, 12]), np.amax(internals[:, 12])
        print 'Phi4', np.amin(internals[:, 13]), np.amax(internals[:, 13])
        print 'Xi4', np.amin(internals[:, 14]), np.amax(internals[:, 14])

        PltHists1D('allH', np.rad2deg(internals[:, 9]), (-360, 360), 'Theta_627', 'doubleCheckInternals/Probability Density', False, symDW)
        PltHists1D('allH', np.rad2deg(internals[:, 10]), (-360, 360), 'Phi_627', 'doubleCheckInternals/Probability Density', False, symDW)
        PltHists1D('allH', np.rad2deg(internals[:, 11]), (-360, 360), 'Xi_627', 'doubleCheckInternals/Probability Density', False, symDW)
        PltHists1D('allH', np.rad2deg(internals[:, 12]), (-360, 360), 'Theta_451', 'doubleCheckInternals/Probability Density', False, symDW)
        PltHists1D('allH', np.rad2deg(internals[:, 13]), (-360, 360), 'Phi_451', 'doubleCheckInternals/Probability Density', False, symDW)
        PltHists1D('allH', np.rad2deg(internals[:, 14]), (-360, 360), 'Xi_451', 'doubleCheckInternals/Probability Density', False, symDW)
        PltHists1D('allH', internals[:, 0], (-2, 2), 'xComp of Shared H9', 'doubleCheckInternals/Probability Density', False,
                   symDW)
        PltHists1D('allH', internals[:, 1], (-2, 2), 'yComp of Shared H', 'doubleCheckInternals/Probability Density', False,
                   symDW)
        PltHists1D('allH', internals[:, 2], (-2, 2), 'zComp of Shared H', 'doubleCheckInternals/Probability Density', False,
                   symDW)
        PltHists1D('allH', internals[:, 3], (-2, 2), 'xComp of Shared H10', 'doubleCheckInternals/Probability Density', False,
                   symDW)
        PltHists1D('allH', internals[:, 4], (-2, 2), 'yComp of Shared H10', 'doubleCheckInternals/Probability Density', False,
                   symDW)
        PltHists1D('allH', internals[:, 5], (-2, 2), 'zComp of Shared H10', 'doubleCheckInternals/Probability Density', False,
                   symDW)
        PltHists1D('allH', internals[:, 6], (1,3), 'rH8', 'doubleCheckInternals/Probability Density', False,
                   symDW)
        PltHists1D('allH', np.rad2deg(internals[:, 7]), (-360,360), 'thH8', 'doubleCheckInternals/Probability Density', False,
                   symDW)
        PltHists1D('allH', np.rad2deg(internals[:, 8]), (-360,360), 'phiH8', 'doubleCheckInternals/Probability Density', False,
                   symDW)
        PltHists1D('allH', internals[:, 15],(1,3), 'rOH41', 'doubleCheckInternals/Probability Density', False,
                   symDW)
        PltHists1D('allH', internals[:, 16], (1,3), 'rOH51', 'doubleCheckInternals/Probability Density', False,
                   symDW)
        PltHists1D('allH', np.rad2deg(internals[:, 17]), (70,180), 'aHOH_451', 'doubleCheckInternals/Probability Density', False,
                   symDW)
        PltHists1D('allH', internals[:, 18], (1,3), 'rOH_26', 'doubleCheckInternals/Probability Density', False,
                   symDW)
        PltHists1D('allH', internals[:, 19], (1,3), 'rOH_27', 'doubleCheckInternals/Probability Density', False,
                   symDW)
        PltHists1D('allH', np.rad2deg(internals[:, 20]), (70,180), 'aHOH_267', 'doubleCheckInternals/Probability Density', False,
                   symDW)
        PltHists1D('allH', internals[:, 21], (1,3), 'rOO_1', 'doubleCheckInternals/Probability Density', False,
                   symDW)
        PltHists1D('allH', internals[:, 22], (1,3), 'rOO_2', 'doubleCheckInternals/Probability Density', False,
                   symDW)
        PltHists1D('allH', np.rad2deg(internals[:, 23]), (70,180), 'aOOO', 'doubleCheckInternals/Probability Density', False,
                   symDW)

        stop

    #print np.shape(internals)
    #print internals[i],'\n',internals[i+nwalkers]
    #print np.average(internals,weights=symDW,axis=0) # get the average of the internal coordinates?)



    #GfileName='TheGMatrix-symmetrized'+str(iwfn)+'-'+molecule+'-stateGround-DWGround-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.gmat'
    #GfileName='TheGMatrix-symmetrized-all-'+molecule+'-stateGround-DWGround-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.gmat'



    HOASpectrum=CalculateSpectrum.HarmonicApproxSpectrum(Wfn,symEckRotCoords,symDW,path,testName) #Initialize class
    #HOASpectrum.calculateG(symEckRotCoords,symDW)
    #print 'dip!' , dip
    print 'gfileName',GfileName
    if 'Eck' in GfileName:
        coordinateSet=coordinateSet+'refGmat'
    fundamentalEnergies,fundamentalIntensities, combinationBandEnergies,combinationBandIntensities=HOASpectrum.calculateSpectrum(symEckRotCoords,symDW,GfileName,pe,dip,coordinateSet,testName,kill)
    fundamentalFile=open(path+'/spectra/Fundamentals_'+coordinateSet+testName+kill,'w')
    for i,(v,intensity) in enumerate(zip(fundamentalEnergies,fundamentalIntensities)):
        fundamentalFile.write(str(i)+"       "+str(v)+"   "+str(intensity)+"\n")
    fundamentalFile.close()
    """combinationFile=open(path+'/spectra/Combinations_'+coordinateSet,'w')
    for i in range(combinationBandEnergies.shape[0]):
        for j in range(i):
            combinationFile.write(str(i)+"  "+str(j)+"       "+str(combinationBandEnergies[i,j])+"    "+str(combinationBandIntensities[i,j])+"\n")
        for i in range(combinationBandEnergies.shape[0]):
            combinationFile.write(str(i)+"  "+str(i)+"       "+str(combinationBandEnergies[i,i])+"    "+str(combinationBandIntensities[i,i])+"\n") #??????????????????????????????
    combinationFile.close()"""


    





#gatheredSymEckRotCoords=np.array(gatheredSymEckRotCoords).reshape(nWalkersTotal,Wfn.nAtoms,3)
#gatheredSymDW=np.array(gatheredSymDW).reshape(nWalkersTotal)
""" LIL MODIFICATION - FOR LOOP JUST GOES OVER ONE THING, WHICH IS WHAT THIS BOTTOM CODE IS SUPPOSED TO DO.
HOASpectrum=CalculateSpectrum.HarmonicApproxSpectrum(Wfn,gatheredSymEckRotCoords,gatheredSymDW,path=path)
#HOASpectrum.calculateG(gatheredSymEckRotCoords,gatheredSymDW)
GfileName='TheGMatrix-symmetrized-all-'+molecule+'-stateGround-DWGround-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.gmat'
fundamentalEnergies,fundamentalIntensities, combinationBandEnergies,combinationBandIntensities=HOASpectrum.calculateSpectrum(gatheredSymEckRotCoords,gatheredSymDW,path+GfileName)
fundamentalFile=open(path+'Fundamentals-from-Wfns-'+fileParameterName+str(nStart)+'-to-'+str(nReps)+'.data','w')
for i,(v,intensity) in enumerate(zip(fundamentalEnergies,fundamentalIntensities)):
    fundamentalFile.write(str(i)+"       "+str(v)+"   "+str(intensity)+"\n")
fundamentalFile.close()
combinationFile=open(path+'Combinations-from-Wfns-'+fileParameterName+str(nStart)+'-to-'+str(nReps)+'.data','w')
for i in range(combinationBandEnergies.shape[0]):
    for j in range(i):
        combinationFile.write(str(i)+"  "+str(j)+"       "+str(combinationBandEnergies[i,j])+"    "+str(combinationBandIntensities[i,j])+"\n")
for i in range(combinationBandEnergies.shape[0]):
    combinationFile.write(str(i)+"  "+str(i)+"       "+str(combinationBandEnergies[i,i])+"    "+str(combinationBandIntensities[i,i])+"\n")
combinationFile.close()

    #print zip(Wfn.molecule.internalName,np.average(internals,weights=symDW,axis=0)*Wfn.molecule.internalConversion)
    
#    Wfn.LoadG(groundPath+GfileName,symEckRotCoords,symDW)"""



