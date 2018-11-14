import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import DMCClusters as dmc
import time
import glob
import Plot
import CalculateSpectrum
import sys, os
angstr = 0.529177

def PltHists1D(cfg, thing, bound, xl, yl, overly, weits):
    theLen, xx = np.histogram(thing, bins=25, range=bound, normed=True, weights=weits)  # WEIGHTS=WEIGHTARRAY
    inin = True
    overlay = False
    bnd = str(bound[0]).replace(".", "").replace("-", "") + str(bound[1]).replace(".", "").replace("-", "")
    print bnd
    mP = Plot.myPlot(cfg, '1d', bnd, bnd, xl, yl, overly, inin, theLen,(xx[1:]+xx[:-1])/2)
    mP.plotIt()

def plotStuff():
    print 'INTERNAL COORDINATES :-O'
    internals = Wfn.molecule.SymInternalsH9O4plus(symEckRotCoords)
    print 'Internal coordinate shape: ', np.shape(internals)
    print 'One attribute shape: ',np.shape(internals[:,0])
    print 'number of dws: ', symDw
    PltHists1D('allH', symEckRotCoords[:, 0, -1], (-2, 2), 'zcomp O1', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 1, -1], (-2, 2), 'zcomp O2', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 2, -1], (-2, 2), 'zcomp O3', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 3, -1], (-2, 2), 'zcomp O4', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 4, -1], (-2, 2), 'zcomp H5', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 5, -1], (-2, 2), 'zcomp H6', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 6, -1], (-2, 2), 'zcomp H7', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 7, -1], (-2, 2), 'zcomp H8', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 8, -1], (-2, 2), 'zcomp H9', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 9, -1], (-2, 2), 'zcomp H10', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 10, -1], (-2, 2), 'zcomp H11', 'tetramerInternals/Probability Denisty',
               False, symDw)

    PltHists1D('allH', internals[:, 0]*angstr, (-2, 2), 'xComp of Shared H11', 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 1]*angstr, (-2, 2), 'yComp of Shared H11', 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 2]*angstr, (-2, 2), 'zComp of Shared H11', 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 3]*angstr, (-2, 2), 'xComp of Shared H12', 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 4]*angstr, (-2, 2), 'yComp of Shared H12', 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 5]*angstr, (-2, 2), 'zComp of Shared H12', 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 6]*angstr, (-5, 5), 'xComp of Shared H13', 'tetramerInternals/Probability Density',
               False,
               symDw)
    PltHists1D('allH', internals[:, 7]*angstr, (-4, 4), 'yComp of Shared H13', 'tetramerInternals/Probability Density',
               False,
               symDw)
    PltHists1D('allH', internals[:, 8]*angstr, (-2, 2), 'zComp of Shared H13', 'tetramerInternals/Probability Density',
               False,
               symDw)

    PltHists1D('allH', np.rad2deg(internals[:, 9]), (-360, 360), 'Theta_651',
               'tetramerInternals/Probability Density', False, symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 10]), (-360, 360), 'Phi_651', 'tetramerInternals/Probability Density',
               False, symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 11]), (-360, 360), 'Chi_651','tetramerInternals/Probability Density', False, symDw)


    PltHists1D('allH', np.rad2deg(internals[:, 12]), (-360, 360), 'Theta_1039',
               'tetramerInternals/Probability Density', False, symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 13]), (-360, 360), 'Phi_1039',
               'tetramerInternals/Probability Density',
               False, symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 14]), (-360, 360), 'Chi_1039',
               'tetramerInternals/Probability Density', False, symDw)

    PltHists1D('allH', np.rad2deg(internals[:, 15]), (-360, 360), 'Theta_728',
               'tetramerInternals/Probability Density', False, symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 16]), (-360, 360), 'Phi_728',
               'tetramerInternals/Probability Density',
               False, symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 17]), (-360, 360), 'Chi_728',
               'tetramerInternals/Probability Density', False, symDw)

    PltHists1D('allH', internals[:, 18],(1,3), 'rO1H5', 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 19], (1,3), 'rO1H6', 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 20]), (70,180), 'aHOH_516', 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 21], (1,3), 'rOH_27', 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 22], (1,3), 'rOH_28', 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 23]), (70,180), 'aHOH_728', 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 24], (1, 3), 'rOH_9', 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 25], (1, 3), 'rOH_10', 'tetramerInternals/Probability Density', False,
               symDw)

    """'HOH516', 'rOH7', 'rOH8', 'HOH728',
                         'rOH9', 'rOH10', 'HOH9310', 'rO4O1', 'rO4O2', 'rO4O3', 'Oumbrella', 'thetaOx', 'thetaOx2']"""

    PltHists1D('allH', np.rad2deg(internals[:, 26]), (70, 180), 'aHOH_1039', 'tetramerInternals/Probability Density',
               False,
               symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 30]), (0,360), 'OUmbrella',
               'tetramerInternals/Probability Density',
               False,
               symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 31]), (0,360), 'thetaOx1',
               'tetramerInternals/Probability Density',
               False,
               symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 32]), (0,360), 'thetaOx2',
               'tetramerInternals/Probability Density',
               False,
               symDw)



# H E R M I T E  P O L Y N O M I A L  A P P R O X I M A T I O N
au2wn=219474.63
nBins=51


starttime=time.time()

stateGround='stateGround'
state='stateGround'
DWstate='DWGround'
molecule='H9O4'
dTau=10

coordinateSet = sys.argv[1]
testName = sys.argv[2]
kill = sys.argv[3]

print 'cd set',coordinateSet
print 'pmz',testName
print 'kill',kill

Wfn=dmc.wavefunction('H9O4+', 1) #define our wavefunction
if 'allD' in coordinateSet:
    Wfn.setIsotope('fullyDeuterated')
if '1He' in coordinateSet:
    Wfn.setIsotope('notDeuteratedOnce_eigen')
if '1Hw' in coordinateSet:
    Wfn.setIsotope('notDeuteratedOnce_fw')

head = '/home/netid.washington.edu/rjdiri/'
#Lindsey's code: Be careful.  Here, we calculate the HOA spectrum on each individual wavefunction, but then also on the overall xyz file.  
#Because of this, in my first step, I really only need to do it for one 'wavefunction' aka my (Anne's) entire simulation

#if 'tet_full' in coordinateSet:
#    cds = head+'Documents/h9o4/freshlyRotatedCoordinates/symmetrizedToFilth'+coordinateSet
#    GfileName = 'allHTesting/spectra/tet_TheGMatrix_'+coordinateSet+'.gmat'
#    dipF = 'eng_dip_'+coordinateSet+'.dat'

cdsPath= '../coordinates/tetramer/'
gPath = '../Gmats/tetramer/'
dipPath='../newDipolesAndPE/tetramer/'
cds =cdsPath+coordinateSet
GfileName = gPath+coordinateSet+'.gmat'
dipF = dipPath+'eng_dip_'+coordinateSet+'.dat'

print cds
print GfileName
print dipF


if 'nz' in kill:
    GfileName = GfileName[:-4]+'nz.gmat'
    #Same as regular
#dipF = "dipoles/"+dipF



if os.path.isfile(cds+'.npy'):
    symCoords=np.load(cds+'.npy')
    symDw = np.load(cds+'_dw.npy')
else:
    symCoords, symDw = Wfn.loadCoords(cds)
    np.save(cds,symCoords)
    np.save(cds+'_dw',symDw)

if os.path.isfile(dipF[:-3]+'npy'):
    pdip=np.load(dipF[:-3]+'npy')
else:
    pdip=np.loadtxt(dipF)
    np.save(dipF[:-3]+'npy',pdip)

print 'PEDIP shape',pdip.shape

print 'Symcoords shape',symCoords.shape
print 'Got symCoords!'
pe = pdip[:,0]
dip = pdip[:,1:]
print 'Shape of dipole: ', np.shape(dip)
print 'NUMBER OF WALKERS IN allH: ',symCoords.shape[0]
symEckRotCoords = symCoords
iwantToPlotStuff=False
path='../spectra/'
if iwantToPlotStuff:
    plotStuff()
else:
    HOASpectrum=CalculateSpectrum.HarmonicApproxSpectrum(Wfn,symEckRotCoords,symDw,path,testName)
    if 'Eck' in GfileName:
        coordinateSet=coordinateSet+'refGmat'
    fundamentalEnergies,fundamentalIntensities, combinationBandEnergies,combinationBandIntensities=HOASpectrum.calculateSpectrum(symEckRotCoords,symDw,GfileName,pe,dip,coordinateSet,testName,kill)
    fundamentalFile=open('../spectra/Fundamentals_'+coordinateSet+testName+kill,'w')
    for i,(v,intensity) in enumerate(zip(fundamentalEnergies,fundamentalIntensities)):
        fundamentalFile.write(str(i)+"       "+str(v)+"   "+str(intensity)+"\n")
    fundamentalFile.close()
