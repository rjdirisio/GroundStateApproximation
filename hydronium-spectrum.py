import numpy as np
import DMCClusters as dmc
import CalculateSpectrum
import sys, os

angstr = 0.529177

# H E R M I T E  P O L Y N O M I A L  A P P R O X I M A T I O N
au2wn = 219474.63
nBins = 51


molecule = 'H3O+'
dTau = 10

coordinateSet = sys.argv[1]
testName = sys.argv[2]
kill = sys.argv[3]

print 'cd set', coordinateSet
print 'pmz', testName
print 'kill', kill

Wfn = dmc.wavefunction(molecule, 1)  # define our wavefunction
head = '/home/netid.washington.edu/rjdiri/'

cdsPath = '../h3oStuff/'
gPath = '../Gmats/h3o/'
dipPath = '../h3oStuff/'
cds = cdsPath + coordinateSet
GfileName = gPath + coordinateSet + '_' + testName + "_" + kill + '.gmat'
dipF = dipPath + 'eng_dip_' + coordinateSet + '.dat'

print cds
print GfileName
print dipF

if 'nz' in kill:
    GfileName = GfileName[:-4] + 'nz.gmat'
    # Same as regular

# lukeOOs = 3
# for f in range(1,lukeOOs+1):
#     a = np.load("dvr_oo_geoms_shift_"+str(f)+".npy")
#     fll = open('lukeGeomz_'+str(f)+'.xyz', 'w+')
#     Wfn.molecule.printCoordsToFile(a,fll)
# stop
# fll = open('cdsWagProblem.xyz','w+')
# a = np.load("../cdsWagProblem.npy")
# # Wfn.molecule.printCoordsToFile(a,fll)
# symDw = np.load("../cdsWagProblem_dw.npy")
# plotStuff(a)
# stop

if os.path.isfile(cds + '.npy'):
    symCoords = np.load(cds + '.npy')
    symDw = np.load(cds + '_dw.npy')
    # symCoords = Wfn.molecule.rotateBackToFrame(symCoords, 3, 2, 1)
    # np.save('../coordinates/trimer/'+'Tfinal_allH'+'.npy', symCoords)
    # np.save('../coordinates/trimer/' + coordinateSet + '_dw.npy', symDw)
    # print 'asdf'
else:
    symCoords, symDw = Wfn.loadCoords(cds)
    np.save('../h3oStuff/' + coordinateSet + '.npy', symCoords)
    np.save('../h3oStuff/' + coordinateSet + '_dw.npy', symDw)

print 'Symcoords shape', symCoords.shape
print 'Got symCoords!'
print 'NUMBER OF WALKERS IN allH: ', symCoords.shape[0]
symEckRotCoords = symCoords
eckt = False
if os.path.isfile(dipPath + 'eng_dip_' + coordinateSet + '_eckart.npy'):
    pdip = np.load(dipPath + 'eng_dip_' + coordinateSet + '_eckart.npy')
    eckt = True
elif os.path.isfile(dipF[:-3] + 'npy'):
    pdip = np.load(dipF[:-3] + 'npy')
else:
    print 'not a npy file'
    pdip = np.loadtxt(dipF)
    np.save(dipF[:-3] + 'npy', pdip)
path = '../spectra/'
print 'PEDIP shape', pdip.shape
pe = pdip[:, 0]
dip = pdip[:, 1:]
print 'Shape of dipole: ', np.shape(dip)
HOASpectrum = CalculateSpectrum.HarmonicApproxSpectrum(Wfn, symEckRotCoords, symDw, path, testName)
# if 'Eck' in GfileName:
#     coordinateSet=coordinateSet+'refGmat'
fundamentalEnergies, fundamentalIntensities, combinationBandEnergies, combinationBandIntensities = HOASpectrum.calculateSpectrum(
    symEckRotCoords, symDw, GfileName, pe, dip, coordinateSet, testName, kill, eckt, dipPath)
fundamentalFile = open('../spectra/Fundamentals_' + coordinateSet + testName + kill, 'w')
for i, (v, intensity) in enumerate(zip(fundamentalEnergies, fundamentalIntensities)):
    fundamentalFile.write(str(i) + "       " + str(v) + "   " + str(intensity) + "\n")
fundamentalFile.close()
