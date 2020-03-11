import numpy as np
import DMCClusters as dmc
import CalculateSpectrum
import sys, os
import Plot
def PltHists1D(cfg, thing, bound, xl, yl, overly, weits):
    theLen, xx = np.histogram(thing, bins=120, range=bound, normed=True, weights=weits)  # WEIGHTS=WEIGHTARRAY
    inin = True
    overlay = False
    bnd = str(bound[0]).replace(".", "").replace("-", "") + str(bound[1]).replace(".", "").replace("-", "")
    print bnd
    mP = Plot.myPlot(cfg, '1d', bnd, bnd, xl, yl, overly, inin, theLen,(xx[1:]+xx[:-1])/2)
    mP.plotIt()
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


if os.path.isfile(cds + '.npy'):
    symCoords = np.load(cds + '.npy')
    symDw = np.load(cds + '_dw.npy')
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

# internals = Wfn.molecule.SymInternalsH3O(symEckRotCoords).T
# internalName = Wfn.molecule.internalName
# unit = ['r','r','r','a','a','a','a']
# for i in range(len(internals)):
#     if unit == 'r':
#         hInt=internals[i]*angstr
#         rng = (0,2)
#     else:
#         hInt=np.degrees(internals[i])
#         rng = (0,160)
#     PltHists1D('h3o',hInt,rng,internalName[i],'H3O/ProbabilityDensity',False,weits=symDw)
# stop
HOASpectrum = CalculateSpectrum.HarmonicApproxSpectrum(Wfn, symEckRotCoords, symDw, path, testName)
# if 'Eck' in GfileName:
#     coordinateSet=coordinateSet+'refGmat'
fundamentalEnergies, fundamentalIntensities, combinationBandEnergies, combinationBandIntensities = HOASpectrum.calculateSpectrum(
    symEckRotCoords, symDw, GfileName, pe, dip, coordinateSet, testName, kill, eckt, dipPath)
fundamentalFile = open('../spectra/Fundamentals_' + coordinateSet + testName + kill, 'w')
for i, (v, intensity) in enumerate(zip(fundamentalEnergies, fundamentalIntensities)):
    fundamentalFile.write(str(i) + "       " + str(v) + "   " + str(intensity) + "\n")
fundamentalFile.close()
