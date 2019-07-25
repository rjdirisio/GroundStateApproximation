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
import numpy.linalg as la
angstr = 0.529177
au2wn=219474.63

def PltHists1D(cfg, thing, bound, xl, yl, overly, weits):
    theLen, xx = np.histogram(thing, bins=25, range=bound, normed=True, weights=weits)  # WEIGHTS=WEIGHTARRAY
    inin = True
    overlay = False
    bnd = str(bound[0]).replace(".", "").replace("-", "") + str(bound[1]).replace(".", "").replace("-", "")
    print bnd
    mP = Plot.myPlot(cfg, '1d', bnd, bnd, xl, yl, overly, inin, theLen,(xx[1:]+xx[:-1])/2)
    mP.plotIt()


def ba(xx, atm1, atm2, atm3):  # left center right
    # Rotation of O1 to xyplane
    atmO = xx[:, atm1, :]
    atmT = xx[:, atm2, :]
    atmH = xx[:, atm3, :]
    left = atmO - atmT
    right = atmH - atmT
    return np.arccos((left * right).sum(axis=1) / (la.norm(left, axis=1) * la.norm(right, axis=1)))

def bL(xx,atm1,atm2):
    #Rotation of O1 to xyplane
    atmO = xx[:,atm1,:]
    atmT = xx[:, atm2, :]
    lens = la.norm(atmO-atmT,axis=1)
    return lens

def anglePlot(cds,wts):
    ang = np.degrees(ba(cds,4-1,13-1,1-1))
    np.savetxt(coordinateSet+"avgAng",[np.average(ang,weights=wts)])
    theLen, xx = np.histogram(ang, bins=25, range=(80,200), density=True, weights=wts)  # WEIGHTS=WEIGHTARRAY
    binz = (xx[1:]+xx[:-1])/2
    np.savetxt(coordinateSet+"ang",zip(binz,theLen))

def bondPlot(cds,wts):
    OO = bL(cds,4-1,1-1)*angstr
    np.savetxt(coordinateSet+"avgooDist",[np.average(OO,weights=wts)])
    theLen, xx = np.histogram(OO, bins=25, range=(2.0, 3.2), density=True, weights=wts)  # WEIGHTS=WEIGHTARRAY
    binz = (xx[1:] + xx[:-1]) / 2
    np.savetxt(coordinateSet + "ooDist", zip(binz, theLen))

def bondAngPlot(cds,wts):
    ang = np.degrees(ba(cds,4-1,13-1,1-1))
    OO = bL(cds,4-1,1-1)*angstr
    theLen, xx,yy = np.histogram2d(OO,ang, bins=25, range=((2.0, 3.2),(80,200)), density=True, weights=wts)  # WEIGHTS=WEIGHTARRAY
    binzx = (xx[1:] + xx[:-1]) / 2
    binzy = (yy[1:] + yy[:-1]) / 2
    np.savetxt(coordinateSet+'2dOOvsOHO',theLen)
    np.savetxt(coordinateSet+'2dOOvsOHO_b',zip(binzx,binzy))

def plotStuff(symEckRotCoords):
    if os.path.isfile("q_" + coordinateSet + ".npy"):
        print 'plotQs'
        q = np.load("q_" + coordinateSet + ".npy")
        for i in range(q.shape[1]):
            PltHists1D('allH', q[:, i], (-100, 100), 'q_' + str(i), 'tetramerInternals/Probability Denisty',
                   False, symDw)
        for i in range(q.shape[1]):
            PltHists1D('allH', np.square(q[:, i]), (-100, 100), 'q2_' + str(i), 'tetramerInternals/Probability Denisty',
                   False, symDw)
        for i in range(q.shape[1]):
            PltHists1D('allH', np.power(q[:, i],4), (-1000, 1000), 'q4_' + str(i), 'tetramerInternals/Probability Denisty',
                   False, symDw)

    print 'INTERNAL COORDINATES :-O'
    internals = Wfn.molecule.SymInternalsH9O4plus(symEckRotCoords)

    nm=Wfn.molecule.internalName
    print internals
    x=np.average(internals,axis=0,weights=symDw)
    print x
    np.savetxt('averageInternalsWithNewEckart_'+coordinateSet,np.average(internals,axis=0,weights=symDw))
    print 'Internal coordinate shape: ', np.shape(internals)
    print 'One attribute shape: ',np.shape(internals[:,0])
    print 'number of dws: ', symDw
    symEckRotCoords*=angstr

    """self.internalName = ['xH11', 'yH11', 'zH11', 'xH12', 'yH12', 'zH12', 'xH13', 'yH13', 'zH13', 'theta651',
                         'phi651', 'Chi651',
                         'theta1039', 'phi1039', 'Chi1039', 'theta728', 'phi728', 'Chi728', 'rOH5', 'rOH6',
                         'HOH516', 'rOH7', 'rOH8', 'HOH728',
                         'rOH9', 'rOH10', 'HOH9310', 'rO1O2', 'rO1O3', 'rO2O3', 'xO4', 'yO4', 'zO4']"""
    #ZComps as a sanity check
    PltHists1D('allH', symEckRotCoords[:, 0, -1], (-2, 2), 'zo1', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 1, -1], (-2, 2), 'zo2', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 2, -1], (-2, 2), 'zo3', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 3, -1], (-2, 2), 'zo4', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 4, -1], (-2, 2), 'zh5', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 5, -1], (-2, 2), 'zh6', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 6, -1], (-2, 2), 'zh7', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 7, -1], (-2, 2), 'zh8', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 8, -1], (-2, 2), 'zh9', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 9, -1], (-2, 2), 'zh10', 'tetramerInternals/Probability Denisty',
               False, symDw)
    PltHists1D('allH', symEckRotCoords[:, 10, -1], (-2, 2), 'zh11', 'tetramerInternals/Probability Denisty',
               False, symDw)
    #rOH1,2,3
    PltHists1D('allH', internals[:, 0]*angstr, (0, 2), nm[0], 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 1]*angstr, (0, 2), nm[1], 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 2]*angstr, (0, 2), nm[2], 'tetramerInternals/Probability Density', False,
               symDw)
    #umbrella and dihedrals

    PltHists1D('allH', np.degrees(internals[:, 3]), (0, 180), nm[3], 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 4]), (0, 180), nm[4], 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 5]), (0, 180), nm[5], 'tetramerInternals/Probability Density', False,
               symDw)

    #Eulers for hydronium
    PltHists1D('allH',np.degrees(internals[:, 6]), (0,180), nm[6], 'tetramerInternals/Probability Density',
               False,
               symDw)
    # PltHists1D('allH', np.degrees(internals[:, 7]), (0,360), nm[7], 'tetramerInternals/Probability Density',
    #            False,
    #            symDw)
    # PltHists1D('allH', np.degrees(internals[:, 8]), (0, 360), nm[8], 'tetramerInternals/Probability Density',
    #            False,
    #            symDw)
    PltHists1D('allH', np.degrees(internals[:, 7]), (-180,180), nm[7], 'tetramerInternals/Probability Density',
               False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 8]), (-180, 180), nm[8], 'tetramerInternals/Probability Density',
               False,
               symDw)


    # PltHists1D('allH', internals[:, 6]*angstr, (-5, 5), nm[6], 'tetramerInternals/Probability Density',
    #            False,
    #            symDw)
    # PltHists1D('allH', internals[:, 7]*angstr, (-4, 4), nm[7], 'tetramerInternals/Probability Density',
    #            False,
    #            symDw)
    # PltHists1D('allH', internals[:, 8]*angstr, (-2, 2), nm[8], 'tetramerInternals/Probability Density',
    #            False,
    #            symDw)

    #eulers
    PltHists1D('allH', np.rad2deg(internals[:, 9]), (0, 360), nm[9],
               'tetramerInternals/Probability Density', False, symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 10]), (-360, 360), nm[10], 'tetramerInternals/Probability Density',
               False, symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 11]), (0, 360), nm[11],'tetramerInternals/Probability Density', False, symDw)


    PltHists1D('allH', np.rad2deg(internals[:, 12]), (0, 360), nm[12],
               'tetramerInternals/Probability Density', False, symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 13]), (-360, 360), nm[13],
               'tetramerInternals/Probability Density',
               False, symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 14]), (0, 360), nm[14],
               'tetramerInternals/Probability Density', False, symDw)

    PltHists1D('allH', np.rad2deg(internals[:, 15]), (0, 360), nm[15],
               'tetramerInternals/Probability Density', False, symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 16]), (-360, 360), nm[16],
               'tetramerInternals/Probability Density',
               False, symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 17]), (0, 360), nm[17],
               'tetramerInternals/Probability Density', False, symDw)



#roh&theta
    PltHists1D('allH', internals[:, 18]*angstr,(1,3), nm[18], 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 19]*angstr, (1,3), nm[19], 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 20]), (70,180), nm[20], 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 21]*angstr, (1,3), nm[21], 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 22]*angstr, (1,3), nm[22], 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.rad2deg(internals[:, 23]), (70,180), nm[23], 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 24]*angstr, (1, 3), nm[24], 'tetramerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 25]*angstr, (1, 3), nm[25], 'tetramerInternals/Probability Density', False,
               symDw)

    PltHists1D('allH', np.rad2deg(internals[:, 26]), (70, 180), nm[26], 'tetramerInternals/Probability Density',
               False,
               symDw)


    PltHists1D('allH', (internals[:, 27])*angstr, (3,6), nm[27], 'tetramerInternals/Probability Density',
               False,
               symDw)
    PltHists1D('allH', (internals[:, 28])*angstr, (3,6), nm[28], 'tetramerInternals/Probability Density',
               False,
               symDw)
    PltHists1D('allH', (internals[:, 29])*angstr, (3,6), nm[29], 'tetramerInternals/Probability Density',
               False,
               symDw)

    # PltHists1D('allH', np.degrees((internals[:, 30])), (0,360),nm[30],
    #            'tetramerInternals/Probability Density',
    #            False,
    #            symDw)
    # PltHists1D('allH', np.degrees((internals[:, 31])), (0,360),nm[31],
    #            'tetramerInternals/Probability Density',
    #            False,
    #            symDw)
    # PltHists1D('allH', np.degrees((internals[:, 32])), (0,360),nm[32],
    #            'tetramerInternals/Probability Density',
    #            False,
    #            symDw)

    PltHists1D('allH', (internals[:, 30])*angstr, (-3,3), nm[30],
               'tetramerInternals/Probability Density',
               False,
               symDw)
    PltHists1D('allH', (internals[:, 31]) * angstr, (-3, 3), nm[31],
               'tetramerInternals/Probability Density',
               False,
               symDw)
    PltHists1D('allH', (internals[:, 32])*angstr, (-3,3), nm[32],
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
    print 'npy file found'
    symCoords=np.load(cds+'.npy')
    if symCoords.shape[0]==13:
        symCoords = np.concatenate((symCoords[np.newaxis,:,:],symCoords[np.newaxis,:,:]))
    #symCoords=Wfn.molecule.rotateBackToFrame(symCoords,2,1,3)
    #np.save("refEck.npy",symCoords)
    symDw = np.load(cds+'_dw.npy')
else:
    symCoords, symDw = Wfn.loadCoords(cds)
    # if 'bad' in cds:
    #     symCoords/=angstr
    # if 'origMin' in cds:
    #     symCoords/=(angstr**2)
    # #symCoords=np.concatenate((symCoords,symCoords))/angstr
    # print symCoords
    # print symCoords.shape
    # if 'min' in cds:
    #     symDw = [1.,1.]
    #     symCoords=Wfn.molecule.rotateBackToFrame(symCoords,2,1,3)
    # # if 'refEck' in cds:
    # #     symDw = [1., 1.]
    # #     symCoords = np.array([Wfn.molecule.pullTetramerRefPos(),Wfn.molecule.pullTetramerRefPos()])*angstr
    # symDw=[1.,1.]
    if symCoords.shape[0]==13:
        symCoords = np.concatenate((symCoords, symCoords))
    symCoords = Wfn.molecule.rotateBackToFrame(symCoords, 2, 1, 3)
    #np.savetxt("eqRotated",symCoords[0])
    #np.save(cds,symCoords)
    #np.save(cds+'_dw',symDw)

# anglePlot(symCoords,symDw)
# bondPlot(symCoords,symDw)
# bondAngPlot(symCoords,symDw)
#
# stop


if 'input' in coordinateSet:
    np.save("rotated"+coordinateSet+".npy",symCoords)
    np.save("rotated"+coordinateSet+"_dw.npy",symDw)
    print 'input file - writing'
    wf = open("../coordinates/tetramer/rotated_"+coordinateSet[-4:],'w+')
    trim = ["O", "O", "O", "O","H","H","H", "H", "H", "H", "H", "H", "H"]
    for wI, walker in enumerate(symCoords):
        wf.write("13\n")
        wf.write("%5.12f 0.0 0.0 0.0 0.0\n" % symDw[wI])
        for aI, atm in enumerate(walker):
            wf.write("%s %5.12f %5.12f %5.12f\n" % (trim[aI], atm[0], atm[1], atm[2]))
        wf.write("\n")
    wf.close()
    stop
elif 'RSwapped' in coordinateSet:
    print 'rotated and symmetrized file - rotatedagain'
    wf = open("../coordinates/tetramer/final_" + coordinateSet[-4:], 'w+')
    trim = ["O", "O", "O", "O","H","H","H", "H", "H", "H", "H", "H", "H"]
    np.save("../coordinates/tetramer/iSymtet"+coordinateSet[-4:],symCoords)
    np.save("../coordinates/tetramer/iSymtet"+coordinateSet[-4:]+"_dw",symDw)
    print 'npy saved finalcds'
    # for wI, walker in enumerate(symCoords):
    #     wf.write("13\n")
    #     wf.write("%5.12f\n" % symDw[wI])
    #     for aI, atm in enumerate(walker):
    #         wf.write("%s %5.12f %5.12f %5.12f\n" % (trim[aI], atm[0], atm[1], atm[2]))
    #     wf.write("\n")
    # wf.close()
    stop

print 'Symcoords shape',symCoords.shape
print 'Got symCoords!'
#print symCoords
print 'NUMBER OF WALKERS IN allH: ',symCoords.shape[0]
symEckRotCoords = symCoords
iwantToPlotStuff=True
path='../spectra/'
if iwantToPlotStuff:
    plotStuff(symEckRotCoords)
else:
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

    print 'PEDIP shape', pdip.shape
    pe = pdip[:, 0]
    dip = pdip[:, 1:]
    print 'Shape of dipole: ', np.shape(dip)
    HOASpectrum=CalculateSpectrum.HarmonicApproxSpectrum(Wfn,symEckRotCoords,symDw,path,testName)
    if 'Eck' in GfileName:
        coordinateSet=coordinateSet+'refGmat'
    fundamentalEnergies,fundamentalIntensities, combinationBandEnergies,combinationBandIntensities=HOASpectrum.calculateSpectrum(symEckRotCoords,symDw,GfileName,pe,dip,coordinateSet,testName,kill,eckt,dipPath)
    fundamentalFile=open('../spectra/Fundamentals_'+coordinateSet+testName+kill,'w')
    for i,(v,intensity) in enumerate(zip(fundamentalEnergies,fundamentalIntensities)):
        fundamentalFile.write(str(i)+"       "+str(v)+"   "+str(intensity)+"\n")
    fundamentalFile.close()
