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

def PltHists2D(cfg,thing1,thing2,bound1,bound2,yl,xl,overly,weits,bins):
    white=False
    theLen, xx,yy = np.histogram2d(thing1,thing2,bins=bins,range=(bound1,bound2),normed=True,weights=weits) #WEIGHTS=WEIGHTARRAY
    inin=True
    overlay=False
    X = (xx[1:] + xx[:-1]) / 2
    Y = (yy[1:] + yy[:-1]) / 2
    bnd1 = str(bound1[0]).replace(".","").replace("-","")+str(bound1[1]).replace(".","").replace("-","")
    bnd2 = str(bound2[0]).replace(".","").replace("-","")+str(bound2[1]).replace(".","").replace("-","")
    mP = Plot.myPlot(cfg,'2d',bnd1,bnd2,xl,yl,overly,inin,theLen,X,Y,white)
    mP.plotIt()

def plotStuff(symEckRotCoords):
    # if os.path.isfile("q_" + coordinateSet + ".npy"):
    #     print 'plotQs'
    #     q = np.load("q_" + coordinateSet + ".npy")
    #     for i in range(q.shape[1]):
    #         PltHists1D('allH', q[:, i], (-100, 100), 'q_' + str(i), 'trimerInternals/Probability Denisty',
    #                False, symDw)
    print 'INTERNAL COORDINATES :-O'

    internals = Wfn.molecule.SymInternalsH7O3plus(symEckRotCoords)
    nm=Wfn.molecule.internalName
    print internals
    x=np.average(internals,axis=0,weights=symDw)
    print x
    np.savetxt('averageInternalsWithNewEckart_'+coordinateSet,np.average(internals,axis=0,weights=symDw))
    print 'Internal coordinate shape: ', np.shape(internals)
    print 'One attribute shape: ',np.shape(internals[:,0])
    print 'number of dws: ', symDw

    # symEckRotCoords*=angstr
    #ZComps as a sanity check
    # PltHists1D('allH', symEckRotCoords[:, 0, -1], (-2, 2), 'zo1', 'trimerInternals/Probability Denisty',
    #            False, symDw)
    # PltHists1D('allH', symEckRotCoords[:, 1, -1], (-2, 2), 'zo2', 'trimerInternals/Probability Denisty',
    #            False, symDw)
    # PltHists1D('allH', symEckRotCoords[:, 2, -1], (-2, 2), 'zo3', 'trimerInternals/Probability Denisty',
    #            False, symDw)
    # PltHists1D('allH', symEckRotCoords[:, 3, -1], (-2, 2), 'zo4', 'trimerInternals/Probability Denisty',
    #            False, symDw)
    # PltHists1D('allH', symEckRotCoords[:, 4, -1], (-2, 2), 'zh5', 'trimerInternals/Probability Denisty',
    #            False, symDw)
    # PltHists1D('allH', symEckRotCoords[:, 5, -1], (-2, 2), 'zh6', 'trimerInternals/Probability Denisty',
    #            False, symDw)
    # PltHists1D('allH', symEckRotCoords[:, 6, -1], (-2, 2), 'zh7', 'trimerInternals/Probability Denisty',
    #            False, symDw)
    # PltHists1D('allH', symEckRotCoords[:, 7, -1], (-2, 2), 'zh8', 'trimerInternals/Probability Denisty',
    #            False, symDw)
    # PltHists1D('allH', symEckRotCoords[:, 8, -1], (-2, 2), 'zh9', 'trimerInternals/Probability Denisty',
    #            False, symDw)
    # PltHists1D('allH', symEckRotCoords[:, 9, -1], (-2, 2), 'zh10', 'trimerInternals/Probability Denisty',
    #            False, symDw)

    # internalName = ['rOH9', 'rOH10', 'spHOH', 'rH8', 'thH8', 'phiH8', 'thH', 'phiH', 'xiH', 'th_627', 'phi_627',
    #                      'xi_627', 'th_514', 'phi_514', 'xi_514', 'rOH_41', 'rOH_51', 'aHOH_451', 'rOH_26',
    #                      'rOH_27', 'aHOH_267',
    #                      'rOO_1', 'rOO_2', 'aOOO']
    #
    # PltHists1D('allH', internals[:, 0]*angstr, (-2, 2), nm[0], 'trimerInternals/Probability Density', False,
    #            symDw)
    # PltHists1D('allH', internals[:, 1]*angstr, (-2, 2), nm[1], 'trimerInternals/Probability Density', False,
    #            symDw)
    # PltHists1D('allH', np.degrees(internals[:, 2]), (0,200), nm[2], 'trimerInternals/Probability Density', False,
    #            symDw)
    #
    # PltHists1D('allH', internals[:, 3]*angstr, (0,2), nm[3], 'trimerInternals/Probability Density', False,
    #            symDw)
    # PltHists1D('allH', np.degrees(internals[:, 4]), (0,180), nm[4], 'trimerInternals/Probability Density', False,
    #            symDw)
    # PltHists1D('allH', np.degrees(internals[:, 5]), (0,360), nm[5], 'trimerInternals/Probability Density', False,
    #            symDw)
    #
    # PltHists1D('allH', np.degrees(internals[:, 6]), (0,180), nm[6], 'trimerInternals/Probability Density',
    #            False,
    #            symDw)
    # print 'maxmin theta', np.degrees(np.amin(internals[:, 6])),np.degrees(np.amax(internals[:, 6]))
    # # PltHists1D('allH', np.rad2deg(internals[:, 7]), (-180,180), nm[7], 'trimerInternals/Probability90180Density',
    # #            False,
    # #            symDw)
    # # PltHists1D('allH', np.rad2deg(internals[:, 8]), (-180,180), nm[8], 'trimerInternals/Probability Density',
    # #            False,
    # #            symDw)
    # PltHists1D('allH', np.rad2deg(internals[:, 7]), (0,360), nm[7], 'trimerInternals/Probability90180Density',
    #            False,
    #            symDw)
    # print 'maxmin phi', np.degrees(np.amin(internals[:, 7])),np.degrees(np.amax(internals[:, 7]))
    #
    # PltHists1D('allH', np.rad2deg(internals[:, 8]), (0,360), nm[8], 'trimerInternals/Probability Density',
    #            False,
    #            symDw)
    # print 'maxmin xi', np.degrees(np.amin(internals[:, 8])),np.degrees(np.amax(internals[:, 8]))
    #
    # PltHists1D('allH', np.rad2deg(internals[:, 9]), (0, 360), nm[9],
    #            'trimerInternals/Probability Density', False, symDw)
    # PltHists1D('allH', np.rad2deg(internals[:, 10]), (-360,360), nm[10], 'trimerInternals/Probability Density',
    #            False, symDw)
    # PltHists1D('allH', np.rad2deg(internals[:, 11]), (0, 360), nm[11],'trimerInternals/Probability Density', False, symDw)
    #
    #
    # PltHists1D('allH', np.rad2deg(internals[:, 12]), (0, 360), nm[12],
    #            'trimerInternals/Probability Density', False, symDw)
    # PltHists1D('allH', np.rad2deg(internals[:, 13]), (-360,360), nm[13],
    #            'trimerInternals/Probability Density',
    #            False, symDw)
    # PltHists1D('allH', np.rad2deg(internals[:, 14]), (-360, 360), nm[14],
    #            'trimerInternals/Probability Density', False, symDw)
    #
    # PltHists1D('allH', np.rad2deg(internals[:, 15]), (0, 360), nm[15],
    #            'trimerInternals/Probability Density', False, symDw)
    # PltHists1D('allH', np.rad2deg(internals[:, 16]), (0, 360), nm[16],
    #            'trimerInternals/Probability Density',
    #            False, symDw)
    # PltHists1D('allH', np.rad2deg(internals[:, 17]), (0, 360), nm[17],
    #            'trimerInternals/Probability Density', False, symDw)
    #
    # PltHists1D('allH', internals[:, 18]*angstr,(0,3), nm[18], 'trimerInternals/Probability Density', False,
    #            symDw)
    # PltHists1D('allH', internals[:, 19]*angstr, (0,3), nm[19], 'trimerInternals/Probability Density', False,
    #            symDw)
    # PltHists1D('allH', np.rad2deg(internals[:, 20]), (70,180), nm[20], 'trimerInternals/Probability Density', False,
    #            symDw)
    # PltHists1D('allH', internals[:, 21]*angstr, (0,3), nm[21], 'trimerInternals/Probability Density', False,
    #            symDw)
    # PltHists1D('allH', internals[:, 22]*angstr, (0,3), nm[22], 'trimerInternals/Probability Density', False,
    #            symDw)
    # PltHists1D('allH', np.rad2deg(internals[:, 23]), (70,180), nm[23], 'trimerInternals/Probability Density', False,
    #            symDw)
    # PltHists2D('allH',np.degrees(internals[:,14]),np.degrees(internals[:,11]),(-180,180),(-180,180),'Xi1','Xi2',False,symDw,30)
    #
    # # PltHists2D('allH',np.degrees(internals[:,7]),np.degrees(internals[:,8]),(-180,180),(-180,180),'PhiH','XiH',False,symDw,30)
    # PltHists2D('allH',np.degrees(internals[:,7]),np.degrees(internals[:,8]),(0,360),(0,360),'PhiH','XiH',
    #            False,symDw,30)
    # PltHists2D('allH', np.degrees(internals[:, 6]), np.degrees(internals[:, 7]), (0, 360), (0, 360), 'ThH', 'PhiH',
    #            False, symDw, 30)
    # PltHists2D('allH', np.degrees(internals[:, 6]), np.degrees(internals[:, 8]), (0, 360), (0, 360), 'ThH', 'XiH',
    #            False, symDw, 30)


    internalName = ['rOH8', 'thH8', 'phiH8','rOH9', 'thH9', 'phiH9','rOH10', 'thH10', 'phiH10',
                             'th_627', 'phi_627','xi_627', 'th_514', 'phi_514', 'xi_514', 'rOH_41',
                             'rOH_51', 'aHOH_451', 'rOH_26','rOH_27', 'aHOH_267','rOO_1', 'rOO_2', 'aOOO']

    PltHists1D('allH', internals[:, 0]*angstr, (0, 3), nm[0], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 1]), (0,180), nm[1], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 2]), (0,360), nm[2], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 3]*angstr, (0, 3), nm[3], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 4]), (0,180), nm[4], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 5]), (-180,180), nm[5], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 6]*angstr, (0, 3), nm[6], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 7]), (0,180), nm[7], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 8]), (-180,180), nm[8], 'trimerInternals/Probability Density', False,
               symDw)




# H E R M I T E  P O L Y N O M I A L  A P P R O X I M A T I O N
au2wn=219474.63
nBins=51


starttime=time.time()

stateGround='stateGround'
state='stateGround'
DWstate='DWGround'
molecule='H7O3'
dTau=10

coordinateSet = sys.argv[1]
testName = sys.argv[2]
kill = sys.argv[3]

print 'cd set',coordinateSet
print 'pmz',testName
print 'kill',kill

Wfn=dmc.wavefunction('H7O3+', 1) #define our wavefunction
if 'allH' in coordinateSet:
    pass
elif 'allD' in coordinateSet:
    Wfn.setIsotope('fullyDeuterated')
elif '1He' in coordinateSet:
    Wfn.setIsotope('notDeuteratedOnce_eigen')
elif '1Hw' in coordinateSet:
    Wfn.setIsotope('notDeuteratedOnce_fw')
elif '1Hh' in coordinateSet:
    Wfn.setIsotope('DeuteratedOnce_hydronium')
elif '1De' in coordinateSet:
    Wfn.setIsotope('DeuteratedOnce_eigen')
elif '1Dw' in coordinateSet:
    Wfn.setIsotope('DeuteratedOnce_fw')
elif '1Dh' in coordinateSet:
    Wfn.setIsotope('DeuteratedOnce_hydronium')
elif 'test' in coordinateSet:
    print 'test,fam'
else:
    raise Exception

head = '/home/netid.washington.edu/rjdiri/'

cdsPath= '../coordinates/trimer/'
gPath = '../Gmats/trimer/'
dipPath='../newDipolesAndPE/trimer/'
cds =cdsPath+coordinateSet
GfileName = gPath+coordinateSet+'_'+testName+"_"+kill+'.gmat'
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
    np.save('../coordinates/trimer/'+coordinateSet+'.npy', symCoords)
    np.save('../coordinates/trimer/'+coordinateSet+'_dw.npy', symDw)
    if symCoords.shape[0] == 0:
        symCoords = Wfn.molecule.rotateBackToFrame(np.array([symCoords,symCoords]), 3, 2, 1)[0]
    else:
        symCoords = Wfn.molecule.rotateBackToFrame(symCoords, 3, 2, 1)

if 'input' in coordinateSet:
    print 'input file - rotating'
    wf = open("../coordinates/trimer/rotated_"+coordinateSet[-4:],'w+')
    trim = ["O","O","O","H","H","H","H","H","H","H"]
    for wI, walker in enumerate(symCoords):
        wf.write("13 0\n")
        wf.write("%5.12f\n" % symDw[wI])
        for aI, atm in enumerate(walker):
            wf.write("%s %5.12f %5.12f %5.12f\n" % (trim[aI], atm[0], atm[1], atm[2]))
        wf.write("\n")
    wf.close()

elif 'Rotated' in coordinateSet:
    print 'rotated and symmetrized file - rotatedagain'
    # wf = open("../coordinates/trimer/final_" + coordinateSet[-4:], 'w+')
    # trim = ["O", "O", "O", "H", "H", "H", "H", "H", "H", "H"]
    np.save("../coordinates/trimer/final_"+coordinateSet[-4:],symCoords)
    np.save("../coordinates/trimer/final_"+coordinateSet[-4:]+"_dw",symDw)
    print 'npy saved finalcds'
else:
    print 'Symcoords shape',symCoords.shape
    print 'Got symCoords!'
    #print symCoords
    print 'NUMBER OF WALKERS IN allH: ',symCoords.shape[0]
    symEckRotCoords = symCoords
    iwantToPlotStuff=False
    if iwantToPlotStuff:
        plotStuff(symEckRotCoords)
    else:
        eckt=False
        if os.path.isfile(dipPath+'eng_dip_'+coordinateSet+'_eckart.npy'):
            pdip = np.load(dipPath+'eng_dip_'+coordinateSet+'_eckart.npy')
            eckt=True
        elif os.path.isfile(dipF[:-3] + 'npy'):
            pdip = np.load(dipF[:-3] + 'npy')
        else:
            print 'not a npy file'
            pdip = np.loadtxt(dipF)
            np.save(dipF[:-3] + 'npy', pdip)
        path='../spectra/'
        print 'PEDIP shape', pdip.shape
        pe = pdip[:, 0]
        dip = pdip[:, 1:]
        print 'Shape of dipole: ', np.shape(dip)
        HOASpectrum=CalculateSpectrum.HarmonicApproxSpectrum(Wfn,symEckRotCoords,symDw,path,testName)
        # if 'Eck' in GfileName:
        #     coordinateSet=coordinateSet+'refGmat'
        fundamentalEnergies,fundamentalIntensities, combinationBandEnergies,combinationBandIntensities=HOASpectrum.calculateSpectrum(symEckRotCoords,symDw,GfileName,pe,dip,coordinateSet,testName,kill,eckt,dipPath)
        fundamentalFile=open('../spectra/Fundamentals_'+coordinateSet+testName+kill,'w')
        for i,(v,intensity) in enumerate(zip(fundamentalEnergies,fundamentalIntensities)):
            fundamentalFile.write(str(i)+"       "+str(v)+"   "+str(intensity)+"\n")
        fundamentalFile.close()
