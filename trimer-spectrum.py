import numpy as np
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import DMCClusters as dmc
import time
import glob
import Plot
import CalculateSpectrum
import sys, os
angstr = 0.529177

def PltHists1D(cfg, thing, bound, xl, yl, overly, weits):
    theLen, xx = np.histogram(thing, bins=100, range=bound, normed=True, weights=weits)  # WEIGHTS=WEIGHTARRAY
    inin = True
    overlay = False
    bnd = str(bound[0]).replace(".", "").replace("-", "") + str(bound[1]).replace(".", "").replace("-", "")
    print bnd
    mP = Plot.myPlot(cfg, '1d', bnd, bnd, xl, yl, overly, inin, theLen,(xx[1:]+xx[:-1])/2)
    mP.plotIt()

def PltHists2D(cfg,thing1,thing2,bound1,bound2,xl,yl,overly,weits,bins):
    white=False
    theLen, xx,yy = np.histogram2d(thing1,thing2,bins=bins,range=(bound1,bound2),normed=True,weights=weits) #WEIGHTS=WEIGHTARRAY
    inin=True
    overlay=False
    X = (xx[1:] + xx[:-1]) / 2
    Y = (yy[1:] + yy[:-1]) / 2
    bnd1 = str(bound1[0]).replace(".","").replace("-","")+str(bound1[1]).replace(".","").replace("-","")
    bnd2 = str(bound2[0]).replace(".","").replace("-","")+str(bound2[1]).replace(".","").replace("-","")
    mP = Plot.myPlot(cfg,'2d',bnd1,bnd2,yl,xl,overly,inin,theLen,X,Y,white)
    mP.plotIt()

def plotStuff(symEckRotCoords):
    import matplotlib.cm as cm

    #waterfall Luke plot, histogram rOH based on Bend angle
    # rOH = Wfn.molecule.bL(symEckRotCoords,9-1,3-1)
    # H3Obend=Wfn.molecule.ba(symEckRotCoords,10-1,3-1,9-1)
    # H3Obend = np.degrees(H3Obend)
    # rOH_eig = Wfn.molecule.bL(symEckRotCoords,10-1,3-1)
    # rOO_eig = Wfn.molecule.bL(symEckRotCoords,2-1,3-1)
    # zcompH = symEckRotCoords[:,9-1,-1]
    # rOH = rOH * angstr
    # edges = np.linspace(70,170,num=11)
    # edges = np.linspace(-1.0,1.0,num=11)
    # edges = np.linspace(0,num=11)
    #
    # edgz=(edges[1:] + edges[:-1])*0.5
    # colors = iter(cm.rainbow(np.linspace(0, 1, 10)))
    # for i in range(len(edges)-1):
    #     qual = (H3Obend > edges[i]) * (H3Obend < edges[i + 1])
    ##     qual = (zcompH > edges[i]) * (zcompH < edges[i + 1])
        #
        # idx = np.where(qual)[0]
        # yar, xx = np.histogram(rOH[idx], bins=15, range=(0.6, 1.5),weights=symDw[idx],
        #                       density=True)  # histogramming OH distances in a range of OO distances
        # yar /= 1e5
        # X = (xx[:-1]+xx[1:]) * 0.5
        # tc = next(colors)
        # plt.fill(X, yar, color=tc, alpha=0.7)
        # plt.plot(X, yar, color='k')
    ## plt.yticks(edgz)
    # plt.savefig("rOH_vsIHBBend.png",dpi=500)
    # plt.close()





    #rotations annd bend and wag
    # roh8, roh9, roh10, thH8,thH9,thH10,phH8,phH9,phH10 = Wfn.molecule.sphericalTrimer(symEckRotCoords)
    # PltHists1D('allH', np.degrees(phH8+phH9+phH10), (-180, 180), "Phi8+Phi9+Phi10", 'pubtrimerInternals/Probability Density',
    #            False,symDw)
    # PltHists1D('allH', np.degrees(-phH8 + phH9 + phH10), (-180, 180), "-Phi8+Phi9+Phi10",
    #            'pubtrimerInternals/Probability Density',
    #            False, symDw)
    # PltHists1D('allH', np.degrees(-phH9 + phH10), (180,360), "-Phi9+Phi10",
    #            'pubtrimerInternals/Probability Density',
    #            False, symDw)

    # strang="q_" + coordinateSet +'_'+testName+"_"+kill+".npy"
    # if os.path.isfile("q_" + coordinateSet +'_'+testName+"_"+kill+".npy"):
    #     print 'plotQs'
    #     q = np.load("q_" + coordinateSet +'_'+testName+"_"+kill+".npy")
        # for i in range(q.shape[1]):
        #     PltHists1D('allH', q[:, i], (-150,150), 'q_' + str(i), 'trimerInternals/Probability Denisty',False, symDw)
        # for j in range(q.shape[1]):
        #     PltHists2D('allH', q[:, 10], q[:, j], (-50, 50), (-50, 50), 'q10', 'q'+str(j),
        #                False, symDw, 30)
        # for k in range(q.shape[1]):
        #     PltHists2D('allH', q[:, 12], q[:, k], (-50, 50), (-50, 50), 'q12', 'q'+str(k),
        #                False, symDw, 30)
    ##
    # qp = np.copy(q)
    # waga = qp[:,10]
    # for f in range(24):
    #     q[:,f]/=np.amax(q[:,f])
    # wag = q[:,10]
    # edges = np.linspace(-1.0,1.0, num=16)
    # # edges = np.append(edgesb, 1.0)
    # for mode in range(q.shape[1]):
    #     choppyList = [[] for x in range(len(edges))]
    #     cLIndex = [[] for x in range(len(edges))]
    #     yar = np.zeros((len(edges), 15))
    #     colors = iter(cm.rainbow(np.linspace(0, 1, 16)))
    #     dx = 0.1
    #     dy = 1
    #     ################
    #     avg = np.zeros(len(edges)-1)
    #     for i in (range(len(edges)-1)):
    #         mde = q[:,mode]
    #         qual = (mde > edges[i]) * (mde < edges[i + 1])
    #         idx = np.where(qual)[0]
    #         if len(idx) != 0:
    #             avg[i] = np.average(waga[idx]**2 * waga[idx]**2 ,weights=symDw[idx])*(np.sum(symDw[idx]))
    #             # avg[i] = np.sum(symDw[idx])
    #         else:
    #             # avg[i] = 0.0
    #             avg[i] = 0.0
    #     X = (edges[:-1]+edges[1:])*0.5
    #     plt.plot(X,avg)
    #     plt.xlabel("q_"+str(mode))
    #     plt.ylabel("<q_wag^2> * weights")
    #     plt.savefig("WagPlot_weightedWeight_Mode_"+str(mode))
    #     plt.close()
    #     ###############
    #     for i in reversed(range(len(edges))):
    #         mde = q[:,mode]
    #         qual = (mde < edges[i]) * (mde > edges[i - 1])
    #         cLIndex[i] = np.where(qual)[0]
    #         choppyList[i] = mde[cLIndex[i]]
    #         yar[i], xx = np.histogram(wag[cLIndex[i]], bins=15, range=(-1.0, 1.0),weights=symDw[cLIndex[i]],
    #                                   density=True)  # histogramming OH distances in a range of OO distances
    #         X = (xx[:-1]+xx[1:]) * 0.5
    #         tc = next(colors)
    #         # plt.fill(X + i * dx, yar[i] + i * dy + 1, color=tc, alpha=0.7)
    #         # plt.plot(X + i * dx, yar[i] + i * dy + 1, color='k')
    #         plt.fill(X, yar[i] + i * dy + 1, color=tc, alpha=0.7)
    #         plt.plot(X, yar[i] + i * dy + 1, color='k')
    #     plt.savefig("cascade_Wag_Mode_"+str(mode))
    #     plt.close()
    print 'INTERNAL COORDINATES :-O'
    # stop
    #theta: align the z axes
    #phi: align the two x axes
    #chi: bring to reference structure(?)
    #Applied Chi -> Theta -> Phi(vector)
    #or Phi -> Theta -> Chi(vector) if applying to xyz as opposed to XYZ
    internals = Wfn.molecule.SymInternals(symEckRotCoords)
    nm=Wfn.molecule.internalName
    nm = [nm[g].lower() for g in range(len(nm))]
    units = []
    print(units)
    for name in nm:
        if 'ph' in name or 'th' in name or 'chi' in name or 'xi' in name or 'hoh' in name:
            units.append('degrees')
        elif 'xh' in name or 'yh' in name or 'zh' in name:
            units.append('pmAngstroms')
        else:
            units.append('angstroms')
    x=np.average(internals,axis=0,weights=symDw)
    print('average',x)
    for qi in range(24):
        if units[qi] == 'degrees':
            PltHists1D('allH', np.degrees(internals[:, qi]), (-180,180), nm[qi], 'pubtrimerInternals/Probability Density', False,
                       symDw)
        elif units[qi] == 'pmAngstroms':
            PltHists1D('allH', internals[:, qi] * angstr, (-2, 2), nm[qi], 'pubtrimerInternals/Probability Density', False,
                       symDw)
        else:
            PltHists1D('allH', internals[:, qi] * angstr, (0, 3), nm[qi], 'pubtrimerInternals/Probability Density', False,
                       symDw)


    PltHists2D('allH', np.degrees(internals[:,1]),np.degrees(internals[:,2]),(0,180),(-90,90), 'Theta', 'Phi',
               False, symDw, 40)
    np.savetxt('averageInternalsWithNewEckart_'+coordinateSet+'_'+testName+"_"+kill,np.average(internals,axis=0,weights=symDw))
    
    print 'Internal coordinate shape: ', np.shape(internals)
    print 'One attribute shape: ',np.shape(internals[:,0])
    print 'number of dws: ', symDw.shape

    # PltHists2D('allH',angstr*(internals[:, 21])+angstr*(internals[:,22]) ,  angstr*(internals[:, 21])-angstr*(internals[:,22]), (4.4,5.6),(-0.6,0.6),'Roo_1+Roo2', 'Roo1-Roo_2',
    #            False, symDw, 40)
    PltHists2D('allH', np.degrees(internals[:, 5]), np.degrees(internals[:, 8]), (-180, 180), (-180, 180), 'PhiH9',
               'PhiH10',False, symDw, 60)

    PltHists2D('allH', np.degrees(internals[:, 5]), np.degrees(internals[:, 4]), (-180, 180), (-180, 180), 'PhiH9',
               'ThH10', False, symDw, 60)

    #            False, symDw, 60)
    PltHists2D('allH', np.degrees(internals[:, 2]), np.degrees(internals[:, 9]), (-180,180), (-180,180), 'phiH8', 'Th627',
               False, symDw, 60)
    PltHists2D('allH', np.degrees(internals[:, 2]), np.degrees(internals[:, 10]), (-180,180), (-180,180), 'phiH8', 'Ph627',
               False, symDw, 60)
    PltHists2D('allH', np.degrees(internals[:, 2]), np.degrees(internals[:, 11]), (-180,180), (-180,180), 'phiH8', 'Chi627',
               False, symDw, 60)
    PltHists2D('allH', np.degrees(internals[:, 2]), np.degrees(internals[:, 12]), (-180,180), (-180,180), 'phiH8', 'Th514',
               False, symDw, 60)
    PltHists2D('allH', np.degrees(internals[:, 2]), np.degrees(internals[:, 13]), (-180,180), (-180,180), 'phiH8', 'Ph514',
               False, symDw, 60)
    PltHists2D('allH', np.degrees(internals[:, 2]), np.degrees(internals[:, 14]), (-180,180), (-180,180), 'phiH8', 'Chi514',
               False, symDw, 60)
    
    PltHists2D('allH', np.degrees(internals[:, 5]), np.degrees(internals[:, 9]), (-180,180), (-180,180), 'phiH9', 'Th627',
               False, symDw, 60)
    PltHists2D('allH', np.degrees(internals[:, 5]), np.degrees(internals[:, 10]), (-180,180), (-180,180), 'phiH9', 'Ph627',
               False, symDw, 60)
    PltHists2D('allH', np.degrees(internals[:, 5]), np.degrees(internals[:, 11]), (-180,180), (-180,180), 'phiH9', 'Chi627',
               False, symDw, 60)
    PltHists2D('allH', np.degrees(internals[:, 5]), np.degrees(internals[:, 12]), (-180,180), (-180,180), 'phiH9', 'Th514',
               False, symDw, 60)
    PltHists2D('allH', np.degrees(internals[:, 5]), np.degrees(internals[:, 13]), (-180,180), (-180,180), 'phiH9', 'Ph514',
               False, symDw, 60)
    PltHists2D('allH', np.degrees(internals[:, 5]), np.degrees(internals[:, 14]), (-180,180), (-180,180), 'phiH9', 'Chi514',
               False, symDw, 60)

    PltHists2D('allH', np.degrees(internals[:, 1]), np.degrees(internals[:, 9]), (-180,180), (-180,180), 'thH8', 'Th627',
               False, symDw, 60)
    PltHists2D('allH', np.degrees(internals[:, 1]), np.degrees(internals[:, 10]), (-180,180), (-180,180), 'thH8', 'Ph627',
               False, symDw, 60)
    PltHists2D('allH', np.degrees(internals[:, 1]), np.degrees(internals[:, 11]), (-180,180), (-180,180), 'thH8', 'Chi627',
               False, symDw, 60)
    PltHists2D('allH', np.degrees(internals[:, 1]), np.degrees(internals[:, 12]), (-180,180), (-180,180), 'thH8', 'Th514',
               False, symDw, 60)
    PltHists2D('allH', np.degrees(internals[:, 1]), np.degrees(internals[:, 13]), (-180,180), (-180,180), 'thH8', 'Ph514',
               False, symDw, 60)
    PltHists2D('allH', np.degrees(internals[:, 1]), np.degrees(internals[:, 14]), (-180,180), (-180,180), 'thH8', 'Chi514',
               False, symDw, 60)
    stop
    internalName = ['rOH8', 'thH8', 'phiH8','rOH9', 'thH9', 'phiH9','rOH10', 'thH10', 'phiH10',
                             'th_627', 'phi_627','xi_627', 'th_514', 'phi_514', 'xi_514', 'rOH_41',
                             'rOH_51', 'aHOH_451', 'rOH_26','rOH_27', 'aHOH_267','rOO_1', 'rOO_2', 'aOOO']

    PltHists1D('allH', internals[:, 0]*angstr, (0, 3), nm[0], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 1]), (0,180), nm[1], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 2]), (-360,360), nm[2], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 3]*angstr, (0, 3), nm[3], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 4]), (0,180), nm[4], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 5]), (-90,90), nm[5], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 6]*angstr, (0, 3), nm[6], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 7]), (0,180), nm[7], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 8]), (-90,90), nm[8], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 9]), (-360,360), nm[9], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 10]), (-180,180), nm[10], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 11]), (-360,360), nm[11], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 12]), (-360,360), nm[12], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 13]), (-180,180), nm[13], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', np.degrees(internals[:, 14]), (-360,360), nm[14], 'trimerInternals/Probability Density', False,
               symDw)
    PltHists1D('allH', internals[:, 15]*angstr, (0.6, 1.4), nm[15], 'trimerInternals/Probability Density', False,
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
    Wfn.setIsotope('notDeuteratedOnce_hydronium')
elif '1De' in coordinateSet:
    Wfn.setIsotope('DeuteratedOnce_eigen')
elif '1Dw' in coordinateSet:
    Wfn.setIsotope('DeuteratedOnce_fw')
elif '1Dh' in coordinateSet:
    Wfn.setIsotope('DeuteratedOnce_hydronium')
elif 'test' in coordinateSet:
    print 'test,fam'
# else:
#     poop

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

# lukeOOs = 3
# for f in range(1,lukeOOs+1):
#     a = np.load("dvr_oo_geoms_shift_"+str(f)+".npy")
#     fll = open('lukeGeomz_'+str(f)+'.xyz', 'w+')
#     Wfn.molecule.printCoordsToFile(a,fll)
# stop
# fll = open('cdsWagProblem.xyz','w+')
# Wfn.molecule.printCoordsToFile(a,fll)
# symDw = np.load("../cdsWagProblem_dw.npy")
# plotStuff(a)
# stop

if os.path.isfile(cds+'.npy'):
    symCoords=np.load(cds+'.npy')
    symDw = np.load(cds+'_dw.npy')
    # symCoords = Wfn.molecule.rotateBackToFrame(symCoords, 3, 2, 1)
    # np.save('../coordinates/trimer/'+'Tfinal_allH'+'.npy', symCoords)
    # np.save('../coordinates/trimer/' + coordinateSet + '_dw.npy', symDw)
    # print 'asdf'
else:
    symCoords, symDw = Wfn.loadCoords(cds)
    np.save('../coordinates/trimer/'+coordinateSet+'.npy', symCoords)
    np.save('../coordinates/trimer/'+coordinateSet+'_dw.npy', symDw)
    # if symCoords.shape[0] == 0:
    #     symCoords = Wfn.molecule.rotateBackToFrame(np.array([symCoords,symCoords]), 3, 2, 1)[0]
    # else:
    #     symCoords = Wfn.molecule.rotateBackToFrame(symCoords, 3, 2, 1)

if 'inpt' in coordinateSet:
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
    np.save("../coordinates/trimer/Tfinal_"+coordinateSet[-4:],symCoords)
    np.save("../coordinates/trimer/Tfinal_"+coordinateSet[-4:]+"_dw",symDw)
    print 'npy saved finalcds'
else:
    print 'Symcoords shape',symCoords.shape
    print 'Got symCoords!'
    #print symCoords
    print 'NUMBER OF WALKERS IN allH: ',symCoords.shape[0]
    symEckRotCoords = symCoords
    # symEckRotCoords = symCoords[:len(symCoords)/2]
    # symDw = symDw[:len(symCoords)/2]
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
