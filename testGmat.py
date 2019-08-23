import numpy as np
import matplotlib.pyplot as plt
import DMCClusters as dmc
import CalculateSpectrum




Wfn=dmc.wavefunction('H7O3+', 1) #define our wavefunction
path = '../spectra'
testName = 'gmatTestStuff'
cds = np.array(
            [
                 [-2.34906009e+00,  4.06869143e+00,  0.00000000e+00],
                 [ 4.69812018e+00,  0.00000000e+00,  0.00000000e+00],
                 [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                 [-2.88862583e+00,  5.00324669e+00,  1.47532198e+00],
                 [-2.88862583e+00,  5.00324669e+00, -1.47532198e+00],
                 [ 5.77725164e+00, -2.46900000e-09,  1.47532198e+00],
                 [ 5.77725164e+00, -2.46900000e-09, -1.47532198e+00],
                 [-9.12352955e-01, -1.58024167e+00,  0.00000000e+00],
                 [-9.76751990e-01,  1.69178407e+00,  0.00000000e+00],
                 [ 1.95350397e+00, -3.53000000e-09,  0.00000000e+00]])
sett = "refGeom_axes"
cds = np.concatenate([cds[np.newaxis],cds[np.newaxis]],axis=0)


# sett = "refGeom_nz"
# cds[:,:,-1]*=-1

# # sett = "test_refGeomSwap"
# # cds[1,:,-1]*=-1

dw = np.ones(2)
# sett = "input_eqH"
# cds = np.load("../coordinates/trimer/"+sett+".npy")
# cds = Wfn.molecule.rotateBackToFrame(cds,3,2,1)
# dw = np.load("../coordinates/trimer/"+sett+"_dw.npy")
HOASpectrum = CalculateSpectrum.HarmonicApproxSpectrum(Wfn, cds, dw, path, testName)
gm = HOASpectrum.calculateG(cds,dw)
np.savetxt("gm_spc_shit"+sett,gm)
print('dun')
# a = np.loadtxt("gmat_onlyZ/gm_input_eqH")
# b = np.loadtxt("gmat_onlyZ/gm_justOinput_eqH")
# c = np.loadtxt("gmat_onlyZ/gm_refGeom")
# d = np.loadtxt("gmat_onlyZ/gm_justOrefGeom")
# e = np.loadtxt("gm_refGeom_nz")
# f = np.loadtxt("gm_refGeom_nz")
# lst = [a,b,c,d,e,f]
# for mat in lst:
#     plt.matshow(mat)
#     plt.colorbar()
# plt.show()