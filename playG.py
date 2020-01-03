import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
cfg = 'ffinal_allH_rn_spc_xfx_NoTweakZIneulers'
M = np.load("mu2ave_"+cfg+".npy")
G = np.loadtxt("../Gmats/trimer/"+cfg+".gmat")
sqrt2o = 1/np.sqrt(2)
sqrt4o = 1/np.sqrt(4)

#########Just shared Protons
# k = np.zeros((2,2))
# np.fill_diagonal(k,1/np.sqrt(2))
# # print(k)
# k[0,1] = sqrt2o
# k[1,0] = -sqrt2o
# print(k)
# Mred = np.copy(M)[np.ix_([3,6],[3,6])]
# print(Mred)
# res = np.matmul(k,np.matmul(Mred,kT))
#########/Just shared protons

#########Including r,theta,phi of shared protons
# k = np.zeros((9,9))
# np.fill_diagonal(k,1/np.sqrt(2))
# k[0,0]=1.0 #rH8
# k[1,1]=1.0 #thH8
# k[2,2]=1.0 #phiH8
# pairs = [(3,6),(4,7),(5,8)]
# pairsrevd = pairs[::-1]
# for pair in pairs:
#     k[pair[0],pair[1]]=sqrt2o
#     k[pair[1],pair[0]]=-1.0*sqrt2o
# Mred = np.copy(M)[:9,:9]
# Gred = np.copy(G)[:9,:9]
# res = np.matmul(k,np.matmul(Mred,k.T))
# resG = np.matmul(k,np.matmul(Gred,k.T))
# print(res)
# plt.matshow(res)
# plt.colorbar()
# plt.show()
#########Including r,theta,phi of shared protons

###Everything except eulers##########
"""
0, rOH8
1, thH8
2, phiH8
3, rOH9
4, thH9
5, phiH9
6, rOH10
7, thH10
8, phiH10

9, R14
10, R15
11, th145
12, R26
13, R27
14, th267

15, Roo1
16, Roo2
17, thOOO
"""
# k = np.zeros((18,18))
# np.fill_diagonal(k,sqrt2o)
# for i in range(10,14):
#     k[i,i] = 1.0
#
# k[0,0]=1.0 #rH8
# k[1,1]=1.0 #thH8
# k[2,2]=1.0 #phiH8
# k[11,11]=sqrt2o
# k[14,14]=sqrt2o
# pairs = [(3,6),(4,7),(5,8),(11,14),(15,16)]
# pairsrevd = pairs[::-1]
# for pair in pairs:
#     k[pair[0],pair[1]]=sqrt2o
#     k[pair[1],pair[0]]=-1.0*sqrt2o
#
# pairsF = [9,10,12,13]
# rowz = np.array([
#           [1,1,1,1],
#           [1,1,-1,-1],
#           [1,-1,1,-1],
#           [1,-1,-1,1]])*sqrt4o
# z=0
# for i in pairsF:
#     k[i,pairsF]=rowz[z]
#     z+=1
#
# Mtest = np.copy(M)
# Mtest2 = np.delete(Mtest, [9,10,11,12,13,14], axis=1)
# Mtest2 = np.delete(Mtest2, [9,10,11,12,13,14], axis=0)
# Mred = np.copy(Mtest2)
# res = np.matmul(k,np.matmul(Mred,k.T))
# print(res)
# plt.matshow(res)
# plt.colorbar()
# plt.show()

#add on eulers
k = np.zeros((24,24))
np.fill_diagonal(k,sqrt2o)
for i in range(10,14):
    k[i,i] = 1.0

k[0,0]=1.0 #rH8
k[1,1]=1.0 #thH8
k[2,2]=1.0 #phiH8
k[11,11]=sqrt2o
k[14,14]=sqrt2o
"""
0, rOH8
1, thH8
2, phiH8
3, rOH9
4, thH9
5, phiH9
6, rOH10
7, thH10
8, phiH10

9, R14
10, R15
11, th145
12, R26
13, R27
14, th267

15, Roo1
16, Roo2
17, thOOO

18 theta1
19 Phi1
20 Chi1
21 theta2
22 Phi2
23 Chi2
"""
pairs = [(3,6),(4,7),(5,8),(11,14),(15,16),(18,21),(19,22),(20,23)]
pairsrevd = pairs[::-1]
for pair in pairs:
    k[pair[0],pair[1]]=sqrt2o
    k[pair[1],pair[0]]=-1.0*sqrt2o

pairsF = [9,10,12,13]
rowz = np.array([
          [1,1,1,1],
          [1,1,-1,-1],
          [1,-1,1,-1],
          [1,-1,-1,1]])*sqrt4o
z=0
for i in pairsF:
    k[i,pairsF]=rowz[z]
    z+=1

Mtest = np.copy(M)
Mtest2 = np.copy(M)
#A[idx,:][:,idx]
idx = np.concatenate((np.arange(9),np.arange(15,24),np.arange(9,15)))
Mred = Mtest2[idx,:][:,idx]
Gred = G[idx,:][:,idx]
# Mtest[[9,10,11,12,13,14]] = Mtest2[[18,19,20,21,22,23]]
# Mtest[[18,19,20,21,22,23]] = M[[9,10,11,12,13,14]]
# Mtest[:,[9,10,11,12,13,14]] = Mtest2[:,[18,19,20,21,22,23]]
# Mtest[:,[18,19,20,21,22,23]] = M[:,[9,10,11,12,13,14]]

# Mtest2 = np.delete(Mtest, [9,10,11,12,13,14], axis=1)
# Mtest2 = np.delete(Mtest2, [9,10,11,12,13,14], axis=0)
# Mred = np.copy(Mtest2)
# res = Mred
# resG = Gred
res = np.matmul(k,np.matmul(Mred,k.T))
resG = np.matmul(k,np.matmul(Gred,k.T))
print(res)
plt.matshow(res,cmap='hot')
plt.title("Transformed M")
plt.colorbar()
plt.matshow(np.abs(res),cmap='hot')
plt.title("Absolute Value of Transformed M")
plt.colorbar()
plt.savefig("AbsM"+cfg,dpi=400)

plt.matshow(resG,cmap='hot')
plt.title("Transformed G")
plt.colorbar()
plt.matshow(np.abs(resG),cmap='hot')
plt.title("Absolute Value of Transformed G")
plt.colorbar()
# plt.show()
plt.savefig("AbsG"+cfg,dpi=400)

resP = res - np.diag(np.diag(res))
resGP = resG - np.diag(np.diag(resG))

plt.matshow(resP,cmap='hot')
plt.title("Transformed M")
plt.colorbar()
plt.matshow(np.abs(resP),cmap='hot')
plt.title("Absolute Value of Transformed M")
plt.colorbar()
plt.savefig("RAbsM"+cfg,dpi=400)

plt.matshow(resGP,cmap='hot')
plt.title("Transformed G")
plt.colorbar()
plt.matshow(np.abs(resGP),cmap='hot')
plt.title("Absolute Value of Transformed G")
plt.colorbar()
# plt.show()
plt.savefig("RAbsG"+cfg,dpi=400)
