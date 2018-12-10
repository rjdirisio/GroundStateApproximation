import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

numAtoms=13

def writeNewWalkers(walkz,title):
    coordAr = ['O','O','O','O','H','H','H','H','H','H','H','H','H']
    v=0
    weight=0
    m=0
    writeFile = open(title,"w")
    for walkn,walker in enumerate(walkz):
        m = 0
        writeFile.write("13\n")
        writeFile.write("0 0 0 0 0\n")
        while m<numAtoms:
            writeFile.write("%s %0.12f  %0.12f  %0.12f\n" % (coordAr[m],walkz[walkn,m,0],walkz[walkn,m,1],walkz[walkn,m,2]))
            m+=1
        writeFile.write("\n")
        v+=3
        weight+=1
    writeFile.close()


#Distance formula r = sqrt((x2-x1) + (y2-y1) + (z2-z1))
angstr = 0.529177

fl=open('../coordinates/tetramer/eqGeometry/eqRot.xyz','r')
linz = fl.readlines()
cd = np.zeros((13,3))
cct=0
for ln,line in enumerate(linz):
    sp = line.split()
    if len(sp)==5:
        wt=sp[0]
    elif len(sp)==1:
        print 'bobobo'
    elif len(sp)==4:
        cd[cct,:]=sp[1:]
        cct+=1
print cd
cd /=angstr #bohr shit
o2 = cd[2-1,:]
o3 = cd[3-1,:]
h7 = cd[7-1,:]
h8 = cd[8-1,:]
h9 = cd[9-1,:]
h10 = cd[10-1,:]
print 'H9 coordinate', h9
#________________

OH7=h7-o2 #vectr
OH8=h8-o2 #vectr
OH9=h9-o3 #vectr
OH10=h10-o3 #vectr
print 'OH9 vector', OH9
#_________________
rOH7 = la.norm(OH7)
rOH8 = la.norm(OH8)
rOH9 = la.norm(OH9)
rOH10 = la.norm(OH10)
print 'ROH9', rOH9
#_________________
uvOH7 = OH7/rOH7
uvOH8= OH8/rOH8
uvOH9= OH9/rOH9
uvOH10= OH10/rOH10
print 'unit vector along OH9', uvOH9
#eq = rOH*uvOH

dx7=21.149006817057
dx8=-21.150920333664
dx9=-20.076217284936
dx10= 20.078104983626
norm = la.norm([dx7,dx8,dx9,dx10])
dx7/=norm
dx8/=norm
dx9/=norm
dx10/=norm

dx7*=0.001 #bohr
dx8*=0.001
dx9*=0.001
dx10*=0.001

nscans = 1000
scannedCds = np.repeat(cd[np.newaxis,:,:],nscans,axis=0)
scannedCds2 = np.repeat(cd[np.newaxis,:,:],nscans,axis=0)

print scannedCds.shape
print scannedCds[0]
for i in range(1,nscans+1):
    newH7=uvOH7*(rOH7+i*dx7)
    newH8=uvOH8*(rOH8+i*dx8)
    newH9=uvOH9*(rOH9+i*dx9)
    newH10=uvOH10*(rOH10+i*dx10)
    #print "newH9 q/coordinate/q",newH9
    newH7+=o2
    newH8+=o2
    newH9+=o3
    newH10+=o3
    scannedCds[i-1,7-1]= newH7
    scannedCds[i-1,8-1]= newH8
    scannedCds[i-1,9-1]= newH9
    scannedCds[i-1,10-1]=newH10
    #print "final coordinate",newH9

for j in range(1,nscans+1):
    newH7 = uvOH7 * (rOH7 + j * -dx7)
    newH8 = uvOH8 * (rOH8 + j * -dx8)
    newH9 = uvOH9 * (rOH9 + j * -dx9)
    newH10 = uvOH10 * (rOH10 + j * -dx10)
    #print "newH9 q/coordinate/q", newH9
    newH7 += o2
    newH8 += o2
    newH9 += o3
    newH10 += o3
    scannedCds2[j - 1, 7 - 1] = newH7
    scannedCds2[j - 1, 8 - 1] = newH8
    scannedCds2[j - 1, 9 - 1] = newH9
    scannedCds2[j - 1, 10 - 1] = newH10
    #print "final coordinate", newH9

print scannedCds.shape
print scannedCds2.shape
print cd.shape
totalScanned = np.concatenate((np.flip(scannedCds2,axis=0),cd[np.newaxis,:,:],scannedCds),axis=0)
numWalkers = len(totalScanned)
writeNewWalkers(totalScanned,'rOHscan')

#rOH8
# -21.150920333664
#rOH7
# 21.149006817057
#rOH10
# 20.078104983626
#rOH9
# -20.076217284936

