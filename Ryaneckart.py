import numpy as np

def eckartRotate(xx):
    com = getCOM()
    xx-=COM[:,np.newaxis,:]
    masses = [massO,massH,...]
    for walker in range(len(xx)):
        #f=multiply mass * eckart atom * walker's atom
	f=3x3 matrix 

xx = np.load("thoseFuckingCoords.npy")
com,allVecs,killList = eckartRotate(xx)
