import numpy as np

def RandomHamming():
    N=1000
    rands=np.random.randint(low=0,high=2**16, size=N)

    hams=[]

    for i in xrange(N):
        for j in xrange(i,N):
            hams.append(bin(rands[i] ^rands[j]).count('1'))
    return hams
