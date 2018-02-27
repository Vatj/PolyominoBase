import numpy as np


def revbits(x,L):
    return int(bin(x)[2:].zfill(L)[::-1], 2)

def RandomHamming(L=16):
     N=50000
     rands=np.random.randint(low=0,high=2**L, size=N*8)
     hams=[]
     bs=[]
     for i in xrange(0,N*8,8):
          for j in xrange(8):
               bs.append(bin(rands[i+j]).count('1'))
               for k in xrange(j,8):
                    hams.append(bin(rands[i+j] ^revbits(rands[i+k],L)).count('1'))
     plt.figure()
     hist,bins=np.histogram(hams,bins=np.linspace(0,L+1,L+2))
     plt.plot(bins[:-1]/L,hist)
     b=binom(L,0.5)
     b2=binom(L/2,0.5)
     plt.plot(np.linspace(0,1,L+1),[50000*(4*7)*b.pmf(i) for i in xrange(L+1)],'kx',markeredgewidth=1)
     plt.plot(np.linspace(0,1,9),[50000*(4*7)*b.pmf(i)+50000*8*b2.pmf(j) for i,j in zip(xrange(0,L+1,2),xrange(0,L/2+1))],'m+',markeredgewidth=1)
     hist,bins=np.histogram(bs,bins=np.linspace(0,L+1,L+2))
     plt.plot(bins[:-1]/L,hist)
     plt.plot(np.linspace(0,1,L+1),[N*b.pmf(i) for i in xrange(L+1)],'ro',markeredgewidth=1)
     plt.yscale('log')
     plt.show(block=False)

from numpy import uint8
def ConjugateInteraction(genotype_str):
    genotype=[int(i) for i in genotype_str.split()]
    for base in genotype:
        print base, int(bin(~uint8(base))[2:].zfill(8)[::-1],2), "({},{})".format(bin(base)[2:].zfill(8),bin(~uint8(base))[2:].zfill(8)), str(int(bin(~uint8(base))[2:].zfill(8)[::-1],2)) in genotype_str
    #return int(bin(~uint8(base))[2:].zfill(8)[::-1],2)
        
