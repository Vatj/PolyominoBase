import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binom
import icy;icy.Use_Seaborn()


def LoadStrengths(T,r_type='',runs=1):
     total_ret=[]
     for r in xrange(runs):
          per_runs=[]
          bulk_data=[{k:v for (k,v) in zip([float(d) for d in line.rstrip().split()[::2]],[int(d) for d in line.rstrip().split()[1::2]])} for line in open("/rscratch/asl47/Bulk_Run/Interfaces/Strengths{}_Run{}.txt".format('_'+r_type if r_type!='' else '',r))]
          for data in bulk_data:
               temp=np.zeros((T+1))
               for k,v in data.iteritems():
                    temp[int(k*T)]=v
               per_runs.append(temp)
          total_ret.append(np.stack(per_runs,axis=1))
     return np.stack(total_ret,axis=2)

def LoadSizes(r_type='',runs=1):
     sizes=[]
     for r in xrange(runs):
          sizes.extend([sum([int(s) for s in line.rstrip().split()][2:])  for line in open("/rscratch/asl47/Bulk_Run/Interfaces/Sizes{}_Run{}.txt".format('_'+r_type if r_type!='' else '',r))])
     return sizes

def LoadFitness(r_type='',runs=1):
     fitnesses=[]
     for r in xrange(runs):
          fitnesses.append(np.loadtxt("/rscratch/asl47/Bulk_Run/Interfaces/Fitness{}_Run{}.txt".format('_'+r_type if r_type!='' else '',r)))
     return np.stack(fitnesses,axis=2)

def PlotInterfaceStrengths(sss,L,T=16):
     fig, axarr = plt.subplots(1, 2, sharey=True,sharex=True)
     pl=plt.cm.plasma
     b_uncorr=binom(T,0.5)
     b_corr=binom(T/2,0.5)
     n=0
     slices=np.zeros((1),dtype=np.int32)
     slices=np.concatenate([slices,np.logspace(0,np.log10(4999),10,dtype=np.int32)])
     pop_size=sss[0][0,0,0]#/(L*(L+1)/2)
     print pop_size
     
     for ss,ax in zip(sss,axarr.reshape(-1)):
     #np.append(data,means[:,-1].reshape(T+1,1))
          means=np.mean(ss,axis=2)
          data=means[:,slices]
          for i,time_slice in enumerate(data.T):
               ax.plot(np.linspace(0,1,T+1),time_slice/pop_size,c=pl(i*1./data.shape[1]),marker='o')
     
          ax.plot(np.linspace(0,1,T+1),[(L-1.)/(L+1)*b_uncorr.pmf(i)+2./(L+1)*b_corr.pmf(i/2)*(1-i%2) for i in xrange(T+1)],'kx',markeredgewidth=1)
     
     axarr[0].set_ylim([10**-5,1.5])
     axarr[0].set_yscale('log',nonposy='mask')
     axarr[0].set_ylabel(r'$\langle f \rangle$')
     axarr[0].set_title(r'Selection')
     axarr[1].set_title(r'Random')
     axarr[1].yaxis.tick_right()

     plt.figtext(0.5,0.01,'Interaction Strength',ha='center',va='bottom')


     fig.subplots_adjust(right=0.8)
     cbar_ax = fig.add_axes([0.895, 0.15, 0.05, 0.7])
     plt.figtext(0.95,0.5,r'increasing generation $\longrightarrow$',rotation=90,ha='left',va='center')
     sm = plt.cm.ScalarMappable(cmap='plasma')
     sm._A = []
     cbar=fig.colorbar(sm,cax=cbar_ax)
     #cbar.outline.set_visible(False)
     cbar.set_ticks([])

     
     plt.show(block=False)
         
def PlotPhenotypeSizes(ss,runs=1):
     fig=plt.figure()
     plot_params={'Selection':{'c':'coral','marker':'o'},'Random':{'c':'royalblue','marker':'D'}}
     poly_sizes=[1, 1, 2, 7, 18, 60, 196, 704, 2500, 9189, 33896, 126759, 476270, 1802312, 6849777, 26152418, 100203194, 385221143, 1485200848, 5741256764, 22245940545, 86383382827, 336093325058]
     for s in ss:
          hist,bins=np.histogram(s[0],bins=np.linspace(0,37,38))
          plt.plot(bins[:-1],hist*1./runs,label=s[1],**plot_params[s[1]])
     plt.plot(range(1,len(poly_sizes)+1),poly_sizes,marker='x',c='k',markeredgewidth=1.25,ls='',label='One-sided Polyominoes')

     plt.ylim((0.8,np.max(hist)*1.5/runs))
     plt.yscale('log',nonposy='mask')
     plt.ylabel(r'$\langle f \rangle$')
     plt.xlabel(r'Phenotype Size')
     plt.legend()
     plt.show(block=False)
c
def PlotFitness(fss):
     fig, axarr = plt.subplots(1, 2, sharey=True,sharex=True)
     log_slices=np.logspace(np.log10(1),np.log10(4999),10,dtype=np.int32)
     #[0,1,2,3,4,5,7,9,12,15,20,50,100,500,1000,2500,4999]
     main_slices=np.logspace(np.log10(1),np.log10(4999),50,dtype=np.int32)
     #range(50)+range(50,5000,50)
     for fs,ax in zip(fss,axarr.reshape(-1)):
          for f in fs.T:
               ax.errorbar(log_slices,f[0][log_slices],yerr=np.sqrt(f[1][log_slices]),alpha=0.25,ls='--')
          ax.errorbar(main_slices,np.mean(fs,axis=2)[main_slices,0],yerr=np.sqrt(np.mean(fs,axis=2)[main_slices,1]),lw=2,c='k')
     #axarr[0].set_xlabel('Generations')
     plt.figtext(0.5,0.01,'Generations',ha='center',va='bottom')
     axarr[0].set_ylabel(r'$\langle \mathcal{F} \rangle$')
     axarr[0].set_ylim((0,1))
     axarr[0].set_xscale('log')
     axarr[0].set_title(r'Selection')
     axarr[1].set_title(r'Random')
     axarr[1].yaxis.tick_right()
     
     plt.show(block=False)


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
          
          #for j in xrange(i,N):
          #     hams.append(bin(rands[i] ^revbits(rands[j],L)).count('1'))

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
