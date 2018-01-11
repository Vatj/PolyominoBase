import matplotlib.pyplot as plt
import numpy as np
import icy;icy.Use_Seaborn()

from scipy.stats import binom
from operator import itemgetter
from collections import defaultdict
from decorator import decorator
from matplotlib.backends.backend_pdf import PdfPages

BASE_FILE_PATH='/scratch/asl47/Data_Runs/Interface/T{2:.6f}/{0}_{1}_T{2:.6f}_Mu{3:.6f}_Gamma{4:.6f}_Run{5}.txt'


def LoadSizes(r_type,temperature,mu,gamma,runs=1):
     sizes=defaultdict(int)
     for r in xrange(250+runs):
          for line in open(BASE_FILE_PATH.format('Sizes',r_type,temperature,mu,gamma,r)):
               (size,count)=[int(i) for i in line.rstrip().split()]
               sizes[size]+=count
     return sizes

def LoadPairSizes(temperature,mu,gamma,runs=1):
     return (LoadSizes('S',temperature,mu,gamma,runs),'Selection',runs),(LoadSizes('R',temperature,mu,gamma,runs),'Random',runs)

def LoadData(d_type,r_type,temperature,mu,gamma,runs=1):
     data=[]
     for r in xrange(250,250+runs):
          data.append(np.loadtxt(BASE_FILE_PATH.format(d_type,r_type,temperature,mu,gamma,r)))
     return np.stack(data,axis=2)

def LoadPairFitness(temperature,mu,gamma,runs=1):
     return LoadData('Fitness','S',temperature,mu,gamma,runs),LoadData('Fitness','R',temperature,mu,gamma,runs)

def LoadPairStrengths(temperature,mu,gamma,runs=1):
     return LoadData('Strengths','S',temperature,mu,gamma,runs),LoadData('Strengths','R',temperature,mu,gamma,runs)

def PlotInterfaceStrengths(sss,L,title_string=''):
     fig, axarr = plt.subplots(1, 2, sharey=True,sharex=True)
     T=sss[0].shape[1]-1
     pl=plt.cm.plasma
     b_uncorr=binom(T,0.5)
     b_corr=binom(T/2,0.5)
     n=0
     slices=np.zeros((1),dtype=np.int32)
     slices=np.concatenate([slices,np.logspace(0,np.log10(sss[0].shape[0]-1),10,dtype=np.int32)])
     pop_size=sss[0][0,0,0]#/(L*(L+1)/2)

     
     for ss,ax in zip(sss,axarr.reshape(-1)):
     #np.append(data,means[:,-1].reshape(T+1,1))
          means=np.mean(ss,axis=2)
          data=means[slices]
          for i,time_slice in enumerate(data):
               ax.plot(np.linspace(0,1,T+1),time_slice*1./pop_size,c=pl(float(i)/data.shape[1]),marker='o')
     
          ax.plot(np.linspace(0,1,T+1),[(L-1.)/(L+1)*b_uncorr.pmf(i)+2./(L+1)*b_corr.pmf(i/2)*(1-i%2) for i in xrange(T+1)],'kx',markeredgewidth=1)
     
     axarr[0].set_ylim([10**-5,1.5])
     axarr[0].set_yscale('log',nonposy='mask')
     axarr[0].set_ylabel(r'$\langle f \rangle$')
     axarr[0].set_title(r'Selection')
     axarr[1].set_title(r'Random')
     #axarr[1].yaxis.tick_right()
     #axarr[1].minorticks_off()
     #plt.setp(axarr[1].get_yticklabels(), visible=False)
     #axarr[1]._off()
     axarr[1].tick_params(axis='y', which='both', left='off') #
     plt.figtext(0.5,0.01,'Interaction Strength',ha='center',va='bottom')


     
     #fig.subplots_adjust(right=0.8)
     #cbar_ax = fig.add_axes([0.895, 0.15, 0.05, 0.7])
     plt.figtext(0.95,0.5,r'increasing generation $\longrightarrow$',rotation=90,ha='left',va='center')
     sm = plt.cm.ScalarMappable(cmap='plasma')
     sm._A = []
     cbar=fig.colorbar(sm,use_gridspec=True)
     #cbar.outline.set_visible(False)
     cbar.set_ticks([])
     fig.set_tight_layout(True)
     fig.suptitle('{}'.format(title_string))
     
     plt.show(block=False)

         
def PlotPhenotypeSizes(ss,title_string=''):
     fig=plt.figure()
     plot_params={'Selection':{'c':'coral','marker':'o'},'Random':{'c':'royalblue','marker':'D'}}
     poly_sizes=[1, 1, 2, 7, 18, 60, 196, 704, 2500, 9189, 33896, 126759, 476270, 1802312, 6849777, 26152418, 100203194, 385221143, 1485200848, 5741256764, 22245940545, 86383382827, 336093325058]

     min_v=-1
     max_v=-1
     
     for s in ss:
          plt.plot(s[0].keys(),[s[0][k]/float(s[2]) for k in s[0].keys()],label=s[1],**plot_params[s[1]])
          if min_v<0:
               min_v=min(s[0].values())/float(s[2])
               max_v=max(s[0].values())/float(s[2])
          else:
               min_v = min(s[0].values())/float(s[2]) if min(s[0].values())/float(s[2])<min_v else min_v
               max_v = max(s[0].values())/float(s[2]) if max(s[0].values())/float(s[2])>max_v else max_v
          print min_v,max_v
               
     plt.plot(range(1,len(poly_sizes)+1),poly_sizes,marker='x',c='k',markeredgewidth=1.25,ls='',label='One-sided Polyominoes')


     plt.ylim((min_v*0.8,max_v*1.5))
     plt.yscale('log',nonposy='mask')
     plt.ylabel(r'$\langle f \rangle$')
     plt.xlabel(r'Phenotype Size')
     plt.legend()
     fig.suptitle('{}'.format(title_string))
     plt.show(block=False)

def PlotFitness(fss,title_string=''):

     fig, axarr = plt.subplots(1, 2, sharey=True,sharex=True)
     log_slices=np.logspace(np.log10(1),np.log10(fss[0].shape[0]-1),8,dtype=np.int32)
     #[0,1,2,3,4,5,7,9,12,15,20,50,100,500,1000,2500,4999]
     main_slices=np.logspace(np.log10(1),np.log10(fss[0].shape[0]-1),30,dtype=np.int32)
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

     #plt.yscale('log',nonposy='clip')
     fig.suptitle('{}'.format(title_string))
     fig.set_tight_layout(True)
     plt.show(block=False)

 

def PlotBindingStrengths(Ts,I_size):
     plt.figure()
     cs=['firebrick','royalblue','darkgreen','coral','orchid','goldenrod']
     for T,c in zip(Ts,cs):
          plt.plot(np.linspace(0,1,101),np.exp(-1*np.linspace(1,0,101)/T),lw=2,ls='--',label='T={}'.format(T),c=c)
          plt.plot(np.linspace(0,1,I_size+1),np.exp(-1*np.linspace(1,0,I_size+1)/T),lw=2,ls='',marker='o',c=c)
          
     plt.legend()
     plt.yscale('log',nonposy='mask')
     plt.show(block=False)
     

     #//DEAD CODE//#

def PlotPhenotypeSizesOld(ss):
     fig=plt.figure()
     plot_params={'Selection':{'c':'coral','marker':'o'},'Random':{'c':'royalblue','marker':'D'}}
     poly_sizes=[1, 1, 2, 7, 18, 60, 196, 704, 2500, 9189, 33896, 126759, 476270, 1802312, 6849777, 26152418, 100203194, 385221143, 1485200848, 5741256764, 22245940545, 86383382827, 336093325058]
     min_v=-1
     max_v=-1
     for s in ss:
          hist,bins=np.histogram(s[0],bins=np.linspace(0,37,38))
          plt.plot(bins[:-1],hist*1./s[2],label=s[1],**plot_params[s[1]])
          if min_v<0:
               min_v=np.min(hist[hist>0])*1./s[2]
               max_v=np.max(hist[hist>0])*1./s[2]
          else:
               min_v = np.min(hist[hist>0])*1./s[2] if np.min(hist[hist>0])*1./s[2]<min_v else min_v
               max_v = np.max(hist[hist>0])*1./s[2] if np.max(hist[hist>0])*1./s[2]>max_v else max_v
               print min_v,max_v    
     plt.plot(range(1,len(poly_sizes)+1),poly_sizes,marker='x',c='k',markeredgewidth=1.25,ls='',label='One-sided Polyominoes')
     plt.ylim((min_v*0.8,max_v*1.5))
     plt.yscale('log',nonposy='mask')
     plt.ylabel(r'$\langle f \rangle$')
     plt.xlabel(r'Phenotype Size')
     plt.legend()
     fig.set_tight_layout(True)
     plt.show(block=False)    

def LoadSizesOld(r_type,temperature,mu,runs=1):
     sizes=[]
     for r in xrange(runs):
          sizes.extend([sum([int(s) for s in line.rstrip().split()][2:])  for line in open("/scratch/asl47/Data_Runs/Interface/T{1:.6f}/Sizes_{0}_T{1:.6f}_Mu{2:.6f}_Run{3}.txt".format(r_type,temperature,mu,r))])
     return sizes

def LoadPairSizesOld(temperature,mu,runs=1):
     return LoadSizesOld('S',temperature,mu,runs),LoadSizes('R',temperature,mu,runs)

def LoadStrengthsOld(r_type,temperature,mu,gamma,runs=1,interface_size=16):
     total_ret=[]
     for r in xrange(runs):
          per_runs=[]
          bulk_data=[{k:v for (k,v) in zip([float(d) for d in line.rstrip().split()[::2]],[int(d) for d in line.rstrip().split()[1::2]])} for line in open(BASE_FILE_PATH.format('Strengths',r_type,temperature,mu,gamma,r))]
          for data in bulk_data:
               temp=np.zeros((interface_size+1))
               for k,v in data.iteritems():
                    temp[int(k*interface_size)]=v
               per_runs.append(temp)
          total_ret.append(np.stack(per_runs,axis=1))
     return np.stack(total_ret,axis=2)
