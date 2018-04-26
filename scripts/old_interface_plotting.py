import matplotlib.pyplot as plt
import numpy as np
import icy;icy.Use_Seaborn()
import seaborn as sns

from scipy.stats import binom,sem,chi2
from operator import itemgetter
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages
from tile_shape_visuals import Visualise_Shape_From_Binary as VSFB

import glob

#BASE_FILE_PATH='/scratch/asl47/Data_Runs/Interface_Cron/{0}_{1}_T{2:.6f}_Mu{3:.6f}_Gamma{4:.6f}_Run{5}.txt'
#BASE_FILE_PATH='/rscratch/asl47/Bulk_Run/Interfaces/{0}_{1}_T{2:.6f}_Mu{3:.6f}_Gamma{4:.6f}_Run{5}.txt'
BASE_FILE_PATH='../{0}_{1}_T{2:.6f}_Mu{3:.6f}_Gamma{4:.6f}_Run{5}.txt'

def GetMaxRun(r_type,temperature,mu,gamma):
     return max([int(s[s.rindex('Run')+3:-4]) for s in glob.glob(BASE_FILE_PATH.format('Sizes',r_type,temperature,mu,gamma,'*'))])
     
def VisualisePhenotypes(r_type,temperature,mu,gamma,run):
     line_count=0
     page_max=24
     shapes_data=sorted([[int(i) for i in line.split()] for line in open(BASE_FILE_PATH.format('Phenotypes',r_type,temperature,mu,gamma,run))],key=lambda z: sum(z[4:]))
     return shapes_data
     for shape in shapes_data:
          if line_count%page_max==0:
               fig, axarr = plt.subplots(6,4,figsize=(8.27,11.69))
     
          VSFB(shape[2:],'',axarr.reshape(-1)[line_count],1,True,'')
          line_count+=1
          line_count%=page_max
     else:
          for l in xrange(line_count,page_max):
               axarr.reshape(-1)[l].set_axis_off()
     
     plt.show(block=False)

def Test():
     un=set()
     for r in xrange(20):
          for v in VisualisePhenotypes('S',0.000001,.5,1,r):
               un.add(tuple(v[2:]))
          print len(un)
     return sorted(list(un),key = lambda z : sum([1 for i in z[2:] if i!=0]))
             
def LoadSizes(r_type,temperature,mu,gamma,runs=0):
     sizes=defaultdict(int)
     phen_count=[]
     for r in xrange(runs):
          phens=0
          for line in open(BASE_FILE_PATH.format('Sizes',r_type,temperature,mu,gamma,r)):
               (size,count)=[int(i) for i in line.rstrip().split()]
               if size and count:
                    sizes[size]+=count
                    phens+=count
          phen_count.append(phens)
     return sizes
     #return phen_count



def LoadPairSizes(temperature,mu,gamma,runs=0):
     if runs==0:
          runs=min(GetMaxRun('S',temperature,mu,gamma),GetMaxRun('R',temperature,mu,0))
          print "Runs evaluated at ",runs
     return (LoadSizes('S',temperature,mu,gamma,runs),'Selection',runs),(LoadSizes('S',temperature,mu,gamma,runs),'Random',runs)

def LoadData(d_type,r_type,temperature,mu,gamma,runs=1):
     data=[]
     for r in xrange(runs):
          data.append(np.loadtxt(BASE_FILE_PATH.format(d_type,r_type,temperature,mu,gamma,r)))
     return np.stack(data,axis=2)

def LoadPairFitness(temperature,mu,gamma,runs=0):
     if runs==0:
          runs=min(GetMaxRun('S',temperature,mu,gamma),GetMaxRun('R',temperature,mu,0))
          print "Runs evaluated at ",runs
     return LoadData('Fitness','S',temperature,mu,gamma,runs),LoadData('Fitness','S',temperature,mu,gamma,runs)

def LoadPairStrengths(temperature,mu,gamma,runs=0):
     if runs==0:
          runs=min(GetMaxRun('S',temperature,mu,gamma),GetMaxRun('R',temperature,mu,0))
          print "Runs evaluated at ",runs
     return LoadData('Strengths','S',temperature,mu,gamma,runs),LoadData('Strengths','S',temperature,mu,gamma,runs)

def PlotStrengthRatios(data_frame,N_tiles,title_string='',temperature=-1):
     fig, axarr = plt.subplots(2, 6, sharey='row',sharex=True)
     interface_size=int((data_frame[0].shape[1]-2)/1.5)
     generations=data_frame[0].shape[0]
     NEGATIVE_INF_LIMIT=-0.5
     

     
     binom_pair=binom(interface_size,0.5)
     binom_self=binom(interface_size/2.,0.5)
     n=0
     slices=np.zeros((1),dtype=np.int32)
     slices=np.concatenate([slices,np.logspace(0,np.log10(generations-1),10,dtype=np.int32)])
     population_size=np.sum(data_frame[0][0,:,0])/((N_tiles*4)*(N_tiles*4+1)/2.)




     pair_strengths=np.linspace(0,1,interface_size+1)
     self_strengths=np.linspace(0,1,interface_size/2+1)
     bar_width=(1./interface_size/2.)

     binding_alpha=0.05
     Tsx=list(np.logspace(np.log10(0.03),-1,6))

     #print pair_strengths,self_strengths,bar_width
     for data,axx in zip(data_frame,axarr.T):
          if temperature>0:
               axx[0].axvline(1+Tsx.pop(0)*np.log(binding_alpha),-1,1)

          run_averaged=np.mean(data,axis=2)
          means=np.mean(run_averaged,axis=0)
          errs=sem(run_averaged,axis=0)
          
          

          #print binom_pair.pmf(pair_strengths*interface_size),binom_self.pmf(self_strengths*interface_size/2)
          #print means[0][:interface_size+1]/population_size,means[0][interface_size+1:]/population_size
          
          pairwise_ratio=(means[:interface_size+1]/population_size/((N_tiles*4-1.)*(N_tiles*4)/2.))/(binom_pair.pmf(pair_strengths*interface_size))

          pairwise_err=(errs[:interface_size+1]/population_size/((N_tiles*4-1.)*(N_tiles*4)/2.))/(binom_pair.pmf(pair_strengths*interface_size))
          pairwise_err_upper=np.log10(pairwise_ratio+pairwise_err)
          pairwise_err_lower=np.log10(pairwise_ratio-pairwise_err)
          
          self_ratio=(means[interface_size+1:]/population_size/(N_tiles*4))/(binom_self.pmf(self_strengths*interface_size/2)) 
          self_err=(errs[interface_size+1:]/population_size/(N_tiles*4))/(binom_self.pmf(self_strengths*interface_size/2))
          self_err_upper=np.log10(self_ratio+self_err)
          self_err_lower=np.log10(self_ratio-self_err)
   
          pairwise_G_stat=2* np.nansum(means[:interface_size+1]*np.log(pairwise_ratio))
          print "G: {}, p-value: {}".format(pairwise_G_stat,chi2.sf(pairwise_G_stat,interface_size))

          self_G_stat=2* np.nansum(means[interface_size+1:]*np.log(self_ratio))
          print "G: {}, p-value: {}".format(self_G_stat,chi2.sf(self_G_stat,interface_size/2))

 
          
          #pair_s=means[-1][:interface_size+1]/population_size/((N_tiles*4-1.)*(N_tiles*4)/2.)
          #pair_b=binom_pair.pmf(pair_strengths*interface_size)
 

                  

          pairwise_ratio[pairwise_ratio < 1] = -1 * np.reciprocal( pairwise_ratio[pairwise_ratio < 1])
          pairwise_ratio[pairwise_ratio>0]=np.log10(pairwise_ratio[pairwise_ratio>0])
          pairwise_ratio[pairwise_ratio<0]=-1*np.log10(-1*pairwise_ratio[pairwise_ratio<0])
          pairwise_colors=['palegreen' if cond else 'darkgreen' for cond in np.isinf(pairwise_ratio)]
          pairwise_ratio[np.isinf(pairwise_ratio)]=NEGATIVE_INF_LIMIT

          self_ratio[self_ratio < 1] = -1 * np.reciprocal( self_ratio[self_ratio < 1])
          self_ratio[self_ratio>0]=np.log10(self_ratio[self_ratio>0])
          self_ratio[self_ratio<0]=-1*np.log10(-1*self_ratio[self_ratio<0])
          self_colors=['lightblue' if cond else 'darkblue' for cond in np.isinf(self_ratio)]
          self_ratio[np.isinf(self_ratio)]=NEGATIVE_INF_LIMIT
          
          #print pairwise_ratio
          #print self_ratio
          

          for strength,(rat,err) in enumerate(zip(pairwise_err_lower,pairwise_err_upper)):
               axx[0].plot([strength/32.]*2,[rat,err],'r--')
          for strength,(rat,err) in enumerate(zip(self_err_lower,self_err_upper)):
               axx[1].plot([strength/16.]*2,[rat,err],'r--')

               
          
          axx[0].bar(pair_strengths,pairwise_ratio,width=bar_width,color=pairwise_colors)
          axx[1].bar(self_strengths,self_ratio,width=bar_width,color=self_colors)
          axx[0].axhline(0,0.02,0.98,c='k',lw=0.75,ls='--')
          axx[1].axhline(0,0.02,0.98,c='k',lw=0.75,ls='--')


     #axarr[0,0].set_title('Selection')
     #axarr[0,1].set_title('Random')
     axarr[0,1].yaxis.tick_right()
     axarr[1,1].yaxis.tick_right()
     fig.text(0.5, 0.04, 'Interface Strengths', ha='center', va='center')
     fig.text(0.03, 0.5, 'log E/O ', ha='center', va='center', rotation='vertical')
     plt.tight_layout(h_pad=0.02, w_pad=0.02)
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
          
               
     plt.plot(range(1,len(poly_sizes)+1),poly_sizes,marker='x',c='k',markeredgewidth=1.25,ls='',label='One-sided Polyominoes')


     plt.ylim((min_v*0.8,max_v*1.5))
     #plt.yscale('log',nonposy='mask')
     plt.ylabel(r'$\langle f \rangle$')
     plt.xlabel(r'Phenotype Size')
     plt.legend()
     fig.suptitle('{}'.format(title_string))
     plt.show(block=False)

def PlotFitness(data_frame,title_string='',g_factor=1):

     fig, axarr = plt.subplots(1, 6, sharey=True,sharex=True)

     #range(50)+range(50,5000,50)
     generations=data_frame[0].shape[0]

     for data,ax in zip(data_frame,axarr.reshape(-1)):
          for fitness in data.T:
               #(gen_vals,fit_vals)=ValueJumps(fitness[0],0.025)
               #gen_vals=np.append(gen_vals,generations)
               #gen_vals=np.insert(gen_vals,0,1)
               #fit_vals=np.append(fit_vals,fit_vals[-1])
               if abs(fitness[0,0]-fitness[0,-1])/float(fitness[0,0])<0.05:
                    ax.plot([1,generations*g_factor],[np.mean(fitness[0])]*2,alpha=0.25,ls='--')
               else:
                    ax.plot(np.linspace(0,generations-1,generations/2)*g_factor+1,np.mean(fitness[0].reshape(-1, 2), axis=1),alpha=0.4,ls='-')
          ax.plot(np.linspace(0,generations-1,generations)*g_factor+1,np.mean(data,axis=2)[:,0],lw=2,c='k')
               #ax.errorbar(range(fitness.shape[1]),fitness[0],yerr=np.sqrt(fitness[1]),alpha=0.25,ls='--')
          #ax.errorbar(range(fitness.shape[1]),np.mean(data,axis=2)[:,0],yerr=np.sqrt(np.mean(data,axis=2)[:,1]),lw=2,c='k')
     #axarr[0].set_xlabel('Generations')
     plt.figtext(0.5,0.01,'Generations',ha='center',va='bottom')
     axarr[0].set_ylabel(r'$\langle \mathcal{F} \rangle$')
     axarr[0].set_ylim((0,5))
     axarr[0].set_xscale('log')
     axarr[0].set_title(r'Selection')
     axarr[1].set_title(r'Random')
     axarr[1].yaxis.tick_right()

     #plt.yscale('log',nonposy='clip')
     fig.suptitle('{}'.format(title_string))
     fig.set_tight_layout(True)
     plt.show(block=False)

          
def PlotInterfaceStrengths(data_frame,N_tiles,title_string=''):
     fig, axarr = plt.subplots(1, 2, sharey=True,sharex=True)
     interface_size=int((data_frame[0].shape[1]-2)/1.5)
     generations=data_frame[0].shape[0]
     pl=plt.cm.plasma
     binom_pair=binom(interface_size,0.5)
     binom_self=binom(interface_size/2.,0.5)
     n=0
     slices=np.zeros((1),dtype=np.int32)
     slices=np.concatenate([slices,np.array([generations-1])])#np.logspace(0,np.log10(generations-1),20,dtype=np.int32)])
     population_size=78000. #float(data_frame[0][0,0,0])#/(L*(L+1)/2)
     print population_size
     
     for data,ax in zip(data_frame,axarr.reshape(-1)):
     #np.append(data,means[:,-1].reshape(T+1,1))
          means=np.mean(data,axis=2)
          data=means[slices]

          for i,time_slice in enumerate(data):
               #ax.plot(np.linspace(0,1,interface_size+1),time_slice[:interface_size+1]/population_size,c=pl(float(i)/data.shape[0]),marker='o')
               
               ax.plot(np.linspace(0,1,interface_size/2+1),time_slice[interface_size+1:]/population_size,c=pl(float(i)/data.shape[0]),marker='s')
     
          #ax.plot(np.linspace(0,1,interface_size+1),[(N_tiles-1.)/(N_tiles+1)*binom_pair.pmf(i)+2./(N_tiles+1)*binom_self.pmf(i/2)*(1-i%2) for i in xrange(interface_size+1)],'kx',markeredgewidth=1)

          #ax.plot(np.linspace(0,1,interface_size+1),[(N_tiles-1.)/(N_tiles+1)*binom_pair.pmf(i) for i in xrange(interface_size+1)],'kx',markeredgewidth=1)
          ax.plot(np.linspace(0,1,interface_size/2+1),[2./(N_tiles+1)*binom_self.pmf(i) for i in xrange(interface_size/2+1)],'k+',markeredgewidth=1)
          
     
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

         
def PlotAll(temp,mu,gamma,runs,N_t):
     PlotFitness(LoadPairFitness(temp,mu,gamma,runs))
     PlotPhenotypeSizes(LoadPairSizes(temp,mu,gamma,runs))
     PlotStrengthRatios(LoadPairStrengths(temp,mu,gamma,runs),N_t)

 

def PlotBindingStrengths(Ts,I_size):
     plt.figure()
     #cs=['firebrick','royalblue','darkgreen','coral','orchid','goldenrod','c','k','y','r','b','g']
     for T in Ts:
          #plt.plot(np.linspace(0,1,I_size+1),np.exp(-1*np.linspace(1,0,I_size+1)/T),lw=2,ls='--',marker='o',label='O={}'.format(T))
          
          plt.plot(np.linspace(0,1,I_size+1),FF(np.linspace(0,1,I_size+1),T),lw=2,marker='s',ls='')#,label='FF={}'.format(T))
          #plt.plot(np.linspace(0,1,I_size+1),np.exp(-1*np.linspace(1,0,I_size+1)/T),lw=2,ls='',marker='o',c=c,markersize=6,markeredgewidth=1)

          #plt.plot(np.linspace(0,1,101),np.linspace(0,1,101)**(5./T),lw=2,ls='--',label='T={}'.format(T),c=c)

     binom_pair=binom(I_size,0.5)
     binom_self=binom(I_size/2.,0.5)

     comb=np.zeros((I_size+1))
     for i in xrange(I_size+1):
          comb[i]+=binom_pair.pmf(i)
          if i%2==0:
               comb[i]+=binom_self.pmf(i/2)
     plt.plot(np.linspace(0,1,I_size+1),comb,c='k',marker='^',ls='--')
     
     plt.legend()
     plt.yscale('log',nonposy='mask')
     plt.ylim((10**(-10),1))
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


def ValueJumps(data, stepsize=1):
     return (np.where(np.diff(data) >= stepsize)[0],np.array([np.mean(cluster) for cluster in np.split(data, np.where(np.diff(data) >= stepsize)[0]+1)]))

import matplotlib.colors as colors
def PlotParamSpace():
     plt.figure()
     param_yield=np.empty((15,6))
     

     j=0
     
     for M in np.logspace(0,np.log10(16),15):
          i=0
          for T in np.logspace(np.log10(0.03),-1,6):
               #print M,T
               param_yield[j,i]=sum(LoadSizes('S',T,M,1,250).values())/250.
               #print param_yield[j,i]
               i+=1
          j+=1

     #print param_yield.min(),param_yield.max()
     plt.pcolormesh(np.logspace(np.log10(0.03),-1,7),np.logspace(0,np.log10(16),16) ,param_yield, cmap='viridis',norm=colors.LogNorm(vmin=param_yield.min(), vmax=param_yield.max()))
     plt.title('Averaged Phenotype Discovery')
     plt.ylabel(r'$\langle \mu \cdot L \cdot I \rangle$')
     plt.xlabel(r'$T$')
     plt.yscale('log')
     plt.xscale('log')
                    
     # set the limits of the plot to the limits of the data
     #plt.axis([x.min(), x.max(), y.min(), y.max()])
     plt.colorbar()
     plt.tight_layout()
     plt.show(block=False)


def HistoryLoad(temperature=0.000001,mu=1,gamma=1,run=0):

     phen_line=True
     phenotype_IDs=[]
     selections=[]
     for line in open(BASE_FILE_PATH.format('PhenotypeHistory','S',temperature,mu,gamma,run)):
          converted=[int(i) for i in line.split()]
          if phen_line:
               phens=[]
               for i in xrange(0,len(converted),2):
                    phens.append((converted[i],converted[i+1]))
               phenotype_IDs.append(phens)
               phen_line=False
          else:
               selections.append(converted)
               phen_line=True
               continue

     
     return phenotype_IDs,selections

def HistoryDiagram(IDs,selections,low=-1,high=-1):
     if low==high:
          low=0
          high=len(IDs)
     plt.figure()
     seen_phens=set()
     seen_phens.add((0,0))
     phen_col={(0,0):[0,0,0]}
     unique_phens= set([item for sublist in IDs[low:high] for item in sublist])
     
     for phen in unique_phens:
          if sum(phen):
               phen_col[phen]=icy.generate_new_color(phen_col.values(),0)
          
     for g,(ID,sel) in enumerate(zip(IDs[low:high],selections[low:high])):
          
          for gen,idp in enumerate(ID):
               #print "hi"
               
               plt.scatter(gen,g,marker='o',c=phen_col[idp],lw=0.5,edgecolor='k',zorder=10)
          if g==0:
               for i,s in enumerate(sel):
                    plt.plot([i,i],[-1,0],'b-',lw=3,zorder=1)
  
          for i,s in enumerate(sel):
                    #print i,s
               plt.plot([s,i],[g,g+1],'k--',lw=0.5,zorder=10)

     for idpx,c in phen_col.iteritems():
          plt.scatter(-10,0,c=c,label=idpx)
     


     ##attempt
     if False:
          rev_sel=selections[::-1]
          fixation=0#rev_sel[0].index(max(set(rev_sel[0]), key=rev_sel[0].count))
          
          printed_paths=set()
          for init in xrange(len(rev_sel[0])):
               fixation=init
               tm=len(selections)-1
               for nex in rev_sel[1:]:
                    next_fix=nex[fixation]
                    path=(fixation,next_fix,tm,tm-1)
                    if path not in printed_paths:
                         printed_paths.add(path)
                         plt.plot([fixation,next_fix],[tm,tm-1],'k-',lw=1.5)
                         fixation=next_fix
                         
                         tm-=1          
          
     plt.legend()
     plt.axis('off')
     plt.tight_layout()
     plt.xlim([-1,len(selections[0])+1])
     plt.show(block=False)

from interface_analysis import LoadT,RandomHistorySampling
from itertools import groupby
def PlotAnalysis(gz,N_samps=50,normed=False):
     
     #return sampled_history_strengths
     #plt.figure()
     for i in xrange(N_samps):
          sampled_history_strengths=RandomHistorySampling(*gz,goback=1999)
          for key,interface_pairing in sampled_history_strengths.iteritems():
               for sequence in interface_pairing:
                    norm_factor=sequence[0] if normed else 1
                    if sequence:
                         plt.plot(range(len(sequence)),np.array(sequence)/norm_factor-1,alpha=0.4)
                         return
                         uniq=[L[0] for L in groupby(sequence)]
                         ups=0
                         downs=0
                         #for (a,b) in zip(uniq,uniq[1:]):
                         #     if b>a:
                         #          ups+=1
                         #     else:
                         #          downs+=1
                         #print ups,downs

     plt.plot([0,2000],[0,0],'k--')
     plt.xlabel('t')
     plt.ylabel(r'$\Delta \hat{S}$')
     sns.despine()
     plt.show(block=False)

from scipy import stats
def FF(xs,t):
     q=t
     t=.05
     mi=stats.norm.cdf(xs[0],q,t)
     ma=stats.norm.cdf(xs[-1],q,t)
     print ma,mi
     scal=1./(ma-mi)
     #print scal,mi
     return stats.norm.cdf((xs-q)/t)#-mi)*scal#*xs**2

def pl_f(L):
     binom_pair=binom(L,0.5)
     binom_self=binom(L/2.,0.5)

     plt.figure()
     plt.plot(np.linspace(0,1,L+1),66./78*binom_pair.pmf(np.linspace(0,L,L+1)),'o')
     plt.plot(np.linspace(0,1,L/2+1),12./78*binom_self.pmf(np.linspace(0,L/2,L/2+1)),'^')
     plt.yscale('log')
     plt.xlabel(r'$\hat{S}$')
     plt.ylabel(r'$Pr$')
     plt.title('Strength Distribution')
     sns.despine()
     plt.show(block=False)
