import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

import icy;icy.Use_Seaborn()
from seaborn import despine

from colorsys import hsv_to_rgb
from random import uniform,choice
from interface_analysis import qBFS, loadManyResults, concatenateResults,RandomWalk,BindingStrength, set_length
from scipy.stats import linregress,binom,scoreatpercentile
from itertools import combinations_with_replacement as cwr
from random import choice

"""RANDOM THEORY SECTION """
def plotRandomTheory(I_size,g_len):
     
     b_asym=binom(I_size,.5)
     b_sym=binom(I_size/2,.5)

     s_hats=np.linspace(0,1,I_size+1)

     plt.figure()
     plt.plot(s_hats,float(g_len-1)/(g_len+1)*b_asym.pmf(s_hats*I_size),c='firebrick',ls='-')
     plt.plot(s_hats[::2],2/float(g_len+1)*b_sym.pmf(s_hats[::2]*I_size/2),c='royalblue',ls='-')

     plt.text(.5,2/float(g_len+1)*b_sym.pmf(I_size/4.)*.9,'sym',ha='center',va='top')

     plt.text(.5,float(g_len-1)/(g_len+1)*b_asym.pmf(I_size/2.)*.9,'asym',ha='center',va='bottom')
     
     plt.yscale('log')
     plt.show(block=False)
     

def plotInterfaceProbability(l_I,l_g,Nsamps=False):

     def SF_sym(T_stars):
          return binom(l_I/2,.5).sf(np.ceil(l_I/2*T_stars)-1)#*(1./(l_g+1))
     def SF_asym(T_stars):
          return binom(l_I,.5).sf(np.ceil(l_I*T_stars)-1)#-sym(T_stars))/2*((l_g-1.)/(l_g+1))

     def sym_factor(A):
          return float(2)/(A+1)
     def asym_factor(A):
          return float(A-1)/(A+1)

     s_hats=np.linspace(0,1,l_I+1)


     fig, ax1 = plt.subplots()
     ax1.plot(s_hats[::2],np.log10(sym_factor(l_g)*SF_sym(s_hats[::2])),ls='',marker='^',c='royalblue')
     ax1.plot(s_hats,np.log10(asym_factor(l_g)*SF_asym(s_hats)),ls='',marker='o',c='firebrick')

    
     ax2 = ax1.twinx()
     
     ratios=np.log10((sym_factor(l_g)*SF_sym(s_hats))/(asym_factor(l_g)*SF_asym(s_hats)))
     ax2.plot(s_hats,ratios,c='darkseagreen')
     crossover=np.where(ratios>0)[0][0]
     ax2.axvline(s_hats[crossover],color='k',ls='--')
     ax2.axhline(color='k',ls='-',lw=0.2)
     Is={8:np.uint8,16:np.uint16,32:np.uint32,64:np.uint64}
     if Nsamps:
          set_length(l_I)
          s_m=np.zeros(l_I+1)
          a_m=np.zeros(l_I+1)
          for _ in xrange(Nsamps):
               indices=choice(list(cwr(range(l_g),2)))
               if indices[0]!=indices[1]:
                    bases=np.random.randint(0,np.iinfo(Is[l_I]).max,dtype=Is[l_I],size=2)
                    
                    a_m[np.where(BindingStrength(*bases)>=s_hats)]+=1
               else:
                    base=np.random.randint(0,np.iinfo(Is[l_I]).max,dtype=Is[l_I])
                    s_m[np.where(BindingStrength(base,base)>=s_hats)]+=1
          s_m2=np.ma.log10(s_m/Nsamps)
          a_m2=np.ma.log10(a_m/Nsamps)
          ax1.plot(s_hats[::2],s_m2[::2],ls='--',c='royalblue')
          ax1.plot(s_hats,a_m2,ls='--',c='firebrick')
     

     crossover_height=np.log10(asym_factor(l_g)*SF_asym(1))/2.
     ax1.text(crossover/float(l_I),crossover_height,'crossover',ha='right',va='center',rotation=90)
     scale_factor=np.log10(asym_factor(l_g)*SF_asym(s_hats))[0]-np.log10(asym_factor(l_g)*SF_asym(s_hats))[-1]
     ax1.text(.2,np.log10(sym_factor(l_g)*SF_sym(.2))-scale_factor*0.05,'sym',va='top')
     ax1.text(.2,np.log10(asym_factor(l_g)*SF_asym(.2)+scale_factor*0.05),'asym',va='bottom')
     
     ax2.text(.1,(ratios[-1]-ratios[0])*.05+ratios[0],'ratio',ha='center',va='bottom')

     ax1.set_ylabel('Interface Probability')
     ax2.set_ylabel('Ratio')
     ax1.set_xlabel(r'$S^*$')

     ax1.spines['top'].set_visible(False)
     ax2.spines['top'].set_visible(False)
     
     plt.show(block=False)     
     

     

def plotInterfaceActivation(I_size,L_size):
     b_asym=binom(I_size,.5)
     b_sym=binom(I_size/2,.5)
     N_interactions=(2**(2*I_size-1))+(2**(I_size-1))

     def sym(T_stars):
          return (2**I_size)*b_sym.sf(np.ceil(I_size/2*T_stars-1))*(1./(L_size+1))
     def asym(T_stars):
          return ((2**(2*I_size))*b_asym.sf(np.ceil(I_size*T_stars-1))-sym(T_stars))/2*((L_size-1.)/(L_size+1))

     s_hats=np.linspace(0,1,I_size*5+1)
     plt.figure()
     plt.plot(s_hats,sym(s_hats)/N_interactions)
     plt.plot(s_hats,asym(s_hats)/N_interactions)
     plt.plot(s_hats,sym(s_hats)/asym(s_hats),'k--')

     plt.text(.6,sym(0.6)/N_interactions*0.9,'sym',ha='right',va='top')
     plt.text(.6,asym(0.6)/N_interactions*.9,'asym',ha='right',va='top')
     plt.text(.6,sym(0.6)/asym(0.6)/.8,'ratio',ha='right',va='bottom')
     

     plt.xlabel(r'$\hat{S}$')
     plt.ylabel('Frac pairwise')
     plt.title(r'$I_{size}=$ %i' % I_size)
     plt.yscale('log')
     plt.tight_layout()
     plt.show(block=False)   
          

"""PHYLOGENCY SECTION """
def plotBulkPhylogeny(data_struct,use_offset=False):
     plt.figure()
     for data in data_struct:
          plotPhylogeny(*data.values(),called=True,use_offset=use_offset)
     plt.ylabel(r'$\hat{S}$')
     plt.xlabel('generations')
     plt.show(block=False)

def plotPhylogeny(mae,mai,ms,called=False,use_offset=False):
     if not called:
          plt.figure()     
     for interface_type,cr in zip([mae,mai,ms],[(.25,.38),(0.58,0.75),(0.91,.08)]):
          for data in interface_type:
               slope, intercept, r_value, p_value, std_err = linregress(xrange(len(data)-1),data[1:])
               h= uniform(*cr) if cr[1]>cr[0] else choice([uniform(cr[0],1),uniform(0,cr[1])])
               s = uniform(0.2, 1)
               v = uniform(0.5, 1)       
               r, g, b = hsv_to_rgb(h, s, v)
               g_offset=data[0] if use_offset else 0
               plt.plot(xrange(g_offset,g_offset+len(data)-1),data[1:],c=(r,g,b),alpha=0.2,zorder=1)
     patches = [plt.plot([],[],c=color,ls='--')[0] for color in [(.21,.8,.16),(0.16,.22,.8),(.8,.16,.16)]]
     plt.legend(patches,['AE','AI','S'], frameon=False)
     if not called:
          plt.ylabel(r'$\hat{S}$')
          plt.xlabel('generations')
          plt.show(block=False)

def bootstrap(data, n_boot=10000, ci=68):
     boot_dist = []
     for i in range(int(n_boot)):
          resampler = np.random.randint(0, data.shape[0], data.shape[0])
          sample = data.take(resampler, axis=0)
          boot_dist.append(np.nanmean(sample, axis=0))
     b= np.array(boot_dist)
     s1 = np.apply_along_axis(scoreatpercentile, 0, b, 50.-ci/2.)
     s2 = np.apply_along_axis(scoreatpercentile, 0, b, 50.+ci/2.)
     return (s1,s2)
    
def tsplotboot(ax,data,title='',**kw):
    x = np.arange(data.shape[1])
    est = np.nanmean(data, axis=0)
    cis = bootstrap(data,10)
    ax.fill_between(x,cis[0],cis[1],alpha=0.3,color='dimgray', **kw)
    
    ax.plot(x,est,c='dimgray',lw=2)
    for i in data:
         ax.plot(x,i,alpha=0.05)
    ax.margins(x=0)
    
    ax.set_ylabel(r'$\langle \hat{S} \rangle$')
    if title!='':
         ax.set_title(title)
    plt.show(block=False)


def plotData(cc,I_size,t_star):
     fig,axs = plt.subplots(3)
     for (k,v),ax in zip(cc.iteritems(),axs):
          print k,len(v)
          if len(v)==0:
               continue
          #tsplotboot(ax,v,k+' {}'.format(len(v)))
          for q in v:
               ax.plot(q[0]+np.linspace(0,len(q)-2,len(q)-1),q[1:])
          if 'A' in k:
               #print k, "plotting here"
               pgs=RandomWalk(64,100,.5,t_star,1,1)
               ax.plot(range(0,1200,12),pgs[:-1],'r--')
               
               ax.plot([1100,1500],[pgs[-1]]*2,'k--')
          else:
               #print k, "plotting there"
               pgs=RandomWalk(32,50,.5,t_star,1,1)
               
               ax.plot(range(0,1200,24),pgs[:-1],'r--')
               ax.plot([1100,1500],[pgs[-1]]*2,'k--')
     plt.xlabel('elapsed generations')
     plt.tight_layout(pad=0)
     fig.suptitle(r'$l_I = %i , S^* = %.2f$' % (I_size,t_star))
     plt.show(block=False)

def plotWs(Ws,I_size,t_star):
     plt.figure()
     for W in Ws:
          lists = sorted(W.items()) # sorted by key, return a list of tuples
          x, y = zip(*lists)
          plt.plot(x,[float(i)/sum(y) for i in y],ls='',marker='o',label=len(W))
     plt.plot(np.linspace(0,1,I_size/2+1),binom(I_size/2,.5).pmf(np.linspace(0,I_size/2,I_size/2+1)),'k--')
     plt.plot(np.linspace(0,1,I_size+1),binom(I_size,.5).pmf(np.linspace(0,I_size,I_size+1)),'k--')
     plt.axvline(t_star,0,1,c='r')

     plt.legend()
     plt.yscale('log')
     plt.ylabel('Prob')
     plt.xlabel(r'$ \hat{S} $')
     plt.title(r'$l_I = %i , S^* = %.2f$' % (I_size,t_star))
     plt.show(block=False)
