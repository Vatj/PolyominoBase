import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

import icy;icy.Use_Seaborn()
from seaborn import despine

from colorsys import hsv_to_rgb
from random import uniform,choice
from interface_analysis import qBFS, loadManyResults, concatenateResults
from scipy.stats import linregress,binom

"""RANDOM THEORY SECTION """
def plotRandomTheory(I_size,g_len):
     
     b_asym=binom(I_size,.5)
     b_sym=binom(I_size/2,.5)

     s_hats=np.linspace(0,1,33)

     plt.figure()
     plt.plot(s_hats,float(g_len-1)/(g_len+1)*b_asym.pmf(s_hats*I_size),c='firebrick',ls='-')
     plt.plot(s_hats[::2],2/float(g_len+1)*b_sym.pmf(s_hats[::2]*I_size/2),c='royalblue',ls='-')

     T_star=0.65
     valid_syms=s_hats[::2][np.where(s_hats[::2]>=T_star)]
     valid_asyms=s_hats[np.where(s_hats>=T_star)]

     prob0_syms=b_sym.cdf(s_hats[::2][np.where(s_hats[::2]<T_star)][-1]*I_size/2)
     prob0_asyms=b_asym.cdf(s_hats[np.where(s_hats<T_star)][-1]*I_size)

     plt.plot(valid_asyms,(b_asym.cdf(valid_asyms*I_size)-prob0_asyms)/(1-prob0_asyms),c='firebrick',ls='-')
     plt.plot(valid_syms,(b_sym.cdf(valid_syms*I_size/2)-prob0_syms)/(1-prob0_syms),c='royalblue',ls='-')

     plt.text(.5,2/float(g_len+1)*b_sym.pmf(I_size/4.)*.9,'sym',ha='center',va='top')

     plt.text(.5,float(g_len-1)/(g_len+1)*b_asym.pmf(I_size/2.)*.9,'asym',ha='center',va='bottom')
     
     plt.yscale('log')
     plt.show(block=False)
     
     

"""PHYLOGENCY SECTION """
def plotBulkPhylogeny(data_struct,use_offset=False):
     plt.figure()
     for data in data_struct:
          plotPhylogeny(*data.values(),called=True,use_offset=use_offset)
     plt.ylabel(r'$\hat{S}$')
     plt.xlabel('generations')
     #plt.text(100,.8,'+={},-={}'.format(ups,downs))
     plt.show(block=False)



def plotPhylogeny(mae,mai,ms,called=False,use_offset=False):
     if not called:
          plt.figure()
     ups=0
     downs=0
     
     for interface_type,cr in zip([mae,mai,ms],[(.25,.38),(0.58,0.75),(0.91,.08)]):
          for data in interface_type:
               slope, intercept, r_value, p_value, std_err = linregress(xrange(len(data)-1),data[1:])
               if slope>0:
                    ups+=1
               else:
                    downs+=1
          
               h= uniform(*cr) if cr[1]>cr[0] else choice([uniform(cr[0],1),uniform(0,cr[1])])
               s = uniform(0.2, 1)
               v = uniform(0.5, 1)       
               r, g, b = hsv_to_rgb(h, s, v)
               g_offset=data[0] if use_offset else 0
               plt.plot(xrange(g_offset,g_offset+len(data)-1),data[1:],c=(r,g,b),alpha=0.2,zorder=1)
               #plt.plot(xrange(g_offset,g_offset+len(data)-1),[slope*x+intercept for x in xrange(len(data)-1)],'--',c=(r,g,b),zorder=10)        

     patches = [plt.plot([],[],c=color,ls='--')[0] for color in [(.21,.8,.16),(0.16,.22,.8),(.8,.16,.16)]]
     plt.legend(patches,['AE','AI','S'], frameon=False)
     
     if not called:
          plt.ylabel(r'$\hat{S}$')
          plt.xlabel('generations')
          plt.show(block=False)
     else:
          return ups,downs


from scipy import stats
def bootstrap(data, n_boot=10000, ci=68):
     boot_dist = []
     for i in range(int(n_boot)):
          resampler = np.random.randint(0, data.shape[0], data.shape[0])
          sample = data.take(resampler, axis=0)
          boot_dist.append(np.nanmean(sample, axis=0))
     b= np.array(boot_dist)
     s1 = np.apply_along_axis(stats.scoreatpercentile, 0, b, 50.-ci/2.)
     s2 = np.apply_along_axis(stats.scoreatpercentile, 0, b, 50.+ci/2.)
     return (s1,s2)
    
def tsplotboot(data,**kw):
    fig,ax = plt.subplots()
    x = np.arange(data.shape[1])
    est = np.nanmean(data, axis=0)
    cis = bootstrap(data)
    print
    ax.fill_between(x,cis[0],cis[1],alpha=0.2, **kw)
    
    ax.plot(x,est,**kw)
    for i in data:
         plt.plot(x,i,alpha=0.2)
    ax.margins(x=0)
    plt.show(block=False)

def form(data):
     length = len(sorted(data,key=len, reverse=True)[0])
     return np.array([xi+[np.nan]*(length-len(xi)) for xi in data])

def stripGen(data):
     return [i[1:] for i in data]
