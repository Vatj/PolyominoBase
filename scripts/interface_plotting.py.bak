import matplotlib.pyplot as plt
from matplotlib.patches import Patch,Rectangle
import numpy as np

import icy;icy.Use_Seaborn()
from tile_shape_visuals import VisualiseSingleShape as VSS
#from seaborn import despine

from colorsys import hsv_to_rgb
from random import uniform,choice
from interface_analysis import qBFS, loadManyResults, concatenateResults,RandomWalk,BindingStrength, set_length
from scipy.stats import linregress,binom,scoreatpercentile
from itertools import combinations_with_replacement as cwr,product
from random import choice
from math import ceil
from copy import deepcopy
from collections import defaultdict

def plotFF(g,s):
    xs=np.linspace(0,1,101)
    plt.figure()
    plt.plot(xs,(xs>s)*((xs-s)**g/(1-s)**g))
    plt.show(block=False)

"""RANDOM THEORY SECTION """
def plotFitnessFunc():
     data=np.loadtxt('../distrs.txt')
     bin_range=(0,10)
     bin_count=500
     hists=np.apply_along_axis(lambda a: np.histogram(a, range=bin_range,bins=bin_count,density=True)[0], 1, data)
     plt.figure()
     for row in hists:
          plt.plot(np.linspace(*bin_range,num=bin_count),row)

     plt.xlabel('Fitness')
     plt.ylabel('PDF')
     plt.show(block=False)


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

     ax2.text(.1,(ratios[-1]-ratios[0])*.015+ratios[0],'ratio',ha='center',va='bottom')

     ax1.set_ylabel(r'$  \log \left( Pr_{\mathrm{interface}} \right) $')
     ax2.set_ylabel(r'$\log \frac{Pr_{\mathrm{sym}}}{Pr_{\mathrm{asym}}}$')
     ax1.set_xlabel(r'$\hat{S}^*$')

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
          plotPhylogeny(*list(data.values()),called=True,use_offset=use_offset)
     plt.ylabel(r'$\hat{S}$')
     plt.xlabel('generations')
     plt.show(block=False)

def plotPhylogeny(mae,mai,ms,called=False,use_offset=False):
     if not called:
          plt.figure()
     for interface_type,cr in zip([mae,mai,ms],[(.25,.38),(0.58,0.75),(0.91,.08)]):
          for data in interface_type:
               slope, intercept, r_value, p_value, std_err = linregress(range(len(data)-1),data[1:])
               h= uniform(*cr) if cr[1]>cr[0] else choice([uniform(cr[0],1),uniform(0,cr[1])])
               s = uniform(0.2, 1)
               v = uniform(0.5, 1)
               r, g, b = hsv_to_rgb(h, s, v)
               g_offset=data[0] if use_offset else 0
               plt.plot(range(g_offset,g_offset+len(data)-1),data[1:],c=(r,g,b),alpha=0.2,zorder=1)
     patches = [plt.plot([],[],c=color,ls='--')[0] for color in [(.21,.8,.16),(0.16,.22,.8),(.8,.16,.16)]]
     plt.legend(patches,['AE','AI','S'], frameon=False)
     if not called:
          plt.ylabel(r'$\hat{S}$')
          plt.xlabel('generations')
          plt.show(block=False)

def bootstrap(data, n_boot=10000, ci=68):
     boot_dist = []
     for i in xrange(int(n_boot)):
          resampler = np.random.randint(0, data.shape[0], data.shape[0])
          sample = data.take(resampler, axis=0)
          sample=sample.astype(float)

          #try:
          boot_dist.append(np.nanmean(sample, axis=0))
          #except:
          #     return sample
          #     return
     b= np.array(boot_dist)
     s1 = np.apply_along_axis(scoreatpercentile, 0, b, 50.-ci/2.)
     s2 = np.apply_along_axis(scoreatpercentile, 0, b, 50.+ci/2.)
     return (s1,s2)

def tsplotboot(ax,data,title='',**kw):
    x = np.arange(data.shape[1])
    est = np.nanmean(data, axis=0)
    cis = bootstrap(data,100)
    ax.fill_between(x,cis[0],cis[1],alpha=0.3,color='dimgray', **kw)

    ax.plot(x,est,c='dimgray',lw=2)
    for i in data:
         ax.plot(x,i,alpha=0.05,rasterized=1)
    ax.margins(x=0)

    ax.set_ylabel(r'$\langle \hat{S} \rangle$')
    if title!='':
         ax.set_title(title)
    plt.show(block=False)


def plotData(cc,I_size,t_star,mu=1,g_size=12):
     fig,axs = plt.subplots(3)
     for (k,v),ax in zip(cc.iteritems(),axs):
          #print k,len(v)
          if v:#.size:
               v=np.array(v)
              # v=v[v[:,0]==t_star]
               #v=v[np.count_nonzero(np.isnan(v),axis=1)<1000]
               #tsplotboot(ax,v,k+' {}'.format(len(v)))


               gen_length=2000#v[0].shape[0]
               co_factor=2 if 'S' in k else 1
               step_length=int(g_size/(2.*mu/co_factor))
               N_steps=int(ceil(gen_length/float(step_length)))
               pgs=RandomWalk(I_size/co_factor,N_steps,.5,t_star,1,1)
               ax.plot(range(0,(N_steps+1)*step_length,step_length),pgs[:-1],'r--')
               ax.plot([0,gen_length],[pgs[-1]]*2,'k--',alpha=0.5)

               step_length=1
               N_steps=int(ceil(gen_length/float(step_length)))
               pgs=RandomWalk(I_size/co_factor,N_steps,float(co_factor)/6,t_star,1,1)
               ax.plot(range(0,(N_steps+1)*step_length,step_length),pgs[:-1],'b--')
               ax.plot([0,gen_length],[pgs[-1]]*2,'k--',alpha=0.5)

     plt.xlabel('elapsed generations')
     plt.tight_layout(pad=0)
     fig.suptitle(r'$l_I = %i , S^* = %.2f$' % (I_size,t_star))
     plt.show(block=False)

def plotWs(Ws,I_size,t_star):
     plt.figure()
     for K,W in Ws.iteritems():
          lists = sorted(W.items()) # sorted by key, return a list of tuples
          x, y = zip(*lists)
          plt.plot(x,[float(i)/sum(y) for i in y],ls='',marker='o',label=K)
     plt.plot(np.linspace(0,1,I_size/2+1),binom(I_size/2,.5).pmf(np.linspace(0,I_size/2,I_size/2+1)),'k--')
     plt.plot(np.linspace(0,1,I_size+1),binom(I_size,.5).pmf(np.linspace(0,I_size,I_size+1)),'k--')
     plt.axvline(t_star,0,1,c='r')

     plt.legend()
     plt.yscale('log')
     plt.ylabel('Prob')
     plt.xlabel(r'$ \hat{S} $')
     plt.title(r'$l_I = %i , S^* = %.2f$' % (I_size,t_star))
     plt.show(block=False)

def plotInterfaceCounts(counts):
     plt.figure()
     means=[]
     for cnt in counts:
          a=0
          for i,j in cnt.iteritems():
               a+=i*j
          means.append(a/np.sum(cnt.values(),dtype=np.float64))
     plt.plot(range(len(counts)),means)
     plt.show(block=False)

def plotFatals(counts):
     plt.figure()
     plt.plot(range(counts.shape[1]),np.mean(counts,axis=0))
     plt.show(block=False)

def PhenotypicTransitions(phen_trans,N=40,crit_factor=0.5):
     print "N set for ",N
     common_transitions=deepcopy(phen_trans)
     for phen_key,trans in phen_trans.iteritems():
          #print "max",phen_key,max(trans.iterkeys(), key=(lambda key: trans[key])),max(trans.values())
          for tran,count in trans.iteritems():
               if count<N*crit_factor:
                    del common_transitions[phen_key][tran]

     for key in common_transitions.keys():
          if not common_transitions[key]:
               del common_transitions[key]
     return common_transitions


def plotTransitions(phen_trans,cdict=None):

     fig =plt.figure()
     ncols=int(ceil(np.sqrt(len(phen_trans))))
     nrows=int(ceil(len(phen_trans)/float(ncols)))
     cdict={None:(0,0,0)} if cdict is None else cdict

     for i,(key,trans) in enumerate(phen_trans.iteritems(),1):
          if len(trans)==0:
               continue
          ax = fig.add_subplot(nrows, ncols, i)
          ax.set_title('T:'+r', '.join(map(str,key)),fontsize=22)
          for phen in trans.keys():
               if phen not in cdict:
                    cdict[phen]=icy.generate_new_color(cdict.values(),0.5)

          cols=[cdict[p] if p!=key else 'k' for p in sorted(trans.keys())]
          wedges,_ = ax.pie([trans[k] for k in sorted(trans.keys())],labels=sorted(trans.keys()),colors=cols,startangle=90,counterclock=False)

          for w in wedges:
               w.set_linewidth(1)
               w.set_edgecolor('k')
          centre_circle = plt.Circle((0,0),0.75,color='black', fc='white',linewidth=1)
          fig = plt.gcf()
          fig.gca().add_artist(centre_circle)
     plt.show(block=False)

     return cdict



def plotTransitionsDetailed(pt):
     scaled_count=defaultdict(dict)
     scaled_count[1][(1,1,1)]=None
     connection_subsets=defaultdict(list)
     for phen in pt.keys():
          if phen==(1,1,1):
               continue
          for sub_phen in pt[phen].keys():
               connection_subsets[phen].append((phen,sub_phen))




     #return connection_subsets
     for k,v in pt.iteritems():
          scaled_count[np.count_nonzero(k[2:])][k]=v
     size_counts={i:len(scaled_count[v]) for i,v in enumerate(sorted(scaled_count.keys()))}
     #return scaled_count,size_counts


     fig,ax = plt.subplots(1)
     connection_dict={}
     phen_dict={}
     for i,c in enumerate(sorted(scaled_count.keys())):
          offset=0 if len(scaled_count[c])%2==1 else .5
          for j,phen in enumerate(sorted(scaled_count[c])):
               phen_dict[phen]=AddPhenotypePatch(ax,phen,(i,len(scaled_count[c])/2 - j -offset))
               if i==0:
                    continue
               total_weight=sum(scaled_count[c][phen].values())
               for connector,weight in scaled_count[c][phen].iteritems():

                    con_x=np.count_nonzero(connector[2:])-1
                    offset2=0 if len(scaled_count[con_x+1])%2==1 else .5
                    con_y=len(scaled_count[con_x+1])/2-sorted(scaled_count[con_x+1]).index(connector)-offset2

                    spline_points=np.array([[con_x+.25,con_y],[con_x+.3,con_y],[i-.3,len(scaled_count[c])/2-j-offset],[i-.25,len(scaled_count[c])/2-j-offset]])
                    dx_f=i-con_x
                    dy_f=spline_points[-1,1]-spline_points[0,1]
                    dx=1
                    dy= np.sign(dy_f) if abs(dy_f)>=1 else 0

                    while dx<dx_f:
                         if int(spline_points[dx,1]*2)%2==size_counts[con_x+dx]%2:
                              bump_factor= 0
                         else:
                              bump_factor=.5 if np.sign(spline_points[-1,1]-spline_points[dx,1])>0 else -.5
                         adjustment_factor=dy+bump_factor
                         if abs(adjustment_factor)>1:
                              adjustment_factor=np.sign(adjustment_factor)*(adjustment_factor%1)
                         spline_points=np.insert(spline_points,dx+1,[con_x+dx,spline_points[dx,1]+adjustment_factor],axis=0)
                         dx+=1
                         dy=dy-(np.sign(dy_f)) if abs(dy)>=1 else 0
                    connection_dict[(phen,connector)]=AddConnectionPatch(ax,spline_points,float(weight)/total_weight)





     ax.set_aspect(1)
     ax.relim()
     ax.autoscale_view()

     ax.grid(False)
     plt.axis('off')
     def onpick(event):
          artis=event.artist
          coords=artis.get_bbox()

          mean_click=np.mean(coords.get_points(),axis=0)
          patch_coord=[int(np.round(mean_click[0])),np.round(mean_click[1]*2)/2.]
          if abs(patch_coord[1])<.25:
               patch_coord[1]=0

          phen_slices=scaled_count[sorted(scaled_count.keys())[patch_coord[0]]]
          j=int(len(phen_slices)/2-(0 if len(phen_slices)%2==1 else .5)-patch_coord[1])

          phen_key=sorted(phen_slices.keys())[j]
          #print phen_key
          phen_set=set()
          phen_set.add(phen_key)
          for phen_pairing, artists in connection_dict.iteritems():
               if phen_key==phen_pairing[0]:
                    phen_set.add(phen_pairing[1])
               alpha_v=1 if phen_key in phen_pairing else .1
               for artist in artists:
                    artist.set_alpha(alpha_v)
               for phen in phen_pairing:
                    for artist in phen_dict[phen]:
                         artist.set_alpha(.1)

          for phen in phen_set:
               for artist in phen_dict[phen]:
                    artist.set_alpha(1)




          fig.canvas.draw()
          return True

     fig.canvas.mpl_connect('pick_event', onpick)
     plt.show(block=False)





def AddPhenotypePatch(ax,shape,xy):
     ar_offsets={0:(0,-.25,0,.25),1:(-.25,0,.25,0),2:(0,.25,0,-.25),3:(.25,0,-.25,0)}
     cols=['darkgreen','royalblue','firebrick','goldenrod','mediumorchid']
     dx=shape[0]
     dy=shape[1]
     artists=[]
     scale=.5/max(dx,dy)
     for i,j in product(xrange(dx),xrange(dy)):
          if(shape[2+i+j*dx]):
               new_x=xy[0]+(i-dx/2.)*scale
               new_y=xy[1]+(dy/2.-j)*scale-(1.*scale)
               artists.append(ax.add_patch(Rectangle((new_x,new_y), scale, scale, facecolor=cols[(shape[2+i+j*dx]-1)/4],edgecolor='slategrey',fill=True,hatch='////',lw=0,picker=10,alpha=0.1)))
               artists.append(ax.add_patch(Rectangle((new_x,new_y), scale, scale,edgecolor='maroon',fill=False,lw=2.5,alpha=0.1)))
               theta=(shape[2+i+j*dx]-1)%4;
               artists.append(ax.arrow(new_x+(.5+ar_offsets[theta][0])*scale,new_y+(.5+ar_offsets[theta][1])*scale, ar_offsets[theta][2]*scale, ar_offsets[theta][3]*scale, head_width=0.075*scale, head_length=0.15*scale, fc='k', ec='k',alpha=0.1))
     return artists


from scipy.interpolate import splprep, splev

def AddConnectionPatch(ax,pts,weight):

     tck, u = splprep(pts.T, u=None, s=0.0,k=3, per=False)
     u_new = np.linspace(u.min(), u.max(), 50)
     x_new, y_new = splev(u_new, tck, der=0)

     ar=ax.arrow(x_new[x_new.shape[0]/2],y_new[y_new.shape[0]/2],x_new[x_new.shape[0]/2+1]-x_new[x_new.shape[0]/2],y_new[y_new.shape[0]/2+1]-y_new[y_new.shape[0]/2], shape='full', lw=0, length_includes_head=True, head_width=.05,alpha=0.1)
     ln=ax.plot(x_new, y_new, 'b--',lw=weight*2,alpha=0.1)[0]
     return (ln,ar)


def getD():
     return {(4, 1, 1, 5, 7, 3): {(1, 1, 1): 93, (2, 1, 1, 3): 123, (2, 1, 1, 5): 96}, (3, 3, 1, 5, 2, 8, 0, 6, 4, 7, 3): {(2, 1, 1, 5): 94}, (4, 3, 0, 0, 1, 0, 5, 9, 11, 7, 0, 3, 0, 0): {(4, 1, 1, 5, 7, 3): 68}, (4, 4, 0, 0, 1, 0, 4, 5, 6, 0, 0, 8, 7, 2, 0, 3, 0, 0): {(1, 1, 1): 63, (2, 1, 1, 5): 85}, (2, 2, 1, 5, 7, 3): {(2, 1, 1, 5): 91}, (3, 2, 0, 1, 5, 7, 3, 0): {(2, 1, 1, 5): 98, (1, 1, 1): 96, (2, 1, 1, 3): 151}, (2, 1, 1, 5): {(1, 1, 1): 206, (2, 1, 1, 3): 130}, (6, 1, 1, 5, 9, 11, 7, 3): {(4, 1, 1, 5, 7, 3): 61}, (2, 2, 0, 1, 4, 5): {(2, 1, 1, 5): 106}, (5, 4, 0, 0, 1, 5, 2, 1, 5, 8, 0, 6, 8, 0, 6, 7, 3, 4, 7, 3, 0, 0): {(4, 1, 1, 5, 7, 3): 62}, (4, 3, 0, 0, 1, 5, 0, 1, 3, 0, 7, 3, 0, 0): {(3, 2, 0, 1, 5, 7, 3, 0): 60}, (3, 2, 0, 1, 5, 9, 3, 0): {(3, 2, 0, 1, 5, 7, 3, 0): 62}, (2, 1, 1, 3): {(1, 1, 1): 382}, (3, 1, 1, 5, 9): {(1, 1, 1): 71, (2, 1, 1, 5): 98}, (4, 1, 1, 5, 7, 9): {(4, 1, 1, 5, 7, 3): 68}, (3, 2, 1, 1, 2, 4, 3, 3): {(2, 1, 1, 3): 132}, (3, 1, 1, 5, 3): {(2, 1, 1, 5): 98}, (4, 3, 0, 0, 0, 1, 5, 9, 11, 7, 3, 0, 0, 0): {(4, 1, 1, 5, 7, 3): 69}, (2, 2, 0, 1, 5, 9): {(1, 1, 1): 88, (2, 1, 1, 5): 106}, (2, 2, 1, 2, 4, 3): {(2, 1, 1, 3): 80, (1, 1, 1): 147}, (2, 2, 1, 2, 4, 5): {(1, 1, 1): 70, (2, 1, 1, 5): 96}, (3, 2, 1, 5, 9, 11, 7, 3): {(3, 2, 0, 1, 5, 7, 3, 0): 62}}
