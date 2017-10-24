import icy
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import sem, linregress

icy.Use_Seaborn()
#import seaborn as sns


from matplotlib.animation import FuncAnimation,ImageMagickWriter
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm

from dynamic_fitness_methods import *


RUNS=500
RUN_TYPE=4

def SetRuns(sr):
    global RUNS
    RUNS=sr
def GetRuns():
    global RUNS
    return RUNS
def SetType(st):
    global RUN_TYPE
    RUN_TYPE=st
def GetType():
    global RUN_TYPE
    return RUN_TYPE
def SetPar(st,sr):
    SetRuns(sr)
    SetType(st)





def PlotMaximalOccupation(mu,I):
    fig, axarr = plt.subplots(2, 2, sharex=True, figsize=(10,8))
    Mu_Sets={32:'0.001563',16:'0.003125',8:'0.006250',4:'0.012500',1:'0.050000'}


    for O,ax in zip([5,25,75,125],axarr.reshape(-1)):
        avg_data=np.empty((RUNS,500))
        for r in xrange(RUNS):
            subfile_name='A2_T20_C200_N10000_Mu{}_O{}_K500_I{}_Run{}'.format(Mu_Sets[mu],O,I,r)
            fitness_import=np.genfromtxt('/rscratch/asl47/Bulk_Run/Modular/{}_Fitness.txt'.format(subfile_name),dtype=np.float64)

            avg_data[r]=fitness_import[:,1]*(fitness_import[:,0]>=1)
            ax.plot([i/(1.*O) for i in xrange(fitness_import.shape[0])],avg_data[r],lw=0.5,alpha=0.6)

            
        ax.plot([i/(1.*O) for i in xrange(0,500)],np.mean(avg_data,axis=0),lw=2,c='k',ls=':')
        
        for t,c in zip(range(3),['r','b','g']):
            slope, intercept, r_value, p_value, std_err = linregress([i/(1.*O) for i in xrange(int(O*.2),int(O*.8))],np.log10(np.mean(avg_data,axis=0)[int(O*t+O*.2):int(O*(t+1)-O*.2)]))
            print t,slope, intercept, r_value, p_value, std_err
            for q in xrange(3):
                ax.plot(np.linspace(q*3+t,q*3+1+t,21),[10**(intercept+slope*i) for i in np.linspace(0,1,21)],lw=2.5,c=c)
    
        ax.set_title(r'$\Omega_{{{0}}}$'.format(O))
        ax.set_yscale('log',nonposy='mask')
        ax.set_xlim([0,8])
        if O ==25 or O==125:
            ax.yaxis.tick_right()

        
    fig.add_subplot(111, frameon=False)
    # hide # TODO: ick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    plt.xlabel("Landscape Changes")
    plt.ylabel("\# of 3 Target fit genotypes")
    fig.suptitle('Initial condition {}'.format(I))
    fig.set_tight_layout(True)
    return fig
    #plt.show(block=False)


def PlotMaximalOccupation2(mu,I):
    fig = plt.figure(figsize=(10,8))
    Mu_Sets={32:'0.001563',16:'0.003125',8:'0.006250',4:'0.012500',1:'0.050000'}



    avg_data=np.empty((RUNS,5000))
    for r in xrange(RUNS):
        subfile_name='A3_T20_C200_N10000_Mu{}_O5000_K5000_I{}_Run{}'.format(Mu_Sets[mu],I,r)
        fitness_import=np.genfromtxt('/rscratch/asl47/Bulk_Run/Modular/{}_Fitness.txt'.format(subfile_name),dtype=np.float64)
        avg_data[r]=fitness_import[:,1]*(fitness_import[:,0]>=1)
        plt.plot(range(5000),avg_data[r],lw=0.5,alpha=0.6)
            
            
    plt.plot(range(5000),np.mean(avg_data,axis=0),lw=2,c='k',ls=':')
    print np.mean(avg_data,axis=0)[:5]
        

                
    plt.title(r'Static (condition {})'.format(I))
    plt.yscale('log',nonposy='mask')
    plt.xlabel('Generations')
    plt.ylabel('Maximally Fit genotypes')
    plt.show(block=False)


def PlotModularityGrid2(A,mu,O):
    fig=plt.figure(figsize=(10,7))
    runs=500.

    modularity_import=np.loadtxt('/rscratch/asl47/Processed/Dynamic/A{}Mu{}O{}_Modular_Data.txt'.format(A,mu,O),dtype=np.uint8).reshape(-1,2)
    


    H, xedges, yedges = np.histogram2d(modularity_import[:,0], modularity_import[:,1], bins=(11, 86),range=((9,20),(14,100)))
    X, Y = np.meshgrid(xedges-0.5, yedges-0.5)
    plt.pcolormesh(X,Y,H.T/runs,norm=LogNorm(10,150),cmap='inferno')#, rasterized=True)
    plt.yscale('log')
    cbar=plt.colorbar()
    cbar.set_label('Relative Frequency', rotation=90)
    
    plt.xlabel('``Active" Genome Size')
    plt.ylabel('Phenotype Size')
    plt.title('Mu {}, {}'.format(mu,'Static' if A==3 else 'Dynamic $\Omega^{{{}}}$'.format(O)))

    fig.set_tight_layout(True)
    #return fig
    plt.show(block=False)
    #return modularity_import


def PlotRobustness(O,R,c_in):

    data_frame=[]
    with open('/rscratch/asl47/Processed/Dynamic/TN_rob_O{}_r{}.txt'.format(O,R)) as data_file:
        for line in data_file:
            if 's' in line:
                continue
            else:
                data_frame.append(np.array([float(i) for i in line.split()],dtype=np.float64))
    plt.errorbar([i/5. for i in xrange(1,len(data_frame)+1)],[np.mean(data_frame[i]) for i in xrange(len(data_frame))],yerr=[sem(data_frame[i]) for i in xrange(len(data_frame))],c=c_in)

    #for i in xrange(8):
    #    plt.plot([1.1+1*i]*2,[0,1],'k')

def PlotAll():
    plt.figure(figsize=(10,7))

    for o,r,c in zip([5]*13+[25]*16+[125]*13,[2,3,4,6,7,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0,1,3,6,7,8,9,10,11,12,13,14,15],['darkgreen']*13+['royalblue']*16+['firebrick']*13):
    #for o,r,c in zip([125]*9+[25]*12+[5]*8,[2,3,4,5,7,10,12,13,15,0,1,3,4,6,7,8,9,10,11,12,15,3,5,6,7,9,12,14,15],['firebrick']*9+['royalblue']*12+['darkgreen']*8):
        PlotRobustness(o,r,c)

    plt.xlabel('Landscape changes')
    plt.ylabel('Averaged Fit Phenotype Deleterious Robustness')
    plt.title(r'$\langle \mu L \rangle$=0.25')
    for c_in,label_in in zip(['firebrick','royalblue','darkgreen'],[r'$\Omega_{125}$',r'$\Omega_{25}$',r'$\Omega_{5}$']):
        plt.errorbar([0],[0],c=c_in,label=label_in)

    plt.legend()
    plt.tight_layout()
        
    
    plt.show(block=False)

    
import scipy as sp
import scipy.optimize


def Expo(x,a,b):
    return a+b*np.log10(x)

def PlotHistogram(data,bins='empty',c_in='hotpink',mark_in='o',label_in='Default'):
    if bins=='empty':
        bins=int(np.sqrt(len(data)+0.5))
    hist,bins=np.histogram(data,bins=np.logspace(2,np.log10(10000),bins))#np.log10(min(data)*0.9)

    xs=np.mean(zip(bins,bins[1:]),axis=1)
    histC=np.cumsum(hist)/(1.*RUNS)


    plt.scatter(xs,histC,c=c_in,marker=mark_in,s=50,label=label_in)
    #return
    
    opt_parms, parm_cov = sp.optimize.curve_fit(Expo, xs[8:],histC[8:], maxfev=1000)
    #plt.plot(xs,[Expo(x,*opt_parms) for x in xs],'r--',label='{:.2f} + {:.2f}*log(x)'.format(*opt_parms))
    plt.plot([10000,10000],[-1,2],'k-',lw=3)
    plt.ylim([-0.05,1.05])

    plt.xlabel(r'Generation')
    plt.ylabel('Solution CDF')
    plt.xscale('log')
    plt.legend()
    plt.show(block=False)
    #return (histC,xs)


def PlotPartialHistogram(data,bins='empty',c_in='hotpink',mark_in='o',label_in='Default',label_on=True):
    if bins=='empty':
        bins=int(np.sqrt(len(data)+0.5))
    hist,bins=np.histogram(data,bins=np.logspace(2,np.log10(25000),bins))#np.log10(min(data)*0.9)
    if label_on:
        return plt.scatter(np.mean(zip(bins,bins[1:]),axis=1),np.cumsum(hist)/(1.*RUNS),c=c_in,marker=mark_in,s=50,label=label_in,zorder=10)
    else:
        return plt.scatter(np.mean(zip(bins,bins[1:]),axis=1),np.cumsum(hist)/(1.*RUNS),c=c_in,marker=mark_in,s=50,zorder=10)

def GetHistCoords(data,bins):
    hist,bin_edges=np.histogram(data,bins=np.logspace(2,np.log10(25000),bins))
    return zip(np.mean(zip(bin_edges,bin_edges[1:]),axis=1),np.cumsum(hist)/(1.*RUNS))


plot_params={
    'A3O25000':{'mark_in':'x','c_in':'darkgreen','label_in':'Static'},
    'A2O5':{'mark_in':'o','c_in':'darkred','label_in':r'Dynamic ($\Omega^2_{5}$)'},
    'A2O25':{'mark_in':'s','c_in':'tomato','label_in':r'Dynamic ($\Omega^2_{25}$ )'},
    'A2O50':{'mark_in':'h','c_in':'gold','label_in':r'Dynamic ($\Omega^2_{50}$)'},
    'A2O75':{'mark_in':'<','c_in':'orchid','label_in':r'Dynamic ($\Omega^2_{75}$)'},
    'A2O100':{'mark_in':'>','c_in':'cornflowerblue','label_in':r'Dynamic ($\Omega^2_{100}$)'},
    'A2O125':{'mark_in':'^','c_in':'navy','label_in':r'Dynamic ($\Omega^2_{125}$)'},
    #'T4O250':{'mark_in':'+','c_in':'olive','label_in':r'Dynamic ($\Omega$ = 250)'},
}

def PlotFilledFractions(mu,needles=[5,50,100,250,500]):
    GEN_LIMIT=25000
    bins=30
    fig=plt.figure(figsize=(10,7))
    data_frame_low={}
    data_frame_high={}
    data_frame_locs=0
    for needle in needles:
        data_dict=LoadExistingData(needle,mu)
        data_dict2=LoadExistingData2(needle,mu)
        sorted_keys=sorted(data_dict.keys())
        sorted_keys[:6]=sorted(sorted_keys[:6],key=lambda x: int(x[3:]))
        #sorted_keys[4:]=sorted(sorted_keys[4:],key=lambda x: int(x[3:]))
        for run_t in sorted_keys:
            data_dict[run_t]+=data_dict2[run_t]
            PlotPartialHistogram(data_dict[run_t],bins,label_on=False,**plot_params[run_t])
            if needle==needles[0]:
                data_frame_locs,data_frame_low[run_t]=zip(*GetHistCoords(data_dict[run_t],bins))

            elif needle==needles[-1]:
                data_frame_locs,data_frame_high[run_t]=zip(*GetHistCoords(data_dict[run_t],bins))
                plt.fill_between(data_frame_locs,data_frame_low[run_t],data_frame_high[run_t],alpha=0.5,label=plot_params[run_t]['label_in'],color=plot_params[run_t]['c_in'],zorder=1)

    plt.plot([25000,25000],[-1,2],'k-',lw=3)
    plt.ylim([-0.05,1.05])
    plt.xlim([100,26000])
    plt.xlabel(r'Generation')
    plt.ylabel('Solution CDF')
    plt.xscale('log')
    plt.title(r'$\langle \mu L \rangle$ = {}, Sustained $g$ = [{}]'.format(4./mu,', '.join(map(str,needles))))
    plt.legend(ncol=2)
    plt.tight_layout()
    plt.show(block=False)
    #return fig

def PlotManyHistograms(needle,mu,bins=35):
    data_dict=LoadExistingData(needle,mu)
    if type(data_dict)==int:
        needle=data_dict
        data_dict=LoadExistingData(needle,mu)
        
    SetRuns(500)
    needle_list=[5,25,50,100,500]#1,10,25,50,75,100,300]

    
    sorted_keys=sorted(data_dict.keys())
    sorted_keys[1:4]=sorted(sorted_keys[1:4],key=lambda x: int(x[3:]))
    sorted_keys[4:]=sorted(sorted_keys[4:],key=lambda x: int(x[3:]))
    fig=plt.figure(figsize=(10,10))
    scatters=[]
    for run_t in sorted_keys:
        print run_t
        scatters.append(PlotPartialHistogram(data_dict[run_t],bins,**plot_params[run_t]))
        
    plt.plot([5000,5000],[-1,2],'k-',lw=3)
    plt.ylim([-0.05,1.05])
    plt.xlabel(r'Generation')
    plt.ylabel('Solution CDF')
    plt.xscale('log')
    plt.title('Mu: {}'.format(mu))
    time_text = plt.text(70, 0.25, 'needle: 1',horizontalalignment='left', verticalalignment='center',)
    plt.legend(ncol=2)
    plt.tight_layout()
    #plt.show(block=False)
    #return

    def update_plot(i,needles,scatter_list):
        data_dict=LoadExistingData(needles[i],mu)
        time_text.set_text('needle: {}'.format(needles[i]))
        for scatter,key in zip(scatter_list,sorted_keys):
            scatter.set_offsets(GetHistCoords(data_dict[key],bins))
        
    def init():
        pass

    anim = FuncAnimation(fig, update_plot,init_func=init,frames=len(needle_list), interval=1500, blit=False,fargs=(needle_list,scatters),repeat=True)
    #writer = ImageMagickWriter(fps=0.75)
    #anim.save('Solutions_Slow.gif', writer=writer)
    
    plt.show(block=False)

    







def PlotFitnessOverTime(x,needle_length=50):
    counts=0
    needle=np.array([1]*needle_length,dtype=np.float64)
    plt.figure()
    c=['r','g','b','k']
    for r in xrange(4*x,4*x+4):
        subfile_name='Modular{}_T20_C200_N1000_Mu0.003125_B5000_Run{}'.format(RUN_TYPE,r)
        fitness_import=np.genfromtxt('/rscratch/asl47/Bulk_Run/Modular/{}_Fitness.txt'.format(subfile_name),dtype=np.float64)
        plt.plot(xrange(1,fitness_import.shape[0]+1),fitness_import[:,0],ls='--',label=r,c=c[r-4*x],alpha=0.6)
        plt.plot(xrange(1,fitness_import.shape[0]+1),fitness_import[:,1],ls='',marker='o',alpha=0.6,c=c[r-4*x])


        haystack=search_sequence_numpy(fitness_import[:,1],needle)
        if haystack:
            print "found on {} at {}".format(r,haystack)
            counts+=1
        else:
            print "not found on {}".format(r)

    plt.legend()
    plt.xlabel('generation')
    plt.ylabel('fitness')
    plt.show(block=False)
    return counts






def roulette_selection(weights,N):
    sorted_indexed_weights = sorted(enumerate(weights), key=operator.itemgetter(1))
    indices, sorted_weights = zip(*sorted_indexed_weights)
    tot_sum=sum(sorted_weights)
    prob = [x*1./tot_sum for x in sorted_weights]
    cum_prob=np.cumsum(prob)
   
    for i in xrange(N):
        random_num=random.random()
        for index_value, cum_prob_value in zip(indices,cum_prob):
            if random_num < cum_prob_value:
                yield weights[index_value]
                break
            
from collections import Counter
import random
import operator
def getr(p,N):
    t=list(roulette_selection(p,N))
    
    print Counter(t)
    t.sort()
    
    return t
    
    
def Iterateit(X,N):
    next_pop=[1]*X+[2]*(N-X)
    fracs=np.empty((10,2))
    fracs[0][0]=next_pop.count(1)
    fracs[0][1]=next_pop.count(2)
    
    for i in xrange(9):
        next_pop2=getr(next_pop,N)
        fracs[i+1][0]=next_pop2.count(1)
        fracs[i+1][1]=next_pop2.count(2)
        next_pop=next_pop2
    plt.figure()
    plt.plot(range(10),fracs[:,0],'r')
    plt.plot(range(10),fracs[:,1],'b')
    plt.show(block=False)
        
def g(n,c):
    return 2./(2+(2./c-2)*2**n)

def PlotModularityGrid(A,mu,O):

    data=np.loadtxt('/rscratch/asl47/Processed/Dynamic/A{}Mu{}O{}_Modular_Data.txt'.format(A,mu,O),dtype=np.uint8).reshape(-1,2)
    ymin, ymax=data[:,1].min(),data[:,1].max()
    xmin, xmax=data[:,0].min(),data[:,0].max()
    #Define grid for subplots
    gs = gridspec.GridSpec(2, 2, width_ratios=[3, 1], height_ratios = [1, 4])
    gs.update(wspace=0.07, hspace=0.07)
    #Create scatter plot
    fig = plt.figure(figsize=(10,7))
    ax = plt.subplot(gs[1, 0],frameon = True)
    H, xedges, yedges = np.histogram2d(data[:,0], data[:,1], bins=(12, 86),range=((9,21),(14,100)),normed=True)
    X, Y = np.meshgrid(xedges, yedges)
    ax.pcolormesh(X,Y,H.T,norm=LogNorm(),cmap='inferno')#, rasterized=True)
    
    #ax.grid(True)
    
    #Create Y-marginal (right)
    axr = plt.subplot(gs[1, 1], sharey=ax, frameon = False,xticks = [] ) #xlim=(0, 1), ylim = (ymin, ymax) xticks=[], yticks=[]
    axr.hist(data[:,1], color = '#5673E0', orientation = 'horizontal', normed = True,bins=86,range=(14,100))
    axr.plot([0,1],[np.mean(data[:,1])]*2,c='r',lw=2,ls=':')
    axr.annotate('{:.2f}'.format(np.mean(data[:,1])), xy=(1,np.mean(data[:,1])),fontsize=12)
    axr.set_xscale('log')
    axr.axis('off')


    #Create X-marginal (top)
    axt = plt.subplot(gs[0,0], sharex=ax,frameon = False,yticks = [])# xticks = [], , ) #xlim = (xmin, xmax), ylim=(0, 1)
    axt.hist(data[:,0], color = '#5673E0', normed = True,bins=12,range=(9,21))
    axt.plot([np.mean(data[:,0])]*2,[0,1],c='r',lw=2,ls=':')
    axt.annotate('{:.2f}'.format(np.mean(data[:,0])), xy=(np.mean(data[:,0]), 1),fontsize=12)
    axt.set_yscale('log')

    axt.axis('off')

    adjust_spines(ax,['top','right'])
    adjust_spines(axt,['bottom'])

    ax.set_xlabel('G Size')
    ax.set_ylabel('P Size')
    locs = ax.xaxis.get_ticklocs()
    ax.set_xticks(locs+0.5)
    ax.set_xticklabels([int(i) for i in locs])
    ax.set_xlim([8.7,21])
    #adjust_spines(axr,['left'])
    #ax.set_yscale('log')
    fig.text(0.7,0.80,'Mu {}\n{}'.format(mu,'Static' if A==3 else 'Dynamic $\Omega^{{{}}}$'.format(O)))

    
    
    try:
        ax.set_title(kwargs['title'])
        ax.set_xlabel(kwargs['xlabel'])
        ax.set_ylabel(kwargs['ylabel'])
    except:
        pass
    
    #Bring the marginals closer to the scatter plot
    #fig.set_tight_layout(True)
    return fig
    #plt.show(block=False)

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))  # outward by 10 points
            spine.set_smart_bounds(True)
    return

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])
