import icy
import matplotlib.pyplot as plt
import numpy as np

icy.Use_Seaborn()
#import seaborn as sns


from matplotlib.animation import FuncAnimation,ImageMagickWriter

RUNS=500
RUN_TYPE=3

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



def AddNewDataToFile(classifier,data_in,needle,mu):
    with open('/rscratch/asl47/Processed/Dynamic/Evolution_Solution_Times_Mu{}_Needle{}.txt'.format(mu,needle), "a") as data_file:
        data_file.write(classifier+': ')
        for value in data_in:
            data_file.write('{} '.format(value))
        data_file.write('\n')

def AddBulkDataToFile(Ts,Os,needle,mu):
    keys=['T{}{}'.format(t,'' if t==3 else 'O{}'.format(o)) for (t,o) in zip(Ts,Os)]
    SetRuns(500)
    
    for key in keys:
        SetType(int(key[1]))
        data=MaxFraction(needle,mu,'' if key[1]=='3' else int(key[3:]))
        AddNewDataToFile(key,data,needle,mu)
        
def LoadExistingData(needle,mu):
    data_dict={key:[int(v) for v in values.split()] for (key,values) in [line.rstrip('\n').split(':') for line in open('/rscratch/asl47/Processed/Dynamic/Evolution_Solution_Times_Mu{}_Needle{}.txt'.format(mu,needle))]}
    return data_dict

    
def MaxFraction(needle_length=30,mu=4,O=3,slicer=0):
    occ=0
    bounce=0
    firsts=[]
    additional=''
    if RUN_TYPE==4 or RUN_TYPE==5:
        additional='_O{}'.format(O)

    Mu_Sets={32:'0.001563',16:'0.003125',8:'0.006250',4:'0.012500',1:'0.050000'}
    #Mu_Files={4:'',2:'Mu2',1:'Mu1'}
    
    
    needle=np.array([1]*needle_length,dtype=np.float64)
    t=0
    for r in xrange(RUNS):
        subfile_name='Modular{}_T20_C200_N500_Mu{}{}_K15000_Run{}'.format(RUN_TYPE,Mu_Sets[mu],additional,r)
        fitness_import=np.genfromtxt('/scratch/asl47/Data_Runs/Dynamic_2/T{}Mu{}{}/{}_Fitness.txt'.format(RUN_TYPE,mu,additional[1:],subfile_name),dtype=np.float64)
        #maxes.append(max(fitness_import[:,1]))
        if max(fitness_import[:,slicer])>=1:
            haystack=search_sequence_numpy(fitness_import[:,slicer],needle)
            if haystack>0:
                occ+=1
                firsts.append(haystack)


    print "counted runs: ",t
    print "seen {} times for a fraction of {}".format(occ,occ*1./RUNS)
    #print "bounced {} times".format(bounce)

    return firsts

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
    hist,bins=np.histogram(data,bins=np.logspace(2,np.log10(15000),bins))#np.log10(min(data)*0.9)
    if label_on:
        return plt.scatter(np.mean(zip(bins,bins[1:]),axis=1),np.cumsum(hist)/(1.*RUNS),c=c_in,marker=mark_in,s=50,label=label_in,zorder=10)
    else:
        return plt.scatter(np.mean(zip(bins,bins[1:]),axis=1),np.cumsum(hist)/(1.*RUNS),c=c_in,marker=mark_in,s=50,zorder=10)

def GetHistCoords(data,bins):
    hist,bin_edges=np.histogram(data,bins=np.logspace(2,np.log10(15000),bins))
    return zip(np.mean(zip(bin_edges,bin_edges[1:]),axis=1),np.cumsum(hist)/(1.*RUNS))


plot_params={
    'T3':{'mark_in':'x','c_in':'darkgreen','label_in':'Static'},
    'T4O5':{'mark_in':'o','c_in':'darkred','label_in':r'Dynamic ($\Omega^2$ = 5)'},
    'T4O25':{'mark_in':'s','c_in':'tomato','label_in':r'Dynamic ($\Omega^2$ = 25)'},
    'T4O125':{'mark_in':'h','c_in':'gold','label_in':r'Dynamic ($\Omega^2$ = 125)'},
    'T5O5':{'mark_in':'<','c_in':'skyblue','label_in':r'Dynamic ($\Omega^1$ = 5)'},
    'T5O25':{'mark_in':'>','c_in':'slateblue','label_in':r'Dynamic ($\Omega^1$ = 25)'},
    'T5O125':{'mark_in':'^','c_in':'navy','label_in':r'Dynamic ($\Omega^1$ = 125)'},
    'T4O250':{'mark_in':'+','c_in':'olive','label_in':r'Dynamic ($\Omega$ = 250)'},
}

def PlotFilledFractions(mu,needles=[5,25,50,100,500]):
    bins=30
    plt.figure()
    data_frame_low={}
    data_frame_high={}
    data_frame_locs=0
    for needle in needles:
        data_dict=LoadExistingData(needle,mu)
        sorted_keys=sorted(data_dict.keys())
        sorted_keys[1:4]=sorted(sorted_keys[1:4],key=lambda x: int(x[3:]))
        sorted_keys[4:]=sorted(sorted_keys[4:],key=lambda x: int(x[3:]))
        for run_t in sorted_keys:
            PlotPartialHistogram(data_dict[run_t],bins,label_on=False,**plot_params[run_t])
            if needle==needles[0]:
                data_frame_locs,data_frame_low[run_t]=zip(*GetHistCoords(data_dict[run_t],bins))

            elif needle==needles[-1]:
                data_frame_locs,data_frame_high[run_t]=zip(*GetHistCoords(data_dict[run_t],bins))
                plt.fill_between(data_frame_locs,data_frame_low[run_t],data_frame_high[run_t],alpha=0.5,label=plot_params[run_t]['label_in'],color=plot_params[run_t]['c_in'],zorder=1)

    plt.plot([15000,15000],[-1,2],'k-',lw=3)
    plt.ylim([-0.05,1.05])
    plt.xlim([100,16000])
    plt.xlabel(r'Generation')
    plt.ylabel('Solution CDF')
    plt.xscale('log')
    plt.title(r'$\langle \mu L \rangle$ = {}, Sustained $g$ = [{}]'.format(4./mu,', '.join(map(str,needles))))
    plt.legend(ncol=2)
    plt.tight_layout()
    plt.show(block=False)
    

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

    





def search_sequence_numpy(arr,seq):
    """ Find sequence in an array using NumPy only.

    Parameters
    ----------    
    arr    : input 1D array
    seq    : input 1D array

    Output
    ------    
    Output : 1D Array of indices in the input array that satisfy the 
    matching of input sequence in the input array.
    In case of no match, empty list is returned.
    """

    # Store sizes of input array and sequence
    Na, Nseq = arr.size, seq.size

    # Range of sequence
    r_seq = np.arange(Nseq)

    # Create 2D array of sliding indices across entire length of input array.
    # Match up with the input sequence & get the matching starting indices.
    M = (arr[np.arange(Na-Nseq+1)[:,None] + r_seq] == seq).all(1)
    # Get the range of those indices as final output
    if M.any>0:
        indices=np.where(np.convolve(M,np.ones((Nseq),dtype=int))>0)[0]
        if indices.shape[0]>0:
            return indices[0]
        else:
            return 0
    else:
        return 0
        return [] 

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


#nice -n 15 ./TileEvolution -Q -T 20 -C 199 -N 1000 -U 1 -K 1 -B 5000 -D 150 -O 25 -M 16 -V 0

