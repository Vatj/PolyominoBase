import icy
icy.Use_Seaborn()
import matplotlib.pyplot as plt
from scipy.stats import binom
from scipy.stats import sem
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as color


################
###ROBUSTNESS###
################

COLOURS=[6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 38, 42, 50, 60, 70, 84, 100, 150, 250, 500]
#COLOURS=[6, 8, 10, 12, 16, 20, 24, 30, 36, 42, 50, 60, 80, 100, 140, 200, 300, 500]
MuLs=np.arange(0.7,1.50001,0.00025)

def Plot_Topology_Robustness(R_T,runs):
    
    data_points=np.zeros((9,len(COLOURS),runs),dtype=np.float64)
    N_runs=-1
    N_mut=-1
    Colours=[]
        
    for run in xrange(runs):
        #lines=[line.rstrip('\n') for line in open('/rscratch/asl47/Topology_Robustness_{}_R{}.txt'.format(R_T,run+1))]
        #lines=[line.rstrip('\n') for line in open('/rscratch/asl47/Topology_Robustness_R{}_2.txt'.format(run+1))]
        lines=[line.rstrip('\n') for line in open('/rscratch/asl47/Topology_Robustness_{}_R{}.txt'.format(R_T,run+1))]
        for line in lines:
            if 'runs' in line:
                N_runs=int(line.split()[-2])
            elif 'N_mutations:' in line:
                N_mut= int(line.split()[-1])
            else:
                for C,R in enumerate(line.split(' ')[3::4]):
                    data_points[N_mut][int(C)][run]=float(R)/N_runs


    pcs = [cm.inferno(x) for x in np.linspace(0, 0.9, len(COLOURS)) ]
    mean_data=np.mean(data_points,axis=2,dtype=np.float64)
    MuLs_local=np.logspace(-2,np.log10(8),1000)
    binom_Set=binom(8,MuLs/8.)
    binom_local=binom(8,MuLs_local/8.)
    c_r={}

    #for C,colour_data in enumerate(mean_data.T):
    #    plt.plot(range(9),colour_data,color=pcs[C],ls='-')
        #plt.plot(MuLs_local,np.sum(binom_local.pmf(T)*(1-mean_data[T][C]) for T in xrange(9)),color=pcs[C],ls=':')
    
    f,(ax,ax2)=plt.subplots(2,1)
    for C,colour_data in enumerate(mean_data.T):
        ax.plot(range(9),colour_data,color=pcs[C])
        ax2.plot(MuLs_local,np.sum(binom_local.pmf(T)*(1-mean_data[T][C]) for T in xrange(9)),color=pcs[C])
        c_r[C]=np.sum(binom_Set.pmf(T)*(1-mean_data[T][C]) for T in xrange(9))

    sm = plt.cm.ScalarMappable(cmap=cm.inferno, norm=plt.Normalize(vmin=0, vmax=0.9))
    sm._A = []

    cax = f.add_axes([0.1, 0.625, 0.5, 0.05])

    cbar=f.colorbar(sm, cax=cax, orientation='horizontal')
    cbar.set_ticks(np.arange(0,0.9,1./(1+len(COLOURS))))
    cbar.set_ticklabels(COLOURS)
    cbar.ax.set_title('Colours')
    ax.set_yscale('log',nonposy='mask')

    ax.set_ylim((10**-8,1))
    ax.set_xlabel(r'$N_{\textrm{mutations}}$',fontsize=18)
    ax.set_ylabel(r'$\chi_{\mathrm{robust}}$',fontsize=18)
    ax.set_title(r'\textbf{{{}}}'.format('Thick Cross' if R_T=='TC' else 'Small Cross'),fontsize=24)

    ax2.set_xlabel(r'$\langle \mu L\rangle$',fontsize=18)
    ax2.set_ylabel(r'$\chi_{\mathrm{del}}$',fontsize=18)
    ax2.set_xscale('log')
    ax2.set_yscale('log')

    

    plt.show(block=False)
    return c_r

def Plot_Other(z_in,N=100,alpha=0.95,mark='none'):
    print "For a population of {}, with alpha {}".format(N,alpha)
    COLOURS=[6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 38, 42, 50, 60, 70, 84, 100, 150, 250, 500]
    predicted_divergence=[]
    for C in z_in.keys():
        binom_Set=binom(N,z_in[C])
        survivals=binom_Set.cdf(N/2-1)
  
        mu_index=np.argmax(survivals<alpha)
        #print mu_index
        predicted_divergence.append(MuLs[mu_index])

 
    plt.plot(COLOURS,predicted_divergence,ls=':',color='royalblue',label=r'Pred $\alpha$={:.2f}'.format(alpha),marker=mark)
    
    return predicted_divergence

def NBD():
    plt.figure()
    for N in [100,1000,10000]:
        binomD=binom(N,.8)
        binomD2=binom(N,.4)
        xs=np.arange(0,N+1,1,dtype=np.float64)
        xp=xs/N
        plt.plot(xp,binomD.cdf(xs),label=N)

    plt.plot([-1,2],[.5,.5],'k--')
    plt.xlim([0,1])
    plt.legend()
    plt.show(block=False)
    
#####################
###DIVERGENCE PLOT###
#####################
def func(x, a, b, c):
    return a *(1- np.exp(-b * x,dtype=np.float64)) + c


from scipy.optimize import curve_fit
def Plot_Divergence_Points(Z_in,R_T):


    #sy=np.linspace(100,1000,s_cut,dtype=np.int64)
    sy=np.logspace(2,4,3,dtype=np.int64)
    #sx=np.array([6, 8, 10, 12, 16, 20, 24, 30, 36, 42, 50, 60, 80, 100, 140, 200, 300, 500],dtype=np.int64)
    sx=np.array([6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 38, 42, 50, 60, 70, 84, 100, 150, 250, 500],dtype=np.int64)
    y, x = np.meshgrid(sx,sy)
    z=np.zeros((sy.shape[0],sx.shape[0]))

    #N_med=np.zeros((sy2.shape[0],sx2.shape[0]))
    #N_err=np.zeros((sy2.shape[0],sx2.shape[0]))

    N_mus=[{} for i in xrange(3)]
    N_med=[{} for i in xrange(3)]
    N_err=[{} for i in xrange(3)]
    N_err2=[{} for i in xrange(3)]
    
    N_runs={10:5000,100:2500,1000:1000,10000:250}


    for i,N in enumerate(sy):
        for j,C in enumerate(sx):
            Fitness_Packed=np.array([[float(value) for value in line.rstrip('\n').split() if ':' not in value] for line in open('/rscratch/asl47/Bulk_Run/Evolution_T2_C{}_N{}_R{}_M4_Targ{}_Thresh10_Reg_0_Run1_Size.txt'.format(C,N,N_runs[N],12 if R_T=='TC' else 5))])

            z[i][j]=Fitness_Packed[np.argmax(np.median(Fitness_Packed[:,1:],axis=1)>0),0]

            
            N_mus[i][j]=Fitness_Packed[:,0]
            N_med[i][j]=np.median(Fitness_Packed[:,1:250],axis=1)
            Fitness_Packed[Fitness_Packed>=10]=np.nan
            N_err[i][j]=sem(Fitness_Packed[:,1:250],axis=1,nan_policy='omit')
            N_err2[i][j]=np.nanstd(Fitness_Packed[:,1:250],axis=1)
            #z_err[i][j]=Fitness_Packed[np.argmax(np.median(Fitness_Packed[:,1:],axis=1)>0),0]

    #for i,row in enumerate(z):
    #    plt.plot(sx,row,label='T={}'.format(targ),marker='D')#label='N={}'.format(1000 if i else 100),alpha=1)
    #return
    #f,ax=plt.subplots(1)
    #im = ax.imshow(z, cmap='plasma',interpolation='none',origin='lower', norm=color.LogNorm(vmin=z.min(), vmax=z.max()))
    #ax.set_xticks(xrange(sx.shape[0]))
    #ax.set_yticks(xrange(sy.shape[0]))
    #ax.set_xticklabels(sx)
    #ax.set_yticklabels(sy)
    #plt.colorbar(im)

    f2,ax2=plt.subplots(1)
    for i,row in enumerate(z):
        plt.plot(sx,row,label='N={}'.format(int(sy[i])),alpha=1)
    

    popt, pcov = curve_fit(func, sx, z[0])
    plt.plot(sx, func(sx, *popt), 'k--',label='Exp fit')
    #print popt

    ax2.set_xlabel(r'Colours',fontsize=18)
    ax2.set_ylabel(r'$\langle \mu L\rangle$',fontsize=18)
    #ax2.set_title(r'\textbf{Divergence Point ({{}})}'.format(R_T),fontsize=24)
    #ax2.set_yscale('log')
    ax2.set_xscale('log')

    q_out=Plot_Other(Z_in,N=100,alpha=0.45,mark='o')
    q_out=Plot_Other(Z_in,N=100,alpha=0.5,mark='h')
    #plt.plot([6, 8, 10, 12, 16, 20, 24, 30, 36, 42, 50, 60, 80, 100, 140, 200, 300, 500],[q+0.025 for q in q_out],color='hotpink',ls=':',label='Pred (offset)',marker='D')
    #Plot_Other(Z_in,alpha=0.92,mark='s')
    #Plot_Other(Z_in,alpha=0.95,N=1000,mark='v')
    plt.legend(ncol=2,loc='lower right')

    plt.show(block=False)
    return "early"
    #f3,ax3=plt.subplots(1)
    cs=['g','b','r','k']
    for j in xrange(20):
        plt.figure()
        for i in xrange(sy.shape[0]):
            plt.scatter(N_mus[i][j],N_err[i][j],c=cs[i])
            plt.scatter(N_mus[i][j],N_err2[i][j],c=cs[i],marker='s')

    
    plt.show(block=False)
    #return z
