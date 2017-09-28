import Degeneracy_Recovery as DR

import matplotlib.colors as colors
import matplotlib.cm as cm
import glob
from scipy.stats import linregress
from math import floor
from scipy.stats import binom

import seaborn as sns
def Use_Seaborn():
    sns.set_context("paper",font_scale=2.2)
    sns.set_style("white",rc={"xtick.major.size": 8, "ytick.major.size": 8,"xtick.minor.size":5, "ytick.minor.size": 5,"axes.linewidth": 2,"axes.edgecolor":"darkgray","font.size":8,"axes.titlesize":8,"axes.labelsize":5})
    

import re
import numpy as np
import matplotlib.pyplot as plt
def plotEvolution(files,labels=False,title=""):
    colours=['firebrick','cornflowerblue','seagreen','orchid','slategrey','black','gold','orange','aqua','mediumblue','lime','deeppink','darksage']*10
    if not labels:
        labels=[""]*len(files)
        
    for n,fileN in enumerate(files):
        
        lines=[line.rstrip('\n') for line in open("Evolution_Runs/{}.txt".format(fileN))]
        pp=[filter(None,re.split('[a-z]|:|_|[A-Z]',line)) for line in lines]
        pp_float=np.array([[float(i) for i in  p] for p in pp])


        #plt.plot(pp_float[:,0],pp_float[:,2],color=colours[n],label=labels[n],lw=2)
        #plt.plot(pp_float[:,0],pp_float[:,1],color=colours[n],ls='--',mfc=None,lw=3)
        #plt.plot(pp_float[:,0],pp_float[:,2],color=colours[n],ls='--',mfc=None,lw=1,alpha=0.6)
        #plt.plot(pp_float[:,0],pp_float[:,3],color=colours[n],ls='--',mfc=None,lw=1,alpha=0.6)
        #plt.fill_between(pp_float[:,0], pp_float[:,2], pp_float[:,3],facecolor='',edgecolor=colours[n],alpha=0.3)#, hatch = '//')
        plt.plot(pp_float[:,0],pp_float[:,4],color=colours[n],ls='-',mfc=None,lw=3,label=labels[n])
        #plt.plot(pp_float[:,0],pp_float[:,5],color=colours[n],ls='-',mfc=None,lw=1,alpha=0.6)
        #plt.plot(pp_float[:,0],pp_float[:,6],color=colours[n],ls='-',mfc=None,lw=1,alpha=0.6)
        #plt.fill_between(pp_float[:,0], pp_float[:,5], pp_float[:,6],facecolor='',edgecolor=colours[n],alpha=0.3)#, hatch = '\\') 


    plt.legend(loc='lower left',ncol=2)
    plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel(r'$\mu L$')
    plt.ylabel(r'Generations')
    plt.title(title)#r'$S_{2,8}$  with $S^{*}=12$')
    sns.despine()
    plt.show()


#################
###LENGTH PLOT###
#################
def Calculate_Genome_Data(Ts):
    for i,T in enumerate(Ts):
        data =np.loadtxt('/scratch/asl47/Data_Runs/Regulation_Lengths_17_06_25/Reg_Lengths_Full_T{}.txt'.format(T),dtype=np.float32)
        mean_Data=np.mean(data,axis=0,dtype=np.float64)
        std_Data=np.std(data,axis=0,dtype=np.float64)
        np.savetxt('/scratch/asl47/Data_Runs/Regulation_Lengths_17_06_25/Reg_Lengths_Full_T{}_Mean.txt'.format(T),mean_Data,delimiter=' ')
        np.savetxt('/scratch/asl47/Data_Runs/Regulation_Lengths_17_06_25/Reg_Lengths_Full_T{}_STD.txt'.format(T),std_Data,delimiter=' ')     
        
def Plot_Genome_Length(Ts):
    Use_Seaborn()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    cs=['dodgerblue','firebrick','darkgreen','darkmagenta','orange']
    lss=['-',':','-.']*3
    alps=[1,0.8,0.5]*3
    L_FIX=1000

    f=plt.figure()
    L_ax = plt.subplot(111)
    
    #F_ax = plt.subplot(111)
    
    for i,T in enumerate(Ts):
        avg=defaultdict(list)
        mfits=defaultdict(list)
        mean_Data =1-np.loadtxt('/scratch/asl47/Data_Runs/Regulation_Lengths_17_06_25/Reg_Lengths_Full_T{}_Mean.txt'.format(T),dtype=np.float64)/T
        std_Data =np.loadtxt('/scratch/asl47/Data_Runs/Regulation_Lengths_17_06_25/Reg_Lengths_Full_T{}_STD.txt'.format(T),dtype=np.float64)/T               

        L_ax.plot([1,2],[0,mean_Data[0]],c=cs[i],zorder=10,lw=3)
        
        #L_ax.plot(xrange(1,L_FIX+1),mean_Data[::2]-std_Data[::2],c=cs[i],zorder=5,lw=0.5)
        #L_ax.plot(xrange(1,L_FIX+1),mean_Data[::2]+std_Data[::2],c=cs[i],zorder=5,lw=0.5)
        
        L_ax.plot(xrange(2,L_FIX+2),mean_Data[::2],c=cs[i],zorder=10,lw=3)
        L_ax.fill_between(xrange(2,L_FIX+2), mean_Data[::2]-std_Data[::2],mean_Data[::2]+std_Data[::2],facecolor=cs[i],alpha=0.3,interpolate=True,zorder=1)
        L_ax.plot(xrange(2,L_FIX+2),mean_Data[::2]-std_Data[::2],c=cs[i],zorder=5,lw=0.5)
        L_ax.plot(xrange(2,L_FIX+2),mean_Data[::2]+std_Data[::2],c=cs[i],zorder=5,lw=0.5)
        
        #F_ax.plot(xrange(1,L_FIX+1),mean_Data[1::2],'-',c=cs[i],zorder=10,lw=3)
        #F_ax.fill_between(xrange(1,L_FIX+1), mean_Data[1::2]-std_Data[1::2],np.minimum(mean_Data[1::2]+std_Data[1::2],[1]),facecolor=cs[i],alpha=0.3,interpolate=True,zorder=1)
        #F_ax.plot(xrange(1,L_FIX+1),mean_Data[1::2]-std_Data[1::2],c=cs[i],zorder=5,lw=0.5)
        #F_ax.plot(xrange(1,L_FIX+1),np.minimum(mean_Data[1::2]+std_Data[1::2],[1]),c=cs[i],zorder=5,lw=0.5)
        
        print "T", np.mean(mean_Data[500::2]), "+-", np.mean(std_Data[500::2])/np.sqrt(100000)

    #F_ax.xaxis.tick_top()
    L_ax.set_xscale('log')
    L_ax.set_xlabel('Generation')
    #L_ax.set_ylabel(r'$\langle N_{T}^{\mathrm{active}} \rangle$')
    L_ax.set_ylabel(r'$\frac{L}{L}^{\small \textrm{off}}$')
    #L_ax.text(1,3.2,r'$N_{T}=4$')
    #L_ax.text(1,4.2,r'$N_{T}=5$')
    #L_ax.text(1,5.2,r'$N_{T}=6$')
    #F_ax.set_ylabel(r'$\langle f_{\mathrm{max}} \rangle$')
    plt.show(block=False)
 
##################
###FRAC OC PLOT###
##################
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
import numpy.ma as ma

def Plot_Frac_Oc(Ts,R):
    Use_Seaborn()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    frac_matrix=np.zeros((10,10))
    frac_matrix2=np.zeros((10,10))
    
    for T in xrange(1,Ts+1):
        lines=[line.rstrip('\n') for line in open('/rscratch/asl47/Frac_Oc/Fractional_Occupation_T{}_Run_{}.txt'.format(T,R))]
        for line in lines:
            frac_matrix[int(line.split(' ')[1])/4 -1][T-1] +=float(line.split(' ')[2])/250000000
            frac_matrix2[int(line.split(' ')[1])/4 -1][T-1] +=float(line.split(' ')[3])/250000000


    masked_frac_M=ma.masked_where(frac_matrix <=0, frac_matrix) 
    masked_frac_M2=ma.masked_where(frac_matrix2 <=0, frac_matrix2)
    
    f, (F_ax,P_ax) = plt.subplots(1,2,sharex=True,sharey=True)
    
    #divider = make_axes_locatable(F_ax)
    #cax = divider.append_axes('bottom', size='5%', pad=0.05)

    minmin=min(np.min(masked_frac_M),np.min(masked_frac_M2))
    maxmax=max(np.max(masked_frac_M),np.max(masked_frac_M2))
    
    f_im = F_ax.pcolor(masked_frac_M, cmap='viridis', edgecolor='black', linestyle=':', lw=1,norm=LogNorm(vmin=minmin/1.1, vmax=maxmax*1.1))
    F_ax.patch.set(hatch='////', edgecolor='black')

    cbar_ax = f.add_axes([0.095, 0.06, 0.885, 0.035])
    
    cb=f.colorbar(f_im,cbar_ax,ticks=np.logspace(int(floor(np.log10(np.min(masked_frac_M)))),0,abs(int(floor(np.log10(np.min(masked_frac_M)))))+1),orientation='horizontal')
    #cb.ax.set_title('Fractional Occupation')
    #divider2 = make_axes_locatable(P_ax)
    #cax2 = divider2.append_axes('bottom', size='5%', pad=0.05)
    p_im = P_ax.pcolor(masked_frac_M2, cmap='viridis', edgecolor='black', linestyle=':', lw=1,norm=LogNorm(vmin=minmin/1.1, vmax=maxmax*1.1))
    P_ax.patch.set(hatch='////', edgecolor='black')
    #cb2=f.colorbar(p_im,cax2,ticks=np.logspace(int(floor(np.log10(np.min(masked_frac_M2)))),0,abs(int(floor(np.log10(np.min(masked_frac_M2)))))+1),orientation='horizontal')
    
 
    F_ax.xaxis.set(ticks=np.arange(0.5, 10), ticklabels=[i for i in xrange(1,11)])
    F_ax.yaxis.set(ticks=np.arange(0.5, 10), ticklabels=[i for i in xrange(4,44,4)])

    P_ax.xaxis.set(ticks=np.arange(0.5, 10), ticklabels=[i for i in xrange(1,11)])
    P_ax.yaxis.set(ticks=np.arange(0.5, 10), ticklabels=[i for i in xrange(4,44,4)])

    F_ax.set_ylabel('Colours')
    f.text(0.54,0.98,r'$N_{T}$',ha='center',va='top')
    
    P_ax.xaxis.tick_top()
    F_ax.xaxis.tick_top()


    plt.figure()
    colors = [ cm.viridis(x) for x in np.linspace(0,1,10)]
    for i in xrange(10):
        plt.plot([j for j in xrange(1,11)], frac_matrix[i,:],c=colors[i])
    #plt.yscale('log',nonposy='mask')
    plt.xlabel(r'$N_{T}$')
    plt.ylabel('B/D Fraction')
    
    plt.show(block=False)
    return frac_matrix,frac_matrix2


            
####################
###EVOLUTION PLOT###
####################
def Plot_Evolution_New(subfolder,files,labels):
    Use_Seaborn()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    colours=['royalblue','firebrick','darkorange','gold','cornflowerblue','darkblue','firebrick','gold','deeppink','seagreen']#['darkcyan','cornflowerblue','mediumblue','silver','slategrey','black','gold','orange','firebrick','orchid','lime','deeppink','darkgreen','seagreen']*2
    
    
    fig, Prob_Ax = plt.subplots()
    Gen_Ax=Prob_Ax.twinx()
    fig2=plt.figure()
    plotss=[]
    for n,fileN in enumerate(files):
        #print n,":",labels[n],fileN
        lines=[line.rstrip('\n') for line in open('/scratch/asl47/Data_Runs/{}/{}.txt'.format(subfolder,fileN))]
        pp=[filter(None,re.split('[a-z]|:|_|[A-Z]',line)) for line in lines]
        pp_float=np.array([[float(i) for i in  p] for p in pp])
        
        
        if 'qqq' in labels[n]:
            plt.plot(pp_float[:,0],pp_float[:,2],color=colours[n],ls='-.',mfc=None,lw=3,label=labels[n])
        else:
            plotss.append(plt.plot(pp_float[:,0],pp_float[:,2],color=colours[n],ls='-',mfc=None,lw=3)[0])
        plotss.append(plt.plot(pp_float[:,0],pp_float[:,1],color=colours[n],ls=':',lw=2,alpha=1)[0])
        Prob_Ax.plot(pp_float[:,0],pp_float[:,3]/1000.,color=colours[n],ls='-',label=labels[n],lw=2)
        Gen_Ax.plot(pp_float[:,0],pp_float[:,4],color=colours[n],ls=':',label=labels[n],lw=0.75)
        
    plotss.append(plt.scatter([1,4,8],[982.438,561.683,553.686],marker='o',s=50,c='goldenrod'))
    plt.xlabel(r'$\langle \mu L \rangle$')
    plt.ylabel(r'Generations')
        
    Gen_Ax.set_yscale('log', nonposy='mask')
    Prob_Ax.set_yscale('log', nonposy='mask')
    Prob_Ax.set_xscale('log', nonposx='mask')
    plt.yscale('log',nonposy='mask')
    plt.xscale('log')

    plt.legend(plotss[:2]+[plotss[-1]]+plotss[2:-1],labels,ncol=2)
    
    Gen_Ax.set_ylabel(r'$\langle g \rangle_{\mathrm{Extinction}}$')
    Prob_Ax.set_ylabel(r'Extinction Probability')
    Prob_Ax.set_xlabel(r'$\langle \mu L \rangle$')
    #plt.title(" Colour, 100 Population Evolution")
    sns.despine(fig2)
    sns.despine(fig,right=False)
    #plt.title("Disjointed Definition")
    plt.show(block=False)

def Plot_Evolution_Reg(subName='Evolution_Regulation_17_06_26'):
    Use_Seaborn()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    fig,ax=plt.subplots(1)
    cs=['dodgerblue','firebrick','darkgreen']
    ms=['o','s','D']
    
    left, bottom, width, height = [0.5, 0.56, 0.35, 0.4]
    Z_ax = fig.add_axes([left, bottom, width, height])
    #ls=[]
    for i,fileN in enumerate(glob.glob('/scratch/asl47/Data_Runs/{}/*.txt'.format(subName))):
        T=int(fileN.split('_')[6][1])
        R=int(fileN.split('_')[14])
        print T,R
        lines=np.array([[float(x) for x in line.rstrip('\n').split()[1:6:2]] for line in open(fileN)])
        #ls.append(lines)
        
        ax.plot(lines[:,0],lines[:,2],markeredgecolor=cs[i],markerfacecolor='none',marker=ms[i],markeredgewidth=1.25,markersize=7,zorder=10,ls='')
        ax.plot(lines[:,0],lines[:,2],c=cs[i],zorder=1,ls='--',alpha=0.6,lw=0.75)

        #ax.plot([0],[0],label=r'$N_T=${}{}'.format(T,', Regulated' if R else ''),markeredgecolor=cs[i],markerfacecolor='none',marker=ms[i],markeredgewidth=1.25,markersize=7,c=cs[i])
        
        Z_ax.plot(lines[:,0],lines[:,1],markeredgecolor=cs[i],markerfacecolor='none',marker=ms[i],markeredgewidth=1.25,markersize=4,zorder=10,ls='')
        
        
        slopeR, interceptR, r_valueR, p_valueR, std_errR = linregress(np.log10(lines[:10,0]),np.log10(lines[:10,1]))
        Z_ax.plot(lines[:,0],[10**(interceptR)*(l**slopeR) for l in lines[:,0]],c=cs[i],zorder=1)
        #print [10**(interceptR)*(l**slopeR) for l in lines[:,1]]
        print "R",slopeR, interceptR, r_valueR**2, std_errR

   
    plt.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    #ax.set_ylim((40,3500))
    #ax.text(.175,195,r'Shape',ha='left',va='bottom',fontsize=18)
    #ax.text(.175,38,r'Size',ha='left',va='bottom',fontsize=18)
    ax.text(.535,190,r"$N_{T}=4'$",ha='left',va='bottom',fontsize=18)
    ax.text(.535,105,r'$N_{T}=5$',ha='left',va='bottom',fontsize=18)
    ax.text(.6,70,r'$N_{T}=4$',ha='right',va='top',fontsize=18)
    #ax.text(2,1300,r'$\langle \tau_{D} \rangle$',ha='right',va='top')
    
    Z_ax.set_xscale('log')
    Z_ax.set_yscale('log')
    Z_ax.text(2,1300,r'$\langle \tau_{D} \rangle$',ha='right',va='top')
    
    ax.set_ylabel('Generations')
    ax.set_xlabel(r'$\langle \mu L \rangle$')
    plt.show(block=False)

    
def Plot_Evolution_Mins(T=2,Cs=[i for i in xrange(6,32,2)]+[i for i in xrange(32,120,4)]+[i for i in xrange(120,200,10)]+[i for i in xrange(200,600,100)] ,N=100,D=10000,M=3,F=12,K=500,R=0):
    Use_Seaborn()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    #plt.figure()

    C_D={6: 2.7435430403264, 10: 2.9776582164672956, 14: 3.169646494364731, 18: 3.228978084086571, 22: 3.364281779315803, 26: 3.3115151379784016, 30: 3.3908287461900177, 40: 3.390571735507829, 50: 3.3814142933461846, 60: 3.4727674134543163, 70: 3.4038281744795076}
    colors = [cm.cubehelix(x) for x in np.linspace(0, 0.8,len(Cs))]
    T_Mins=[]
    D_Mins=[]
    R_xs=[]
    R_Mins=[]
    for i,C in enumerate(Cs):
        lines=np.array([[float(x) for x in line.rstrip('\n').split()[1:6:2]] for line in open('/scratch/asl47/Data_Runs/Evolution_2T_C_17_06_15/Evolution_T{}_C{}_N{}_R{}_M{}_Targ{}_Thresh{}_Reg_{}_Run1_Size.txt'.format(T,C,N,D,M,F,K,R))])
        T_Mins.append(np.min(lines[:,2]))
        D_Mins.append(np.min(lines[:,1]))
        #return lines
        #plt.plot(lines[:,0],lines[:,2],lw=1.5,c=colors[i])
        #plt.plot(lines[:,0],lines[:,2],lw=1.5,marker='*',c=colors[i])

    for fileN in glob.glob('/scratch/asl47/Data_Runs/Random_Walks/*.txt'):
        R_xs.append(int(fileN.split('_')[6]))
        try:
            R_Mins.append([[float(x.split()[0]) for x in line.rstrip('\n').split(':')[2:]] for line in open(fileN)][2])
        except:
            R_Mins.append([[float(x.split()[0]) for x in line.rstrip('\n').split(':')[2:]] for line in open(fileN)][0])
        #print R_Mins[-1][0]
        if R_xs[-1] in C_D:
            R_Mins[-1][0]*=C_D[R_xs[-1]]
            R_Mins[-1][1]*=C_D[R_xs[-1]]*30
        else:
            R_Mins[-1][0]*=C_D[70]
            R_Mins[-1][1]*=C_D[70]*30
    #########
    #Lin fig#
    #########
    R_MinsArr=np.array(R_Mins) 

    fig, T_ax = plt.subplots()
    
    T_ax.plot(Cs,T_Mins,lw=0.75,ls='',marker='s',markersize=7,markeredgecolor='firebrick',label=r'$\langle \tau \rangle$',zorder=10,fillstyle='full',markeredgewidth=1.25,markerfacecolor='white')
    
    T_ax.plot(Cs,D_Mins,lw=0.75,ls='',marker='o',markersize=7,markeredgecolor='dodgerblue',label=r'$\langle \tau_{D} \rangle$',zorder=10,fillstyle='full',markeredgewidth=1.25,markerfacecolor='white')
    #T_ax.plot(Cs,D_Mins,lw=0.75,ls='--',c='dodgerblue',zorder=1,alpha=0.6)
    #S_ax = T_ax.twinx()
    sorted_Walk_X,sorted_Walk_Y=zip(*sorted(zip(R_xs,R_MinsArr[:,0])))
    #print sorted_Walk_Y
    T_ax.plot(*zip(*sorted(zip(R_xs,R_MinsArr[:,0]))),ls='',marker='D',markersize=7,markeredgecolor='darkgreen',label='Random Walk',zorder=10,fillstyle='full',markeredgewidth=1.25,markerfacecolor='white')
    #T_ax.plot(*zip(*sorted(zip(R_xs,R_MinsArr[:,0]))),ls=':',zorder=1,alpha=0.6,c='darkgreen')
    T_ax.errorbar(*zip(*sorted(zip(R_xs,R_MinsArr[:,0]))),yerr=zip(*sorted(zip(R_xs,R_MinsArr[:,1])))[1],fmt='',ls='',ecolor='darkgreen')

    slopeA, interceptA, r_valueA, p_valueA, std_errA = linregress(Cs[-5:],T_Mins[-5:])
    slopeD, interceptD, r_valueD, p_valueD, std_errD = linregress(Cs[-5:],D_Mins[-5:])
    slopeR, interceptR, r_valueR, p_valueR, std_errR = linregress(np.log10(sorted_Walk_X[-5:]),np.log10(sorted_Walk_Y[-5:]))
    print slopeA,slopeD,slopeR
    print interceptA,interceptD,interceptR

    print "A",slopeA, interceptA, r_valueA**2, std_errA
    print "D",slopeD, interceptD, r_valueD**2, std_errD
    print "R",slopeR, interceptR, r_valueR**2, std_errR
        
    T_ax.plot(Cs[-35:],[(C**slopeR)*(10**interceptR) for C in Cs[-35:]],lw=1,ls='--',c='darkgreen',alpha=1,zorder=1)
   
    left, bottom, width, height = [0.09, 0.62, 0.35, 0.32]
    Z_ax = fig.add_axes([left, bottom, width, height])
    Z_ax.plot(Cs,T_Mins,lw=0.75,ls='',marker='s',markersize=6,markeredgecolor='firebrick',label=r'$\langle \tau \rangle$',zorder=10,fillstyle='full',markeredgewidth=1.25,markerfacecolor='white')
    
    Z_ax.plot(Cs,D_Mins,lw=0.75,ls='',marker='o',markersize=6,markeredgecolor='dodgerblue',label=r'$\langle \tau_{D} \rangle$',zorder=10,fillstyle='full',markeredgewidth=1.25,markerfacecolor='white')


    Z_ax.plot(Cs,[(C*slopeA)+interceptA for C in Cs],lw=1,ls='--',c='firebrick',alpha=1,zorder=1)
    Z_ax.plot(Cs,[(C*slopeD)+interceptD for C in Cs],lw=1,ls='--',c='dodgerblue',alpha=1,zorder=1)
    
    Z_ax.yaxis.tick_right()

    T_ax.set_xscale('log')
    T_ax.set_yscale('log')
    #S_ax.set_yscale('log')

    T_ax.set_xlabel('Colours')
    T_ax.set_ylabel('Generations')
    
    T_ax.text(160,160,r'$\langle \tau \rangle$',va='bottom',ha='right')
    T_ax.text(160,25,r'$\langle \tau_{D} \rangle$',va='top',ha='right')
    T_ax.text(81.5,73000,r'Random Walk',va='top',ha='left')

    
    plt.show(block=False)

def temp():
    R_xs=[]
    R_ys=[]
    zz={}
    for fileN in glob.glob('/scratch/asl47/Data_Runs/Random_Walks/*.txt'):#'/rscratch/asl47/Bulk_Walks/*.txt'):
        try:
            n=int(fileN.split('_')[6])
            temx=[[float(x.split()[0]) for x in line.rstrip('\n').split(':')[2:]] for line in open(fileN)]
            print temx[2]
            zz[n]=temx[2][0]
            R_ys.append(temx[1][0]/temx[2][0])
            R_xs.append(int(fileN.split('_')[6]))
            #return temx
            #R_Mins.append([[float(x.split()[0]) for x in line.rstrip('\n').split(':')[2:]] for line in open(fileN)])
        except:
            continue
    plt.plot(*zip(*sorted(zip(R_xs,R_ys))))
    plt.show(block=False)
    return zz
from math import floor,ceil
def Plot_Max_Size():
    xs=range(1,11)
    plt.figure()
    plt.plot(xs,[4*(x**2) for x in xs],ls='--',c='k',alpha=.7,zorder=1)
    #plt.plot(xs,[(2*ceil(x/2.)+floor(x/2.))**2 for x in xs],ls='--',c='c',alpha=.7,zorder=1)
    plt.plot(xs,[4*(x**2) for x in xs],ls='',marker='s',markeredgewidth=1.5,markerfacecolor='none',markeredgecolor='k',zorder=10,markersize=7)
    S_star=np.array([[4,12,20,36,56,72,104,132,152,192],[4,12,20,36,56,72,152,160,176,180]])
    #plt.plot([i for i in xrange(1,1+len(S_star))],S_star,ls='-.',marker='H',c='fuchsia')
    #plt.plot([i for i in xrange(1,1+len(S_star2))],S_star2,ls=':',marker='X',c='coral')
    plt.plot(xs,np.amax(S_star,axis=0),ls='--',c='orangered',alpha=.7,zorder=1)
    plt.plot(xs,np.amax(S_star,axis=0),ls='',marker='D',markeredgewidth=1.5,markerfacecolor='none',markeredgecolor='orangered',zorder=10,markersize=7)


    plt.plot(xs,[4*(x+x*(x/2)-(x/2)**2) for x in xs],ls='--',c='darkmagenta',alpha=.7,zorder=1)
    plt.plot(xs,[4*(x+x*(x/2)-(x/2)**2) for x in xs],ls='',marker='o',markeredgewidth=1.5,markerfacecolor='none',markeredgecolor='darkmagenta',zorder=10,markersize=7)
        
           
        #plt.plot([L for L in xrange(1,x+1)],[4*(x+x*L-L**2) for L in xrange(1,x+1)],ls=':',c=colors[x-1])

    plt.xlabel(r'$N_{T}$')
    plt.ylabel('Size')
    plt.yscale('log')
    plt.xscale('log')
    plt.text(11,215,'Evolved',va='bottom',ha='right',fontsize=14)
    plt.text(6,165,'Upper Bound',va='bottom',ha='right',fontsize=14)
    plt.text(6,50,'Optimal Loop',va='top',ha='left',fontsize=14)
    plt.show(block=False)

    

def Get_Max(T,C,Mu,R):
    maxs=[]
    for r in xrange(R):
        lines=np.array([[float(x) for x in line.rstrip('\n').split()[1:6:2]] for line in open('/rscratch/asl47/Bulk_Run/Fitness/T_{}C_{}N_{}Mu_{}_Osc_{}_Gen_{}_Run_{}.txt'.format(T,C,500,Mu,20000,20000,r))])
        maxs.append(np.max(lines[:,2]))
        if np.max(lines[:,2])==152.0:
            print r
            del maxs[-1]
    print max(maxs)




    
    


######################
###FITNESS EVO PLOT###
######################
def Plot_Fitness_Evolution(T,C,Mus,Osc_Rate,labels,runs,Gen_Limit,maxes=True):
    Use_Seaborn()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    Smooth_Factor=10
    Low_Range=2500
    High_Range=2500
    colours={0:'royalblue',13:'firebrick',15:'darkgreen',17:'goldenrod'}
    singular_alpha=1
    Averaged_Data={}

    Adaptations={}
    Peaks={}

    
    fig = plt.figure()
    Avg_Ax_Low = fig.add_subplot(2, 2, 1)
    Avg_Ax_High = fig.add_subplot(2, 2, 2,sharey = Avg_Ax_Low)
    Max_Ax_Low = fig.add_subplot(2, 2, 3,sharex = Avg_Ax_Low)
    Max_Ax_High = fig.add_subplot(2, 2, 4,sharey = Max_Ax_Low,sharex = Avg_Ax_High)

    fig2, Adapt_Ax = plt.subplots()
    Peak_Ax=Adapt_Ax.twinx()

    Avg_Ax_Low.spines['right'].set_visible(False)
    Avg_Ax_High.spines['left'].set_visible(False)
    Max_Ax_Low.spines['right'].set_visible(False)
    Max_Ax_High.spines['left'].set_visible(False)
    line_caps=[]

    temp_cases=[0]+[n*Osc_Rate for n in xrange(1,Gen_Limit/Osc_Rate)]
    Alt= 2 if maxes else 3
    
    ms={0:'o', 13:'s', 17:'^', 15:'D'}
    zo={0:10, 13:1, 17:5, 15:4}
    for n,Mu in enumerate(Mus):
        Averaged_Data[Mu]=np.zeros([Gen_Limit,4])
        Adaptations[Mu]=[]
        Peaks[Mu]=[]
        #AD={}
        ii=0
        runss=[0,15,17,13]
        for run in [0,13,17]:
  
            Fitness_Packed=np.array([[float(value) for value in line.rstrip('\n').split()[1::2]] for line in open('/rscratch/asl47/Bulk_Run/Fitness/T_{}C_{}N_10000Mu_{}_Osc_{}_Gen_{}_Run_{}.txt'.format(T,C,Mu,Osc_Rate,Gen_Limit,run))])
            #AD[run]=Fitness_Packed
            Averaged_Data[Mu][:,1]+=(Fitness_Packed[:,1]/runs)
            Averaged_Data[Mu][:,Alt]+=(Fitness_Packed[:,Alt]/runs)

            Avg_Ax_Low.plot(Fitness_Packed[:,0][:Low_Range:5],Fitness_Packed[:,1][:Low_Range:5],ls='-',marker=ms[run-ii],lw=.5,alpha=singular_alpha,c=colours[run-ii],zorder=zo[run])
            Max_Ax_Low.plot(Fitness_Packed[:,0][:Low_Range:5],Fitness_Packed[:,Alt][:Low_Range:5],ls='-',lw=.5,marker=ms[run-ii],alpha=singular_alpha,c=colours[run-ii],zorder=zo[run])
            Avg_Ax_High.plot(Fitness_Packed[:,0][-High_Range::5],Fitness_Packed[:,1][-High_Range::5],ls='-',lw=.5,marker=ms[run-ii],alpha=singular_alpha,c=colours[run-ii],zorder=zo[run])
            Max_Ax_High.plot(Fitness_Packed[:,0][-High_Range::5],Fitness_Packed[:,Alt][-High_Range::5],ls='-',lw=.5,marker=ms[run-ii],alpha=singular_alpha,c=colours[run-ii],zorder=zo[run])
             

        line_result,=Avg_Ax_Low.plot(Fitness_Packed[:,0][:Low_Range],Averaged_Data[Mu][:,1][:Low_Range],ls='-',lw=1.75,c=colours[n],alpha=0)
        line_caps.append(line_result)
        #Max_Ax_Low.plot(Fitness_Packed[:,0][:Low_Range],Averaged_Data[Mu][:,Alt][:Low_Range],ls='-',lw=1.75,c=colours[n])
       
        #Avg_Ax_High.plot(Fitness_Packed[:,0][-High_Range:],Averaged_Data[Mu][:,1][-High_Range:],ls='-',lw=1.75,c=colours[n])
        #Max_Ax_High.plot(Fitness_Packed[:,0][-High_Range:],Averaged_Data[Mu][:,Alt][-High_Range:],ls='-',lw=1.75,c=colours[n])

        
        #Adaptation_Times(AD,Osc_Rate,Gen_Limit)
        if Mu=='0.500000' or Mu=='0.050000':
            continue
        for oscillation_L,oscillation_H in zip(temp_cases,temp_cases[1:]):
            Oscillation_Slice=Averaged_Data[Mu][:,Alt][oscillation_L:oscillation_H]
            index=np.argmax(Oscillation_Slice>0.5)
            if index==0 and Oscillation_Slice[index]<0.5:
                Adaptations[Mu].append(Osc_Rate)
                continue
            else:
                Adaptations[Mu].append(index)
                Peaks[Mu].append(max(Oscillation_Slice))
        
        Adapt_Ax.plot(Adaptations[Mu],label=r'N$\mu$L={}'.format(labels[n]),c=colours[n])
        #Peak_Ax.plot(Peaks[Mu],'--',c=colours[n])

    Adapt_Ax.set_xlabel('Oscillation Cycle')
    Adapt_Ax.set_ylabel('Adapt Time')
    #Peak_Ax.set_ylabel('Adapt Fraction')
    
    Adapt_Ax.legend()
    plt.show(block=False)

    
    #Make and plot the Bounded Boxes
    Target_1=([1]*Osc_Rate+[0.331946]*Osc_Rate)*(Gen_Limit/(2*Osc_Rate))
    Target_2=Target_1[Osc_Rate:]+Target_1[:Osc_Rate]
    #Avg_Ax.plot(Fitness_Packed[:,0],Target_1,ls='--',lw=0.75,c='lawngreen',alpha=1)
    #Avg_Ax.plot(Fitness_Packed[:,0],Target_2,ls='--',lw=0.75,c='aqua',alpha=1)

    d = .015 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=Avg_Ax_Low.transAxes, color='k', clip_on=False,zorder=999)
    d_F=3
    Avg_Ax_Low.plot((1-d/d_F,1+d/d_F), (-d,+d), **kwargs)
    Avg_Ax_Low.plot((1-d/d_F,1+d/d_F),(1-d,1+d), **kwargs)
    
    kwargs.update(transform=Avg_Ax_High.transAxes)  # switch to the bottom axes
    Avg_Ax_High.plot((-d/d_F,+d/d_F), (1-d,1+d), **kwargs)
    Avg_Ax_High.plot((-d/d_F,+d/d_F), (-d,+d), **kwargs)
    
    kwargs.update(transform=Max_Ax_Low.transAxes)  # switch to the bottom axes
    Max_Ax_Low.plot((1-d/d_F,1+d/d_F), (-d,+d), **kwargs)
    Max_Ax_Low.plot((1-d/d_F,1+d/d_F),(1-d,1+d), **kwargs)
    
    kwargs.update(transform=Max_Ax_High.transAxes)  # switch to the bottom axes
    Max_Ax_High.plot((-d/d_F,+d/d_F), (1-d,1+d), **kwargs)
    Max_Ax_High.plot((-d/d_F,+d/d_F), (-d,+d), **kwargs)

    

    Max_Ax_Low.axhline(.5,c='k',lw=1.5,ls='--')
    Max_Ax_High.axhline(.5,c='k',lw=1.5,ls='--')

    for jj in xrange(1,51):
        Avg_Ax_Low.plot([jj*75]*2,[0,1],c='gray',lw=1.5,ls='--',alpha=0.6)
        Avg_Ax_High.plot([22500+jj*75]*2,[0,1],c='gray',lw=1.5,ls='--',alpha=0.6)
        Max_Ax_Low.plot([jj*75]*2,[0.0000001,1],c='gray',lw=1.5,ls='--',alpha=0.6)
        Max_Ax_High.plot([22500+jj*75]*2,[0.0000001,1],c='gray',lw=1.5,ls='--',alpha=0.6)


    #fig.suptitle('18 Colours, 4 Tiles with Regulation, N=1000')
    
    fig.text(0.525, 0.04, 'Generation', ha='center', va='center')
    fig.text(0.02, 0.75, r'$\langle \mathrm{f} \rangle$', ha='center', va='center')
    fig.text(0.02, 0.30, r'$n_{\mathrm{f^{*}}}$', ha='center', va='center')
    #fig.text(0.94, 0.30, r'$\langle \rho_{\mathrm{F^{0}}} \rangle$', ha='center', va='center',)
    
    #Avg_Ax.set_ylabel(r'$\bar{\mathrm{F}}$')
    #Max_Ax.set_xlabel('Generations')
    #Max_Ax.set_ylabel('Maximally Fit Fraction')

    Avg_Ax_High.yaxis.tick_right()
    Max_Ax_High.yaxis.tick_right()
    plt.setp(Avg_Ax_Low.get_xticklabels(), visible=False)
    plt.setp(Avg_Ax_High.get_xticklabels(), visible=False)
    plt.setp(Max_Ax_High.get_yticklabels(), visible=False)
    plt.setp(Avg_Ax_High.get_yticklabels(), visible=False)
    plt.setp(Max_Ax_High.get_yticklines(),visible=False)
    Max_Ax_High.yaxis.set_ticks_position('none')
    Avg_Ax_High.yaxis.set_ticks_position('none') 

    
    #Max_Ax_High.yaxis.set_tick_params(size=0)
    #plt.setp(ax2.get_yticklines(),visible=False)
    Max_Ax_Low.set_yscale('log',nonposy='mask')
    Max_Ax_High.set_yscale('log',nonposy='mask')
    
    plt.figure(fig.number)
    #plt.figlegend(line_caps, [r'N$\mu$L={}'.format(lab) for lab in labels], loc = 'upper center', ncol=4, labelspacing=0. )
    #Avg_Ax_Low.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=6, mode="expand", borderaxespad=0.)
    plt.show(block=False)
    #return Averaged_Data

def Adaptation_Times(Averaged_DataD,Osc_Rate,Gen_Limit):
    Adaptations={}
    Peaks={}

    fig, Adapt_Ax = plt.subplots()
    #Peak_Ax=Adapt_Ax.twinx()
    for Mu in Averaged_DataD.keys():
        print "on mu: ",Mu
        Adaptations[Mu]=[]
        Peaks[Mu]=[]
        temp_cases=[0]+[n*Osc_Rate for n in xrange(1,Gen_Limit/Osc_Rate)]
        for oscillation_L,oscillation_H in zip(temp_cases,temp_cases[1:]):
            Oscillation_Slice=Averaged_DataD[Mu][:,2][oscillation_L:oscillation_H]
            
            index=np.argmax(Oscillation_Slice>0.5)
            if index==0 and Oscillation_Slice[index]<0.5:
                Adaptations[Mu].append(Osc_Rate)
                continue
            else:
                Adaptations[Mu].append(index)
                Peaks[Mu].append(max(Oscillation_Slice))
                
        Adapt_Ax.plot([Osc_Rate*x for x in xrange(len(Adaptations[Mu])) ],Adaptations[Mu],label=Mu,lw=0.5)
        #Peak_Ax.plot(Peaks[Mu],'--')

    Adapt_Ax.set_xlabel('Generations (Cycle={})'.format(Osc_Rate))
    Adapt_Ax.set_ylabel('Adapt Time')
    Adapt_Ax.set_ylim([-5,45])
    #Peak_Ax.set_ylabel('Adapt Fraction')
    #Adapt_Ax.legend()
    plt.show(block=False)

def Plot_CC_Size(Mus):
    Use_Seaborn()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.figure()
    colours=['royalblue','firebrick','darkgreen','goldenrod','deeppink','navy','lime']*2
    xl=[]
    for Mu in Mus:
        for run in xrange(10):
            Fitness_Packed=np.array([[float(value) for value in line.rstrip('\n').split()[1::2]] for line in open('/rscratch/asl47/Bulk_Run/Fitness/Mu_{}_Osc_{}_Gen_{}_Run_{}_3.txt'.format(Mu,15,1000,run))])
            print Fitness_Packed.shape
            plt.plot(Fitness_Packed[:,4][~np.isnan(Fitness_Packed[:,4])],ls='-',c=colours[run])
            #print "Run ",run, " average is ", np.nanmean(Fitness_Packed[:,4][~np.isnan(Fitness_Packed[:,4])])
            xl.append(np.nanmean(Fitness_Packed[:,4][~np.isnan(Fitness_Packed[:,4])]))
            plt.scatter([-10],[np.nanmean(Fitness_Packed[:,4][~np.isnan(Fitness_Packed[:,4])])],marker='x',s=80,c=colours[run])
    print np.nanmean(np.array(xl))
    #plt.title(Mus[0])
    plt.show(block=False)
    
                
######################
###HAMMING DISTANCE###
######################
def Hamming_Load():
    distances=np.loadtxt('/scratch/asl47/Hamming_Distance.txt',dtype=np.dtype('i1'))
    #distances=np.loadtxt('/scratch/asl47/Test2.txt',dtype=np.dtype('i1'))
    max_d=np.amax(distances)
    mean_d=np.mean(distances)
    min_d=np.amin(distances)
    std_d=np.std(distances)
    IQ1=np.percentile(distances,25)
    IQ2=np.percentile(distances,75)
    counter_d=dict(zip(*np.unique(distances, return_counts=True)))

    with open('/scratch/asl47/Hamming_Details.txt','w') as f:
        f.write(str(max_d)+'\n')
        f.write(str(mean_d)+'\n')
        f.write(str(min_d)+'\n')
        f.write(str(std_d)+'\n')
        f.write(str(IQ1)+'\n')
        f.write(str(IQ2)+'\n')
        f.write(str(counter_d))
    
def Hamming_Plot():
    Use_Seaborn()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    lines=[line.rstrip('\n') for line in open('/scratch/asl47/Hamming_Details.txt')]
    max_d=float(lines[0])
    mean_d=float(lines[1])
    min_d=float(lines[2])
    std_d=float(lines[3])
    IQ1=float(lines[4])
    IQ2=float(lines[5])
    counter_d=eval(lines[6])
    ax=plt.subplot(1,1,1)
    plt.scatter(counter_d.keys(),[counter_d[key]*1./sum(counter_d.values()) for key in counter_d.keys()],marker='o',s=70,c='darkred',label='Distribution')
    plt.plot(counter_d.keys(),[counter_d[key]*1./sum(counter_d.values()) for key in counter_d.keys()],ls='--',c='darkred',lw=1,alpha=0.75)
    plt.yscale('log')
    #####FIXED POSITION#####
    ymin,ymax=ax.get_ylim()
    print ymin,ymax
    
    ymean=10**(0.1*(3*np.log10(ymin)+np.log10(ymax)))
    (_,caps,_)=plt.errorbar([mean_d],[ymean], xerr=std_d, fmt='',lw=2,c='forestgreen', capsize=10)
    print mean_d,std_d
    for cap in caps:
        cap.set_markeredgewidth(2)
        
    #plt.plot([min_d,max_d],[ymean]*2,c='forestgreen',ls='-',lw=0.5,alpha=0.7)
    plt.scatter([mean_d],[ymean],marker='x',s=70,c='forestgreen')
    #plt.scatter([IQ1,IQ2],[ymean]*2,marker='D',c='dodgerblue',s=70)
    #plt.scatter([min_d,max_d],[ymean]*2,marker='s',c='dodgerblue',s=50)

    #plt.title(r'Hollow Square vs Cross $S_{2,8}$')
    
    plt.xlabel('Hamming Distance')
    plt.ylabel('Normalised Frequency')
    plt.show(block=False)       


####################
###PHENOTYPE HELP###
####################
from collections import defaultdict
def Get_Topology_Phenotype_Map():
    TP_Map=defaultdict(list)
    
    lines=[line.rstrip('\n') for line in open('/rscratch/asl47/Phenotype_Map_C10.txt')]
    for line in lines:
        TP_Map[int(line.split()[-1])].append([int(i) for i in line.split(':')[0][:-1] if i !=' '])

    return TP_Map

def Get_Corrective_Map():
    C_Map={}
    lineF=[line.rstrip('\n') for line in open('/rscratch/asl47/Topologies_Multi_C10.txt')]
    Counts=[int(ln.split()[-1]) for ln in lineF]
    linePT=[line.rstrip('\n') for line in open('/rscratch/asl47/Topologies_RAW_C10_V2.txt')]
    lineFP=[[int(ite) for ite in ln.split()] for ln in linePT]
    for i,line in enumerate(lineFP):
        C_Map[tuple(line)]=1.*Counts[i]/DR.Find_Overall_Degeneracy([line[:4],line[4:]],10)
    return C_Map
    
####################
###PHENOTYPE PLOT###
####################
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.colors as color
import numpy.ma as ma
from matplotlib.ticker import MultipleLocator
def Plot_Phenotypes_Colour(mapped=True):
    Use_Seaborn()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    minorLocator = MultipleLocator(5)
    
    TP_Map=Get_Topology_Phenotype_Map()
    C_Map=Get_Corrective_Map()
    Manual_Ordering={-1:-1,10:0,9:1,8:2,7:3,0:6,4:4,2:5,1:7,5:8,6:9,3:10}

    sy=np.linspace(-1,11,13,dtype=np.int64)
    sx=np.linspace(2,20,10,dtype=np.int64)
    sx=np.append(sx,np.logspace(np.log10(24),np.log10(502),15,dtype=np.int64))
    y, x = np.meshgrid(sx,sy)#np.mgrid[slice(2, 502, 2),slice(-1, 12 , 1)]
    y=y.T
    x=x.T
    
    z=np.zeros([y.shape[0], x.shape[1]])
    for c,colour in enumerate(y[:,0]):
        #print "C:",colour
        for phenotype in x[0,:]:
            for topology in TP_Map[phenotype]:
                if max(topology)+2>colour:
                    continue
                else:
                    z[c,Manual_Ordering[phenotype]+1]+=C_Map[tuple(topology)]*DR.Find_Overall_Degeneracy([topology[:4],topology[4:]],colour)

    row_sums = z.sum(axis=1, keepdims=True)
    new_matrix = z / row_sums
    zm=ma.masked_where(new_matrix == 0, new_matrix)
    zz=ma.masked_where(new_matrix != 0, new_matrix)


    
    labels=['UND','Small Square','Hollow Square','Small Cross','Thick Cross','T Shape','L Catherine Wheel','R Catherine Wheel','L Shape','I Shape','Domino','Monomino']
    labels_Ordered=['']*len(labels)
    for k,v in Manual_Ordering.iteritems():
        labels_Ordered[v+1]=labels[k+1]
    if not mapped:
        cols=['black','firebrick','dodgerblue','royalblue','steelblue','skyblue','cornflowerblue','goldenrod','lime','seagreen','darkgreen','darkolivegreen']
        styles=['-','-','-.','-.','-.','-.','-.','-.','-','-','-','-.']
        markers=['X','D','o','s','^','v','H','>','o','s','D','<']
        fig=plt.figure()
        for phenotype in x[0,:-1]:
            #print phenotype
            plt.plot(y[:,0],zm[:,phenotype+1],label=labels_Ordered[phenotype+1],c=cols[phenotype+1],ls=styles[phenotype+1],marker=markers[phenotype+1],markersize=8)
            slopeR, interceptR, r_valueR, p_valueR, std_errR = linregress(np.log10(y[-5:,0]),np.log10(zm[-5:,phenotype+1]))
            print labels_Ordered[phenotype+1],"R",round(slopeR)#, interceptR, r_valueR**2, std_errR

            
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'C')
        plt.ylabel(r'$\rho_{\mathrm{p}}$',fontsize=22)
        plt.legend(ncol=3)#fontsize='large')
        #sns.despine()
        #plt.title(r'Two Tile Asymmetric Interactions')
        plt.show(block=False)
        

        
        return
    cmap = plt.get_cmap('cubehelix_r')

    

    
    
    fig, (ax, ax2) = plt.subplots(1, 2, sharey=True,gridspec_kw = {'width_ratios':[3.5, 1]})
    im = ax.pcolormesh(x, y, zm, cmap=cmap, norm=color.LogNorm(vmin=zm[:,:11].min(), vmax=zm[:,:11].max()))#,edgecolors='k')
    blackMap=color.ListedColormap(['k'])
    ax.pcolor(x, y, zz,alpha=0,hatch='//')
    #ax.set_yscale('log')
    
    ax.set_xticks(x[0,:]+0.5)
    ax.set_xticklabels(labels_Ordered, rotation=45)
    
    y_ticks=np.linspace(2,20,10)+1
    cb_range=np.logspace(round(np.log10(zm[:,:11].min())),0,abs(round(np.log10(zm[:,:11].min())))+1)

    cb = fig.colorbar(im, ax = ax,ticks=cb_range)
    cb.ax.set_yticklabels([r'$10^{{{}}}$'.format(int(np.log10(cb_val))) for cb_val in cb_range])

    cb.ax.get_yaxis().labelpad = 15
    cb.ax.set_ylabel(r'Normalised $\rho_{\mathrm{p}}$', rotation=90)
    ax.set_ylabel(r'Number of Colours')
    
    for x_l in xrange(-1,12,1):
        ax.axvline(x=x_l,color='snow',linestyle='-',linewidth=1)
    ax.set_xlim([-1,11])

    
    
    #plt.show(block=False)
    cmap2 = plt.get_cmap('magma_r')
    #fig2, ax2 = plt.subplots()

    im2 = ax2.pcolormesh(x[:,:4], y[:,:4], zm[:,:4], cmap=cmap2, norm=color.LogNorm(vmin=zm[:,:3].min(), vmax=zm[:,:3].max()))#,edgecolors='k')
    ax2.set_xticks(x[0,:3]+0.5)
    ax2.set_xticklabels(labels_Ordered[:3], rotation=45)
    for x_l in xrange(-1,3,1):
        ax2.axvline(x=x_l,color='snow',linestyle='-',linewidth=1)
    fig.colorbar(im2, ax=ax2)

    plt.show(block=False)
    return zm

    
def Plot_Phenotypes_Colour_2():
    Use_Seaborn()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    minorLocator = MultipleLocator(5)
    
    TP_Map=Get_Topology_Phenotype_Map()
    C_Map=Get_Corrective_Map()
    Manual_Ordering={-1:-1,10:0,9:1,8:2,7:3,0:4,4:5,2:6,1:7,5:8,6:9,3:10}    
    
        
from scipy.stats import linregress
from numpy import log10
def Find_Log_Slope(fileI):
    lines=[line.rstrip('\n') for line in open("Output/{}.txt".format(fileI))]
    pp=[filter(None,re.split('[a-z]|:|_|[A-Z]',line)) for line in lines]
    pp_float=np.array([[float(i) for i in  p] for p in pp])

    round_Factor=0.75
    index_min_D=np.argmin(pp_float[:,1])
    index_start_D=0
    while(pp_float[:,1][index_start_D]>=20000):
        index_start_D+=1
    print index_start_D
    adjusted_min_D=int(round(round_Factor*index_min_D))
    slope_D,intercept_D, R_D,P_D,err_D=linregress(log10(pp_float[:,0][index_start_D:adjusted_min_D]),log10(pp_float[:,1][index_start_D:adjusted_min_D]))
    

    index_min_A=np.argmin(pp_float[:,4])
    index_start_A=0
    while(pp_float[:,4][index_start_A]>=20000):
        index_start_A+=1
    print index_start_A
    adjusted_min_A=int(round(round_Factor*index_min_A))
    slope_A,intercept_A, R_A,P_A,err_A=linregress(log10(pp_float[:,0][index_start_A:adjusted_min_A]),log10(pp_float[:,4][index_start_A:adjusted_min_A]))

    
    print "Discovery slope: ",slope_D
    print "Adaptation slope: ",slope_A

    
def makeFileName(T,C,N,R,M,Targ,Thresh,Reg,Run):
    return 'Evolution_T{}_C{}_N{}_R{}_M{}_Targ{}_Thresh{}_Reg_{}_Run{}'.format(T,C,N,R,M,Targ,Thresh,Reg,Run)
    if Reg:
        if k_on:
            return "Evolution_T{}_C{}_N{}_R{}_M{}_Regulated_Targ{}_Thresh{}_Run{}".format(T,C,N,R,M,Targ,K,Run)
        else:
            return "Evolution_T{}_C{}_N{}_R{}_M{}_Regulated_Targ{}_Run{}".format(T,C,N,R,M,Targ,Run)
    else:
        if k_on:
            return "Evolution_T{}_C{}_N{}_R{}_M{}_Targ{}_Thresh{}_Run{}".format(T,C,N,R,M,Targ,K,Run)
        else:
            return "Evolution_T{}_C{}_N{}_R{}_M{}_Targ{}_Run{}".format(T,C,N,R,M,Targ,Run)


def Plot_Fit_Funcs():
    plt.plot(np.linspace(0,1,100),np.linspace(0,1,100),'g',label='linear')
    plt.plot(np.linspace(0,1,100),[x**3 for x in np.linspace(0,1,100)],'b',label='cubic')
    
    plt.plot(np.linspace(0,1,100),[(np.exp(x)-1.)/(np.exp(1)-1) for x in np.linspace(0,1,100)],c='peru',label='exp')
    plt.plot(np.linspace(0,1,100),[(np.exp(12.5*x)-1.)/(np.exp(12.5)-1) for x in np.linspace(0,1,100)],'violet',label='exp 12.5')
    plt.plot(np.linspace(0,1,100),[(np.exp(25*x)-1.)/(np.exp(25)-1) for x in np.linspace(0,1,100)],'darkmagenta',label='exp 25')
    plt.plot(np.linspace(0,1,21),[((np.exp(12.5*x)-1.)/(np.exp(12.5)-1)+1/12.5*x)/(1+1/12.5) for x in np.linspace(0,1,21)],'salmon',label=r'$\xi$',marker='o')
    plt.plot(np.linspace(0,1,51),[((np.exp(50*x)-1.)/(np.exp(50)-1)+1/4.*x)/(1+1/4.) for x in np.linspace(0,1,51)],'lime',label=r'$\eta$',marker='o')
    plt.plot(np.linspace(0,1,51),[(1/2.*x)/(1+1/2.) for x in np.linspace(0,1,51)],'orchid',label=r'$\zeta$',marker='x')
    plt.xscale('log')
    plt.yscale('log')
   
    #plt.plot(np.linspace(0,1,26),[((np.exp(6.25*x)-1.)/(np.exp(6.25)-1)+1/6.25*x)/(1+1/6.25) for x in np.linspace(0,1,26)],'cadetblue',label=r'$\xi^{\prime}$',marker='D')
    #plt.plot(np.linspace(0,1,100),[(np.exp(x)-1)/(np.exp(1)-1)+0.25 for x in np.linspace(0,1,100)],c='peru',ls=':')
    #plt.plot(np.linspace(0,1,100),[(np.exp(10*x)-1)/(np.exp(10)-1)+0.25 for x in np.linspace(0,1,100)],'orange',ls=':')
    #plt.plot(np.linspace(0,1,100),[(np.exp(12.5*x)-1)/(np.exp(12.5)-1)+0.25 for x in np.linspace(0,1,100)],'violet',ls=':')
    #plt.plot(np.linspace(0,1,100),[(np.exp(25*x)-1)/(np.exp(25)-1)+0.25 for x in np.linspace(0,1,100)],'darkmagenta',ls=':')
    #plt.xscale('log')
    #plt.yscale('log')
    plt.legend()

    plt.show(block=False)
def Make_Colour_Min_Plot(files,xes):
    minsA=[]
    minsA_err_Upper=[]
    minsA_err_Lower=[]
    minsD=[]
    minsD_err_Upper=[]
    minsD_err_Lower=[]
    for fil in files:
        lines=[line.rstrip('\n') for line in open("Evolution_Runs/{}.txt".format(fil))]
        pp=[filter(None,re.split('[a-z]|:|_|[A-Z]',line)) for line in lines]
        pp_float=np.array([[float(i) for i in  p] for p in pp])
        minsA.append(min(pp_float[:,4]))
        minsA_Ind=np.argmin(pp_float[:,4])
        minsA_err_Upper.append(pp_float[:,6][minsA_Ind]-minsA[-1])
        minsA_err_Lower.append(minsA[-1]-pp_float[:,5][minsA_Ind])
        print minsA[-1],minsA_err_Upper[-1],minsA_err_Lower[-1]
        minsD.append(min(pp_float[:,1]))
        minsD_Ind=np.argmin(pp_float[:,1])
        minsD_err_Upper.append(pp_float[:,3][minsD_Ind]-minsD[-1])
        minsD_err_Lower.append(minsD[-1]-pp_float[:,2][minsD_Ind])
        print minsD[-1],minsD_err_Upper[-1],minsD_err_Lower[-1]
    Use_Seaborn()
    plt.errorbar(xes,minsA,yerr=[minsA_err_Lower,minsA_err_Upper],marker='',ls='',c='cornflowerblue')
    plt.errorbar(xes,minsD,yerr=[minsD_err_Lower,minsD_err_Upper],marker='',ls='',c='orangered')
    plt.plot(xes,minsA,marker='o',ls='',lw=2.5,c='cornflowerblue',label='Adaptation')
    plt.plot(xes,minsD,marker='o',ls='',lw=2.5,c='orangered',label='Discovery')
    slope_A,intercept_A, R_A,P_A,err_A=linregress(xes[2:],minsA[2:])
    slope_D,intercept_D, R_D,P_D,err_D=linregress(xes[2:],minsD[2:])
    plt.plot(xes[2:],[slope_A*x+intercept_A for x in xes[2:]],c='cornflowerblue',ls='--',label='Linear Fit, m={:.2f}'.format(slope_A))
    plt.plot(xes[2:],[slope_D*x+intercept_D for x in xes[2:]],c='orangered',ls='--',label='linear Fit, m={:.2f}'.format(slope_D))
    plt.ylabel('Generations')
    plt.xlabel('Colours')
    plt.title('2 Tile, 12 Fitness Target, Unregulated')
    plt.legend(loc='upper left')
    sns.despine()
    plt.show()





######################
###RANDOM WALK PLOT###
######################
from fractions import Fraction
from cycler import cycler
def Plot_Random_Walks(Input_Colours):
    Use_Seaborn()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    ERROR_MAGNIFICATION_FACTOR=1
    Tile_Kit=[[0,0,0,1],[2,2,3,4]]
    Mutation_Steps_Avg=defaultdict(list)
    Mutation_Steps_dev=defaultdict(list)
    
    Logspaced_Fractions=[C**(len(Tile_Kit)*4.)/DR.Find_Overall_Degeneracy(Tile_Kit,C) for C in Input_Colours]
    
    for Col in Input_Colours:
        lines_Raw=[line.rstrip('\n') for line in open('/scratch/asl47/Data_Runs/Random_Walks/Random_Walk_L8_C_{}_W_1000000_Z_0.txt'.format(Col))]
        Filtered_lines=[filter(None,re.split('[a-z]|:|_|[A-Z]',line)) for line in lines_Raw]
        for run in Filtered_lines:
            Mutation_Steps_Avg[float(run[0])].append(float(run[1]))
            Mutation_Steps_dev[float(run[0])].append(float(run[2])*ERROR_MAGNIFICATION_FACTOR)


    fig, Walk_Ax = plt.subplots()
    Walk_Ax.set_prop_cycle(cycler('color',['midnightblue','midnightblue','slategrey','slategrey','darkgreen','darkgreen','firebrick','firebrick'])+cycler('marker',['^','^','s','s','D','D','o','o']))
    
    
    for Mu in Mutation_Steps_Avg.keys():
        Walk_Ax.errorbar(Input_Colours,Mutation_Steps_Avg[Mu],yerr=Mutation_Steps_dev[Mu],markeredgecolor='none',ls='',lw=2)
        Walk_Ax.plot(Input_Colours,Mutation_Steps_Avg[Mu],ls=':',lw=1.5,label=r'$\mu={:}$'.format(Fraction(Mu)),markeredgewidth=2,markersize=7,markerfacecolor='white')
    Walk_Ax.plot(Input_Colours,Logspaced_Fractions,marker='P',markeredgewidth=1.2,markeredgecolor='darkorchid',markersize=8,markerfacecolor='white',c='darkorchid',ls='--',lw=1.5,zorder=10,label=r'$1/\rho_{p}$')
    Walk_Ax.set_yscale('log')
    #Frac_Ax.set_yscale('log')
    Walk_Ax.set_xscale('log')
    Walk_Ax.legend()
    Walk_Ax.set_xlabel(r'Number of Colors')
    Walk_Ax.set_ylabel(r'$\langle \tau_{D} \rangle$')
    Walk_Ax.set_title(r'Random Mutating Walk on Two Tiles')
    sns.despine(fig)
    
    plt.show(block=False)

    fig2, Walk_Ax2 = plt.subplots()
    for Mu in Mutation_Steps_Avg.keys():
        Walk_Ax2.plot(Input_Colours,[a*1./b for (a,b) in zip(Mutation_Steps_Avg[Mu],Logspaced_Fractions)],ls=':',lw=1.5,label=r'$\mu={:}$'.format(Fraction(Mu)),markeredgewidth=2,markersize=7,markerfacecolor='white')
    Walk_Ax2.set_xscale('log')
    Walk_Ax2.set_yscale('log')
    plt.show(block=False)

###################
###MUTATION PLOT###
###################

def Mu_L_Probability():
    Use_Seaborn()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    N_f=80
    MuLs=np.logspace(-2,np.log10(N_f),500)
    mutations=defaultdict(list)

    for mu_L in MuLs:
        binomD=binom(N_f,mu_L/(1.*N_f))
        for i in xrange(N_f+1):
            mutations[i].append(binomD.pmf(i))

    ax.set_prop_cycle(cycler('linestyle',[':','--','-.','-',':','--','-.','-',':'])+cycler('color',['midnightblue','lightslategrey','darkgreen','darkred','orange','lime','pink','m','olive']))
    #Colors=
    min_higher=3
    for i in xrange(min_higher):     
        ax.plot(MuLs,mutations[i],lw=2.5,label='{} Mutation{plural}'.format(i,plural='s' if i!=1 else ''))
        higher_order=np.array(mutations[min_higher])
    for i in xrange(min_higher+1,N_f+1):
        higher_order+=mutations[i]
    ax.plot(MuLs,higher_order,lw=2.5,label='Higher Order Mutations')
    
    ax.set_xscale('log')
    ax.set_yscale('log',nonposy='mask')
    sns.despine(fig)
    ax.legend(loc='center left')
    fig.tight_layout()
    ax.set_ylabel(r'Probability')
    ax.set_xlabel(r'$\langle \mu L\rangle$')
    ax.set_ylim([0,1.01])
    ax.set_title(r'Expected Mutations on a Two Tile System')
    plt.show(block=0)
    
#####################
###DEGENERACY PLOT###
#####################
import matplotlib.ticker as plticker 
from numpy import log10
from numpy import logspace
def Plot_Degeneracy_Rate_Overlayed(Tile_Kit=[[0,0,0,1],[2,2,3,4]],C_Min=6):
    
    Use_Seaborn()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig, Frac_Ax = plt.subplots()
    Search_Ax=Frac_Ax.twinx()
    Frac_Ax.set_xscale('log')
    Frac_Ax.set_yscale('log')
    Search_Ax.set_yscale('log')
    Search_Ax.set_xscale('log')
    Search_Ax.set_ylabel('Size of Sequence Space')
    Frac_Ax.set_xlabel(r'Number of Colours')

    #Cs=[6,8,10,12,14,16,20,24,28,32,40,46,50,60,75]
    Cs=[i for i in xrange(C_Min,100,2)]
    Colours=[i for i in xrange(C_Min,24,2)]+[int(val) for val in logspace(log10(22),2,10)][1:]
    Fit_Cut=len(Colours)
    Colours+=[int(val) for val in logspace(2.13,3.55,10)]
    P_fit=[int(val) for val in logspace(2.13,3.55,10)]

    
    Logspaced_Degeneracies=[DR.Find_Overall_Degeneracy(Tile_Kit,C) for C in Colours]
    Logspaced_Fractions=[Degen*1./C**(len(Tile_Kit)*4) for (Degen,C) in zip(Logspaced_Degeneracies,Colours)]
    
    T_line,=Search_Ax.plot(Colours,Logspaced_Degeneracies,ls='',marker='D',markeredgewidth=2,markeredgecolor='cornflowerblue',markersize=6,markerfacecolor='white',zorder=10)
    S_line,=Search_Ax.plot(Colours,[C**(len(Tile_Kit)*4) for C in Colours],ls='',marker='s',markeredgewidth=2,markeredgecolor='seagreen',markersize=6,markerfacecolor='white',zorder=10)
    F_line,=Frac_Ax.plot(Colours,Logspaced_Fractions,ls='',marker='P',markeredgewidth=1.2,markeredgecolor='darkorchid',markersize=7,markerfacecolor='white',zorder=10)

    
    

    #Frac_Ax.set_ylim([0,0.005])
    #Search_Ax.set_xlim([5,75])
    Frac_Ax.set_ylabel(r'Fractional Occupation')
    Frac_Ax.tick_params('y',colors='darkorchid')
    Frac_Ax.spines['left'].set_color('darkorchid')
    Search_Ax.spines['left'].set_color('darkorchid')
    for label in Search_Ax.yaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    for label in Frac_Ax.yaxis.get_ticklabels()[::2]:
        label.set_visible(False)

    
    #slope_D,intercept_D, R_D,P_D,err_D=linregress(Cs_Fit,[Find_Overall_Degeneracy(Tile_Kit,C)*1./C**(len(Tile_Kit)*4) for C in Cs_Fit])

    P_Fit_Logged=log10(P_fit)
    P_slope,P_intercept, P_R,P_P,P_err=(1,2,3,4,5)#linregress(P_Fit_Logged,log10(Logspaced_Degeneracies[Fit_Cut:]))
    
    P_slope_F,P_intercept_F, P_R_F,P_P_F,P_err_F=linregress(P_Fit_Logged,log10(Logspaced_Fractions[Fit_Cut:]))

    
    print "Multiplicty ~ O(N^{:.2f})".format(P_slope)
    print "Fraction ~ O(N^{:.2f})".format(P_slope_F)
    #P_slope_F,P_intercept_F, P_R_F,P_P_F,P_err_F
    
    T_Fit_line,=Search_Ax.plot(Colours,[10**P_intercept*C**P_slope for C in Colours],ls='--',lw=1.5,c='cornflowerblue',zorder=1,alpha=0.7)
    F_Fit_line,=Frac_Ax.plot(Colours,[10**P_intercept_F*C**P_slope_F for C in Colours],ls='--',lw=1.5,c='darkorchid',zorder=1,alpha=0.7)
    
    sns.despine(fig=fig,top=True,right=False)
    fig.legend([S_line,T_line,F_line],[r'$\Phi$',r'$\phi_{\mathrm{p}}$',r'$\rho_{\mathrm{p}}$'],ncol=1,loc='upper center',fontsize='large')
    fig.tight_layout()
    #fig.suptitle('{} {} Topology'.format(*Tile_Kit),y=0.1)
    plt.show(block=0)




    
from collections import Counter
def GP_Robustness(C):
    Use_Seaborn()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    lines_Raw=[line.rstrip('\n') for line in open('/rscratch/asl47/Topologies_Multi_C8.txt')]
    Counts=[int(ln.split()[-1]) for ln in lines_Raw]
    lines_Raw2=[line.rstrip('\n') for line in open('/rscratch/asl47/Topologies_Multi_C10.txt')]
    Counts2=[int(ln.split()[-1]) for ln in lines_Raw2]

    lines=[line.rstrip('\n') for line in open('/rscratch/asl47/Phenotype_Map_C10.txt')]
    colour_mask=[0 if int(line.split()[-1])==-1 else 1 for line in lines]

    BD_Phen=[i for i, x in enumerate(colour_mask) if x == 1]
    
    
    counted_GPs=Counter(Counts)
    counted_GPs2=Counter(Counts2)
    ba_x,ba_y = log_binning(counted_GPs,10)
    ba_x2,ba_y2 = log_binning(counted_GPs2,10)
    D_min=min(Counts)
    D_max=max(Counts)

    #plt.plot(ba_x,ba_y,ls='--',c='orangered',marker='D',label='C=8')
    #plt.plot(counted_GPs.keys(),[cnt*1./sum(counted_GPs.values()) for cnt in counted_GPs.values()],'ro',alpha=0.3)
    plt.plot([key*1./sum(Counts2) for key in counted_GPs2.keys()],[cnt*1./sum(counted_GPs2.values()) for cnt in counted_GPs2.values()],'o',c='firebrick',alpha=0.8,label='Topologies')
    
    BDs=[[i for i, x in enumerate(counted_GPs2.keys()) if x == Counts2[phen]] for phen in BD_Phen]
    BD_Flat=[item for sublist in BDs for item in sublist]
    print BD_Flat
    print [key*1./sum(Counts2) for key in [counted_GPs2.keys()[phen] for phen in BD_Flat]]
    print [cnt*1./sum(counted_GPs2.values()) for cnt in [counted_GPs2.values()[phen] for phen in BD_Flat]]
    
    plt.plot([key*1./sum(Counts2) for key in [counted_GPs2.keys()[phen] for phen in BD_Flat]],[cnt*1./sum(counted_GPs2.values()) for cnt in [counted_GPs2.values()[phen] for phen in BD_Flat]],'o',c='slateblue',alpha=0.8,label='B/D Topologies')
    
    plt.plot([b2*1./sum(Counts2) for b2 in ba_x2],ba_y2,ls='--',lw=2,c='firebrick',alpha=0.8,marker='',label='Fit')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='upper left')
    plt.xlabel(r'$\rho_{\mathrm{T}}$')
    plt.ylabel('Normalised Frequency')
    plt.title(r'Two Tile Asymmetric Topologies')
    
    sns.despine()
    plt.show(block=False)
    
    return counted_GPs2

def drop_zeros(a_list):
    return [i for i in a_list if i>0]

def log_binning(counter_dict,bin_count=35):

    max_x = log10(max(counter_dict.keys()))
    max_y = log10(max(counter_dict.values()))
    max_base = max([max_x,max_y])

    min_x = log10(min(drop_zeros(counter_dict.keys())))

    bins = np.logspace(min_x,max_base,num=bin_count)

    # Based off of: http://stackoverflow.com/questions/6163334/binning-data-in-python-with-scipy-numpy
    bin_means_y = (np.histogram(counter_dict.keys(),bins,weights=counter_dict.values())[0] / np.histogram(counter_dict.keys(),bins)[0])*1./sum(counter_dict.values())
    bin_means_x = (np.histogram(counter_dict.keys(),bins,weights=counter_dict.keys())[0] / np.histogram(counter_dict.keys(),bins)[0])

    return bin_means_x,bin_means_y




from math import atan2,cos,sin,degrees
import numpy as np

#Label line with line2D label data
def labelLine(line,x,label=None,align=True,**kwargs):

    ax = line.get_axes()
    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if (x < xdata[0]) or (x > xdata[-1]):
        print('x label location is outside data range!')
        return

    #Find corresponding y co-ordinate and angle of the
    ip = 1
    for i in range(len(xdata)):
        if x < xdata[i]:
            ip = i
            break

    y = ydata[ip-1] + (ydata[ip]-ydata[ip-1])*(x-xdata[ip-1])/(xdata[ip]-xdata[ip-1])

    if not label:
        label = line.get_label()

    if align:
        #Compute the slope
        dx = xdata[ip] - xdata[ip-1]
        dy = ydata[ip] - ydata[ip-1]
        ang = degrees(atan2(dy,dx))

        #Transform to screen co-ordinates
        pt = np.array([x,y]).reshape((1,2))
        trans_angle = ax.transData.transform_angles(np.array((ang,)),pt)[0]

    else:
        trans_angle = 0

    #Set a bunch of keyword arguments
    if 'color' not in kwargs:
        kwargs['color'] = line.get_color()

    if ('horizontalalignment' not in kwargs) and ('ha' not in kwargs):
        kwargs['ha'] = 'center'

    if ('verticalalignment' not in kwargs) and ('va' not in kwargs):
        kwargs['va'] = 'center'

    if 'backgroundcolor' not in kwargs:
        kwargs['backgroundcolor'] = ax.get_axis_bgcolor()

    if 'clip_on' not in kwargs:
        kwargs['clip_on'] = True

    if 'zorder' not in kwargs:
        kwargs['zorder'] = 2.5
    #if 'offset' not in kwargs:
    #    offset_F=15
    offset_F=15
    xprime=x-offset_F*cos(atan2(dy,dx))
    yprime=y+offset_F*sin(atan2(dy,dx))
    
    ax.text(xprime,yprime,label,rotation=trans_angle,**kwargs)

def labelLines(lines,align=True,xvals=None,**kwargs):

    ax = lines[0].get_axes()
    labLines = []
    labels = []

    #Take only the lines which have labels other than the default ones
    for line in lines:
        label = line.get_label()
        if "_line" not in label:
            labLines.append(line)
            labels.append(label)

    if xvals is None:
        xmin,xmax = ax.get_xlim()
        xvals = np.linspace(xmin,xmax,len(labLines)+2)[1:-1]

    for line,x,label in zip(labLines,xvals,labels):
        labelLine(line,x,label,align,**kwargs)
