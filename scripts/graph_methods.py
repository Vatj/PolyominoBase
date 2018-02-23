import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict,deque,Counter
import itertools
import seaborn as sns
import numpy as np

def Transform_Graph_From_List(tile_kit):
    graph_kit=nx.MultiDiGraph()
    graph_kit.add_nodes_from(xrange(len(tile_kit)))

    ## Add edges for internal structure in clockwise orientation 
    for internal_edge in xrange(len(tile_kit)/4):
        graph_kit.add_edge(internal_edge*4+0,internal_edge*4+1)#,color='k')
        graph_kit.add_edge(internal_edge*4+1,internal_edge*4+2)#,color='k')
        graph_kit.add_edge(internal_edge*4+2,internal_edge*4+3)#,color='k')
        graph_kit.add_edge(internal_edge*4+3,internal_edge*4+0)#,color='k')

    ## Add edges for graph connections
    for index,face in enumerate(tile_kit):
        search_offset=0
        if face==0:
            continue
        else:
            while(Interaction_Matrix(face) in tile_kit[index+search_offset:]):
                search_offset+=tile_kit[index+search_offset:].index(Interaction_Matrix(face))
                graph_kit.add_edge(index,index+search_offset)#,color='r')
                graph_kit.add_edge(index+search_offset,index)#,color='r')
                search_offset+=1                               
                
    return graph_kit

def Interaction_Matrix(colour):
    return colour if colour <=0 else  (1-colour%2)*(colour-1)+(colour%2)*(colour+1)


def Draw_Graph(graph,kit):
    layout=nx.spring_layout(graph)
    nx.draw_networkx_labels(graph,layout,{j:str(i) for j,i in enumerate(kit)})
    nx.draw_networkx_nodes(graph,layout)
    nx.draw_networkx_edges(graph,layout)
    plt.show(block=False)


def StripIsomorphisms(file_in):
    tile_kits=[[int(i) for i in line.rstrip('\n').split()] for line in open('/rscratch/asl47/Bulk_Run/Modular/{}.txt'.format(file_in))]
    assembly_graphs=zip(range(len(tile_kits)),[Transform_Graph_From_List(tile_kit) for tile_kit in tile_kits])

    unique_assembly_graphs_indices=[]
    
    while len(assembly_graphs)>1:
        #print len(assembly_graphs)
        unique_assembly_graphs_indices.append(assembly_graphs[0][0])
        assembly_graphs[:]=[assembly_graph for assembly_graph in assembly_graphs[1:] if not nx.is_isomorphic(assembly_graph[1],assembly_graphs[0][1])]

    with open('/rscratch/asl47/Bulk_Run/Modular/{}_stripped.txt'.format(file_in), 'w') as outfile:
        for index in unique_assembly_graphs_indices:
            genotype= ' '.join(map(str,tile_kits[index]))+'\n'
            outfile.write(genotype)

from multiprocessing import Pool
def StripInParallel(runs):
    pool = Pool(processes=4)
    pool.map_async(StripIsomorphisms,["A2_T20_C200_N500_Mu0.003125_O25_K15000_Run{}_Genotype".format(r) for r in runs])
    pool.close()
    

    


def Load_Tiles(fileN):
    Raw_Topologies_Input = [line.rstrip('\n') for line in open('/rscratch/asl47/{}.txt'.format(fileN))]
    Raw_Topologies = [[int(face) for face in line.split()] for line in Raw_Topologies_Input]
    topology_dict=defaultdict(list)
    
    for tile_kit in Raw_Topologies:
        zeroes=Order_By_Zeroes(tile_kit)
        topology_dict[zeroes].append(tile_kit)
        
    return topology_dict

def temp(kit):
    n=0
    newk=[]
    while(n<20):
        if kit[n]==0:
            n+=5
        else:
            for i in xrange(1,5):
                newk.append(kit[n+i])
            n+=5
    return [0 if kk==5 else kk for kk in newk]


def Trim_Topologies():
    TD=Load_Tiles('Randomized_Topologies')
    unique_C=0
    for key in reversed(sorted(TD.keys())):
        del_count=0
        print "on ",key, len(TD[key])
        for i,comp1 in enumerate(TD[key]):
            G1=Transform_Graph_From_List(comp1)
            del_count=0
            for j,comp2 in enumerate(TD[key][i+1:]):
                G2=Transform_Graph_From_List(comp2)
                if nx.is_isomorphic(G1,G2):
                    #print i+1+j-del_count,len(value),i,j,del_count
                    del TD[key][i+1+j-del_count]
                    del_count+=1
            unique_C+=1
    print unique_C

    
        
                
def Order_By_Zeroes(tile_kit):
    z_Map=defaultdict(int)
    for i in xrange(len(tile_kit)/4):
        tile=tile_kit[i*4:i*4+4]
        if 0 in tile:
            z_Map[tile.count(0)]+=1
        else:
            z_Map[0]+=1
                       
    return (z_Map[0],z_Map[1],z_Map[2],z_Map[3],z_Map[4])


    while(len(z_Map)>0):
        new_tile_kit.extend(tile_kit[max(z_Map, key=z_Map.get)*4:max(z_Map, key=z_Map.get)*4+4])
        del z_Map[max(z_Map, key=z_Map.get)]

    return new_tile_kit


def Enumerate_Topology(tile_kit):
    queue_kit=[deque(tile) for tile in tile_kit]

    
    for tile_rot in xrange(4**len(tile_kit)):
        for tile_index in xrange(len(tile_kit)-1,-1,-1):
            if tile_rot%(4**tile_index)==0:
                queue_kit[tile_index].rotate(-1)
                yield itertools.permutations(queue_kit)
                break
            

def cycle_tile(tile):
    for rotation in xrange(4):
        tile.rotate(-1)
        yield tile

def Hamming_Distance(s1,s2):
    """Calculate the Hamming distance between two bit strings"""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def Hamming_Distance_Of_Topologies(T1,T2):
    Flat_T1=[item for sublist in T1 for item in sublist]
    Possibles=[]
    for i in Enumerate_Topology(T2):
        for j in i:
            Possibles.append([item for sublist in j for item in sublist])
            
    return  min([Hamming_Distance(Flat_T1,trial) for trial in Possibles])

def dump(C):
    for i in xrange(C):
        for j in xrange(C):
            for k in xrange(C):
                for l in xrange(C):
                    for m in xrange(C):
                        for n in xrange(C):
                            for o in xrange(C):
                                for p in xrange(C):
                                    yield [i,j,k,l,m,n,o,p]

def dump2():
    G1=Transform_Graph_From_List([0,1,0,3,4,0,0,2])
    #G1=Transform_Graph_From_List([0,0,0,1,2,2,2,2])
    with open('/rscratch/asl47/Loop_Genotypes.txt','w') as f:
        for n,T2 in enumerate(dump(8)):
            if n%50000==0:
                print n
            G2=Transform_Graph_From_List(T2)
            if nx.is_isomorphic(G1,G2):
                x=''
                for i in T2:
                    x+=str(i)+' '
                f.write(x[:-1]+'\n')
                

def Phenotype_Hamming_Distances():
    BP_Genotypes_Input = [line.rstrip('\n') for line in open('/rscratch/asl47/BP_Genotypes.txt')]
    BP_Genotypes = [[int(face) for face in line.split()] for line in BP_Genotypes_Input]
    Loop_Genotypes_Input = [line.rstrip('\n') for line in open('/rscratch/asl47/Loop_Genotypes.txt')]
    Loop_Genotypes = [[int(face) for face in line.split()] for line in Loop_Genotypes_Input]
    print "loops",len(Loop_Genotypes)
    print "BPs",len(BP_Genotypes)
    H_D=[]
    with open('/scratch/asl47/Hamming_Distance.txt','w') as f:
        for n,BP in enumerate(BP_Genotypes):
            if n%100==0:
                print "Now on {}".format(n)
            for LP in Loop_Genotypes:
                f.write(str(Hamming_Distance_Of_Topologies([BP[x:x+4] for x in xrange(0,len(BP),4)],[LP[x:x+4] for x in xrange(0,len(LP),4)]))+'\n')
    print "all finished"
    return True
                        
    Counts=Counter(H_D)
    np_H=np.array(H_D)
    print Counts
    print min(H_D)
    print max(H_D)
    print np.percentile(np_H,25)
    print np.percentile(np_H,75)
    print np.std(np_H)
    ax = sns.violinplot(x=H_D)
    plt.show(block=False)
         
     
if __name__ == "__main__":
    #Trim_Topologies()
    #dump2()
    Phenotype_Hamming_Distances()


############
##PLOTTING##
############

def Use_Seaborn():
    sns.set_context("paper",font_scale=2.2)
    sns.set_style("white",rc={"xtick.major.size": 8, "ytick.major.size": 8,"xtick.minor.size":5, "ytick.minor.size": 5,"axes.linewidth": 2,"axes.edgecolor":"darkgray","font.size":8,"axes.titlesize":8,"axes.labelsize":5})
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')


    

def Plot_Method_Testing(Ts,R):
    Use_Seaborn()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    #timing_list=[defaultdict(list) for i in xrange(3)]
    #accuracy_list=[defaultdict(list) for i in xrange(3)]
    
    timing_list=[{K:np.zeros((Ts,R)) for K in [0,5,10,20]} for i in xrange(3)]
    accuracy_list=[{K:np.zeros((Ts,R)) for K in [0,5,10,20]} for i in xrange(3)]
    
    
    
    for T in xrange(1,Ts+1):
        lines=[]
        lines2=[]
        if T<=10:
            lines=[line.rstrip('\n') for line in open('/rscratch/asl47/Method_Analysis_V2_T{}.txt'.format(T))]
            lines2=[line.rstrip('\n') for line in open('/rscratch/asl47/Method_Analysis_V3_T{}.txt'.format(T))]
        if T>10:
            lines=[line.rstrip('\n') for line in open('/rscratch/asl47/Method_Analysis_V3_T{}.txt'.format(T))]
            
        load_time=0
        condition=0
        Run_Total=0
        Run_Quantity=0
        RUN=0


        for line in lines+lines2:
            #print line
            if 'condition' in line:
                Run_Quantity=int(line.split()[1])
                condition=int(line.split()[4])
                
            if 'STARTING RUN' in line:
                RUN=int(line.split()[2])
            if RUN>=R:
                continue
                
            if 'Load time' in line:
                load_time=float(line.split()[-1])
                
            if 'B/D' in line:
                Run_Total=int(line.split()[1])+int(line.split()[-1])
            
            if 'Runtime' in line:
                timing_list[condition][int(line.split()[3])][T-1][RUN]=(float(line.split()[7])-load_time)*(1000000./Run_Quantity)
            
                
            if 'Over' in line:
                accuracy_list[condition][int(line.split()[1])][T-1][RUN]=(float(line.split()[3])*(1000000./Run_Quantity))
            
    figT, Ax_Arr_T = plt.subplots(3, 2)#, sharex=True, sharey='row')

    column_titles=['Runtime (s) per million','Errors per million']
    row_titles=['B/D','CC','None']
    k_colours={0:'dodgerblue',5:'darkgreen',20:'darkorange',10:'firebrick'}
    K_ms={0:'o',5:'D',20:'X',10:'s'}
    z_orders={0:10,5:9,10:8,20:7}
    Ks=[0,10]
    
    #return timing_list
    #return accuracy_list
    for n,ax in enumerate(Ax_Arr_T.reshape(-1)):
        #print "on ax ",n
        if n<2:
            ax.set_title(column_titles[n])                 
        if n%2==0:
            ax.set_ylabel(row_titles[n/2])
        
        if n%2==0: #Timing plots
            results=timing_list[n/2]
            
        else: #Accuracy plots
            results=accuracy_list[n/2]
        for K in Ks:
            ax.plot([T*4 for T in xrange(1,Ts+1)],np.mean(results[K],axis=1),markeredgecolor=k_colours[K],ls='',markeredgewidth=1.25,markerfacecolor='none',marker=K_ms[K],markersize=7,zorder=z_orders[K])
                
            ax.plot([T*4 for T in xrange(1,Ts+1)],np.mean(results[K],axis=1),c=k_colours[K],ls='--',lw=0.75,alpha=0.75,zorder=z_orders[K])
            ax.errorbar([T*4 for T in xrange(1,Ts+1)],np.mean(results[K],axis=1),yerr=np.std(results[K],axis=1,ddof=1)/np.sqrt(R),c=k_colours[K],ls='--',zorder=z_orders[K],lw=2)
            if K==20 and n==5:
                print results[K][3]
                print np.std(results[K],axis=1,ddof=1,dtype=np.float64)/np.sqrt(R)
                print np.mean(results[K],axis=1,dtype=np.float64)
                ax.plot([16],[0.02666667],markeredgecolor=k_colours[K],ls='',markeredgewidth=1.25,markerfacecolor='none',marker=K_ms[K],markersize=7,zorder=z_orders[K])
            #if n==0 and K==0:
            #    return results[K],np.std(results[K],axis=1,ddof=1)/np.sqrt(R)

                
        if n%2==0:
            ax.set_yscale('log',nonposy='mask')
            #ax.set_xscale('log',nonposy='mask')
        if n%2==1:
            ax.set_yscale('log',nonposy='mask')
        
        
    for n,ax in enumerate(Ax_Arr_T.reshape(-1)):
        if n<4:
            ax.set_xticks([])
        if n%2==1:
            ax.yaxis.tick_right()
            #ax.set_yticks([])
            #ax.minorticks_on()
            #continue
        if n>=4:
            ax.xaxis.set(ticks=[T*4 for T in xrange(1,Ts+1)],ticklabels=[T for T in xrange(1,Ts+1)])
            
    ax0=Ax_Arr_T.reshape(-1)[0]
    ax0.text(30,37.5,'Graph-like',ha='left',va='center',fontsize=15)
    ax0.text(30,250,'K=5',ha='left',va='center',fontsize=15)
    ax0.text(30,500,'K=10',ha='left',va='center',fontsize=15)
    ax0.text(30,1000,'K=20',ha='left',va='center',fontsize=15)

    Ax_Arr_T.reshape(-1)[1].set_ylim(Ax_Arr_T.reshape(-1)[3].get_ylim())
    figT.text(0.5, 0.02, r'$N_{T}$', ha='center', va='center')


    fig2, (S_ax,A_ax) = plt.subplots(2,1,sharex=True)
    results=timing_list[2]
    for K in Ks:
        S_ax.plot([T*4 for T in xrange(1,Ts+1)],np.mean(results[K],axis=1),markeredgecolor=k_colours[K],ls='',markeredgewidth=1.25,markerfacecolor='none',marker=K_ms[K],markersize=7,zorder=z_orders[K])
        
        S_ax.plot([T*4 for T in xrange(1,Ts+1)],np.mean(results[K],axis=1),c=k_colours[K],ls='',lw=0.75,alpha=0.75,zorder=z_orders[K])
        S_ax.errorbar([T*4 for T in xrange(1,Ts+1)],np.mean(results[K],axis=1),yerr=np.std(results[K],axis=1,ddof=1)/np.sqrt(R),c=k_colours[K],ls='',zorder=z_orders[K],lw=2)

    results=accuracy_list[2]
    for K in Ks:
        A_ax.plot([T*4 for T in xrange(1,Ts+1)],np.mean(results[K],axis=1),markeredgecolor=k_colours[K],ls='',markeredgewidth=1.25,markerfacecolor='none',marker=K_ms[K],markersize=7,zorder=z_orders[K])
        
        A_ax.plot([T*4 for T in xrange(1,Ts+1)],np.mean(results[K],axis=1),c=k_colours[K],ls='',lw=0.75,alpha=0.75,zorder=z_orders[K])
        A_ax.errorbar([T*4 for T in xrange(1,Ts+1)],np.mean(results[K],axis=1),yerr=np.std(results[K],axis=1,ddof=1)/np.sqrt(R),c=k_colours[K],ls='',zorder=z_orders[K],lw=2)

    S_ax.xaxis.set(ticks=[T*4 for T in xrange(1,Ts+1)],ticklabels=[T for T in xrange(1,Ts+1)])
    A_ax.xaxis.set(ticks=[T*4 for T in xrange(1,Ts+1)],ticklabels=[T for T in xrange(1,Ts+1)])
    S_ax.set_yscale('log',nonposy='mask')
    A_ax.set_yscale('log',nonposy='mask')

    A_ax.set_xlabel(r'$N$')
    A_ax.set_ylabel(r'Misclassifications')
    S_ax.set_ylabel(r'CPU Time (s)')

    S_ax.text(10,65,'K=10',ha='left',va='bottom',fontsize=15)
    S_ax.text(10,6.1,'Graph',ha='left',va='top',fontsize=15)

    plt.figure()
    results=timing_list[2]
    rats=np.mean(results[10],axis=1)/np.mean(results[0],axis=1)
    plt.plot(range(1,Ts+1),rats)
    
    
    

    plt.show(block=False)
    #return accuracy_list[2]

import re
from scipy import stats
def LoadNewtimings(sbst):
    timings=[]
    tempL=[]
    for line in open('/rscratch/asl47/binaries/timingS{}.txt'.format(sbst)):
        if 'user' not in line:
            continue
        else:
            tempL.append(getSeconds(re.split('(\d+)', line.split()[1])))
        if len(tempL)==3:
            timings.append(tempL)
            tempL=[]
            
    raw=np.array(timings)
    raw[:,1]-=raw[:,0]
    raw[:,2]-=raw[:,0]
    ratios=raw[:,1]/raw[:,2]
    var=np.mean(ratios)*np.sqrt( (stats.sem(raw[:,1])/np.mean(raw[:,1]))**2 + (stats.sem(raw[:,2])/np.mean(raw[:,2]))**2 )
    return np.mean(ratios),var

def LoadNewtimingserr():
    errs=0
    tot=0
    rats=[]
    for line in open('/rscratch/asl47/Bulk_Run/Regulation/Evolution_T20_C100_N500_K500_M1_R0_I2.txt'):
        if 'Mc:' not in line:
            tot=sum([int(i) for i in line[line.index('A:')+2:].split()])
            continue
        else:
            errs=int(line.split()[-1])
            rats.append(float(errs)/((500.*tot)/1000000.))
    raw=np.array(rats[12:])
    return np.mean(rats),stats.sem(rats)

def getSeconds(timein):
    return int(timein[1])*60+int(timein[3])+int(timein[5])/1000.

def LoadNewGCs(T):
    raw_time=[]
    raw_err=[]
    for i in xrange(10):
        timings=[]
        errors=[]
        for line in open('/rscratch/asl47/binaries/GC_Space_{}_{}.txt'.format(T,i)):
            if 'Runtime' in line:
                timings.append(float(line.split()[2]))
            elif 'Over' in line:
                errors.append(int(line.split()[3])+int(line.split()[5]))
            if len(timings)==2:
                raw_time.append(timings)
                timings=[]
            elif len(errors)==2:
                raw_err.append(errors)
                errors=[]
                
    raw=np.array(raw_time)
    ratios=raw[:,1]/raw[:,0]
    var=np.mean(ratios)*np.sqrt( (stats.sem(raw[:,0])/np.mean(raw[:,0]))**2 + (stats.sem(raw[:,1])/np.mean(raw[:,1]))**2 )
    rawe=np.array(raw_err)
    
    
    return np.mean(ratios),var,np.mean(rawe[:,1]),stats.sem(rawe[:,1])
from matplotlib import gridspec
def PlotNewComp():
    gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1])
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1],sharey=ax)

    for i in range(1,6)+[7,9,11]:
        data=LoadNewGCs(i)
        ax.errorbar(1+data[2],data[0],yerr=data[1],xerr=data[3],fmt='o')
        ax.annotate(r'$\mathrm{{GP}}_{{{}}}$'.format(i),xy=(1+data[2], data[0]), xytext=(-10, 10),textcoords='offset points', ha='right', va='bottom',arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
        ax2.errorbar(1+data[2],data[0],yerr=data[1],xerr=data[3],fmt='o')
        ax2.annotate(r'$\mathrm{{GP}}_{{{}}}$'.format(i),xy=(1+data[2], data[0]), xytext=(-10, 10),textcoords='offset points', ha='right', va='bottom',arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))


    for i in [(2,8),(3,10),(20,100)]:
        t=LoadNewtimings(''.join(map(str,i)))
        er=[]
        if i[0]!=20:
            er=(0,0)
        else:
            er=LoadNewtimingserr()
        ax.errorbar(1+er[0],t[0],yerr=t[1],xerr=er[1],fmt='o')
        ax.annotate(r'$S_{{{},{}}}$'.format(*i),xy=(1+er[0],t[0]), xytext=(-10, 10),textcoords='offset points', ha='right', va='bottom',arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
        ax2.errorbar(1+er[0],t[0],yerr=t[1],xerr=er[1],fmt='o')
        ax2.annotate(r'$S_{{{},{}}}$'.format(*i),xy=(1+er[0],t[0]), xytext=(-10, 10),textcoords='offset points', ha='right', va='bottom',arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

    ax.set_xlabel('Accuracy Improvement')
    ax.set_ylabel('Speed Improvement')
    ax.set_yscale('log')
    ax.set_xscale('log')

    ax2.set_yscale('log')
    ax2.set_xscale('log')
    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax.yaxis.tick_left()
    ax.tick_params(labelright='off')  # don't put tick labels at the top
    ax2.yaxis.tick_right()

    d = .015  # how big to make the diagonal lines in axes coordinates
  
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((1 - d, 1 + d),(1 - d, 1 + d), **kwargs)        # top-left diagonal
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d*4, +d*4), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((-d*4, +d*4),(-d, +d), **kwargs) 


    plt.show(block=False)
    
