import networkx as nx

def Transform_Graph_From_List(tile_kit):
    graph_kit=nx.MultiDiGraph()
    graph_kit.add_nodes_from(xrange(len(tile_kit)))

    #Add edges for internal structure in clockwise orientation 
    for internal_edge in xrange(len(tile_kit)/4):
        graph_kit.add_edge(internal_edge*4+0,internal_edge*4+1)#,color='k')
        graph_kit.add_edge(internal_edge*4+1,internal_edge*4+2)#,color='k')
        graph_kit.add_edge(internal_edge*4+2,internal_edge*4+3)#,color='k')
        graph_kit.add_edge(internal_edge*4+3,internal_edge*4+0)#,color='k')

    #Add edges for graph connections
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

import matplotlib.pyplot as plt
def Draw_Graph(graph,kit):
    layout=nx.spring_layout(graph)
    nx.draw_networkx_labels(graph,layout,{j:str(i) for j,i in enumerate(kit)})
    nx.draw_networkx_nodes(graph,layout)
    nx.draw_networkx_edges(graph,layout)
    plt.show(block=False)


from collections import defaultdict

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

from collections import deque
import itertools
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
                
from collections import Counter
import seaborn as sns
import numpy as np
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
import seaborn as sns
def Use_Seaborn():
    sns.set_context("paper",font_scale=2.2)
    sns.set_style("white",rc={"xtick.major.size": 8, "ytick.major.size": 8,"xtick.minor.size":5, "ytick.minor.size": 5,"axes.linewidth": 2,"axes.edgecolor":"darkgray","font.size":8,"axes.titlesize":8,"axes.labelsize":5})
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')


    
from collections import defaultdict
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
        
        S_ax.plot([T*4 for T in xrange(1,Ts+1)],np.mean(results[K],axis=1),c=k_colours[K],ls='--',lw=0.75,alpha=0.75,zorder=z_orders[K])
        S_ax.errorbar([T*4 for T in xrange(1,Ts+1)],np.mean(results[K],axis=1),yerr=np.std(results[K],axis=1,ddof=1)/np.sqrt(R),c=k_colours[K],ls='--',zorder=z_orders[K],lw=2)

    results=accuracy_list[2]
    for K in Ks:
        A_ax.plot([T*4 for T in xrange(1,Ts+1)],np.mean(results[K],axis=1),markeredgecolor=k_colours[K],ls='',markeredgewidth=1.25,markerfacecolor='none',marker=K_ms[K],markersize=7,zorder=z_orders[K])
        
        A_ax.plot([T*4 for T in xrange(1,Ts+1)],np.mean(results[K],axis=1),c=k_colours[K],ls='--',lw=0.75,alpha=0.75,zorder=z_orders[K])
        A_ax.errorbar([T*4 for T in xrange(1,Ts+1)],np.mean(results[K],axis=1),yerr=np.std(results[K],axis=1,ddof=1)/np.sqrt(R),c=k_colours[K],ls='--',zorder=z_orders[K],lw=2)

    S_ax.xaxis.set(ticks=[T*4 for T in xrange(1,Ts+1)],ticklabels=[T for T in xrange(1,Ts+1)])
    A_ax.xaxis.set(ticks=[T*4 for T in xrange(1,Ts+1)],ticklabels=[T for T in xrange(1,Ts+1)])
    S_ax.set_yscale('log',nonposy='mask')
    A_ax.set_yscale('log',nonposy='mask')

    A_ax.set_xlabel(r'$N_{T}$')
    A_ax.set_ylabel(r'Misclassifications')
    S_ax.set_ylabel(r'CPU Time (s)')

    S_ax.text(8,65,'K=10',ha='left',va='bottom',fontsize=15)
    S_ax.text(8,6.1,'Graph',ha='left',va='top',fontsize=15)

    plt.figure()
    results=timing_list[2]
    rats=np.mean(results[10],axis=1)/np.mean(results[0],axis=1)
    plt.plot(range(1,Ts+1),rats)
    
    
    

    plt.show(block=False)
    return accuracy_list[2]

