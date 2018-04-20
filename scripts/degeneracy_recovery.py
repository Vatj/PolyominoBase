###############################
###DEGENERACY RECOVERY CODE####
###  Alexander S Leonard   ####
###############################

import matplotlib.pyplot as plt

from math import factorial
from operator import mul
from itertools import combinations_with_replacement

def choose(n, k):
    ntok = 1
    ktok = 1
    for t in xrange(1, min(k, n - k) + 1):
        ntok *= n
        ktok *= t
        n -= 1
    return ntok /  ktok


def product_of(list_iterable):
    return reduce(mul, list_iterable, 1L)

def Find_Neutral_Degeneracy(F,N,P):
    print [[choose(F,U)*N**(F-U)]+[[product_of([P-2*Upp for Upp in xrange(0,Up-1+1)])]+[sum(product_of(factors) for factors in combinations_with_replacement(xrange(1,Up+1),U-Up))] for Up in xrange(0,U+1)] for U in xrange(0,F+1)]
    return sum([choose(F,U)*N**(F-U)*sum([product_of([P-2*Upp for Upp in xrange(0,Up-1+1)])*sum(product_of(factors) for factors in combinations_with_replacement(xrange(1,Up+1),U-Up)) for Up in xrange(0,U+1)]) for U in xrange(0,F+1)])
   
def Find_Interacting_Degeneracy(N_Interacting_Colours,N_Interacting_Tuples):
    return product_of([N_Interacting_Colours-2*x for x in xrange(0,N_Interacting_Tuples)])

def Find_Tile_Set_Interchange_Degeneracy(Tiles):
    return float(factorial(len(Tiles)))/product_of([factorial(Tiles.count(list(unique_tile))) for unique_tile in set(tuple(tile) for tile in Tiles)])


def Find_Tile_Rotational_Degeneracy(Tile):
    zero_count=Tile.count(0)
    if zero_count==2:
        if Tile[0]==0 and Tile[1]==0:
            return 4
        else: # Tile[0]==0 and Tile[2]==0:
            if Tile[1]==Tile[3]:
                return 2
            else:
                return 4;
    elif zero_count==0:
        return Find_Tile_Rotational_Degeneracy_Without_Neutrals(Tile)
    else:
        return 1 if zero_count==4 else 4
        
        
def Find_Tile_Rotational_Degeneracy_Without_Neutrals(Tile):
    unique_Faces=len(set(tuple(Tile)))
    if unique_Faces==2:
        if Tile.count(Tile[0])==2:
            if Tile[0]==Tile[1] and Tile[2]==Tile[3]: #1,1,2,2 type topology
                return 4
            elif Tile[0]==Tile[2] and Tile[1]==Tile[3]: #1,2,1,2 type topology
                return 2
        else: #1,1,1,2 type topology
            return 4
    else:
        return 1 if unique_Faces==1 else 4
        

def Find_Overall_Degeneracy(Tile_Set,N_Colours,verbose=False):
    Combined_Degeneracy=1
    #Neutal_Space
    Tile_SetF=[face for tile in Tile_Set for face in tile]
    if Tile_SetF.count(0)>=1:
        Combined_Degeneracy*=Find_Neutral_Degeneracy(Tile_SetF.count(0),2,N_Colours-max(Tile_SetF)-2)
        if verbose:
            print "Neutral: ",Find_Neutral_Degeneracy(Tile_SetF.count(0),2,N_Colours-max(Tile_SetF)-2)
        
    #Tile Interchange
    Combined_Degeneracy*=Find_Tile_Set_Interchange_Degeneracy(Tile_Set)
    if verbose:
        print "interchange: ",Find_Tile_Set_Interchange_Degeneracy(Tile_Set)
    #Individual Tile Symmetry

    #for i,Tile in enumerate(list(set([tuple(tt) for tt in Tile_Set]))):
    
    for i,Tile in enumerate(Tile_Set):
        Combined_Degeneracy*=Find_Tile_Rotational_Degeneracy(Tile)
        if verbose:
            print "Tile: ",Find_Tile_Rotational_Degeneracy(Tile)
    #Interacting Pair Degeneracy
    Combined_Degeneracy*=Find_Interacting_Degeneracy(N_Colours-2,max(Tile_SetF)/2)
    if verbose:
        print "Interacting: ",Find_Interacting_Degeneracy(N_Colours-2,max(Tile_SetF)/2)

    #Combined_Degeneracy/=Find_Relabelling_Collisions(Tile_Set)
    if verbose:
        print "Relabels: ",Find_Relabelling_Collisions(Tile_Set)
    return int(Combined_Degeneracy)


def Find_Relabelling_Collisions(Tile_Set):
    duplicate_Factor=1
    flat_Tiles=[item for sublist in Tile_Set for item in sublist]
    
    for i in xrange(1,max(flat_Tiles)+2,2):
        if i in flat_Tiles:
            if flat_Tiles.count(i)==flat_Tiles.count(i+1):
                duplicate_Factor*=2

    return duplicate_Factor
    
#def Effective_Pair_Labelling(Tile    


def Test_Degen():
    lineF=[line.rstrip('\n') for line in open('/rscratch/asl47/Topologies_Multi_C10.txt')]
    Counts=[int(ln.split()[-1]) for ln in lineF]

    linePT=[line.rstrip('\n') for line in open('/rscratch/asl47/Topologies_RAW_C10_V2.txt')]
    lineFP=[[int(ite) for ite in ln.split()] for ln in linePT]
    misM=0
    misS=0
    for i,line in enumerate(lineFP):
        misS+=Find_Overall_Degeneracy([line[:4],line[4:]],10)
        if Find_Overall_Degeneracy([line[:4],line[4:]],10)!=Counts[i]:
            misM+=1
            if Find_Overall_Degeneracy([line[:4],line[4:]],10)<Counts[i]:
                print "Not ideal"
            else:
                print line,"N:",i,"  off by factor ",Find_Overall_Degeneracy([line[:4],line[4:]],10)*1./Counts[i]
            #print "fixed attempt is ", Find_Overall_Degeneracy([line[:4],line[4:]],10)*1./Find_Relabelling_Collisions([line[:4],line[4:]])/Counts[i]
            
    print "Mismatch total was ",misM
    print "total recovery was ",misS




def Compare_Methods():
    F_Max=16
    C_Max=50
    
    z=[[0 for x in xrange(1,F_Max+1)] for y in xrange(0,C_Max+1,2)]
    print z[0]
    for F in xrange(1,F_Max+1):
        print "on F ",F
        for C in xrange(0,C_Max+1,2):
            print "on C ",C
            #rint "on ",F," and ",C
            #print "filling indexes ",C/2," ",F-1
            z[C/2][F-1]=Find_Neutral_Degeneracy(F,C,0,[],0)-Find_Neutral_Degeneracy_Experimental_2(F,2,C)
            

            

    #print z
    plt.imshow(z,cmap='viridis', interpolation='none')
    plt.show(block=False)

def Plot_Degens():
    from matplotlib.colors import LogNorm
    F_Max=12
    C_Max=40
    z=[[0 for x in xrange(1,F_Max+1)] for y in xrange(0,C_Max+1,2)]
    #print z[0]
    for F in xrange(1,F_Max+1):
        #print "on F ",F
        for C in xrange(0,C_Max+1,2):
            #print "on C ",C
            #rint "on ",F," and ",C
            #print "filling indexes ",C/2," ",F-1
            z[C/2][F-1]=Find_Neutral_Degeneracy(F,2,C)
    #print z
    plt.imshow(z,cmap='viridis', interpolation='nearest',norm=LogNorm())
    plt.show(block=False)
            
#######################################
#######################################
########TOPOLOGY GENERATOR CODE########
#######################################
#######################################
    
##############
##TODO###
########

#Several places hardcoded, particularly for the 3 and above faces, need to generalise
#Some sort of sanity check that count of conjugates isn't that much higher than count of faces
#Same for count of higher pairs isn't greater than count of lower pairs
#symmetry handling for full tile


def Topology_Generator(Input_Stubs):
    Overall_Topology=[]
    
    for stub in Input_Stubs:
        if stub==0:
            Overall_Topology.append(0)
        elif 1 not in Overall_Topology:
            Overall_Topology.append(1)
        else:
            Overall_Topology.append('*')

    if '*' not in Overall_Topology:
        return Overall_Topology
    else:
        #Start the recursion
        for launch in xrange(Overall_Topology.count('*')):
            replaceIndex=Overall_Topology.index('*')
            for iterate in xrange(launch):
                replaceIndex+=Overall_Topology[replaceIndex+1:].index('*')+1
                
            New_Topology=Overall_Topology[:]
            New_Topology[replaceIndex]=2           
            Conj_Loc={2:replaceIndex}
            Recursive_Topology_Step(New_Topology,Conj_Loc)
                             

def Recursive_Topology_Step(Recursive_Topology,Conjugate_Locations):
    if '*' not in Recursive_Topology:
        for sliceRange in xrange(0,len(Recursive_Topology)-4,4):
            if Check_Smaller_Tile(Recursive_Topology[sliceRange:sliceRange+4],Recursive_Topology[sliceRange+4:sliceRange+8]):
                #print "REJECTED: ",Recursive_Topology
                break
            elif Check_More_Fundemental_Topology_Tile(Recursive_Topology[sliceRange:sliceRange+4],Recursive_Topology[sliceRange+4:sliceRange+8]):
                #print "REJECTED2: ",Recursive_Topology
                break
        else:
            if Recursive_Topology[:4]==[0,1,1,2]:
                print Recursive_Topology
    else:
        nextIndex=Recursive_Topology.index('*')

        New_Topology=Recursive_Topology[:]
        New_Topology[nextIndex]=1
        Conj_Loc=dict(Conjugate_Locations)
        Recursive_Topology_Step(New_Topology,Conj_Loc)
        if nextIndex > Conjugate_Locations[2]:
            New_Topology=Recursive_Topology[:]
            New_Topology[nextIndex]=2
            Recursive_Topology_Step(New_Topology,Conj_Loc)

        if Recursive_Topology.count('*') >=2:
            New_Topology=Recursive_Topology[:]
            New_Topology[nextIndex]=3
            nextIndex=New_Topology.index('*')
            New_Topology[nextIndex]=4
            Conj_Loc=dict(Conjugate_Locations)
            Conj_Loc[4]=nextIndex
            Recursive_Topology_Step(New_Topology,Conj_Loc)


def Check_Smaller_Tile(First_Tile,Second_Tile):
    """Check if the second tile is smaller than the first, and thus reject it.

    Return:
    True if second is larger
    """
    for face in xrange(4):
        if First_Tile[face]<Second_Tile[face]:
            return False
        elif First_Tile[face]==Second_Tile[face]:
            continue
        else:
            return True
def Check_More_Fundemental_Topology_Tile(First_Tile,Second_Tile):
    """Check if the second tile is more fundmental than the first, and thus reject it.

    Return:
    True if second is more fundemental
    """
    First_Clone=First_Tile[:]
    Second_Clone=Second_Tile[:]
    Relabel_Tile(First_Clone)
    Relabel_Tile(Second_Clone)
    return Check_Smaller_Tile(First_Clone,Second_Clone)

def Relabel_Tile(Tile):
    Base_Replacements={}
    active_Relabelling=1
    for face in Tile:
        if face==0:
            continue
        else:
            if face not in Base_Replacements:
                conjugate_Face=(1-face%2)*(face-1)+(face%2)*(face+1)
                Base_Replacements[face]=active_Relabelling
	        Base_Replacements[conjugate_Face]=active_Relabelling+1
	        active_Relabelling+=2
                
    for i,face in enumerate(Tile):
        if face==0:
            continue
        Tile[i]=Base_Replacements[face]

