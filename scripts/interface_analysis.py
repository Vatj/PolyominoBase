import numpy as np
from copy import deepcopy

import sys

import pickle
from multiprocessing import Pool
from functools import partial
from collections import defaultdict,Counter
from itertools import combinations
from numpy import linalg as LA




def setBasePath(path):
     global BASE_FILE_PATH
     if path=='scratch':
          BASE_FILE_PATH='/scratch/asl47/Data_Runs/Bulk_Data/{0}_I{1}_T{2:.6f}_Mu{3:.6f}_Gamma{4:.6f}_Run{5}.txt'
     elif path=='rscratch':
          BASE_FILE_PATH='/rscratch/asl47/Bulk_Run/Interfaces/{0}_I{1}_T{2:.6f}_Mu{3:.6f}_Gamma{4:.6f}_Run{5}.txt'
     else:
          BASE_FILE_PATH='/home/icyhawaiian/Documents/Data/{0}_I{1}_T{2:.6f}_Mu{3:.6f}_Gamma{4:.6f}_Run{5}.txt'

setBasePath('')
          
interface_length=64
interface_type={8:np.uint8,16:np.uint16,32:np.uint32,64:np.uint64}[interface_length]

def set_length(length):
     global interface_length
     interface_length=length
     global interface_type
     interface_type={8:np.uint8,16:np.uint16,32:np.uint32,64:np.uint64}[interface_length]
def convint(x):
     return interface_type(x)


 
def LoadEvolutionHistory(temperature=0.000001,mu=1,gamma=1,run=0):
     phen_line=True
     phenotype_IDs=[]
     selections=[]
     for line in open(BASE_FILE_PATH.format('PhenotypeHistory',interface_length,temperature,mu,gamma,run)):
          converted=[int(i) for i in line.split()]
          if phen_line:
               phens=[]
               for i in xrange(0,len(converted),2):
                    phens.append((converted[i],converted[i+1]))
               phenotype_IDs.append(phens)
               phen_line=False
          else:
               selections.append(converted)
               phen_line=True
               continue
     
     return phenotype_IDs,np.array(selections,np.uint16)

def LoadGenotypeHistory(n_tiles,temperature=0.000001,mu=1,gamma=1,run=0):
     genotypes=[]
     for line in open(BASE_FILE_PATH.format('GenotypeHistory',interface_length,temperature,mu,gamma,run)):
          converted=[int(i) for i in line.split()]
          genotypes.append([[int(i) for i in j] for j in [converted[i:i + 4*n_tiles] for i in xrange(0, len(converted), 4*n_tiles)]])

     return np.array(genotypes,dtype=interface_type)

def LoadStrengthHistory(temperature=0.000001,mu=1,gamma=1,run=0):
    strengths=[]
    for line in open(BASE_FILE_PATH.format('Strengths',interface_length,temperature,mu,gamma,run)):
         strengths.append([[tuple(np.uint8(i) for i in j.split()) for j in tmp.split(',')[:-1]] for tmp in line.split('.')[:-1]])
    return strengths

def LoadStrengthHistory(temperature=0.000001,mu=1,gamma=1,run=0):
    strengths=[]
    for line in open(BASE_FILE_PATH.format('Strengths',interface_length,temperature,mu,gamma,run)):
         strengths.append([[tuple(np.uint8(i) for i in j.split()) for j in tmp.split(',')[:-1]] for tmp in line.split('.')[:-1]])
    return strengths

def LoadPhenotypeTable(temperature=0.000001,mu=1,gamma=1,run=0):
     phenotype_table= sorted([[int(i) for i in line.split()] for line in open(BASE_FILE_PATH.format('Phenotypes',64,temperature,mu,gamma,run))],key=lambda z: z[0])
     return {tuple(px[:2]): px[2:] for px in phenotype_table}

def LoadT(mu=1,t=0.35,run=0):
     g=LoadGenotypeHistory(3,mu=mu,temperature=t,run=run)
     st=LoadStrengthHistory(mu=mu,temperature=t,run=run)
     p,s=LoadEvolutionHistory(mu=mu,temperature=t,run=run)
     return (g,s,p,st)

def BindingStrength(base1,base2):
     return 1-bin(np.bitwise_xor(convint(base1),revbits(base2))).count("1")/float(interface_length)

def revbits(x):
     return interface_type(int(bin(~convint(x))[2:].zfill(interface_length)[::-1], 2))


""" DRIFT SECTION """
def RandomWalk(I_size=32,n_steps=1000,phi=0.5,T_star=0.6,renorm=False,return_prog=False):
     s_hats=np.linspace(0,1,I_size+1)
     N=int(I_size*(1-T_star))+1
     states=np.zeros(N)
     states[0]=1
     progressive_states=[]

     for i in xrange(n_steps):
          states=UpdateStates(states,s_hats[I_size-N+1:],phi,renorm)
          progressive_states.append(np.sum(s_hats[I_size-N+1:]*states))
     if not renorm:
          states/=np.sum(states)
     analytic_states=getSteadyStates(mmatrix(N,phi,s_hats[I_size-N+1:]))[1]
     #print analytic_states
     return analytic_states,s_hats[I_size-N+1:]
     if return_prog:
          progressive_states.append(np.sum(s_hats[I_size-N+1:]*analytic_states))
          return progressive_states
     
     fig, ax1 = plt.subplots()
     msize=12
     ax1.plot(np.linspace(T_star,1,N),states,marker='x',ls='',mew=1,ms=14,mec='orangered',mfc='None')
     ax1.plot(np.linspace(T_star,1,N),analytic_states,marker='o',ls='',ms=14,mew=1,mfc='None',mec='royalblue')
     #ax1.plot([s_hats[int(T_star*I_size)-1]]+[s_hats[int(T_star*I_size)]]*2,[0,0,1],c='royalblue',mec='royalblue',mfc='None',mew=1,ls='--',lw=0.75,marker='s',ms=msize)

     plt.xlabel(r'$\hat{S}$')
     #ax1.text(0.4,0.05,r'$Pr_{\mathrm{binding}}$')
     #ax1.text(0.8,0.1,r'genetic drift')
     #ax1.axes.get_yaxis().set_visible(False)
     #ax1.set_xlim((.3,1.05))
     plt.show(block=False)
     


def UpdateStates(states,val,phi=0.25,renorm=False):
     states_updating=deepcopy(states)
     for i in xrange(states.shape[0]):
          states_updating[i]-=states[i]*phi
          if i!=0:
               states_updating[i]+=states[i-1]*phi*(1-val[i-1])
          if i!=states.shape[0]-1:
               states_updating[i]+=states[i+1]*phi*val[i+1]
     return states_updating if not renorm else states_updating/np.sum(states_updating)

def mmatrix(N_states,mu,val):
     rows=[[1-mu,mu*(1-val[0])]+[0]*(N_states-2)]
     for i in xrange(1,N_states-1):
          rows.append([0]*(i-1)+[mu*val[i],1-mu,mu*(1-val[i])]+[0]*(N_states-2-i))
     rows.append([0]*(N_states-2)+[mu*val[-1],1-mu])
     return np.vstack(rows).T

def getSteadyStates(matrix):
     eigval,eigvec=LA.eig(matrix)
     for va,ve in zip(eigval,eigvec.T):
          if np.all(np.sign(ve)==-1) or np.all(np.sign(ve)==1):
               return va,ve/sum(ve)          



""" PHYLOGENCY SECTION """     

def writeResults(I,M,t,runs,offset=0):
     global interface_length
     interface_length=I
     global interface_type
     interface_type={8:np.uint8,16:np.uint16,32:np.uint32,64:np.uint64}[interface_length]
     setBasePath('scratch')
     pool = Pool()
     data_struct=pool.map(partial(AnalysePhylogeneticStrengths, mu=M,t=t), xrange(offset,offset+runs)) 
     pool.close()
     with open('I{}Mu{}T{}O{}.pkl'.format(I,M,t,offset), 'wb') as f:
          pickle.dump(data_struct, f)

def loadResults(I,M,t,offset=0):
     with open('/scratch/asl47/Data_Runs/Bulk_Data/I{}Mu{}T{}O{}.pkl'.format(I,M,t,offset), 'rb') as f:
          return pickle.load(f)

def loadManyResults(I,M,t,runs):
     bulk_results=[]
     for run in runs:
          bulk_results.extend(loadResults(I,M,t,run))
     return bulk_results

def concatenateResults(data_struct,trim_gen=True):
     conc_data=defaultdict(list)
     Wa=Counter()
     Ws=Counter()
     active_interface_count=[Counter() for _ in xrange(len(data_struct[0]['Cnt']))]
     fatal_phens=np.empty((len(data_struct),len(data_struct[0]['Fatal'])))
     phen_trans=defaultdict(lambda: defaultdict(int))
     slice_start=1 if trim_gen else 0
     for outer_index,data in enumerate(data_struct):
          for k,v in data.iteritems():
               if k=='Wa':
                    Wa+=Counter(v)
               elif k=='Ws':
                    Ws+=Counter(v)
               elif k=='Cnt':
                    for i,cn in enumerate(v):
                         active_interface_count[i]+=Counter(cn)
               elif k=='Fatal':
                    fatal_phens[outer_index]=np.array(v)
               elif k=='PhenTrans':
                    for sub_k,sub_v in v.iteritems():
                         for phen in sub_v:
                              phen_trans[sub_k][phen]+=1
               else:
                    conc_data[k].extend([i[slice_start:] for i in v])

     for k,v in conc_data.iteritems():
          if len(v) and 'W' not in k:
               length = len(sorted(v,key=len, reverse=True)[0])
               conc_data[k]=np.array([xi+[np.nan]*(length-len(xi)) for xi in v])
          
     return conc_data,(Wa,Ws),active_interface_count,fatal_phens,phen_trans
          
def AnalysePhylogeneticStrengths(r,mu,t):
     g,s,p,st=LoadT(mu=mu,t=t,run=r)
     pt=LoadPhenotypeTable(temperature=t,mu=mu,run=r)
     trans_probs=PhenOrder(g,pt,p,s,st)

     mae,mai,ms,wa,ws,cnt_ai,fp=qBFS(deepcopy(g),deepcopy(p),deepcopy(s),deepcopy(st))
     return {'AsymE':mae,'AsymI':mai,'Sym':ms,'Wa':wa,'Ws':ws,'Cnt':cnt_ai,'Fatal':fp,'PhenTrans':trans_probs}
     
     
def qBFS(genotypes,phenotypes,selections,strengths):
     #ANALYSIS PARAMETERS
     WEAKNESS_GOBACK=50
     #RUNTIME PARAMETERS
     gen_limit=len(genotypes)
     pop_size=len(genotypes[0])
     MSE_ae=[]
     MSE_ai=[]
     MSE_s=[]
     fatal_phens=[]
     active_interfaces=[]

     W_a=Counter({K:0 for K in np.linspace(0,1,interface_length+1)})
     W_s=Counter({K:0 for K in np.linspace(0,1,interface_length/2+1)})
     for generation in xrange(gen_limit):
          temp_count=defaultdict(int)
          for index in xrange(pop_size):
               if phenotypes[generation][index][0]!=0:
                    temp_count[len(strengths[generation][index])]+=1
          active_interfaces.append(temp_count)
               
               
     for generation in xrange(gen_limit-1):
          fatal_phens.append(0)
          for index in xrange(pop_size):
               for new_bond in strengths[generation][index]:
                    stren_tree=qDFS(index,new_bond,genotypes[generation:],selections[generation:],strengths[generation:])
                    if stren_tree:
                         stren_tree.insert(0,generation)
                         if new_bond[0]==new_bond[1]:
                              MSE_s.append(stren_tree)
                         elif new_bond[0]/4==new_bond[1]/4:
                              MSE_ai.append(stren_tree)
                         else:
                              MSE_ae.append(stren_tree)
               if phenotypes[generation][index][0]==0:
                    fatal_phens[-1]+=1
               elif generation>gen_limit-WEAKNESS_GOBACK:
                    weak_a,weak_s=qWeak(strengths[generation][index],genotypes[generation][index])
                    W_a+=weak_a
                    W_s+=weak_s
          
     return filter(None,MSE_ae),filter(None,MSE_ai),filter(None,MSE_s),dict(W_a),dict(W_s),active_interfaces,fatal_phens

def qWeak(strength,genotype):
     weaknesses_asym=Counter({K:0 for K in np.linspace(0,1,interface_length+1)})
     weaknesses_sym=Counter({K:0 for K in np.linspace(0,1,interface_length/2+1)})
     g_length=12
     used_interfaces=set([item for sublist in strength for item in sublist])
     non_interactings=set(xrange(g_length))-used_interfaces
     for pair in combinations(non_interactings,2):
          weaknesses_asym[BindingStrength(*genotype[[i for i in pair]])]+=1
     for face in non_interactings:
          weaknesses_sym[BindingStrength(*genotype[[face]*2])]+=1
     return weaknesses_asym,weaknesses_sym

def qDFS(index,new_bond,genotypes,selections,strengths):
     #ANALYSIS PARAMETERS
     MIN_LEN=50
     #RUNTIME PARAMETERS
     gen_limit=len(genotypes)
     pop_size=len(genotypes[0])

     mean_str=[{BindingStrength(*genotypes[0][index][[i for i in new_bond]]):1}]
     strengths[0][index].remove(new_bond)
     descendents=list(np.where(selections[0]==index)[0])

     gen_index=1
     while gen_index<gen_limit and descendents:
          str_temp=defaultdict(float)
          descendents_temp=[]
          for descendent in descendents:
               if new_bond in strengths[gen_index][descendent]:
                    str_temp[BindingStrength(*genotypes[gen_index][descendent][[i for i in new_bond]])]+=1
                    descendents_temp.extend(np.where(selections[gen_index]==descendent)[0])
                    strengths[gen_index][descendent].remove(new_bond)
          descendents=descendents_temp
          gen_index+=1

          mean_str.append(str_temp)
          
     return mean_str if len(mean_str)>MIN_LEN else None


""" Phenotype Matching """
def PhenOrder(genotypes,phen_dict,phenotype_IDs,selections,strengths):
     gen_limit=len(genotypes)
     pop_size=len(genotypes[0])
     

     SEARCH_PHENS=[[4,1,1,5,7,3],[4,4,0,0,1,0,4,5,6,0,0,8,7,2,0,3,0,0]]
     PHEN_KEYS={}

     transition_probabilities=defaultdict(lambda: defaultdict(int))
     
     for key,phen in phen_dict.iteritems():
          for i,phen_key in enumerate(SEARCH_PHENS):
               if phen_key == phen:
                    PHEN_KEYS[tuple(key)]=phen
 	            
     for generation in xrange(gen_limit-1):
          for index in xrange(pop_size):
               for phen_key,phen_value in PHEN_KEYS.iteritems():
                    if phen_key==phenotype_IDs[generation][index]:
                         sel=selections[generation-1][index]
                         #if phen_key[0]==8:
                         #     print generation,index
                         #     print genotypes[generation][index]
                         #     print genotypes[generation-1][sel]
                         #     print phen_dict[phenotype_IDs[generation-1][sel]]
                         if phenotype_IDs[generation-1][sel]!=phen_key:
                              transition_probabilities[tuple(phen_value)][tuple(phen_dict[phenotype_IDs[generation-1][sel]])]+=1
                         #PurgeDescendents(generation,index,selections,phenotype_IDs)

     TP_dict=dict(transition_probabilities)
     for key in TP_dict:
          TP_dict[key]=dict(TP_dict[key])
     return TP_dict

def PurgeDescendents(generation,index,selections,phenotype_IDs):
     descendents=list(np.where(selections[generation]==index)[0])
     gen_index=generation+1
     while gen_index<len(selections) and descendents:
          descendents_temp=[]
          
          for descendent in descendents:
               phenotype_IDs[gen_index][descendent]=(0,0)
               descendents_temp.extend(np.where(selections[gen_index]==descendent)[0])

          descendents=descendents_temp
          gen_index+=1

def PhenotypicTransitions(phen_trans):
     return "hi"
     
     
     
def main(argv):
     writeResults(int(argv[1]),*(float(i) for i in argv[2:4]),runs=int(argv[4]),offset=int(argv[5]))
     return
if __name__ == "__main__":
    main(sys.argv)






    

"""from random import randint
from collections import defaultdict
def RandomHistorySampling(genotypes,selections,phenotypes,strengths,goback=10):
     printer=False
     rg=randint(genotypes.shape[0]/2,genotypes.shape[0]-1)
     rp=randint(0,genotypes.shape[1]-1)
     strength_tracker=defaultdict(list)
     if(goback>rg):
          rg=goback
     assert rg>=goback, "going back too far"
     for stren in strengths[rg][rp]:
          strength_tracker[stren].append([BindingStrength(*[genotypes[rg][rp][j] for j in stren])])    
     for bg in xrange(1,goback+1):
          rp=selections[rg-bg][rp]
          for stren in strengths[rg-bg][rp]:
               if stren in strength_tracker:
                    strength_tracker[stren][-1].append(BindingStrength(*[genotypes[rg-bg][rp][j] for j in stren]))
               else:
                    strength_tracker[stren].append([BindingStrength(*[genotypes[rg-bg][rp][j] for j in stren])])
          for stren in strength_tracker.keys():
               if stren not in strengths[rg-bg][rp] and len(strength_tracker[stren][-1])!=0 :
                    break
                    strength_tracker[stren].append([])                   
     return strength_tracker



def xx(mu,runs):
     plt.figure()
     ups=0
     downs=0
     for r in xrange(runs):
          g,s,p,st =LoadT(mu=mu,run=r)
          print r
          ma,ms=qBFS(g,s,st)
          #return ma,ms
          u,d=pp(ma,ms,True)
          ups+=u
          downs+=d

     plt.ylabel(r'$\hat{S}$')
     plt.xlabel('generations')
     plt.text(100,.8,'+={},-={}'.format(ups,downs))
     plt.show(block=False)
"""
