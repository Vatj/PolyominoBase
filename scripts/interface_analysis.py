import numpy as np
from copy import deepcopy

import sys

import pickle
from multiprocessing import Pool
from functools import partial
from collections import defaultdict,Counter
from itertools import combinations
from numpy import linalg as LA


def GG(N):
     a=[]
     while len(a)<N:
	  r1=random.getrandbits(64)
 	  r2=random.getrandbits(64)
 	  if BindingStrength(r1,r2)>0.7:
               a.append((r1,r2))
     return a
def setBasePath(path):
     global BASE_FILE_PATH
     if path=='scratch':
          BASE_FILE_PATH='/scratch/asl47/Data_Runs/Bulk_Data/{0}_I{1}_T{2:.6f}_Mu{3:.6f}_Gamma{4:.6f}_Run{5}.txt'
     elif path=='rscratch':
          BASE_FILE_PATH='/rscratch/asl47/Pickles/{0}_I{1}_T{2:.6f}_Mu{3:.6f}_Gamma{4:.6f}_Run{5}.txt'
     else:
          BASE_FILE_PATH='/home/icyhawaiian/Documents/Data/{0}_I{1}_T{2:.6f}_Mu{3:.6f}_Gamma{4:.6f}_Run{5}.txt'

setBasePath('rscratch')
          
interface_length=64
interface_type={8:np.uint8,16:np.uint16,32:np.uint32,64:np.uint64}[interface_length]

def set_length(length):
     global interface_length
     interface_length=length
     global interface_type
     interface_type={8:np.uint8,16:np.uint16,32:np.uint32,64:np.uint64}[interface_length]
     
def convint(x):
     return interface_type(x)

def BindingStrength(base1,base2):
     return 1-bin(np.bitwise_xor(convint(base1),revbits(base2))).count("1")/float(interface_length)

def revbits(x):
     return interface_type(int(bin(~convint(x))[2:].zfill(interface_length)[::-1], 2))
 
def LoadEvolutionHistory(temperature=0.000001,mu=1,gamma=1,run=0):
     phen_line=True
     phenotype_IDs=[]
     selections=[]
     for line in open(BASE_FILE_PATH.format('PhenotypeHistory',interface_length,temperature,mu,gamma,run)):
          converted=[int(i) for i in line.split()]
          if phen_line:
               phenotype_IDs.append(zip(*(iter(converted),) * 2))
               phen_line=False
          else:
               selections.append(converted)
               phen_line=True     
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

def LoadPhenotypeTable(temperature=0.000001,mu=1,gamma=1,run=0):
     phenotype_table= sorted([[int(i) for i in line.split()] for line in open(BASE_FILE_PATH.format('Phenotypes',64,temperature,mu,gamma,run))],key=lambda z: z[0])
     return {tuple(px[:2]): tuple(px[2:]) for px in phenotype_table}

def LoadT(mu=1,t=0.35,run=0):
     g=LoadGenotypeHistory(3,mu=mu,temperature=t,run=run)
     st=LoadStrengthHistory(mu=mu,temperature=t,run=run)
     p,s=LoadEvolutionHistory(mu=mu,temperature=t,run=run)
     return (g,s,p,st)



""" DRIFT SECTION """
def RandomWalk(I_size=64,n_steps=1000,phi=0.5,T_star=0.6,renorm=False,return_prog=False):
     if return_prog:
          renorm=True
     s_hats=np.linspace(0,1,I_size+1)
     N=int(I_size*(1-T_star))+1
     analytic_states=getSteadyStates(mmatrix(N,phi,s_hats[-N:]))[1]
     steady_state=np.sum(s_hats[-N:]*analytic_states)
     print steady_state
     states=np.zeros(N)
     states[0]=1
     
     progressive_states=[np.sum(s_hats[-N:]*states)]

     for i in xrange(n_steps):
          states=UpdateStates(states,s_hats[-N:],phi,renorm)
          progressive_states.append(np.sum(s_hats[-N:]*states))
          #print progressive_states
          if progressive_states[-1]>=steady_state*.995:
               print i
               return i
          
     if not renorm:
          states/=np.sum(states)
     

     return mmatrix(N,phi,s_hats[-N:])
     #print analytic_states
     
     if return_prog:
          progressive_states.append(np.sum(s_hats[-N:]*analytic_states))
          return progressive_states
     return analytic_states,s_hats[-N:]
     
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
          bulk_results.extend(loadResults(I,float(M),t,run))
     return bulk_results

def concatenateResults(data_struct,trim_gen=True):
     ##params##
     count_thresh=10
     strengths={K:[] for K in ['E','I','S']}
     neutrals={K:Counter() for K in ['Sym','Asym']}

     #bindings=[Counter() for _ in xrange(len(data_struct[0]['N_binds']))]
     #fatal_phens=np.empty((len(data_struct),len(data_struct[0]['Unbounds'])))
     phen_trans=defaultdict(lambda: defaultdict(int))
     first_trans=defaultdict(lambda: defaultdict(int))
     slice_start=1 if trim_gen else 0
     for outer_index,run_data in enumerate(data_struct):
          for k,v in run_data.iteritems():
               if 'Neutral' in k:
                    neutrals[k.split('_')[-1]]+=Counter(v)
               elif 'Stren' in k:
                    strengths[k.split('_')[-1]].extend([i[slice_start:]for i in v])
               elif k=='N_binds':
                    for i,cn in enumerate(v):
                         bindings[i]+=Counter(cn)
               elif k=='Unbounds':
                    fatal_phens[outer_index]=np.array(v)
               elif k=='PhenTrans':
                    for sub_k,sub_v in v.iteritems():
                         for phen in sub_v:
                              phen_trans[sub_k][phen]+=1
               elif k=='PhenFirsts':
                    for sub_k,sub_v in v.iteritems():
                         first_trans[sub_k][sub_v]+=1
              

     for stren_type,stren_matrix in strengths.iteritems():
          for k in xrange(len(stren_matrix)):
               stren_matrix[k]=[trim for trim in stren_matrix[k] if sum(trim.values())>count_thresh]
          if stren_matrix:
               length = len(sorted(stren_matrix,key=len, reverse=True)[0])
               for i,strens in enumerate(stren_matrix):
                    strengths[stren_type][i]=np.array([np.sum([a*b for a,b in value.iteritems()])/np.sum(value.values()) for value in strens]+[np.nan]*(length-len(strens)),dtype=np.double)
               smarr=np.array(stren_matrix)     
               strengths[stren_type]=smarr[~np.all(np.isnan(smarr),axis=1)]
                         
     return strengths,neutrals,{k:dict(v) for k,v in phen_trans.iteritems()},{k:dict(v) for k,v in first_trans.iteritems()}
#dict(phen_trans),dict(first_trans)#,bindings,fatal_phens

          
def AnalysePhylogeneticStrengths(r,mu,t):
     g,s,p,st=LoadT(mu=mu,t=t,run=r)
     phen_table=LoadPhenotypeTable(temperature=t,mu=mu,run=r)

     trans_probs,first_trans= PhenOrderTail(g,phen_table,p,s,st)
     mae,mai,ms,wa,ws,cnt_ai,fp=qBFS(g,p,s,st)                   
                      
     return {'Stren_E':mae,'Stren_I':mai,'Stren_S':ms,'Neutral_Asym':wa,'Neutral_Sym':ws,'PhenTrans':trans_probs,'PhenFirsts':first_trans}#,'N_binds':cnt_ai,'Unbounds':fp
     
     
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

def qDFS(index,new_bond,genotypes,selections,strengths):
     detailed_distr=True
     #ANALYSIS PARAMETERS
     MIN_LEN=50
     #RUNTIME PARAMETERS
     gen_limit=len(genotypes)
     pop_size=len(genotypes[0])

     mean_str=[{BindingStrength(*genotypes[0][index][[i for i in new_bond]]):1} if detailed_distr else BindingStrength(*genotypes[0][index][[i for i in new_bond]])]
     strengths[0][index].remove(new_bond)
     descendents=list(np.where(selections[0]==index)[0])

     gen_index=1
     while gen_index<gen_limit and descendents:
          str_temp=defaultdict(float) if detailed_distr else []
          descendents_temp=[]
          for descendent in descendents:
               if new_bond in strengths[gen_index][descendent]:
                    new_strength=BindingStrength(*genotypes[gen_index][descendent][[i for i in new_bond]])
                    if detailed_distr:
                         str_temp[new_strength]+=1
                    else:
                         str_temp.append(new_strength)
                    descendents_temp.extend(np.where(selections[gen_index]==descendent)[0])
                    strengths[gen_index][descendent].remove(new_bond)
          descendents=descendents_temp
          gen_index+=1

          if detailed_distr:
               mean_str.append(dict(str_temp))
          elif str_temp:
               mean_str.append(np.mean(str_temp))
          
     return mean_str if len(mean_str)>MIN_LEN else None

def qWeak(strength,genotype):
     weak_asym=Counter({K:0 for K in np.linspace(0,1,interface_length+1)})
     weak_sym=Counter({K:0 for K in np.linspace(0,1,interface_length/2+1)})
     g_length=len(genotype)
     used_interfaces=set([item for sublist in strength for item in sublist])
     non_interactings=set(xrange(g_length))-used_interfaces
     for pair in combinations(non_interactings,2):
          weak_asym[BindingStrength(*genotype[[i for i in pair]])]+=1
     for face in non_interactings:
          weak_sym[BindingStrength(*genotype[[face]*2])]+=1
     return weak_asym,weak_sym

""" Phenotype Matching """
def PhenOrderTail(genotypes,phen_table,phenotype_IDs,selections,strengths):    
     gen_limit=len(genotypes)
     pop_size=len(genotypes[0])        

     transition_probabilities=defaultdict(lambda:defaultdict(int))
     tree_found_states=[(0,p,[(1,0)]) for p in xrange(pop_size)]
     first_transitions={(1,1,1):None}

     while tree_found_states:
          generation,index,seen_phens=tree_found_states.pop()
          children=np.where(selections[generation]==index)[0]
          for child in children:
               states_local=deepcopy(seen_phens)
               pid=phenotype_IDs[generation+1][child]
               if pid!=(0,0):
                    if phen_table[pid] not in first_transitions:
                         first_transitions[phen_table[pid]]=phen_table[phenotype_IDs[generation][index]]
                    if pid not in states_local:
                         transition_probabilities[pid][phenotype_IDs[generation][index]]+=1
                         states_local.append(pid)
                    
               if generation<(gen_limit-2):
                    tree_found_states.append((generation+1,child,states_local))

     transition_probabilities_PHENOTYPES={}
     for key in transition_probabilities:
          transition_probabilities_PHENOTYPES[phen_table[key]]={phen_table[phen_ID]:count for phen_ID,count in transition_probabilities[key].iteritems()}
     del first_transitions[(1,1,1)]
     return transition_probabilities_PHENOTYPES,first_transitions

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
     
     
     
def main(argv):
     writeResults(int(argv[1]),*(float(i) for i in argv[2:4]),runs=int(argv[4]),offset=int(argv[5]))
     return
if __name__ == "__main__":
    main(sys.argv)

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

