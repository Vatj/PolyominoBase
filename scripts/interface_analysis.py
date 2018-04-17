import numpy as np
from scipy import stats

#BASE_FILE_PATH='/scratch/asl47/Data_Runs/Interface_Cron/{0}_{1}_T{2:.6f}_Mu{3:.6f}_Gamma{4:.6f}_Run{5}.txt'
#BASE_FILE_PATH='/rscratch/asl47/Bulk_Run/Interfaces/{0}_{1}_T{2:.6f}_Mu{3:.6f}_Gamma{4:.6f}_Run{5}.txt'
BASE_FILE_PATH='../output/{0}_{1}_T{2:.6f}_Mu{3:.6f}_Gamma{4:.6f}_Run{5}.txt'


def LoadEvolutionHistory(temperature=0.000001,mu=1,gamma=1,run=0):

     phen_line=True
     phenotype_IDs=[]
     selections=[]
     for line in open(BASE_FILE_PATH.format('PhenotypeHistory','S',temperature,mu,gamma,run)):
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
     for line in open(BASE_FILE_PATH.format('GenotypeHistory','S',temperature,mu,gamma,run)):
          converted=[int(i) for i in line.split()]
          genotypes.append([[int(i) for i in j] for j in [converted[i:i + 4*n_tiles] for i in xrange(0, len(converted), 4*n_tiles)]])
          #print [[int(i) for i in j] for j in [converted[i:i + 4*n_tiles] for i in xrange(0, len(converted), 4*n_tiles)]]
          #break
     return np.array(genotypes,dtype=np.uint32)

def LoadStrengthHistory(temperature=0.000001,mu=1,gamma=1,run=0):
    strengths=[]
    for line in open(BASE_FILE_PATH.format('Strengths','S',temperature,mu,gamma,run)):
         strengths.append([[tuple(np.uint8(i) for i in j.split()) for j in tmp.split(',')[:-1]] for tmp in line.split('.')[:-1]])


    return strengths

def LoadT(mu=1,t=0.35,run=0):
     g=LoadGenotypeHistory(3,mu=mu,temperature=t,run=run)
     st=LoadStrengthHistory(mu=mu,temperature=t,run=run)
     p,s=LoadEvolutionHistory(mu=mu,temperature=t,run=run)
     return (g,s,p,st)


from random import randint
from collections import defaultdict
def RandomHistorySampling(genotypes,selections,phenotypes,strengths,goback=10):
     printer=False
     rg=randint(genotypes.shape[0]/2,genotypes.shape[0]-1)
     rp=randint(0,genotypes.shape[1]-1)

     #rg=150


     strength_tracker=defaultdict(list)

     if(goback>rg):
          rg=goback
     assert rg>=goback, "going back too far"
     if printer:
          print "random sampling from g: ",rg," and p: ",rp
          print "started from ",genotypes[rg][rp],phenotypes[rg][rp]
          print "staring strengths", strengths[rg][rp]
          
     for stren in strengths[rg][rp]:
          strength_tracker[stren].append([BindingStrength(*[genotypes[rg][rp][j] for j in stren])])

     

     for bg in xrange(1,goback+1):
          if printer:
               print "selected from: ",selections[rg-bg][rp]
         
          rp=selections[rg-bg][rp]
          if printer:
               print "now ",genotypes[rg-bg][rp],phenotypes[rg-bg][rp]
               print "new inters", strengths[rg-bg][rp]
               print "strs ",[BindingStrength(*[genotypes[rg-bg][rp][j] for j in i]) for i in strengths[rg-bg][rp]]

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



        
interface_length=32
def convint(x):
     return np.uint32(x)

def SeqDiff(genotype1,genotype2):
     return [i for i in xrange(len(genotype1)) if genotype1[i] != genotype2[i]]

def BindingStrength(base1,base2):
     
     #assert type(base1)==np.uint8 and type(base2)==np.uint8 , "wrong type"
     return 1-bin(np.bitwise_xor(convint(base1),revbits(base2))).count("1")/float(interface_length)

def revbits(x):
     return int(bin(~convint(x))[2:].zfill(interface_length)[::-1], 2)



def SampleStrengths():
     return 2

import matplotlib.pyplot as plt
import icy;icy.Use_Seaborn()
import seaborn as sns
def RandomWalk(I_size=32,n_steps=1000,phi=0.5,scale_factor=-1):
     
     T_star=0.75
     s_hats=np.linspace(0,1,I_size+1)
     msize=12

     N=int(I_size*(1-T_star))+1
     if scale_factor<0:
          scale_factor=N
     states=np.zeros(N+1)
     #states[1:]=1./N
     states[1]=1
     #states[2]=.5
     #states[1]=.5
     #states[46]=.5
     #states[-1]=1
     #print states
     #plt.plot(range(N+1),states,'b--')

     fig, ax1 = plt.subplots()
     
     for i in xrange(n_steps):
          states=UpdateStates(states,phi)
          
     #plt.plot(range(0,N*8+1,8),states/8.,'rh')
     scaf=.9/states[-3]
     #print states
     ax1.plot(np.linspace(T_star,1,N),scaf*states[1:]/(scale_factor*1./N),marker='^',ls='',ms=14,mfc='orangered')
     
     ff=2-states[1]/(1-states[1]*phi/2.)
     
     #print "ff",ff
     for i in xrange(2,N):
          ax1.scatter(T_star+(i-1)*1./I_size,scaf*(states[i-1]*ff-states[i-2]),c='k',marker='x',s=500)

     ax1.scatter(1,scaf*states[-2]/ff,c='k',marker='x',s=500)
     
     #ax2=ax1.twinx()
     
     ax1.plot(s_hats[:int(T_star*I_size)],[0 for i in xrange(int(T_star*I_size))],c='royalblue',ls='--',lw=0.75,marker='o',ms=msize)
     ax1.plot(s_hats[int(T_star*I_size):],[1 for i in xrange(int(T_star*I_size),I_size+1)],c='royalblue',ls='--',lw=0.75,marker='o',ms=msize)

     ax1.plot([s_hats[int(T_star*I_size)-1]]+[s_hats[int(T_star*I_size)]]*2,[0,0,1],c='royalblue',mec='royalblue',mfc='None',mew=1,ls='--',lw=0.75,marker='o',ms=msize)

     plt.xlabel(r'$\hat{S}$')
     ax1.text(0.4,0.05,r'$Pr_{\textrm{binding}}$')
     ax1.text(0.8,0.1,r'genetic drift')
     ax1.text(.35,.5,r'$P_i^{t+1}=(1-\mu)P_i^t$''\n'r'$+\frac{\mu}{2}(P_{i-1}^t+P_{i+1}^t)$''\n'r'$ + \frac{(\frac{\mu}{2} P_1^t)P_i^t}{1-\frac{\mu}{2} P_1^t}$')
     ax1.axes.get_yaxis().set_visible(False)
     sns.despine(left=1,top=1,right=1)
     ax1.set_xlim((.3,1.05))
     plt.show(block=False)
     #return states
     
from copy import deepcopy
def UpdateStates(states,phi=0.25):
     states_updating=deepcopy(states)
     for i in xrange(1,states.shape[0]):
          if i!=states.shape[0]-1:
               states_updating[i+1]+=states[i]*phi/2.
               states_updating[i-1]+=states[i]*phi/2.
               states_updating[i]-=states[i]*phi
          else:
               states_updating[i-1]+=states[i]*phi
               states_updating[i]-=states[i]*phi
          
          
     
     #print states
     #print "redistr",states_updating
     birth_factor=states_updating[0]
     survive_factor=np.sum(states_updating[1:])
     if(birth_factor!=0):
          states_updating+= birth_factor*states_updating/survive_factor
          states_updating[0]=0
     #print "reborn ",states_updating, sum(states_updating)
     return states_updating


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

from colorsys import hsv_to_rgb
from random import uniform,choice

def pp(mae,mai,ms,called=False):
     if not called:
          plt.figure()
     ups=0
     downs=0
     
     for interface_type,color_range in zip([mae,mai,ms],[(.25,.38),(0.58,0.75),(0.91,.08)]):
          for data in interface_type:
          
               slope, intercept, r_value, p_value, std_err = stats.linregress(xrange(len(data)-1),data[1:])
               if slope>0:
                    ups+=1
               else:
                    downs+=1
          
               if color_range[1]>color_range[0]:
                    h = uniform(*color_range)
               else:
                    h=choice([uniform(color_range[0],1),uniform(0,color_range[1])])
               s = uniform(0.2, 1)
               v = uniform(0.5, 1)
                                 
               r, g, b = hsv_to_rgb(h, s, v)
               plt.plot(xrange(data[0],data[0]+len(data)-1),data[1:],c=(r,g,b),alpha=0.2,zorder=1)
          
               plt.plot(xrange(data[0],data[0]+len(data)-1),[slope*x+intercept for x in xrange(len(data)-1)],'--',c=(r,g,b),zorder=10)
          
          
     if not called:
          plt.show(block=False)
     else:
          return ups,downs
     
import pickle
def writeResults(I,M,t,runs):
     datastruct={}
     for run in xrange(runs):
          datastruct[run]=AnalysePhylogeneticStrengths(M,t,run)
          
     with open('I{}Mu{}T{}.pkl'.format(I,M,t), 'wb') as f:
          pickle.dump(datastruct, f)

def AnalysePhylogeneticStrengths(mu,t,r):
     g,s,p,st=LoadT(mu=mu,t=t,run=r)
     mae,mai,ms=qBFS(g,s,st)
     
     return {'AsymE':mae,'AsymI':mai,'Sym':ms}
     
     
     
def qBFS(genotypes,selections,strengths):
     gen_limit=len(genotypes)
     pop_size=len(genotypes[0])
     
     MSE_ae=[]
     MSE_ai=[]
     MSE_s=[]
     
     first_bond=False
     for generation in xrange(gen_limit-1):
          for index in xrange(pop_size):
               if generation==0:
                    for new_bond in strengths[generation][index]:
                         res=[generation]+qDFS(generation,index,new_bond,genotypes,selections,strengths)
                         if new_bond[0]==new_bond[1]:
                              MSE_s.append(res)
                         else:
                              if new_bond[0]/4==new_bond[1]/4:
                                   MSE_ai.append(res)
                              else:
                                   MSE_ae.append(res)
                           
               else:
                    diff=set(strengths[generation][index])-set(strengths[generation-1][selections[generation][index]])
                    if len(diff)>0:
                         for new_bond in diff:
                              res=qDFS(generation,index,new_bond,genotypes,selections,strengths)
                              if res is not None:
                                   res.insert(0,generation)
                                   
                                   if new_bond[0]==new_bond[1]:
                                        MSE_s.append(res)
                                   else:
                                        if new_bond[0]/4==new_bond[1]/4:
                                             MSE_ai.append(res)
                                        else:
                                             MSE_ae.append(res)
                              

     return filter(None,MSE_ae),filter(None,MSE_ai),filter(None,MSE_s)


def qDFS(gen0,index,new_bond,genotypes,selections,strengths):
     #ANALYSIS PARAMETERS
     MIN_LEN=50
     #RUNTIME PARAMETERS
     gen_limit=len(genotypes)
     pop_size=len(genotypes[0])

     mean_str=[BindingStrength(*genotypes[gen0][index][[i for i in new_bond]])]
     strengths[gen0][index].remove(new_bond)
     descendents=np.where(selections[gen0]==index)[0]

     gen_index=gen0+1
     while True:
          str_temp=[]
          descendents_temp=[]
          for descendent in descendents:
               if new_bond in strengths[gen_index][descendent]:
                    str_temp.append(BindingStrength(*genotypes[gen_index][descendent][[i for i in new_bond]]))
                    descendents_temp.extend(np.where(selections[gen_index]==descendent)[0])
                    strengths[gen_index][descendent].remove(new_bond)
               
          descendents=descendents_temp
          gen_index+=1
          if len(str_temp)==0 or gen_index==gen_limit:
               break
          mean_str.append(np.mean(str_temp))
          
     return mean_str if len(mean_str)>MIN_LEN else None
