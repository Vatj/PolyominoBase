import numpy as np

#BASE_FILE_PATH='/scratch/asl47/Data_Runs/Interface_Cron/{0}_{1}_T{2:.6f}_Mu{3:.6f}_Gamma{4:.6f}_Run{5}.txt'
#BASE_FILE_PATH='/rscratch/asl47/Bulk_Run/Interfaces/{0}_{1}_T{2:.6f}_Mu{3:.6f}_Gamma{4:.6f}_Run{5}.txt'
BASE_FILE_PATH='../{0}_{1}_T{2:.6f}_Mu{3:.6f}_Gamma{4:.6f}_Run{5}.txt'


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
     
     return np.array(phenotype_IDs,dtype=np.uint8),np.array(selections,np.uint16)

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

def LoadT(mu=1,t=0.35):
     g=LoadGenotypeHistory(2,mu=mu,temperature=t)
     st=LoadStrengthHistory(mu=mu,temperature=t)
     p,s=LoadEvolutionHistory(mu=mu,temperature=t)
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
