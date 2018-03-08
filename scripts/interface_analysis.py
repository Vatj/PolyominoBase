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
    print len(genotypes)
    return np.array(genotypes,dtype=np.uint8)

from random import randint
def RandomHistorySampling(genotypes,selections,goback=20):
    rg=randint(genotypes.shape[0]/2,genotypes.shape[0])
    rp=randint(0,genotypes.shape[1])
    rg=20
    rp=7
    assert rg>=goback, "going back too far"
    print "random sampling from g: ",rg," and p: ",rp
    print genotypes[rg][rp]
    
    for bg in xrange(1,goback+1):
        print "selected from: ",selections[rg-bg][rp]
        rp=selections[rg-bg][rp]
        print genotypes[rg-bg][rp]

        
        
def SeqDiff(genotype1,genotype2):
    return [i for i in xrange(len(genotype1)) if genotype1[i] != genotype2[i]]

def SammingDistance(base1,base2):
    assert type(base1)==np.uint8 and type(base2)==np.uint8 , "wrong type"
    return bin(np.bitwise_xor(base1,revbits(base2))).count("1")

def revbits(x):
    return int(bin(~np.uint8(x))[2:].zfill(8)[::-1], 2)

def SampleStrengths():
    return 2
