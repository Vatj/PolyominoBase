import icy
import matplotlib.pyplot as plt
import numpy as np

icy.Use_Seaborn()

def PlotEvolution():
    T=4
    C=40
    N=500
    K=5000
    M=13
    R=0
    I=1
    D=20
    muLs= [float(line.rstrip('\n')) for line in open('/rscratch/asl47/Bulk_Run/Configs/MuL_Values.txt')]
  
    discovery_data=np.empty([len(muLs),D])
    adaptation_data=np.empty([len(muLs),D])
    mu_count=-1
    with open('/scratch/asl47/Data_Runs/Regulation/Evolution_T{}_C{}_N{}_K{}_M{}_R{}_I{}.txt'.format(T,C,N,K,M,R,I),'r') as input_file:
        for line in input_file:
            if 'muL' in line:
                mu_count+=1
            elif 'D' in line:
                discovery_data[mu_count,:]=[int(i) for i in line.split()[1:]]
            elif 'A' in line:
                adaptation_data[mu_count,:]=[int(i) for i in line.split()[1:]]
            
    plt.plot(muLs,np.mean(adaptation_data,axis=1))
    plt.plot(muLs,np.mean(discovery_data,axis=1),'r--')
    plt.show(block=False)
