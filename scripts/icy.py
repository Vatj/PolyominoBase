import matplotlib.pyplot as plt
import seaborn as sns

def Use_Seaborn():
    sns.set_context("paper",font_scale=2.2)
    sns.set_style("white",rc={"xtick.major.size": 8, "ytick.major.size": 8,"xtick.minor.size":5, "ytick.minor.size": 5,"axes.linewidth": 2,"axes.edgecolor":"darkgray","font.size":8,"axes.titlesize":8,"axes.labelsize":5})
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')


import glob
import numpy as np
import os
#BASE_FILE_PATH='/scratch/asl47/Data_Runs/Interface/T{2:.6f}/{0}_{1}_T{2:.6f}_Mu{3:.6f}_Gamma{4:.6f}_Run{5}.txt'
BASE_FILE_PATH='/rscratch/asl47/Bulk_Run/Interfaces/{0}_{1}_T{2:.6f}_Mu{3:.6f}_Gamma{4:.6f}_Run{5}.txt'

    
def DeleteEmptyFiles(r_type,temperature,mu,gamma):
       runs=sorted([int(s[s.rindex('Run')+3:-4]) for s in glob.glob(BASE_FILE_PATH.format('Sizes',r_type,temperature,mu,gamma,'*'))])
       for run in runs:
           if os.path.getsize(BASE_FILE_PATH.format('Sizes',r_type,temperature,mu,gamma,run))==0:
               for f in ['Sizes','Phenotypes','Strengths','Fitness']:
                   os.remove(BASE_FILE_PATH.format(f,r_type,temperature,mu,gamma,run))
               
def AlignFiles(r_type,temperature,mu,gamma):
    
    runs=sorted([int(s[s.rindex('Run')+3:-4]) for s in glob.glob(BASE_FILE_PATH.format('Sizes',r_type,temperature,mu,gamma,'*'))])
    print len(runs)
    run_splits= np.split(runs, np.where(np.diff(runs) > 1)[0]+1)
    print run_splits
    #return run_splits
    file_deltas=np.cumsum([run_splits[i+1][0]-run_splits[i][-1]-1 for i in xrange(len(run_splits)-1)])
    for splits,delta in zip(run_splits[1:],file_deltas):
        for run in splits:
            #print BASE_FILE_PATH.format('Sizes',r_type,temperature,mu,gamma,str(run-delta))
            #print "ren",BASE_FILE_PATH.format(f_type,r_type,temperature,mu,gamma,str(run)),BASE_FILE_PATH.format(f_type,r_type,temperature,mu,gamma,str(run-delta))
            for f in ['Sizes','Phenotypes','Strengths','Fitness']:
                os.rename(BASE_FILE_PATH.format(f,r_type,temperature,mu,gamma,str(run)),BASE_FILE_PATH.format(f,r_type,temperature,mu,gamma,str(run-delta)))
