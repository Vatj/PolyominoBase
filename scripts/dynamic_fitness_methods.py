import numpy as np
import subprocess

def AddNewDataToFile(classifier,data_in,needle,mu):
    with open('/rscratch/asl47/Processed/Dynamic/Evolution_Solution_Times_q2_Mu{}_Needle{}.txt'.format(mu,needle), "a") as data_file:
        data_file.write(classifier+': ')
        for value in data_in:
            data_file.write('{} '.format(value))
        data_file.write('\n')

def AddBulkDataToFile(As,Os,needle,mu):    
    for (A,O) in zip(As,Os):
        data=SolutionCDF(needle,A,mu,O,500)
        AddNewDataToFile('A{}O{}'.format(A,O),data,needle,mu)
        
def LoadExistingData(needle,mu):
    data_dict={key:[int(v) for v in values.split()] for (key,values) in [line.rstrip('\n').split(':') for line in open('/rscratch/asl47/Processed/Dynamic/Evolution_Solution_Times_q2_Mu{}_Needle{}.txt'.format(mu,needle))]}
    return data_dict

def LoadExistingData2(needle,mu):
    data_dict={key:[int(v) for v in values.split()] for (key,values) in [line.rstrip('\n').split(':') for line in open('/rscratch/asl47/Processed/Dynamic/Evolution_Solution_Times3_Mu{}_Needle{}.txt'.format(mu,needle))]}
    return data_dict

def LoadPartials(N=100,mu=16):
    d={}
    for A,O in [(3,50000),(2,5),(2,25),(2,75),(1,5),(1,15),(1,25)]:
        print 'loading for A,O'
        d['A{}O{}'.format(A,O)]=SolutionCDF_Partials(500,A,mu,O,N)
    return d
def SavePartial(d,mu):
        np.savez('/rscratch/asl47/Processed/Dynamic/fitness_increments_Mu{}_K5000'.format(mu),d.values(),args=d.keys(),allow_pickle=True)
def ReadPartial(mu):
    loaded= np.load('/rscratch/asl47/Processed/Dynamic/fitness_increments_Mu{}_K5000.npz'.format(mu),allow_pickle=True)
    data={}
    for i,k in enumerate(loaded['args']):
        data[k]=loaded['arr_0'][i]
    return data
        
def SolutionCDF(needle_length=30,A=2,mu=4,O=3,runs=1):
    occ=0
    firsts=[]



    Mu_Sets={32:'0.001563',16:'0.003125',8:'0.006250',4:'0.012500',1:'0.050000'}
     
    
    needle=np.ones(needle_length,dtype=np.float64)
    
    for r in xrange(runs):
        subfile_name='A{}_T20_C200_N500_Mu{}_O{}_K25000_I0_Run{}'.format(A,Mu_Sets[mu],O,r)
        fitness_import=np.loadtxt('/scratch/asl47/Data_Runs/Dynamic_3/A{}Mu{}O{}/{}_Fitness.txt'.format(A,mu,O,subfile_name))
        #/A{}Mu{}O{}/
        #fitness_import[fitness_import >0]=1
        if max(fitness_import[:,0])>=1:
            haystack=search_sequence_numpy(fitness_import[:,0],needle)
            
            if haystack>0:
                occ+=1
                firsts.append(haystack)

    print "seen {} times for a fraction of {}".format(occ,occ*1./runs)
    return firsts
x=[]
def SolutionCDF_Partials(pop_amount=500,A=2,mu=4,O=3,runs=1):
    firsts=np.zeros((runs,4))
    Mu_Sets={32:'0.001563',16:'0.003125',8:'0.006250',4:'0.012500',1:'0.050000'}
     
    
    for r in xrange(runs):

        subfile_name='A{}_T20_C200_N1000_Mu{}_O{}_K50000_I0_Run{}'.format(A,Mu_Sets[mu],O,r)
        fitness_import=np.loadtxt('/scratch/asl47/Data_Runs/Dynamic_T2/Mu{}/{}_Fitness.txt'.format(mu,subfile_name))
        
        #return fitness_import
        try:
            
        
            triple_target=np.argmax(fitness_import[:,0]>=pop_amount)
            if triple_target==0 and fitness_import.shape[0]!=50000:
                raise Exception


            
            double_targets= (fitness_import[:,1:4]>=pop_amount).argmax(axis=0)
            if np.any(double_targets):
                double_target=double_targets[np.nonzero(double_targets)].min()
            else:
                double_target=triple_target
            
                
            single_targets= (fitness_import[:,4:7]>=pop_amount).argmax(axis=0)
            if np.any(single_targets):
                single_target=single_targets[np.nonzero(single_targets)].min()
            else:
                single_target=double_target

            
            
            if triple_target!=0:
                firsts[r]=np.array([0,single_target,double_target,triple_target])
        
            
        except Exception as e:
            print r,"excepted",fitness_import.shape
        


            
    return firsts[~np.all(firsts == 0, axis=1)]

def DifferentialDiscovery(target_times):
    return np.diff(target_times,axis=1)

def search_sequence_numpy(arr,seq):
    """ Find sequence in an array using NumPy only.

    Parameters
    ----------    
    arr    : input 1D array
    seq    : input 1D array

    Output
    ------    
    Output : 1D Array of indices in the input array that satisfy the 
    matching of input sequence in the input array.
    In case of no match, empty list is returned.
    """

    # Store sizes of input array and sequence
    Na, Nseq = arr.size, seq.size

    # Range of sequence
    r_seq = np.arange(Nseq)

    # Create 2D array of sliding indices across entire length of input array.
    # Match up with the input sequence & get the matching starting indices.
    M = (arr[np.arange(Na-Nseq+1)[:,None] + r_seq] == seq).all(1)
    # Get the range of those indices as final output
    if M.any>0:
        indices=np.where(np.convolve(M,np.ones((Nseq),dtype=int))>0)[0]
        if indices.shape[0]>0:
            return indices[0]
        else:
            return 0
    else:
        return 0
        return [] 

   

def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, 
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')
