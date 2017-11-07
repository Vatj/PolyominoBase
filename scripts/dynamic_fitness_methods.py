import numpy as np
import subprocess

def AddNewDataToFile(classifier,data_in,needle,mu):
    with open('/rscratch/asl47/Processed/Dynamic/Evolution_Solution_Times_new__Mu{}_Needle{}.txt'.format(mu,needle), "a") as data_file:
        data_file.write(classifier+': ')
        for value in data_in:
            data_file.write('{} '.format(value))
        data_file.write('\n')

def AddBulkDataToFile(As,Os,needle,mu):    
    for (A,O) in zip(As,Os):
        data=SolutionCDF(needle,A,mu,O,275)
        AddNewDataToFile('A{}O{}'.format(A,O),data,needle,mu)
        
def LoadExistingData(needle,mu):
    data_dict={key:[int(v) for v in values.split()] for (key,values) in [line.rstrip('\n').split(':') for line in open('/rscratch/asl47/Processed/Dynamic/Evolution_Solution_Times_new_Mu{}_Needle{}.txt'.format(mu,needle))]}
    return data_dict

def LoadExistingData2(needle,mu):
    data_dict={key:[int(v) for v in values.split()] for (key,values) in [line.rstrip('\n').split(':') for line in open('/rscratch/asl47/Processed/Dynamic/Evolution_Solution_Times3_Mu{}_Needle{}.txt'.format(mu,needle))]}
    return data_dict

    
def SolutionCDF(needle_length=30,A=2,mu=4,O=3,runs=1):
    occ=0
    firsts=[]



    Mu_Sets={32:'0.001563',16:'0.003125',8:'0.006250',4:'0.012500',1:'0.050000'}
     
    
    needle=np.ones(needle_length,dtype=np.float64)
    
    for r in xrange(runs):
        subfile_name='A{}_T20_C200_N1000_Mu{}_O{}_K15000_I0_Run{}'.format(A,Mu_Sets[mu],O,r)
        fitness_import=np.loadtxt('/rscratch/asl47/Bulk_Run/Modular/{}_Fitness.txt'.format(subfile_name))
        #/A{}Mu{}O{}/
        fitness_import[fitness_import >0]=1
        if max(fitness_import[:,1])>=1:
            haystack=search_sequence_numpy(fitness_import[:,1],needle)
            
            if haystack>0:
                occ+=1
                firsts.append(haystack)

    print "seen {} times for a fraction of {}".format(occ,occ*1./runs)
    return firsts

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
