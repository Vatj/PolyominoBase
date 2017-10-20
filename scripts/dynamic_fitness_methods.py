import numpy as np

def AddNewDataToFile(classifier,data_in,needle,mu):
    with open('/rscratch/asl47/Processed/Dynamic/Evolution_Solution_Times_Mu{}_Needle{}.txt'.format(mu,needle), "a") as data_file:
        data_file.write(classifier+': ')
        for value in data_in:
            data_file.write('{} '.format(value))
        data_file.write('\n')

def AddBulkDataToFile(Ts,Os,needle,mu):
    keys=['T{}{}'.format(t,'' if t==3 else 'O{}'.format(o)) for (t,o) in zip(Ts,Os)]
    SetRuns(500)
    
    for key in keys:
        SetType(int(key[1]))
        data=SolutionCDF(needle,mu,'' if key[1]=='3' else int(key[3:]))
        AddNewDataToFile(key,data,needle,mu)
        
def LoadExistingData(needle,mu):
    data_dict={key:[int(v) for v in values.split()] for (key,values) in [line.rstrip('\n').split(':') for line in open('/rscratch/asl47/Processed/Dynamic/Evolution_Solution_Times_Mu{}_Needle{}.txt'.format(mu,needle))]}
    return data_dict

    
def SolutionCDF(needle_length=30,mu=4,O=3,slicer=0):
    occ=0
    firsts=[]
    additional=''
    if RUN_TYPE==4 or RUN_TYPE==5:
        additional='_O{}'.format(O)

    Mu_Sets={32:'0.001563',16:'0.003125',8:'0.006250',4:'0.012500',1:'0.050000'}
    #Mu_Files={4:'',2:'Mu2',1:'Mu1'}
    
    
    needle=np.array([1]*needle_length,dtype=np.float64)
    for r in xrange(RUNS):
        subfile_name='Modular{}_T20_C200_N500_Mu{}{}_K15000_Run{}'.format(RUN_TYPE,Mu_Sets[mu],additional,r)
        fitness_import=np.genfromtxt('/scratch/asl47/Data_Runs/Dynamic_2/T{}Mu{}{}/{}_Fitness.txt'.format(RUN_TYPE,mu,addition))
        if max(fitness_import[:,slicer])>=1:
            haystack=search_sequence_numpy(fitness_import[:,slicer],needle)
            if haystack>0:
                occ+=1
                firsts.append(haystack)

    print "seen {} times for a fraction of {}".format(occ,occ*1./RUNS)
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
