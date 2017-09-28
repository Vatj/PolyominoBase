#IMPORTS
import seaborn as sns
import matplotlib.pyplot as plt
import re
import numpy as np

#SETTINGS
def Use_Seaborn():
    sns.set_context("paper",font_scale=2.2)
    sns.set_style("white",rc={"xtick.major.size": 8, "ytick.major.size": 8,"xtick.minor.size":5, "ytick.minor.size": 5,"axes.linewidth": 2,"axes.edgecolor":"darkgray","font.size":8,"axes.titlesize":8,"axes.labelsize":5})

#PLOTTING METHODS

def Plot_Brute_Force_vs_Graph_Comparison():
    classification_List_Raw=[line.rstrip('\n') for line in open('/rscratch/asl47/Method_Comparisons_Sample_Run1.txt')]
    classification_List_Ks=[line for line in classification_List_Raw if line[0]=='K']
    classification_List_Formatted=[filter(None,re.split('[a-z]|[A-Z]|:|_|>|<|=',line)) for line in classification_List_Ks]
    classification_List_Integers=[[int(element) for element in line] for line in classification_List_Formatted]
    classification_List_Numpy=np.array(classification_List_Integers)

    plt.plot(classification_List_Numpy.T[0],1-classification_List_Numpy[:,2]*1./classification_List_Numpy[:,1])
    #plt.yscale('log')
    #plt.xscale('log')
    plt.ylim([1-1.05*max(classification_List_Numpy[:,2]*1./classification_List_Numpy[:,1]),1])

    plt.xlabel('Maximum Repeated Checks')
    plt.ylabel('Underclassification Error')
    
    
    plt.show()
    
    return classification_List_Integers

    
    

########################
####GRAPH METHOD PRO####
########################

def Speed_Plot():
    Use_Seaborn()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')


    #plt.plot([3,6,9,12],[278.4
    upp=101
    plt.plot(xrange(1,upp),[0.875**x for x in xrange(1,upp)])
    plt.plot(xrange(1,upp),[0.01]*(upp-1))
    plt.yscale('log')
    plt.show(block=False)
