def Produce_Log_Space_Configs(lower,upper,nums):
    import numpy as np
    lowerLog=np.log10(lower)
    upperLog=np.log10(upper)
    if type(nums)==list:
        full=np.array([])
        delta=(upperLog-lowerLog)*1./len(nums)
        for i,n in enumerate(nums):
            if i!=len(nums)-1:
                logs=np.logspace(lowerLog+i*delta,lowerLog+(i+1)*delta,n)[:-1]
            else:
                logs=np.logspace(lowerLog+i*delta,lowerLog+(i+1)*delta,n)
            full=np.append(full,logs)
        with open('Logspace_Config_{}_{}_2.txt'.format(len(nums),sum(nums)), 'wb') as the_file:
            np.savetxt(the_file,full,fmt='%.8f')
    elif type(nums)==int:
        log_spaced=np.logspace(np.log10(lower),np.log10(upper),nums)
        with open('Logspace_Config_{}.txt'.format(nums), 'wb') as the_file:
            np.savetxt(the_file,log_spaced,fmt='%.8f')
        
import numpy as np
def Partial_Log_Space(lower,upper,nums,nparray=np.array([])):
    lowerLog=np.log10(lower)
    upperLog=np.log10(upper)
    pls=np.logspace(lowerLog,upperLog,nums)
    if nparray.shape[0]==0:
        nparray= np.append(nparray,pls)
    else:
        nparray= np.append(nparray,pls[1:])
    print(pls)
    return nparray
def Write_Log_File(fName,nparray):
    with open('Logspace_Config_{}.txt'.format(fName), 'w') as the_file:
        np.savetxt(the_file,nparray,fmt='%.8f')
