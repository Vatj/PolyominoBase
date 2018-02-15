#!/usr/bin/env python
import ctypes
import os
from sys import platform
import numpy as np



Poly_Lib=ctypes.cdll.LoadLibrary('./AGF.so')
#np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='C_CONTIGUOUS')

Poly_Lib.WrappedGetPhenotypeID.restype=None
Poly_Lib.WrappedGetPhenotypeID.argtypes=[ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.POINTER(ctypes.c_int)]

genotype=[1,1,1,1,2,0,0,0]


def Do(genotype,phens):
    g_P=(ctypes.c_int*len(genotype))(*genotype)
    buf_len=(len(genotype)**2)/4+2
    flat_phens=[item for sublist in phens for item in sublist]
    shapes=flat_phens+[0]*buf_len
    s_P=(ctypes.c_int*len(shapes))(*shapes)
    s_s=ctypes.c_int(len(s_P))
    n_s=ctypes.c_int(len(phens))
    ID_p=(ctypes.c_int*int(len(genotype)/4))()
    Poly_Lib.WrappedGetPhenotypeID(len(genotype),g_P,s_s,s_P,n_s,ID_p)
    IDs=[i for i in ID_p if i!=-2]
    
    new_phen=[i for i in s_P if i!=-1][len(flat_phens):]
    if new_phen:
        sub_phens=[]
        n_ind=0
        while n_ind <len(new_phen):
            sub_phens.append(new_phen[n_ind:n_ind+new_phen[n_ind]*new_phen[n_ind+1]+2])
            n_ind+=new_phen[n_ind]*new_phen[n_ind+1]+2
        phens.extend(sub_phens)
    return IDs,phens
        
    
    

    

    

Phen_list=[[1,1,1], [3,1,1,1,1], [2,1,1,1], [2,2,1,1,1,0]]
gen=[1,1,1,1, 2,0,0,0, 3,4,0,0, 5,0,0,0, 6,0,0,0]
