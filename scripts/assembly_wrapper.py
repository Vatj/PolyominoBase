#!/usr/bin/env python
import ctypes
import os
from sys import platform
import numpy as np



Poly_Lib=ctypes.cdll.LoadLibrary('./AGF.so')
#np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='C_CONTIGUOUS')

Poly_Lib.WrappedGetPhenotypeID.restype=ctypes.c_int
Poly_Lib.WrappedGetPhenotypeID.argtypes=[ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.c_int]

genotype=[1,1,1,1,2,0,0,0]


def Do(genotype,phens):
    g_P=(ctypes.c_int*len(genotype))(*genotype)
    buf_len=(len(genotype)**2)/4+2
    flat_phens=[item for sublist in phens for item in sublist]
    shapes=flat_phens+[0]*buf_len
    s_P=(ctypes.c_int*len(shapes))(*shapes)
    s_s=ctypes.c_int(len(s_P))
    n_s=ctypes.c_int(len(phens))
    ID= int(Poly_Lib.WrappedGetPhenotypeID(len(genotype),g_P,s_s,s_P,n_s))
    new_phen=[i for i in s_P if i!=-1]

    phens.append(new_phen[len(flat_phens):])
    return ID,phens
        
    
    

    

    

Phen_list=[[1,1,1], [3,1,1,1,1], [2,1,1,1], [2,2,1,1,1,0]]
gen=[1,1,1,1, 2,0,0,0]
