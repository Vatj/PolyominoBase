#!/usr/bin/env python
import ctypes
import os
import os.path
from sys import platform
import numpy as np



Poly_Lib=ctypes.cdll.LoadLibrary('./AGF.so')
#np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='C_CONTIGUOUS')

def NewMethod(fileName,out_path='/',file_of_genotypes=True,cols=0,ngenes=0):
    Poly_Lib.WrappedGetPhenotypesID.restype=None
    Poly_Lib.WrappedGetPhenotypesID.argtypes=[ctypes.POINTER(ctypes.c_char),ctypes.POINTER(ctypes.c_char),ctypes.c_bool,ctypes.c_int,ctypes.c_int]
    c_=ctypes.c_buffer(fileName+".txt")
    o_=ctypes.c_buffer(out_path)
    assert os.path.isfile(fileName+".txt"), "not a valid file"
    
    Poly_Lib.WrappedGetPhenotypesID(c_,o_,file_of_genotypes,cols,ngenes)

def WriteMethod(fileName,file_of_genotypes=True,cols=0,ngenes=0):
    Poly_Lib.GGenerator.restype=None
    Poly_Lib.GGenerator.argtypes=[ctypes.POINTER(ctypes.c_char),ctypes.c_bool,ctypes.c_uint8,ctypes.c_uint8]
    c_=ctypes.c_buffer(fileName+".txt")

    Poly_Lib.GGenerator(c_,file_of_genotypes,cols,ngenes)
    

def G2I(genotype,colours,n_genes):
    Poly_Lib.genotype_to_index.restype=ctypes.c_uint64
    Poly_Lib.genotype_to_index.argtypes=[ctypes.POINTER(ctypes.c_uint8),ctypes.c_uint8,ctypes.c_uint8]
    g_P=(ctypes.c_uint8*len(genotype))(*genotype)
    print g_P
    return Poly_Lib.genotype_to_index(g_P,colours,n_genes)

def I2G(index,colours,n_genes):
    Poly_Lib.index_to_genotype.restype=None
    Poly_Lib.index_to_genotype.argtypes=[ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int]
    g_P=(ctypes.c_int*(n_genes*4))()
    Poly_Lib.index_to_genotype(index,g_P,colours,n_genes)
    return list(g_P)





def OldMethod():
    Poly_Lib.WrappedGetPhenotypeID.restype=None
    Poly_Lib.WrappedGetPhenotypeID.argtypes=[ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.POINTER(ctypes.c_int)]

    genotype=[1,1,1,1,2,0,0,0]
    Phen_list=[[1,1,1], [3,1,1,1,1], [2,1,1,1], [2,2,1,1,1,0]]
    gen=[1,1,1,1, 2,0,0,0, 3,4,0,0, 5,0,0,0, 6,0,0,0]

def Do(genotype,phens):
    g_P=(ctypes.c_int*len(genotype))(*genotype)
    buf_len=(len(genotype)**2)/4+2
    flat_phens=[item for sublist in phens for item in sublist]
    shapes=flat_phens+[0]*buf_len
    s_P=(ctypes.c_int*len(shapes))(*shapes)
    s_s=ctypes.c_int(len(s_P))
    n_s=ctypes.c_int(len(phens)-1)
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
        
import tarfile    
def LoadGenotypes(part):
    with tarfile.open("/scratch/asl47/GenL{}".format(part), "r:gz") as tar:
        for member in tar.getmembers():
            if  member.isfile():
                filex=tar.extractfile(member)
                return np.fromstring(filex.read(), dtype=NP.uint8)
                


def CompilePhenotypes(parts):

    phen_list=[]
    phen_xfer=[dict() for _ in xrange(len(parts)-1)]
    phen_list = [[int(i) for i in line.rstrip('\n').split()] for line in open('/scratch/asl47/Phenotype_Split_{}.txt'.format(parts[0]))]
    
    for p,part in enumerate(parts[1:]):
        phens = [[int(i) for i in line.rstrip('\n').split()] for line in open('/scratch/asl47/Phenotype_Split_{}.txt'.format(part))]
        for ind,phen in enumerate(phens):
            if phen in phen_list:
                phen_xfer[p][phen_list.index(phen)]=ind

    return phen_xfer



    

    


