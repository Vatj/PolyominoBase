#!/usr/bin/env python
import ctypes
import os
import os.path
from sys import platform
import numpy as np



Poly_Lib=ctypes.cdll.LoadLibrary('./AGF.so')

def NewMethod(fileName,out_path='/',file_of_genotypes=True,cols=0,ngenes=0):
    Poly_Lib.WrappedGetPhenotypesID.restype=None
    Poly_Lib.WrappedGetPhenotypesID.argtypes=[ctypes.POINTER(ctypes.c_char),ctypes.POINTER(ctypes.c_char),ctypes.c_bool,ctypes.c_uint8,ctypes.c_uint8]
    c_=ctypes.c_buffer(fileName+".txt")
    o_=ctypes.c_buffer(out_path)
    assert os.path.isfile(fileName+".txt"), "not a valid file"
    
    Poly_Lib.WrappedGetPhenotypesID(c_,o_,file_of_genotypes,cols,ngenes)

def WriteMethod(fileName,file_of_genotypes=True,cols=0,ngenes=0):
    Poly_Lib.GGenerator.restype=None
    Poly_Lib.GGenerator.argtypes=[ctypes.POINTER(ctypes.c_char),ctypes.c_bool,ctypes.c_uint8,ctypes.c_uint8]
    c_=ctypes.c_buffer(fileName+".txt")

    Poly_Lib.GGenerator(c_,file_of_genotypes,cols,ngenes)

def SampleMinimalMethod_wrapper(file_path,ngenes,cols,samples,dups,g_or_i):
    Poly_Lib.SampleMinimalGenotypes.restype=None
    Poly_Lib.SampleMinimalGenotypes.argtypes=[ctypes.POINTER(ctypes.c_char),ctypes.c_uint8,ctypes.c_uint8,ctypes.c_uint32,ctypes.c_bool,ctypes.c_bool]
    c_=ctypes.c_buffer(file_path)
    Poly_Lib.SampleMinimalGenotypes(c_,ngenes,cols,samples,dups,g_or_i)

def GPMap_wrapper(file_path,ngenes,rcols,cols,g_or_i):
    Poly_Lib.GP_MapSampler.restype=None
    Poly_Lib.GP_MapSampler.argtypes=[ctypes.POINTER(ctypes.c_char),ctypes.c_uint8,ctypes.c_uint8,ctypes.c_uint8,ctypes.c_bool]
    c_=ctypes.c_buffer(file_path)
    Poly_Lib.GP_MapSampler(c_,ngenes,rcols,cols,g_or_i)
    


