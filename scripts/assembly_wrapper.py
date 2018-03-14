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

def GenMethod(ngenes,cols,samples):
    Poly_Lib.SampleMinimalGenotypes.restype=None
    Poly_Lib.SampleMinimalGenotypes.argtypes=[ctypes.c_uint8,ctypes.c_uint8,ctypes.c_uint32,ctypes.c_bool]
    print "here"
    Poly_Lib.GGenerator(ngenes,cols,samples,False)

    


