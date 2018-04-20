#!/usr/bin/env python
import ctypes
import os
import os.path
from sys import platform
import numpy as np



Poly_Lib=ctypes.cdll.LoadLibrary('./AGF.so')

def GetPhenotypeIDs_wrapper(file_path,file_name,ngenes,cols,g_or_i): 
    Poly_Lib.GetPhenotypeIDs.restype=None
    Poly_Lib.GetPhenotypeIDs.argtypes=[ctypes.POINTER(ctypes.c_char),ctypes.POINTER(ctypes.c_char),ctypes.c_uint8,ctypes.c_uint8,ctypes.c_bool]
    c_=ctypes.c_buffer(file_path)
    o_=ctypes.c_buffer(file_name)
    assert os.path.isfile(file_path+file_name), "not a valid file"
    
    Poly_Lib.GetPhenotypeIDs(c_,o_,ngenes,cols,g_or_i)

def ExhaustiveMinimalMethod_wrapper(file_path,ngenes,cols,g_or_i):
    Poly_Lib.ExhaustiveMinimalGenotypes.restype=None
    Poly_Lib.ExhaustiveMinimalGenotypes.argtypes=[ctypes.POINTER(ctypes.c_char),ctypes.c_uint8,ctypes.c_uint8,ctypes.c_bool]
    c_=ctypes.c_buffer(file_path)
    Poly_Lib.ExhaustiveMinimalGenotypes(c_,ngenes,cols,g_or_i)

def SampleMinimalMethod_wrapper(file_path,ngenes,cols,samples,dups,g_or_i):
    #file path is *ONLY* the path of the directory to output the file in, the file name is fixed with suffix _ngenes_colous
    Poly_Lib.SampleMinimalGenotypes.restype=None
    Poly_Lib.SampleMinimalGenotypes.argtypes=[ctypes.POINTER(ctypes.c_char),ctypes.c_uint8,ctypes.c_uint8,ctypes.c_uint32,ctypes.c_bool,ctypes.c_bool]
    c_=ctypes.c_buffer(file_path)
    Poly_Lib.SampleMinimalGenotypes(c_,ngenes,cols,samples,dups,g_or_i)

def GPMap_wrapper(file_path,ngenes,rcols,cols,g_or_i):
    Poly_Lib.GP_MapSampler.restype=None
    Poly_Lib.GP_MapSampler.argtypes=[ctypes.POINTER(ctypes.c_char),ctypes.c_uint8,ctypes.c_uint8,ctypes.c_uint8,ctypes.c_bool]
    c_=ctypes.c_buffer(file_path)
    Poly_Lib.GP_MapSampler(c_,ngenes,rcols,cols,g_or_i)
    


