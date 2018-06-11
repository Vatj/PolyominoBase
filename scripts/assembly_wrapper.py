#!/usr/bin/env python
import ctypes
import os
import os.path
from sys import platform
import numpy as np

from graph_methods import Trim_Topologies

Poly_Lib = ctypes.cdll.LoadLibrary('/rscratch/vatj2/Polyominoes/PolyominoBase/scripts/AGF.so')


def GetPhenotypesIDs_wrapper(file_path, file_name, ngenes, colours):
    Poly_Lib.GetPhenotypesIDs.restype = None
    Poly_Lib.GetPhenotypesIDs.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.POINTER(ctypes.c_char), ctypes.c_uint8, ctypes.c_uint8]
    c_ = ctypes.c_buffer(file_path)
    o_ = ctypes.c_buffer(file_name)

    Poly_Lib.GetPhenotypesIDs(c_, o_, ngenes, colours)

def PreProcessWrite_wrapper(file_path, file_name, ngenes, colours):
    Poly_Lib.PreProcessWrite.restype = None
    Poly_Lib.PreProcessWrite.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.POINTER(ctypes.c_char), ctypes.c_uint8, ctypes.c_uint8]
    c_ = ctypes.c_buffer(file_path)
    o_ = ctypes.c_buffer(file_name)

    Poly_Lib.PreProcessWrite(c_, o_, ngenes, colours)


def ExhaustiveMinimalMethod_wrapper(file_path, ngenes, colours):
    Poly_Lib.ExhaustiveMinimalGenotypes.restype = None
    Poly_Lib.ExhaustiveMinimalGenotypes.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_uint8, ctypes.c_uint8]
    c_ = ctypes.c_buffer(file_path)
    Poly_Lib.ExhaustiveMinimalGenotypes(c_, ngenes, colours)


def SampleMinimalMethod_wrapper(file_path, ngenes, colours, samples, dups):
    Poly_Lib.SampleMinimalGenotypes.restype = None
    Poly_Lib.SampleMinimalGenotypes.argtypes = [ctypes.POINTER(ctypes.c_char),ctypes.c_uint8,ctypes.c_uint8,ctypes.c_uint32,ctypes.c_bool]
    c_ = ctypes.c_buffer(file_path)
    Poly_Lib.SampleMinimalGenotypes(c_, ngenes, colours, samples, dups)


def GPMap_wrapper(file_path, ngenes, rcolours, colours):
    Poly_Lib.GP_MapSampler.restype = None
    Poly_Lib.GP_MapSampler.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_uint8, ctypes.c_uint8, ctypes.c_uint8]
    c_ = ctypes.c_buffer(file_path)
    Poly_Lib.GP_MapSampler(c_, ngenes, rcolours, colours)


def GenerateGenotypes(file_path, ngenes, colours, samples=-1):
    if samples == -1:
        ExhaustiveMinimalMethod_wrapper(file_path, ngenes, colours)
        PreProcessWrite_wrapper(file_path, "ExhaustiveGenotypes".encode("utf-8"), ngenes, colours)
    else:
        SampleMinimalMethod_wrapper(file_path, ngenes, colours, samples, True)
