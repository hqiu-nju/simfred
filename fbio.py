#!/usr/bin/env python2
"""
Writes numpy files into filterbank, needs sigproc.py which runs in python2
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
# import bilby
# from astropy import units as u
import sigproc3 as sgp

MIN_FLOAT = sys.float_info[3]
__author__ = "Harry Qiu"



class makefilterbank:
    def __init__(self,filename,header={'az_start': 0.0,
    'barycentric': None,
    'data_type': 1,
    'fch1':1464.0,
    'fchannel': None,
    'foff': -1.0,
    'machine_id': 0,
    'nbits': 8,
    'nchans': 336,
    'nifs': 1,
    'nsamples': None,
    'period': None,
    'pulsarcentric': None,
    'rawdatafile': None,
    'refdm': None,
    'source_name': 'None',
    #'src_dej': -32441.084833752,
    #'src_raj': 215344.63079648,
    'src_raj':174540.1662,
    'src_dej':-290029.896,
    'telescope_id': 7,
    'tsamp': 0.00126646875,
    'tstart': 57946.52703893818,
    'za_start': 0.0}):
        # if header== "Empty":
        #     print("using default header with fch1=1464, nchan=336, foff=-1")
        self.header=header
        self.fbank=sgp.SigprocFile(filename,'w',header)
        self.fbank.seek_data()

    def writeblock(self,input):
        input.T.tofile(self.fbank.fin)
    def writenoise(self,nsamp,std,base):
        noise=(np.random.randn(self.header['nchans'], nsamp)*std + base).astype(np.uint8)
        noise.T.tofile(self.fbank.fin)
    def closefile(self):
        self.fbank.fin.flush()
        self.fbank.fin.close()
