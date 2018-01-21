#!/usr/bin/env python3
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pylab import *
from matplotlib.ticker import  MultipleLocator
from matplotlib.ticker import  FormatStrFormatter
from scipy import stats

'''
dev:python3

python2 will have int/float bugs in this code, plz see comments below
'''
__author__ = "CRAFT Harry Qiu <hqiu0129@physics.usyd.edu.au>"

#def _main():
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
parser.add_argument('-s','--show', action='store_true', help='Show')
parser.add_argument('-x','--xaxis',type=int,default=0)
parser.add_argument('-y','--yaxis',type=int,default=0)
parser.add_argument('-o','--output',type=str,default='freddacheck')
parser.add_argument('-d','--set',type=str,default='testset_')
parser.add_argument('--sncut',type=float,default=50.0)
parser.add_argument('--scatter', action='store_true', help='Show')
parser.add_argument('--line', action='store_true', help='Show')
parser.add_argument('--errornone', action='store_true', help='Show')
parser.add_argument('--errorbar', default='std',type=str, help='Show')
parser.add_argument('-l','--label',type=int,default=1,help=' 1 for dm label, 2 for fluence label, 3 for width label')
parser.add_argument('--binmode', type=str,default='mean',help='Show')
parser.add_argument(dest='files', nargs='+')
parser.set_defaults(verbose=False)
values = parser.parse_args()
if values.verbose:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

pscat=values.scatter
pguide=values.line
limit=2000/1.2
dlim=100
dbin=values.binmode
derr=values.errorbar

###number links just in case you forget
'''
sn=0
sampno=1
time=2
boxcar=3
idt=4
dm=5
beamno=6
'''
units={}
units[0]='S/N'
units[1]='Sample No.'
units[2]='Time from Start (s)'
units[3]='Boxcar Width (tsamps)'
units[4]='idt (tsamps)'
units[5]='DM (pc cm$^{-3}$)'
units[6]='Beam No.'
units[8]='Intrinsic Width (ms)'
units[9]='Insert Offset'
name={}
name[0]='_sn'
name[1]='_Samp'
name[2]='_time'
name[3]='_boxcar'
name[4]='_idt'
name[5]='_dm'
name[6]='_beam'
name[8]='_iwd'
name[9]='_off'
label={}
label[1]=' (pc cm$^{-3}$)'
label[2]=""##' (Jy)'
label[3]=' (ms)'
label[4]=' Samples'
mark={}
mark[0]='v'
mark[1]='o'
mark[2]='d'
mark[3]='s'
mark[4]='o'
col={}
col[0]='red'
col[1]='black'
col[2]='#3e82fc'
col[3]='orange' ###dodger blue
col[4]='green'
col[5]='#9a0eea'  ###violet
col[6]='#d90166'  ###dark hot pink
col[7]='black'
ol=open(values.output+'outlier.txt','w')
histo=open(values.output+"histodata.txt",'w')
histo.write("#### time error, s/n truth, s/n fredda, dm, dm_fredda, width_intrinsic, boxcar_fredda \n")
x=values.xaxis
y=values.yaxis
ident=values.set
upbound=values.sncut



ol.close()
histo.close()
