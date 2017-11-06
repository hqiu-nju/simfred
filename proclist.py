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
import sigproc as sgp
import time
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

__author__ = "Harry Qiu"

def _main():
    date=time.strftime('%Y_%m_%d',time.localtime())

    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-d', '--dmax',type=float, default=2000.)
    parser.add_argument('-s', '--stepdm',type=float, default=10)
    parser.add_argument('-w', '--width',type=float, default=0.1)
    parser.add_argument('--scranwidth',type=int, default=0, help='0 for standard, 1 for list, 2 for random scramble')
    parser.add_argument('-f', '--fmax',type=float, default=10.)
    parser.add_argument('--fchange',action='store_true')
    parser.add_argument('--dmin',type=float, default=100)
    parser.add_argument('--fmin',type=float, default=0.5)
    parser.add_argument('--stepflu',type=float, default=0.5)
    parser.add_argument('-n','--name',type=str, default='testset_')
    #parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    dmax=values.dmax
    dmin=values.dmin
    dstep=values.stepdm
    wswitch=values.scranwidth
    fmax=values.fmax
    fmin=values.fmin
    fswitch=values.fchange
    width=values.width
    ident=values.name
    fstep=values.stepflu
    flu=fmax
    if wswitch ==0:
        if fswitch:
            for j in range(int((fmax-fmin)/fstep)+1):
                flu=fmin+j*fstep
                for i in range(int((dmax-dmin)/dstep)+1):
                    dm=i*dstep+dmin
                    name=ident+str(int(dm))+'_'+str(width)+'_'+str(flu)+'_fixed.fil'
                    print ('width, fluence, dm')
                    print (width,flu,dm)
                    print ('process file:'+name+' \n')
                    os.system('python filsim.py -d '+str(dm)+' -w '+str(width)+' -f '+str(flu)+' -o '+name)
        else:
            for i in xrange(int((dmax-dmin)/dstep)+1):
                dm=i*dstep+dmin
                name=ident+str(int(dm))+'_'+str((width))+'_fixed.fil'
                print ('width, fluence, dm')
                print (width,flu,dm)
                print ('process file:'+name+' \n')
                os.system('python filsim.py -d '+str(dm)+' -w '+str(width)+' -f '+str(fmax)+' -o '+name)

    #os.system('tar -cvzf '+ident+"*.candlist candlist.tar.gz")



if __name__ == '__main__':
    _main()
