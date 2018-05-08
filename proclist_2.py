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

def run_filsim(ident,dm,width,flu,x,index):
    name=ident+"{0:04}".format(int(dm))+'_'+"{0:03}".format(width)+'_'+"{0:03}".format(flu)+'_fixed'
    print ('\n---------')
    print ('width, fluence, dm')
    print (width,flu,dm)
    print ('process file:'+name+' \n')
    if index is not None:
        print("spectral index generated")
        os.system('python filsimv3_compare.py -d '+str(dm)+' -w '+str(width)+' -f '+str(flu)+' -o '+name+' -x '+str(x)+' -i -I '+str(index))
    else:
        os.system('python filsimv3_compare.py -d '+str(dm)+' -w '+str(width)+' -f '+str(flu)+' -o '+name+' -x '+str(x))

def sn_sim(ident,dm,width,sn,x,index):
    name=ident+"{0:04}".format(int(dm))+'_'+"{0:03}".format(width)+'_'+"{0:03}".format(sn)+'_fixed'
    print ('\n---------')
    print ('width, S/N, dm')
    print (width,sn,dm)
    print ('process file:'+name+' \n')
    if index is not None:
        print("spectral index generated")
        os.system('python filsimv3_compare.py -d '+str(dm)+' -w '+str(width)+' -a -A '+str(sn)+' -o '+name+' -x '+str(x)+' -i -I '+str(index))
    else:
        os.system('python filsimv3_compare.py -d '+str(dm)+' -w '+str(width)+' -a -A '+str(sn)+' -o '+name+' -x '+str(x))

def _main():
    date=time.strftime('%Y_%m_%d',time.localtime())

    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-d', '--dmax',type=float, default=2000.)
    parser.add_argument('-s', '--stepdm',type=float, default=10)
    parser.add_argument('--dmin',type=float, default=100)
    parser.add_argument('-w', '--width',type=float, default=0.1)
    parser.add_argument('--wmax',type=float, default=10.0)
    parser.add_argument('--stepw',type=float, default=0.5)
    parser.add_argument('--fchange',action='store_true')
    parser.add_argument('--wchange',action='store_true')
    parser.add_argument('-l','--lmode',action='store_true',help='list command dm,width,flu')
    parser.add_argument('-a','--snmode',action='store_true',help='fix sn')
    parser.add_argument('--list',type=str,default='list')
    parser.add_argument('-f', '--fmax',type=float, default=10.)
    parser.add_argument('--fmin',type=float, default=0.5)
    parser.add_argument('--stepflu',type=float, default=0.5)
    parser.add_argument('-n','--name',type=str, default='testset_')
    parser.add_argument('--offset',type = float, default =0.5)
    parser.add_argument('--index',type = float)
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
    wswitch=values.wchange
    wstep=values.stepw
    wmin=values.width
    wmax=values.wmax
    fmax=values.fmax
    fmin=values.fmin
    fswitch=values.fchange
    ident=values.name
    fstep=values.stepflu
    flu=fmax
    width=wmin
    off=values.offset
    index=values.index
    if values.lmode:
        print('run list')
        exc=np.loadtxt(values.list,delimiter=',')
        if values.snmode:
            for i in range(len(exc)):
                sn_sim(ident,exc[i][0],exc[i][1],exc[i][2],off,index)
        else:
            for i in range(len(exc)):
                run_filsim(ident,exc[i][0],exc[i][1],exc[i][2],off,index)
    elif values.snmode:
        ###flu is sn
        print ("snr fix mode")
        if wswitch != True:
            print("fix width")
            if fswitch:
                for j in range(int((fmax-fmin)/fstep)+1):
                    flu=fmin+j*fstep
                    for i in range(int((dmax-dmin)/dstep)+1):
                        dm=i*dstep+dmin
                        sn_sim(ident,dm,width,flu,off,index)
            else:
                for i in xrange(int((dmax-dmin)/dstep)+1):
                    dm=i*dstep+dmin
                    sn_sim(ident,dm,width,flu,off,index)
        else:
            for k in range(int((wmax-wmin)/wstep)+1):
                width=wmin+k*wstep
                if fswitch:
                    for j in range(int((fmax-fmin)/fstep)+1):
                        flu=fmin+j*fstep
                        for i in range(int((dmax-dmin)/dstep)+1):
                            dm=i*dstep+dmin
                            sn_sim(ident,dm,width,flu,off,index)
                else:
                    for i in xrange(int((dmax-dmin)/dstep)+1):
                        dm=i*dstep+dmin
                        sn_sim(ident,dm,width,flu,off,index)



    else:
        if wswitch != True:
            print("fix width")
            if fswitch:
                for j in range(int((fmax-fmin)/fstep)+1):
                    flu=fmin+j*fstep
                    for i in range(int((dmax-dmin)/dstep)+1):
                        dm=i*dstep+dmin
                        run_filsim(ident,dm,width,flu,off,index)
            else:
                for i in xrange(int((dmax-dmin)/dstep)+1):
                    dm=i*dstep+dmin
                    run_filsim(ident,dm,width,flu,off,index)
        else:
            for k in range(int((wmax-wmin)/wstep)+1):
                width=wmin+k*wstep
                if fswitch:
                    for j in range(int((fmax-fmin)/fstep)+1):
                        flu=fmin+j*fstep
                        for i in range(int((dmax-dmin)/dstep)+1):
                            dm=i*dstep+dmin
                            run_filsim(ident,dm,width,flu,off,index)
                else:
                    for i in xrange(int((dmax-dmin)/dstep)+1):
                        dm=i*dstep+dmin
                        run_filsim(ident,dm,width,flu,off,index)
    #os.system('tar -cvzf '+ident+"*.candlist candlist.tar.gz")


if __name__ == '__main__':
    _main()
