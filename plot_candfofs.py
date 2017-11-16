#!/usr/bin/env python
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


'''
dev:python3

python2 will have int/float bugs in this code, plz see comments below
'''
__author__ = "CRAFT Harry Qiu <hqiu0129@physics.usyd.edu.au>"

#def _main():
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
parser.add_argument('-t','--truth',type=str,default='candlist')
parser.add_argument('-f','--file',type=str,default='fredda')
parser.add_argument('-m','--mode',type=int,default=0,help='0 parameter comparison draw(2 fof/candlists), 1 pd fa draw (must be list of files), 2 offset draw(list of files)')
parser.add_argument('-s','--show', action='store_true', help='Show')
parser.add_argument('-x','--xaxis',type=int,default=0)
parser.add_argument('-y','--yaxis',type=int,default=0)
parser.add_argument('-o','--output',type=str,default='freddacheck')
#parser.add_argument(dest='files', nargs='+')
parser.set_defaults(verbose=False)
values = parser.parse_args()
if values.verbose:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)
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
name={}
name[0]='_sn'
name[1]='_Samp'
name[2]='_time'
name[3]='_boxcar'
name[4]='_idt'
name[5]='_dm'
name[6]='_beam'
x=values.xaxis
y=values.yaxis
if values.mode==0:  ######### for simply comparing truth and fredda for one parameter
    fig, ax = plt.subplots(figsize=(14,9))
    tru=np.loadtxt(values.truth,dtype=float)
    fred=np.loadtxt(values.file,dtype=float)
    y=x
    ax.set_xlabel('Truth '+units[x])
    ax.set_ylabel('Fredda '+units[y])
    divider = make_axes_locatable(ax)
    axHistx = divider.append_axes("top", 0.8, pad=0.1, sharex=ax)
    axHisty = divider.append_axes("right", 0.8, pad=0.1, sharey=ax)
    ### mathching with 2 time and 5 dm
    lt=len(tru.T[0])
    lf=len(fred.T[0])
    pd=np.zeros((len(tru)),dtype=float)
    bpd= np.zeros((len(tru)),dtype=bool) ### probability of detection
    fa=np.zeros((lf),dtype=float)
    bfa= np.zeros((lf),dtype=bool) ### false acquistion
    for i in range(lt):
        dcheck=tru.T[5][i]
        tcheck=tru.T[2][i]
        dpos=np.intersect1d(np.where(fred.T[5]<dcheck+10),np.where(fred.T[5]>dcheck-10))
        tpos=np.intersect1d(np.where(fred.T[2]<tcheck+2.5),np.where(fred.T[2]>tcheck-2.5))
        pos=np.intersect1d(dpos,tpos)
        if len(pos) >0:
            sigcheck=np.min(abs(fred.T[0][pos]-tru.T[0][i]))
            if sigcheck <= 4:
                match=np.where(abs(fred.T[0][pos]-tru.T[0][i])==sigcheck)
                pd[i]= fred.T[x][pos][match][0]
                bpd[i]=1
    for i in range(lf):
        dcheck=fred.T[5][i]
        tcheck=fred.T[2][i]
        dpos=np.intersect1d(np.where(tru.T[5]<dcheck+10),np.where(tru.T[5]>dcheck-10))
        tpos=np.intersect1d(np.where(tru.T[2]<tcheck+2.5),np.where(tru.T[2]>tcheck-2.5))
        pos=np.intersect1d(dpos,tpos)
        if len(pos) >0:
            sigcheck=np.min(abs(tru.T[0][pos]-fred.T[0][i]))
            if sigcheck <= 4:
                match=np.where(abs(tru.T[0][pos]-fred.T[0][i])==sigcheck)
                fa[i]= tru.T[x][pos][match][0]
                bfa[i]=1
    print ('detection rate')
    print (str('samp'),i,float(sum(bpd))/len(bpd))
    #print ('false acquistion')
    #print (1-float(sum(bfa))/len(bfa))
    #### note that the dividing here will become wrong in python2
    #print (pd*bpd)
    meanie=np.arange(int(np.max(fa))+1)
    print(int(np.max(fa))+1)
    ax.plot(meanie,meanie,color='orange')
    ax.scatter(tru.T[x],pd,color='darkblue')
    ax.scatter(fa,fred.T[x],color='darkblue')

    axHistx.hist(fa, bins=1000)

    axHisty.hist(pd, bins=1000, orientation='horizontal')
    if values.show:
        plt.show()

    plt.savefig(values.output+name[x]+name[y]+".png")
    plt.close()
if values.mode==1:   ### drawing the probability of detected versus false acquistion
    plt.xlabel("False Acquistion Rate")
    plt.ylabel("Detection Rate")
    plt.xlim(-0.01,1.01)
    plt.ylim(-0.01,1.01)
    pdy=[]
    pdx=[]
    ###
    t=open(values.truth,'r')
    t.seek(0)
    a=t.readlines() ### truth
    ###detection Rate
    plrange=len(a)
    pylab.figure()
    pylab.xlabel('Truth '+units[x])
    pylab.ylabel('Fredda '+units[y])
    xamax=10
    yamax=10
    for i in range(plrange):
        tru=np.loadtxt(a[i][:-1],dtype=float)
        lt=len(tru.T[0])
        pd=np.zeros((len(tru)),dtype=float)
        bpd= np.zeros((len(tru)),dtype=bool) ### probability of detection
        fname=a[i][:-9]+'fil.cand.fof'
        #print(fname,a[i])
        if os.path.exists(fname):
            fred=np.loadtxt(fname,dtype=float)
            lf=len(fred.flatten())/12
            #print(fred)
            #print(lf)
            fa=np.zeros((lf),dtype=float)
            bfa= np.zeros((lf),dtype=bool) ### false acquistion
            if lf >1:
                for i in range(lt):
                    dcheck=tru.T[5][i]
                    tcheck=tru.T[2][i]
                    dpos=np.intersect1d(np.where(fred.T[5]<dcheck+10),np.where(fred.T[5]>dcheck-10))
                    tpos=np.intersect1d(np.where(fred.T[2]<tcheck+2.5),np.where(fred.T[2]>tcheck-2.5))
                    pos=np.intersect1d(dpos,tpos)
                    if len(pos) >0:
                        sigcheck=np.min(abs(fred.T[0][pos]-tru.T[0][i]))
                        match=np.where(abs(fred.T[0][pos]-tru.T[0][i])==sigcheck)
                        pd[i]= fred.T[x][pos][match][0]
                        bpd[i]=1
                for i in range(lf):
                    dcheck=fred.T[5][i]
                    tcheck=fred.T[2][i]
                    dpos=np.intersect1d(np.where(tru.T[5]<dcheck+10),np.where(tru.T[5]>dcheck-10))
                    tpos=np.intersect1d(np.where(tru.T[2]<tcheck+2.5),np.where(tru.T[2]>tcheck-2.5))
                    pos=np.intersect1d(dpos,tpos)
                    if len(pos) >0:
                        sigcheck=np.min(abs(tru.T[0][pos]-fred.T[0][i]))
                        match=np.where(abs(tru.T[0][pos]-fred.T[0][i])==sigcheck)
                        fa[i]= tru.T[x][pos][match][0]
                        bfa[i]=1
                print (bfa,bpd)
                pdx.append(1.-float(sum(bfa))/len(bfa))
                pdy.append(float(sum(bpd))/len(bpd))
            else:
                #print('single fredda')
                for i in range(lt):
                    dcheck=tru.T[5][i]
                    tcheck=tru.T[2][i]
                    dpos=np.intersect1d(np.where(fred.T[5]<dcheck+10),np.where(fred.T[5]>dcheck-10))
                    tpos=np.intersect1d(np.where(fred.T[2]<tcheck+2.5),np.where(fred.T[2]>tcheck-2.5))
                    pos=np.intersect1d(dpos,tpos)
                    if len(pos) > 0:
                        pd[i]= fred.T[x]
                        bpd[i]=1

                dcheck=fred.T[5]
                tcheck=fred.T[2]
                dpos=np.intersect1d(np.where(tru.T[5]<dcheck+10),np.where(tru.T[5]>dcheck-10))
                tpos=np.intersect1d(np.where(tru.T[2]<tcheck+2.5),np.where(tru.T[2]>tcheck-2.5))
                pos=np.intersect1d(dpos,tpos)
                if len(pos) >0:
                    sigcheck=np.min(abs(tru.T[0][pos]-fred.T[0]))
                    match=np.where(abs(tru.T[0][pos]-fred.T[0])==sigcheck)
                    fa[i]= tru.T[x][pos][match][0]
                    bfa[i]=1
                print (bfa,bpd)
                pdx.append(1.-float(sum(bfa))/len(bfa))
                pdy.append(float(sum(bpd))/len(bpd))
            meanie=np.arange(int(np.max(fa))+1)
            pylab.plot(meanie,meanie,color='orange')
            if xamax < int(np.max(fa))+1:
                xamax=int(np.max(fa))+5
            if yamax < int(np.max(pd))+1:
                yamax=int(np.max(pd))+5
            pylab.scatter(tru.T[x],pd,color='darkblue')
            pylab.scatter(fa,fred.T[x],color='darkblue')
        else:
            #print('no match file')
            pdx.append(0.)
            pdy.append(0.)
            pylab.scatter(tru.T[x],pd,color='darkblue')
    pylab.xlim(-0.1,xamax)
    pylab.ylim(-0.1,yamax)
    pylab.savefig(values.output+"compare.png")
    if values.show:
        pylab.show()
    plt.xlabel("False Acquistion Rate")
    plt.ylabel("Detection Rate")
    plt.xlim(-0.01,1.01)
    plt.ylim(-0.01,1.01)
    plt.scatter(pdx,pdy)
    if values.show:
        plt.show()

    plt.savefig(values.output+"pdpfa.png")
    plt.close()
