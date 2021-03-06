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
parser.add_argument('-m','--mode',type=int,default=1,help='1 default, 0 listmode, 2 allplots')
parser.add_argument('-s','--show', action='store_true', help='Show')
#parser.add_argument('--tlist',type=str,default='candlist')
#parser.add_argument('--flist',type=str,default='fredda')
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
psn=0
psampno=1
ptime=2
pboxcar=3
pidt=4
pdm=5
pbeamno=6
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
if values.mode == 1:
    fig, ax = plt.subplots(figsize=(14,9))
    ax.set_xlabel('Truth '+units[x])
    ax.set_ylabel('Fredda '+units[y])
    divider = make_axes_locatable(ax)
    axHistx = divider.append_axes("top", 0.8, pad=0.1, sharex=ax)
    axHisty = divider.append_axes("right", 0.8, pad=0.1, sharey=ax)
    tru=np.loadtxt(values.truth,dtype=float)
    fred=np.loadtxt(values.file,dtype=float)
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
        tpos=np.intersect1d(np.where(fred.T[2]<tcheck+2),np.where(fred.T[2]>tcheck-2))
        pos=np.intersect1d(dpos,tpos)
        if len(pos) >0:
            sigcheck=np.min(abs(fred.T[0][pos]-tru.T[0][i]))
            if sigcheck <= 4:
                match=np.where(abs(fred.T[0][pos]-tru.T[0][i])==sigcheck)
                pd[i]= fred.T[x][pos][match]
                bpd[i]=1
    for i in range(lf):
        dcheck=fred.T[5][i]
        tcheck=fred.T[2][i]
        dpos=np.intersect1d(np.where(tru.T[5]<dcheck+10),np.where(tru.T[5]>dcheck-10))
        tpos=np.intersect1d(np.where(tru.T[2]<tcheck+2),np.where(tru.T[2]>tcheck-2))
        pos=np.intersect1d(dpos,tpos)
        if len(pos) >0:
            sigcheck=np.min(abs(tru.T[0][pos]-fred.T[0][i]))
            if sigcheck <= 4:
                match=np.where(abs(tru.T[0][pos]-fred.T[0][i])==sigcheck)
                fa[i]= tru.T[x][pos][match[0]]
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

elif values.mode == 0:
    t=open(values.truth,'r')
    f=open(values.file,'r')
    t.seek(0)
    f.seek(0)
    a=t.readlines() ### truth
    b=f.readlines() #### fredda
    plrange=len(a)
    print len(a),len(b)
    assert len(a)== len(b)
    fig, ax = plt.subplots(figsize=(14,9))
    ax.set_xlabel('Truth '+units[x])
    ax.set_ylabel('Fredda '+units[y])
    divider = make_axes_locatable(ax)
    axHistx = divider.append_axes("top", 0.8, pad=0.1, sharex=ax)
    axHisty = divider.append_axes("right", 0.8, pad=0.1, sharey=ax)
    #tru=np.loadtxt(a[0][:-1],dtype=float)
    #fred=np.loadtxt(b[0][:-1],dtype=float)
    for i in range(0,plrange):
        tru=np.loadtxt(a[i][:-1],dtype=float)
        fred=np.loadtxt(b[i][:-1],dtype=float)
        lt=len(tru)
        if len(fred)==12:
            lf=1
        else:
            lf=len(fred)
        #print tru
        #print fred

        pd=np.zeros((len(tru)),dtype=float)
        bpd= np.zeros((len(tru)),dtype=bool) ### probability of detection
        fa=np.zeros((lf),dtype=float)
        bfa= np.zeros((lf),dtype=bool) ### false acquistion
        if lf !=1:
            for i in range(lt):
                check=tru.T[2][i]
                if len(np.intersect1d(np.where(fred.T[2]<check+2),np.where(fred.T[2]>check-2))):
                    box=fred.T[2][np.intersect1d(np.where(fred.T[2]<check+2),np.where(fred.T[2]>check-2))] ### time values
                    #print (box)
                    if len(box) !=1:
                        match=box[np.where((box-check)==np.min(box-check))] ##### find
                    else:
                        match=box
                    #print (type(match))
                    #pd[i]= match
                    pd[i]= fred.T[x][np.where(fred.T[2]==float(match))]
                    bpd[i]=1
            for j in range(lf):
                check=fred.T[2][j]
                if len(np.intersect1d(np.where(tru.T[2]<check+2),np.where(tru.T[2]>check-2))):
                    box=tru.T[2][np.intersect1d(np.where(tru.T[2]<check+2),np.where(tru.T[2]>check-2))]  ### time values
                    if len(box) !=1:
                        match=box[np.where((box-check)==np.min(box-check))] ##### find
                    else:
                        match=box
                    fa[j]= tru.T[x][np.where(tru.T[2]==float(match))]
                    bfa[j]=1
        elif lf ==1:
            for i in range(lt):
                check = fred.T[2]
                if check <= tru.T[2][i]+2 and check >= tru.T[2][i]-2:
                    fa[0]= tru.T[x][i]
                    bfa[0]=1
                    pd[i]= fred.T[x]
                    bpd[i]=1
                    break

        #print ('results found detected/accurate')
        #print (sum(bpd),len(bpd),sum(bfa))
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
elif values.mode ==2:
    t=open(values.truth,'r')
    t.seek(0)
    a=t.readlines() ### truth
    ###detection Rate
    plrange=len(a)
    #assert len(a)== len(b)
    fig, ax = plt.subplots(figsize=(14,9))
    ax.set_xlabel('Truth '+units[x])
    ax.set_ylabel('Fredda '+units[y])
    divider = make_axes_locatable(ax)
    axHistx = divider.append_axes("top", 0.8, pad=0.1, sharex=ax)
    axHisty = divider.append_axes("right", 0.8, pad=0.1, sharey=ax)
    #tru=np.loadtxt(a[0][:-1],dtype=float)
    #fred=np.loadtxt(b[0][:-1],dtype=float)
    for i in range(0,plrange):
        tru=np.loadtxt(a[i][:-1],dtype=float)
        lt=len(tru)

        pd=np.zeros((len(tru)),dtype=float)
        bpd= np.zeros((len(tru)),dtype=bool) ### probability of detection
        fname=a[i][:-9]+'fil.cand.fof'
        #print(fname,a[i])
        if os.path.exists(fname):
            fred=np.loadtxt(fname,dtype=float)
            lf=len(fred.flatten())/12
            fa=np.zeros((lf),dtype=float)
            bfa= np.zeros((lf),dtype=bool)
            #print(fred)
            if lf >1:
                #print(lf)
                fa=np.zeros((lf),dtype=float)
                bfa= np.zeros((lf),dtype=bool)
                for i in range(len(tru.T[2])):
                    check=tru.T[2][i]
                    if len(np.intersect1d(np.where(fred.T[2]<check+2),np.where(fred.T[2]>check-2))):
                        box=fred.T[2][np.intersect1d(np.where(fred.T[2]<check+2),np.where(fred.T[2]>check-2))] ### time values
                        #print (box)
                        if len(box) !=1:
                            match=box[np.where((box-check)==np.min(box-check))] ##### find
                        else:
                            match=box
                        #print (type(match))
                        #pd[i]= match
                        pd[i]= fred.T[x][np.where(fred.T[2]==float(match))]
                        bpd[i]=1
                for j in range(lf):
                    check=fred.T[2][j]
                    if len(np.intersect1d(np.where(tru.T[2]<check+2),np.where(tru.T[2]>check-2))):
                        box=tru.T[2][np.intersect1d(np.where(tru.T[2]<check+2),np.where(tru.T[2]>check-2))]  ### time values
                        if len(box) !=1:
                            match=box[np.where((box-check)==np.min(box-check))] ##### find
                        else:
                            match=box
                        fa[j]= tru.T[x][np.where(tru.T[2]==float(match))]
                        bfa[j]=1
            elif lf ==1:
                for i in range(lt):
                    check = fred.T[2]
                    if check <= tru.T[2][i]+2 and check >= tru.T[2][i]-2:
                        fa[0]= tru.T[x][i]
                        bfa[0]=1
                        pd[i]= fred.T[x]
                        bpd[i]=1
                        break

        #print ('results found detected/accurate')
        #print (sum(bpd),len(bpd),sum(bfa))
        print ('detection rate')
        print (str('samp'),float(sum(bpd))/len(bpd))
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










#if __name__ == '__main__':
#    _main()
