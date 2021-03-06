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


'''
dev:python3

python2 will have int/float bugs in this code, please see comments below
'''
__author__ = "CRAFT Harry Qiu <hqiu0129@uni.sydney.edu.au>"

#def _main():
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
parser.add_argument('-t','--truth',type=str,default='candlist')
parser.add_argument('-f','--file',type=str,default='fredda')
parser.add_argument('-m','--mode',type=int,default=0,help='0 parameter comparison draw(2 fof/candlists), 1 pd fa draw (must be list of files), 2 dm label pd fa, 3 fluence label pd fa, 4 width label pd fa')
parser.add_argument('-s','--show', action='store_true', help='Show')
parser.add_argument('-x','--xaxis',type=int,default=0)
parser.add_argument('-y','--yaxis',type=int,default=0)
parser.add_argument('-o','--output',type=str,default='freddacheck')
parser.add_argument('-d','--set',type=str,default='testset_')
parser.add_argument('--sncut',type=float,default=50.0)
#parser.add_argument('-l','--label',type=int,default='',help='this is for mode 2, 1 for dm label, 2 for fluence label, 3 for width label')
#parser.add_argument(dest='files', nargs='+')
parser.set_defaults(verbose=False)
values = parser.parse_args()
if values.verbose:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

limit=10
dlim=5
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
label={}
label[1]=' (pc cm$^{-3}$)'
label[2]=' (Jy)'
label[3]=' (ms)'
mark={}
mark[0]='^'
mark[1]='o'
mark[2]='d'
col={}
col[0]='red'
col[1]='#3e82fc'
col[2]='orange'
col[3]='green' ###dodger blue
col[4]='black'
col[5]='#9a0eea'  ###violet
col[6]='#d90166'  ###dark hot pink
col[7]='black'
ol=open(values.output+'outlier.txt','w')
histo=open(values.output+"histodata.txt",'w')
histo.write("#### time error, s/n truth, s/n fredda, dm, dm_fredda, width_intrinsic, boxcar_fredda \n")
x=values.xaxis
y=values.yaxis
ident=values.set
#plabel=values.label
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
        tcheck=tru.T[1][i]
        dpos=np.intersect1d(np.where(fred.T[5]<dcheck+dlim),np.where(fred.T[5]>dcheck-dlim))
        tpos=np.intersect1d(np.where(fred.T[1]<tcheck+limit),np.where(fred.T[1]>tcheck-limit))
        pos=np.intersect1d(dpos,tpos)
        if len(pos) >0:
            sigcheck=np.min(abs(fred.T[5][pos]-tru.T[5][i]))
            if sigcheck <= 4:
                match=np.where(abs(fred.T[5][pos]-tru.T[5][i])==sigcheck)
                pd[i]= fred.T[x][pos][match][0]
                bpd[i]=1
    for i in range(lf):
        dcheck=fred.T[5][i]
        tcheck=fred.T[1][i]
        dpos=np.intersect1d(np.where(tru.T[5]<dcheck+10),np.where(tru.T[5]>dcheck-10))
        tpos=np.intersect1d(np.where(tru.T[1]<tcheck+limit),np.where(tru.T[1]>tcheck-limit))
        pos=np.intersect1d(dpos,tpos)
        if len(pos) >0:
            sigcheck=np.min(abs(tru.T[0][pos]-fred.T[0][i]))
            if sigcheck <= 4:
                match=np.where(abs(tru.T[5][pos]-fred.T[5][i])==sigcheck)
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
    dmlegend=0
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
                    tcheck=tru.T[1][i]
                    dpos=np.intersect1d(np.where(fred.T[5]<dcheck+dlim),np.where(fred.T[5]>dcheck-dlim))
                    tpos=np.intersect1d(np.where(fred.T[1]<tcheck+limit),np.where(fred.T[1]>tcheck-limit))
                    pos=np.intersect1d(dpos,tpos)
                    if len(pos) >0:
                        sigcheck=np.min(abs(fred.T[5][pos]-tru.T[5][i]))
                        match=np.where(abs(fred.T[5][pos]-tru.T[5][i])==sigcheck)
                        pd[i]= fred.T[y][pos][match][0]
                        bpd[i]=1
                for i in range(lf):
                    dcheck=fred.T[5][i]
                    tcheck=fred.T[1][i]
                    dpos=np.intersect1d(np.where(tru.T[5]<dcheck+10),np.where(tru.T[5]>dcheck-10))
                    tpos=np.intersect1d(np.where(tru.T[1]<tcheck+limit),np.where(tru.T[1]>tcheck-limit))
                    pos=np.intersect1d(dpos,tpos)
                    if len(pos) >0:
                        sigcheck=np.min(abs(tru.T[0][pos]-fred.T[0][i]))
                        match=np.where(abs(tru.T[5][pos]-fred.T[5][i])==sigcheck)
                        fa[i]= tru.T[x][pos][match][0]
                        bfa[i]=1
                print (bfa,bpd)
                pdx.append(1.-float(sum(bfa))/len(bfa))
                pdy.append(float(sum(bpd))/len(bpd))
            else:
                #print('single fredda')
                for i in range(lt):
                    dcheck=tru.T[5][i]
                    tcheck=tru.T[1][i]
                    dpos=np.intersect1d(np.where(fred.T[5]<dcheck+dlim),np.where(fred.T[5]>dcheck-dlim))
                    tpos=np.intersect1d(np.where(fred.T[1]<tcheck+limit),np.where(fred.T[1]>tcheck-limit))
                    pos=np.intersect1d(dpos,tpos)
                    if len(pos) > 0:
                        pd[i]= fred.T[y]
                        bpd[i]=1

                dcheck=fred.T[5]
                tcheck=fred.T[1]
                dpos=np.intersect1d(np.where(tru.T[5]<dcheck+10),np.where(tru.T[5]>dcheck-10))
                tpos=np.intersect1d(np.where(tru.T[1]<tcheck+limit),np.where(tru.T[1]>tcheck-limit))
                pos=np.intersect1d(dpos,tpos)
                if len(pos) >0:
                    sigcheck=np.min(abs(tru.T[0][pos]-fred.T[0]))
                    match=np.where(abs(tru.T[5][pos]-fred.T[5])==sigcheck)
                    fa[0]= tru.T[x][pos][match][0]
                    bfa[0]=1
                print (bfa,bpd)
                if bpd < 0.9:
                    ol.write("pd "+fredfile+" "+str(dp)+" "+str(wp)+" "+str(fp)+" "+str(bpd)+"\n")
                if bfa > 0.1:
                    ol.write("fa "+fredfile+" "+str(dp)+" "+str(wp)+" "+str(fp)+" "+str(bfa)+"\n")
                pdx.append(1.-float(sum(bfa))/len(bfa))
                pdy.append(float(sum(bpd))/len(bpd))
            if xamax < int(np.max(fa))+1:
                xamax=int(np.max(fa))+5
            if yamax < int(np.max(pd))+1:
                yamax=int(np.max(pd))+5

            pylab.scatter(tru.T[x],pd,color='darkblue')
            pylab.scatter(fa,fred.T[y],color='darkblue')
        else:
            #print('no match file')
            pdx.append(0.)
            pdy.append(0.)
            pylab.scatter(tru.T[x],pd,color='darkblue')
    meanie=np.arange(xamax)
    if x==y:
        pylab.plot(meanie,meanie,color='orange')
    plt.legend(loc=0)
    pylab.xlim(-0.1,xamax)
    pylab.ylim(-0.1,yamax)
    if values.show:
        pylab.show()
    pylab.savefig(values.output+"compare.png")
    pylab.close()
    plt.xlabel("False Acquistion Rate")
    plt.ylabel("Detection Rate")
    plt.xlim(-0.01,1.01)
    plt.ylim(-0.01,1.01)
    plt.scatter(pdx,pdy)
    if values.show:
        plt.show()

    plt.savefig(values.output+"pdpfa.png")
    plt.close()

if values.mode==2:  ### mixed data open fof files with different flu dm and width as indicator to match each set dm divide
    plt.figure(1,figsize=(12, 9))
    plt.xlabel("False Acquistion Rate",fontsize=15)
    plt.ylabel("Detection Rate",fontsize=15)
    plt.xlim(-0.01,1.01)
    plt.ylim(-0.01,1.01)
    ###
    #ol=open(ident+".ol.dat","w")
    ####
    plt.figure(2,figsize=(12, 9))
    plt.xlabel('Truth '+units[x],fontsize=15)
    plt.ylabel('Fredda '+units[y],fontsize=15)
    #####
    tom=np.loadtxt(values.truth,dtype=float)
    dil=np.where(tom.T[0]<=values.sncut)
    t=tom[dil]
    dmlist=np.unique(t.T[5])
    flulist=np.unique(t.T[7])
    widthlist=np.unique(t.T[8])
    ####boxcar= t.T[3] idt=
    a=len(dmlist)
    b=len(flulist)
    c=len(widthlist)
    yamax=10.
    xamax=10.
    dmax = (np.max(dmlist))
    #ddf = int(dmax)/3 +1
    pluck = 0
    cluck = 1
    kluck = 0
    for d in range(a):
        dp=dmlist[d]
        kk=kluck
        if dp <= dmax/3:
            pluck = 0
            if dmlist[d+1] > dmax/3:
                kluck = -1
        elif dp <= dmax/3*2:
            pluck = 1
            if dmlist[d+1] > dmax/3*2:
                kluck = -1
        else:
            pluck = 2
            if dmlist[d] == dmax:
                kluck = -1
        print (pluck,dp)
        #if
        for j in range(b):
            for k in range(c):
                fp=flulist[j]
                wp=widthlist[k]
                depo_dm=np.where(t.T[5]==dp)
                depo_fl=np.where(t.T[7]==fp)
                depo_wd=np.where(t.T[8]==wp)
                mix=np.intersect1d(depo_dm,depo_fl)
                mix=np.intersect1d(mix,depo_wd)
                limit=wp/1.2
                if len(mix) == 0:
                    continue
                tru=t[mix]
                fredfile=ident+"{0:04}".format(int(dp))+'_'+"{0:03}".format(wp)+'_'+"{0:03}".format(fp)+'_fixed.fil.cand.fof'
                if os.path.exists(fredfile):
                    fred=np.loadtxt(fredfile,dtype=float)
                    lt=len(mix)
                    lf=int(len(fred.flatten())/12)
                    pd=np.zeros(lt,dtype=float)
                    fa=np.zeros(lf,dtype=float)
                    bpd=np.zeros(lt,dtype=bool)
                    bfa=np.zeros(lf,dtype=bool)
                    if lf >1:
                        for i in range(lt):
                            dcheck=tru.T[5][i]
                            tcheck=tru.T[1][i]
                            dpos=np.intersect1d(np.where(fred.T[5]<dcheck+dlim),np.where(fred.T[5]>dcheck-dlim))
                            tpos=np.intersect1d(np.where(fred.T[1]<tcheck+limit),np.where(fred.T[1]>tcheck-limit))
                            pos=np.intersect1d(tpos,dpos)
                            if len(pos) >0:
                                sigcheck=np.min(abs(fred.T[0][pos]-tru.T[0][i]))
                                match=np.where(abs(fred.T[0][pos]-tru.T[0][i])==sigcheck)
                                pd[i]= fred.T[y][pos][match][0]
                                bpd[i]=1
                        for i in range(lf):
                            tcheck=fred.T[1][i]
                            dcheck=fred.T[5][i]
                            tpos=np.intersect1d(np.where(tru.T[1]<tcheck+limit),np.where(tru.T[1]>tcheck-limit))
                            dpos=np.intersect1d(np.where(tru.T[5]<dcheck+10),np.where(tru.T[5]>dcheck-10))
                            pos=np.intersect1d(tpos,dpos)
                            if len(pos) >0:
                                sigcheck=np.min(abs(tru.T[0][pos]-fred.T[0][i]))
                                match=np.where(abs(tru.T[0][pos]-fred.T[0][i])==sigcheck)
                                fa[i]= tru.T[x][pos][match][0]
                                bfa[i]=1
                                #print(tcheck,dcheck,tru.T[5][pos][match][0])
                                histo.write("%d %f %f %f %f %f %d\n"%((fred.T[1][i]-tru.T[1][pos][match][0]),tru.T[0][pos][match][0],fred.T[0][i],tru.T[5][pos][match][0],fred.T[5][i],tru.T[8][pos][match][0],fred.T[3][i]))
                        pdx=(1.-float(sum(bfa))/len(bfa))
                        pdy=(float(sum(bpd))/len(bpd))
                    else:
                        #print('single fredda')
                        for i in range(lt):
                            tcheck=tru.T[1][i]
                            dcheck=tru.T[5][i]
                            dpos=np.intersect1d(np.where(fred.T[5]<dcheck+dlim),np.where(fred.T[5]>dcheck-dlim))
                            tpos=np.intersect1d(np.where(fred.T[1]<tcheck+limit),np.where(fred.T[1]>tcheck-limit))
                            pos=np.intersect1d(tpos,dpos)
                            if len(pos) > 0:
                                sigcheck=abs(tru.T[0][i]-fred.T[0])
                                pd[i]= fred.T[y]
                                bpd[i]=1
                                #if fred.T[0] < 15 and (tru.T[0][i]>25):
                                #     ol.write(fredfile+" "+str(fred.T[2])+" "+str(fred.T[0])+" "+str(tru.T[0][i])+"\n")
                        tcheck=fred.T[1]
                        dcheck=fred.T[5]
                        tpos=np.intersect1d(np.where(tru.T[1]<tcheck+limit),np.where(tru.T[1]>tcheck-limit))
                        dpos=np.intersect1d(np.where(tru.T[5]<dcheck+10),np.where(tru.T[5]>dcheck-10))
                        pos=np.intersect1d(tpos,dpos)
                        if len(pos) >0:
                            sigcheck=np.min(abs(tru.T[0][pos]-fred.T[0]))
                            match=np.where(abs(tru.T[0][pos]-fred.T[0])==sigcheck)
                            fa[0]= tru.T[x][pos][match][0]
                            bfa[0]=1
                            #print(tcheck,dcheck,tru.T[5][pos][match][0])
                            #histo.write("#### time error, s/n truth, s/n fredda, dm, dm fredda, boxcar, boxcar fredda \n")
                            histo.write("%d %f %f %f %f %d %d\n"%((fred.T[1]-tru.T[1][pos][match][0]),tru.T[0][pos][match][0],fred.T[0],tru.T[5][pos][match][0],fred.T[5],tru.T[3][pos][match][0],fred.T[3]))
                        pdx=(1.-float(sum(bfa))/len(bfa))
                        pdy=(float(sum(bpd))/len(bpd))
                    if xamax < int(np.max(fa))+1:
                        xamax=int(np.max(fa))+5
                    if yamax < int(np.max(pd))+1:
                        yamax=int(np.max(pd))+5
                    #if (sum((tru.T[x]-pd)>20)) > 0:
                        #print(tru.T[x]-pd)
                        #print(dp,wp,fp)

                    #print(fa-fred.T[y])
                    plt.scatter(tru.T[x],pd,color=col[kk],marker=mark[pluck],alpha=0.5,s=15)
                    plt.scatter(fa,fred.T[y],color=col[kk],marker=mark[pluck],alpha=0.5,s=15)
                else:
                    pdx=0.0
                    pdy=0.0
                ol.write("fa "+fredfile+" "+str(dp)+" "+str(wp)+" "+str(fp)+" "+str(pdx)+"\n")
                ol.write("pd "+fredfile+" "+str(dp)+" "+str(wp)+" "+str(fp)+" "+str(pdy)+"\n")
                plt.figure(1)
                plt.scatter(pdx,pdy,color=col[kk],marker=mark[pluck],alpha=0.5,s=15)
                plt.figure(2)
        plt.scatter(tru.T[x],pd,color=col[kk],marker=mark[pluck],label='DM ='+str(dp)+label[1],alpha=0.5,s=7)
        plt.figure(1)
        plt.scatter(pdx,pdy,color=col[kk],marker=mark[pluck],label='DM ='+str(dp)+label[1],alpha=0.5,s=7)
        plt.figure(2)
        kluck+=1
    meanie=np.arange(xamax)
    if x==y:
        plt.plot(meanie,meanie,color='orange')
    plt.legend(loc=0,fontsize=15)
    plt.xlim(-0.1,xamax)
    plt.ylim(-0.1,yamax)
    if values.show:
        plt.show()
    plt.savefig(values.output+"compare.png")
    plt.close()
    plt.figure(1)
    plt.legend(loc=0,fontsize=15)


    plt.savefig(values.output+"pdpfa.png")
    plt.close()

ol.close()
histo.close()
