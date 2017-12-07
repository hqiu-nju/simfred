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
parser.add_argument('-l','--label',type=int,default=1,help=' 1 for dm label, 2 for fluence label, 3 for width label')
parser.add_argument(dest='files', nargs='+')
parser.set_defaults(verbose=False)
values = parser.parse_args()
if values.verbose:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)
pscat=values.scatter
pguide=values.line
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
tom=np.loadtxt(values.files[0],dtype=float)
dil=np.where(tom.T[0]<=values.sncut)
t=tom[dil]
if values.label == 1:
    alist=np.unique(t.T[5]) #dm
    blist=np.unique(t.T[7]) #fluence
    clist=np.unique(t.T[8]) #width
if values.label == 2:
    alist=np.unique(t.T[7])
    blist=np.unique(t.T[5])
    clist=np.unique(t.T[8])
if values.label == 3:
    alist=np.unique(t.T[8])
    blist=np.unique(t.T[5])
    clist=np.unique(t.T[7])
ltag=values.label
####boxcar= t.T[3] idt=
a=len(alist)
b=len(blist)
c=len(clist)
yamax=10.
xamax=10.
dmax = (np.max(alist))
#ddf = int(dmax)/3 +1
pluck = 0

pdx_array=np.array([])
pdy_array=np.array([])
fax_array=np.array([])
fay_array=np.array([])
for d in range(a):
    if ltag== 1:
        plabel = 'DM = '
        punit=label[1]
    if ltag== 2:
        plabel = 'Fluence = '
        punit=label[2]
    if ltag== 3:
        plabel = 'Width = '
        punit=label[3]
    dp=alist[d]

    if dp <= dmax/3:
        pluck = 0

    elif dp <= dmax/3*2:
        pluck = 1

    else:
        pluck = 2
    kk=pluck
    print (pluck,dp)
    #if
    for j in range(b):
        for k in range(c):
            fp=blist[j]
            wp=clist[k]
            ######## dp- alist fp - blist wp -clist
            if ltag ==1 :  ########  a -dp- dm  b-fp-fl c -wp-wd
                depo_dm=np.where(t.T[5]==dp)
                depo_fl=np.where(t.T[7]==fp)
                depo_wd=np.where(t.T[8]==wp)
            if ltag ==2 : ########  a -dp- fl  b-fp-dm c -wp-wd
                depo_dm=np.where(t.T[5]==fp)
                depo_fl=np.where(t.T[7]==dp)
                depo_wd=np.where(t.T[8]==wp)
            if ltag ==3 : ########  a -dp- wd  b-fp-dm c -wp-fl
                depo_dm=np.where(t.T[5]==fp)
                depo_fl=np.where(t.T[7]==wp)
                depo_wd=np.where(t.T[8]==dp)
            mix=np.intersect1d(depo_dm,depo_fl)
            mix=np.intersect1d(mix,depo_wd)

            if len(mix) == 0:
                continue
            tru=t[mix]
            dmp=t.T[5][depo_dm][0]
            flp=t.T[7][depo_fl][0]
            wdp=t.T[8][depo_wd][0]
            limit=wdp/1.2
            fredfile=ident+"{0:04}".format(int(dmp))+'_'+"{0:03}".format(wdp)+'_'+"{0:03}".format(flp)+'_fixed.fil.cand.fof'
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
                #plt.scatter(tru.T[x],pd,color=col[kk],marker=mark[pluck],alpha=0.5,s=15)
                #plt.scatter(fa,fred.T[y],color=col[kk],marker=mark[pluck],alpha=0.5,s=15)
                #print(dmp,flp,wdp,len(pd),len(tru.T[x]),len(pdx_array),len(pdy_array))
                pdx_array=np.append(pdx_array,tru.T[x])
                pdy_array=np.append(pdy_array,pd)
                fax_array=np.append(fax_array,fa[fa==0])
                fay_array=np.append(fay_array,fred.T[y][fa==0])
                if len(pdx_array)!= len(pdy_array):
                    print('-----bad')
                    print(dmp,flp,wdp,len(pd),len(tru.T[x]),len(pdx_array),len(pdy_array))
            else:
                pdx=0.0
                pdy=0.0
            ol.write("fa "+fredfile+" "+str(dp)+" "+str(wp)+" "+str(fp)+" "+str(pdx)+"\n")
            ol.write("pd "+fredfile+" "+str(dp)+" "+str(wp)+" "+str(fp)+" "+str(pdy)+"\n")
            plt.figure(1)
            plt.scatter(pdx,pdy,color=col[kk],marker=mark[pluck],alpha=0.5,s=15)
            plt.figure(2)
    pd_std, bin_edges, binnumber=stats.binned_statistic(pdx_array[pdy_array>0],pdy_array[pdy_array>0], statistic='std', bins=int(len(pdx_array)/100)+1)
    pd_mean =stats.binned_statistic(pdx_array[pdy_array>0],pdy_array[pdy_array>0], statistic='mean', bins=(len(pdx_array)/100)+1)[0]
    bin_width = (bin_edges[1] - bin_edges[0])
    xbinned=bin_edges[:-1]+ bin_width
    ybinned=pd_mean
    if pscat:
        plt.scatter(pdx_array,pdy_array,color=col[kk],marker=mark[pluck],alpha=0.5,s=15)
    else:
        plt.scatter(pdx_array[pdy_array==0],pdy_array[pdy_array==0],color=col[kk],marker=mark[pluck],alpha=0.5,s=15)
    plt.scatter(fax_array,fay_array,color=col[kk],marker=mark[pluck],alpha=0.5,s=15)
    if pguide:
        binmark=mark[pluck]+"--"
    else:
        binmark=mark[pluck]
    plt.errorbar(xbinned,ybinned,yerr=pd_std,color=col[kk],fmt=binmark,alpha=0.5,ms=15,label=plabel+str(dp)+punit)
    plt.figure(1)
    plt.scatter(pdx,pdy,color=col[kk],marker=mark[pluck],label=plabel+str(dp)+punit,alpha=0.5,s=15)
    plt.figure(2)
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
