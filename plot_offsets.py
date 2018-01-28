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


def scatter_plotter(x,y,tag,bmode=False):
    plt.figure(2)
    if bmode:
        pd_mean, bin_edges, binnumber=stats.binned_statistic(x[y>0],y[y>0], statistic='mean')
    else:
        plt.scatter(x,y,label=tag)

def box_plotter(x,y,tag,derr='std',dbin='mean',perr=True):
    markersize=15
    plt.figure(2)
    pd_mean, bin_edges, binnumber=stats.binned_statistic(x[y>0],y[y>0], statistic=dbin, bins=int(len(x)/25))
    bin_width = (bin_edges[1] - bin_edges[0])
    xbinned=bin_edges[:-1]+ bin_width/2
    ybinned=pd_mean
    if perr:
        pd_std=0
    else:
        pd_std =stats.binned_statistic(x[y>0],y[y>0], statistic=derr, bins=(len(x)/25))[0]
    plt.errorbar(xbinned,ybinned,yerr=pd_std,alpha=0.5,ms=markersize,label=tag)



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
parser.add_argument('--sncut',type=float,default=50.0)
parser.add_argument('--scatter', action='store_true', help='Show')
parser.add_argument('--line', action='store_true', help='Show')
parser.add_argument('--errornone', action='store_true', help='Show')
parser.add_argument('--errorbar', default='std',type=str, help='Show')
parser.add_argument('-l','--label',type=int,default=5,help=' 5 for dm label, 1/7 for fluence/sn label, 3/8 for width label')
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
limit=1000/1.2  ### 2seconds worth of samples for arrival time box
dlim=100
dbin=values.binmode
derr=values.errorbar
x=values.xaxis
y=values.yaxis
upbound=values.sncut

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
units[1]='Incident Time(Samples)'
units[2]='Incident Time(Seconds)'
units[3]='Boxcar Width'
units[4]='idt'
units[5]='DM'
units[6]='Beam No.'
units[7]="Fluence"
units[8]='Intrinsic Width'
units[9]='Offset'
name={}
name[0]='_sn'
name[1]='_Samp'
name[2]='_time'
name[3]='_boxcar'
name[4]='_idt'
name[5]='_dm'
name[6]='_beam'
name[7]="_flu"
name[8]='_iwd'
name[9]='_off'
label={}
label[1]=' (pc cm$^{-3}$)'
label[2]=" Jy"##' (Jy)'
label[3]=' (ms)'
label[4]=' Samples'
mark={}
mark[0]='v'
mark[1]='o'
mark[2]='d'
mark[3]='s'
mark[4]='o'
ol=open(values.output+'outlier.txt','w')
histo=open(values.output+"histodata.txt",'w')
histo.write("#### time error, s/n truth, s/n fredda, dm, dm_fredda, width_intrinsic, boxcar_fredda \n")

tag=values.label

plt.figure(1,figsize=(12, 9))
plt.xlabel("False Acquistion Rate",fontsize=15)
plt.ylabel("Detection Rate",fontsize=15)
plt.xlim(-0.01,1.01)
plt.ylim(-0.01,1.01)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
###
####
plt.figure(2,figsize=(12, 9))
plt.xlabel('Truth '+units[x],fontsize=15)
plt.ylabel('Fredda '+units[y],fontsize=15)
#####
copies=open(values.files[0],'r')
filelist=copies.readlines()


###create array

x_axis=np.array([])
y_axis=np.array([])
label_axis=np.array([])

prob_pd=[]
prob_fa=[]


### for all items in filelist remember to -1 for \n in string, -9 for .candlist  so that is -10 in total.
for cd in filelist:
    cand_name= cd[:-1]
    fof_name=cd[:-9]+"fil.cand.fof"
    if os.path.exists(fof_name):
        cand_array=np.loadtxt(cand_name,dtype=float)
        fof_array=np.loadtxt(fof_name,dtype=float)
        if (len(fof_array.flatten())/12)>1: #### if fredda found more than 1 candidate, fredda fofs have 12 columns
            cycle=len(cand_array.flatten())/10  #### candlist files have 10 columns
            for i in range(int(cycle)):
                mk_time=cand_array[i][1]  #### time on column 2
                mk_dm=cand_array[i][5]  #### dm column 6
                time_check=np.intersect1d(np.where(fof_array.T[1]>mk_time-limit),np.where(fof_array.T[1]<mk_time+limit))
                ### locating the time of the burst to match in the true results
                if time_check.size: ##### if they found something within the time range
                    if len(time_check) > 1:  ### more than 1 candidate in true file
                        tc=np.where(fof_array.T[0]==fof_array.T[0][time_check].max())[0][0]
                    else:
                        tc=time_check[0]
                    dmcheck=(fof_array[tc][5]>mk_dm-dlim)*(fof_array[tc][5]<mk_dm+dlim) ### check if the dm matches
                    if dmcheck:  ### assign true detection information
                        x_axis=np.append(x_axis,cand_array[i][x])
                        y_axis=np.append(y_axis,fof_array[tc][y])
                        label_axis=np.append(label_axis,cand_array[i][tag])
                        prob_pd.append(1)
                    else: ### assign no detection information because of dm no match
                        x_axis=np.append(x_axis,cand_array[i][x])
                        y_axis=np.append(y_axis,0)
                        label_axis=np.append(label_axis,cand_array[i][tag])
                        prob_pd.append(0)
                else:  ### assign no detection because of empty within time range
                    x_axis=np.append(x_axis,cand_array[i][x])
                    y_axis=np.append(y_axis,0)
                    label_axis=np.append(label_axis,cand_array[i][tag])
                    prob_pd.append(0)
        ####false Acquistion check
            cycle=len(fof_array.flatten())/12
            for i in range(int(cycle)):
                dmcheck=0
                mk_time=fof_array[i][1]  #### time on column 2
                mk_dm=fof_array[i][5]  #### dm column 6
                time_check=np.intersect1d(np.where(cand_array.T[1]>mk_time-limit),np.where(cand_array.T[1]<mk_time+limit))
                if time_check.size:
                    tc=time_check[0]
                    dmcheck=(cand_array[tc][5]>mk_dm-dlim)*(cand_array[tc][5]<mk_dm+dlim)
                    #print(dmcheck)
                    if dmcheck:
                        prob_fa.append(1)
                    else:
                        prob_fa.append(0)
                else:
                    prob_fa.append(0)
    else: ### if no fredda result
        cand_array=np.loadtxt(cand_name,dtype=float)
        cycle=len(cand_array.flatten())/12
        x_axis=np.append(x_axis,cand_array.T[x])
        y_axis=np.append(y_axis,np.zeros(cycle))
        label_axis=np.append(label_axis,cand_array.T[tag])
        for i in range(cycle):
            prob_pd.append(0)


labels=np.unique(label_axis)






ol.close()
histo.close()
