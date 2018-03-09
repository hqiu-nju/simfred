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



def box_plotter(x_axis,y_axis,label_axis,labels,tag1='DM',tag2='(pc cm$^{-3}$)',derr='std'):
    for i in labels:
        print (i)
        pd_std, bin_edges, binnumber=stats.binned_statistic(x_axis[label_axis==i],y_axis[label_axis==i], statistic=derr,bins=400)
        pd_mean, bin_edges, binnumber=stats.binned_statistic(x_axis[label_axis==i],y_axis[label_axis==i], statistic='mean',bins=400)
        xbinned=bin_edges[:-1]+((bin_edges[1] - bin_edges[0])/2)
        ybinned=pd_mean
        plt.errorbar(xbinned,ybinned,yerr=pd_std,alpha=0.5,markersize=15,label=tag1+'= '+str(i)+tag2,fmt='s')
    plt.legend(loc=0,fontsize=10)
    plt.tight_layout()




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
parser.add_argument('-l','--label',type=int,default=5,help=' 5 for dm label, 1/7 for fluence/sn label, 3/8 for width label')
parser.add_argument('-u','--units',type=int,default=0,help=" 0=(pc cm$^{-3}$) 1=jy 2=ms 3=samples 4=none ")
parser.add_argument('--prefix', type=str,default='',help='fof file prefix')
parser.add_argument(dest='files', nargs='+')
parser.set_defaults(verbose=False)
values = parser.parse_args()
if values.verbose:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)


limit=1000/1.2  ### 2seconds worth of samples for arrival time box
dlim=100
x=values.xaxis
y=values.yaxis


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
label={}
label[0]='S/N'
label[1]='Incident Time(Samples)'
label[2]='Incident Time(Seconds)'
label[3]='Boxcar Width'
label[4]='idt'
label[5]='DM'
label[6]='Beam No.'
label[7]="Fluence"
label[8]='Intrinsic Width'
label[9]='Offset'
label[10]='Spectral Index'
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
name[10]="_specinx"
units={}
units[0]=' (pc cm$^{-3}$)'
units[1]=" Jy"##' (Jy)'
units[2]=' (ms)'
units[3]=' Samples'
units[4]=' '

#ol=open(values.output+'outlier.txt','w')
#histo=open(values.output+"histodata.txt",'w')
#histo.write("#### time error, s/n truth, s/n fredda, dm, dm_fredda, width_intrinsic, boxcar_fredda \n")

tag=values.label
p_label=label[tag]
name_label=name[tag]
p_unit=units[values.units]
'''
plt.figure(1,figsize=(8, 6))
plt.xlabel("False Acquistion Rate",fontsize=15)
plt.ylabel("Detection Rate",fontsize=15)
plt.xlim(-0.01,1.01)
plt.ylim(-0.01,1.01)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
'''
###
####
plt.figure(figsize=(8, 6))
plt.xlabel('Truth '+label[x],fontsize=15)
plt.ylabel('Fredda '+label[y],fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#####
copies=open(values.files[0],'r')
filelist=copies.readlines()
copies.close()

###create array

x_axis=np.array([])
y_axis=np.array([])
label_axis=np.array([])

prob_pd=[]
prob_fa=[]

fof_pref=values.prefix

### for all items in filelist remember to -1 for \n in string, -9 for .candlist  so that is -10 in total.
for cd in filelist:
    cand_name= cd[:-1]
    fof_name=fof_pref+cd[:-9]+"fil.cand.fof"
    if os.path.exists(fof_name):
        cand_array=np.loadtxt(cand_name)
        fof_array=np.loadtxt(fof_name)
        if (len(fof_array.flatten())/12)>1: #### if fredda found more than 1 candidate, fredda fofs have 12 columns
            cycle=len(cand_array.flatten())/11  #### candlist files have 10 columns
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
        cand_array=np.loadtxt(cand_name)
        cycle=len(cand_array.flatten())/11
        x_axis=np.append(x_axis,cand_array.T[x])
        y_axis=np.append(y_axis,np.zeros(int(cycle)))
        label_axis=np.append(label_axis,cand_array.T[tag])
        for i in range(int(cycle)):
            prob_pd.append(0)


labels=np.unique(label_axis)
box_plotter(x_axis,y_axis,label_axis,labels,tag1=p_label,tag2=p_unit)
'''
In [12]: for i in labels:
    ...:     print (i)
    ...:     pd_std, bin_edges, binnumber=stats.binned_statistic(x_axis[label_axis==i],y_axis[label_axis==i], statistic='std',bins=20)
    ...:     pd_mean, bin_edges, binnumber=stats.binned_statistic(x_axis[label_axis==i],y_axis[label_axis==i], statistic='mean',bins=20)
    ...:     xbinned=bin_edges[:-1]+((bin_edges[1] - bin_edges[0])/2)
    ...:     ybinned=pd_mean
    ...:     plt.errorbar(xbinned,ybinned,yerr=pd_std,alpha=0.5,markersize=15,label="Width="+str(i)+" ms",fmt='s')
'''
#plt.show()
plt.savefig(values.output+name_label+'.pdf')
plt.close()



#ol.close()
#histo.close()
