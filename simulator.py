
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from numpy import convolve
from matplotlib.gridspec import GridSpec
import bilby
from astropy import units as u

MIN_FLOAT = sys.float_info[3]
__author__ = "Harry Qiu"



def dmdelays(dm,vi,ftop): ## dispersion time delay offset
    beta=2
    #ftop=1464/1000 ## MHz--->GHz
    #fbot=1128/1000 ## MHz--->GHz
    ### ftop is in GHz
    ti=4.15*dm*np.abs(ftop**(-beta)-vi**(-beta)) ### ms
    return ti ### ms
def delaypos(f,bwchan,fch1,dm):
    ftop=fch1/1000
    dw=bwchan/1000
    print(bwchan)
    times=[]
    for i in f:
        fchan=i/1000
        #print(ftop,i,dw)
        print(fchan,dmdelays(dm,fchan,ftop)*(-1*dw/np.abs(dw)))
        times.append(dmdelays(dm,fchan,ftop)*(-1*dw/np.abs(dw)))
    print ("times generated, max delay is",times[-1])
    return times
def freq_splitter_idx(n,skip,end,bwchan,fch1):
    dw=(end-skip)/n
    print(dw,bwchan)
    vi=np.arange(n)*dw*bwchan+0.5*dw*bwchan
    base=fch1+skip*bwchan
    vi=base+vi
    chan_idx=np.arange(n)*dw+skip
    chan_idx=np.append(chan_idx,end)
    chan_idx=chan_idx.astype(np.int)
    return vi,chan_idx


### basic functions here ###
### Basic functions should be general purpose use for all data formats,
### you should be able to copy them to other versions directly, please do not change the units in these functions

def gaus_func(sigi,t0,t,ti): ### gaussian function
    sit=1/np.sqrt(np.pi*2*(sigi**2))*np.exp(-(t-t0-ti)**2/sigi**2) ### model 0 in ravi 2018
    return sit

### adjust dm
def tidm(dmerr,vi): ## dispersion time delay offset
    vi=vi/1000 ## MHz---> GHz
    dmerr=dmerr
    beta=2
    ### ftop is in GHz
    ti=4.15*dmerr*np.abs(ftop**(-beta)-vi**(-beta)) ### ms
    return ti ### ms

### scattering
def scat(t,t0,tau1,alpha,v):
    ###tau=tau1/1000 ## ms
    flux=np.zeros(len(t)) + MIN_FLOAT
    #flux[t>=t0]=np.exp(-(t[t>=t0]-t0)/(tau1*(v/fcentre)**(-alpha)))
    flux[t>=t0]=np.exp(-(t[t>=t0]-t0)/(tau1*(v/fcentre)**(-alpha)))

    return flux

### dm smearing ## note that this should not be necessary for htr/xtr data


def delta_t(dm,v): ### calculate dm smearing
    v=v/1000 ###MHz ---> GHz
    B=bwchan ###1MHz channels in fly's eye change if needed
    inverse_v=1/v #### 1/GHz
    #print(v)
    dt=8.3*dm*(inverse_v**3)*B/2 #### unit:us, 2 sigma---> 1 sigma
    return dt/1000 ###us ---> ms



# Design functions here

def single_pulse(t,t0,dm,dmerr,sigma,a,vi):
    ### vi is MHz

    #dmerr=dmerr ### smaller
    ti=tidm(dmerr,vi) ##ms
    smear=delta_t(dmerr,vi) ##ms
    width=np.sqrt(sigma**2+smear**2)
    pulse=gaus_func(width,t0,t,ti) ## create pulse
    flux=pulse
    flux/=np.max(flux) ### normalise
    return a*flux
def single_pulse_smear(t,t0,dm,dmerr,sigma,a,vi):
    ### vi is MHz
    #dmerr=dmerr ### smaller
    ti=tidm(dmerr,vi) ##ms
    smear=delta_t(dm+dmerr,vi) ##ms
    width=np.sqrt(sigma**2+smear**2)
    pulse=gaus_func(width,t0,t,ti) ## create pulse
    flux=pulse
    flux/=np.max(flux) ### normalise
    return a*flux
def scat_pulse(t,t0,tau1,dm,dmerr,sigma,alpha,a,vi):
    ### vi is MHz
    vi=vi
    #dmerr=dmerr ## smaller
    ti=tidm(dmerr,vi)##ms
    smear=delta_t(dmerr,vi) ##ms
    width=np.sqrt(sigma**2+smear**2)
    gt0=np.mean(t)
    pulse=gaus_func(width,t0,t,ti) ## create pulse
    scat_corr=scat(t,gt0,tau1,alpha,vi) ## create scatter kernel
    flux=convolve(scat_corr,pulse,'same')
    flux/=np.max(flux) ### normalise
    return a*flux
def scat_pulse_smear(t,t0,tau1,dm,dmerr,sigma,alpha,a,vi):
    ### vi is MHz
    vi=vi
    #dmerr=dmerr ## smaller
    ti=tidm(dmerr,vi)##ms
    smear=delta_t(dm+dmerr,vi) ##ms
    width=np.sqrt(sigma**2+smear**2)
    gt0=np.mean(t)
    pulse=gaus_func(width,t0,t,ti) ## create pulse
    scat_corr=scat(t,gt0,tau1,alpha,vi) ## create scatter kernel
    flux=convolve(scat_corr,pulse,'same')
    flux/=np.max(flux) ### normalise
    return a*flux



from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
parser.add_argument('-d', '--dm',type=float, default=200)
parser.add_argument('--dedisperse',type=float, default=0)
parser.add_argument('-w', '--width',type=float, default=0.1,help='millisecond width')
parser.add_argument('-o', '--output',type=str, default='testfilterbank')
parser.add_argument('-s','--show', action='store_true', help='Show')
parser.add_argument('--nchan',type=int,default=336)
parser.add_argument('--tsamp',type=float,default=1,help='millisecond tsamp')
parser.add_argument('-A', '--snfac',type=float, default=10)
parser.add_argument('-t', '--tau',type=float, default=0)
parser.add_argument('-a', '--alpha',type=float, default=0)
parser.add_argument('-I', '--spectralindex',type=float, default=0)
parser.add_argument('-x','--offset',type=float,default=0.5, help='Offset within sample')
parser.add_argument('-N','--nfrb',type=int,default=1, help='how many FRBs')
parser.add_argument('-L','--nsamp',type=int,default=1024, help='data length')
#parser.add_argument(dest='files', nargs='+')
parser.set_defaults(verbose=False)
values = parser.parse_args()
dm0=values.dm
dmerr=values.dedisperse
sigma=values.width
N=values.nfrb
dm=dm0
tsamp=values.tsamp ##ms
nsamp=values.nsamp
nchan=values.nchan
fch1=1100
ftop=fch1/1000  # GHz
fcentre=fch1+336/2 ## MHz
bwchan=336/nchan
o='simtest'
tau1=values.tau
alpha=values.alpha
width=values.width
amp=values.snfac
##create base datasample
print(bwchan)
#np.random.seed(25)
base = np.random.randn(nchan, nsamp)
time=np.arange(nsamp)*tsamp
vif,chan_idx=freq_splitter_idx(336,0,336,bwchan,fch1)
toas=np.array(delaypos(vif,bwchan,fch1,dm))



for i in range(nchan):
    t0=nsamp//2*tsamp+toas[i]
    ampx=amp*np.random.rand(1)
    #print (ampx)
    #scat_pulse(t,t0,tau1,dm,dmerr,sigma,alpha,a,vi)
    if tau1 !=0:
        base[i]+=scat_pulse_smear(time,t0,tau1,dm,0,width,alpha,ampx,vif[i])
    else: #single_pulse_smear(t,t0,dm,dmerr,sigma,a,vi)
        base[i]+=single_pulse_smear(time,t0,dm,0,width,ampx,vif[i])
if values.show:
    plt.imshow(base,aspect='auto')
    plt.show()
np.save(arr=base,file=values.output)
if dmerr !=0:
    base = np.random.randn(nchan, nsamp)
    time=np.arange(nsamp)*tsamp
    toas=np.array(delaypos(vif,bwchan,fch1,dm-dmerr))
    for i in range(nchan):
        t0=nsamp//2*tsamp+toas[i]
        ampx=amp*np.random.rand(1)
        #print (ampx)
        #scat_pulse(t,t0,tau1,dm,dmerr,sigma,alpha,a,vi)
        if tau1 !=0:
            base[i]+=scat_pulse_smear(time,t0,tau1,dm,0,width,alpha,ampx,vif[i])
        else: #single_pulse_smear(t,t0,dm,dmerr,sigma,a,vi)
            base[i]+=single_pulse_smear(time,t0,dm,0,width,ampx,vif[i])
    if values.show:
        plt.figure()
        plt.imshow(base,aspect='auto')
        plt.show()
    np.save(arr=base,file=values.output+"_shifted")
