
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
# import bilby
# from astropy import units as u
import fbio
import sigproc3 as sgp
MIN_FLOAT = sys.float_info[3]
__author__ = "Harry Qiu"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-d', '--dm',type=float, default=200,help='DM of pulse')
    parser.add_argument('--dedisp',type=float, default=0,help='How much is dedispersed, leave as zero for original pulse')
    parser.add_argument('-w', '--width',type=float, default=0.1,help='millisecond pulse width')
    parser.add_argument('-o', '--output',type=str, default='testfilterbank',help='output filename')
    parser.add_argument('-s','--show', action='store_true', help='Show')
    parser.add_argument('--nchan',type=int,default=336,help='number of channels')
    parser.add_argument('--bw',type=float,default=336,help='bandwidth MHz')
    parser.add_argument('--ftop',type=int,default=1100.5,help='fch1 frequency (MHz)')
    parser.add_argument('--tsamp',type=float,default=1.26,help='millisecond tsamp')
    parser.add_argument('-A', '--snfac',type=float, default=20)
    parser.add_argument('-t', '--tau',type=float, default=0)
    parser.add_argument('-a', '--alpha',type=float, default=0)
    parser.add_argument('-I', '--spectralindex',type=float, default=0)
    parser.add_argument('-x','--offset',type=float,default=0.5, help='Offset within sample')
    parser.add_argument("--noise",type=int,default=1,help='noise level adjustment')
    parser.add_argument('--inject',action='store_true',help='make ASKAP filterbanks, only works on py27 at the moment')
    parser.add_argument('-N','--nfrb',type=int,default=1, help='how many FRBs to inject')
    parser.add_argument("--fbstd",type=int,default=18,help='filterbank units')
    parser.add_argument("--fbbase",type=int,default=128,help='filterbank baseline')
    parser.add_argument('-L','--nsamp',type=int,default=1024, help='data length')
    #parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    dm0=values.dm
    dmerr=values.dedisp
    sigma=values.width
    N=values.nfrb
    dm=dm0
    global ftop,fch1,bwchan
    noise=values.noise
    tsamp=values.tsamp ##ms
    nsamp=values.nsamp
    nchan=values.nchan
    fch1=values.ftop
    ftop=fch1/1000  # GHz
    bandwidth=values.bw
    show=values.show
    output=values.output
    fcentre=fch1+bandwidth/2 ## MHz
    bwchan=bandwidth/nchan
    tau1=values.tau
    alpha=values.alpha
    width=values.width
    amp=values.snfac
    fbstd=values.fbstd
    base=values.fbbase
    offset=values.offset
    ##create base datasample
    #print(bwchan)
    #np.random.seed(25)
    if values.inject :
        mockheader=makeheader(fch1,bwchan,nchan,nsamp,dmerr)
        filterbank=fbio.makefilterbank(output+".fil",header=mockheader)
        # filterbank=sgp.SigprocFile(output+'.fil','w',mockheader)
        # print filterbank.header
        filterbank.writenoise(int(10000//tsamp),fbstd*noise,base)
        # noise=(np.random.randn(nchan, nsamp)*fbstd + fbbase).astype(np.uint8)
        # noise.T.tofile(filterbank.fin)
        burst0=dispersion_waterfall(nchan,nsamp,0,tsamp,bwchan,fch1,dm,amp,tau1,alpha,width,dmerr,offset,show=False)
        for i in range(N):
            np.random.seed(i)
            background=(np.random.randn(nchan, nsamp)*fbstd*noise + base).astype(np.uint8)
            burst=(burst0*fbstd).astype(np.uint8)+background
            if show:
                plt.imshow(burst0,aspect='auto')
                plt.show()
                plt.imshow(burst,aspect='auto')
                plt.show()
            filterbank.writeblock(burst)
            # burst.T.tofile(filterbank.fin)
        filterbank.writenoise(int(10000//tsamp),fbstd*noise,base)
        filterbank.closefile()
        # filterbank.fin.flush()
        # filterbank.fin.close()
    else:
        burst=dispersion_waterfall(nchan,nsamp,noise,tsamp,bwchan,fch1,dm,amp,tau1,alpha,width,dmerr,offset,show)
        np.save(arr=burst,file=output)

def makeheader(freqaskap,bwchan,nchan,nsamp,dmerr):
    header={'az_start': 0.0,
    'barycentric': None,
    'data_type': 1,
    'fch1':freqaskap-0.5*bwchan ,
    'fchannel': None,
    'foff': bwchan,
    'machine_id': 0,
    'nbits': 8,
    'nchans': nchan,
    'nifs': 1,
    'nsamples': None,
    'period': None,
    'pulsarcentric': None,
    'rawdatafile': None,
    'refdm': None,
    'source_name': 'Fake FRB',
    #'src_dej': -32441.084833752,
    #'src_raj': 215344.63079648,
    'src_raj':174540.1662,
    'src_dej':-290029.896,
    'telescope_id': 7,
    'tsamp': 0.00126,
    'tstart': 57946.52703893818,
    'za_start': 0.0}
    return header
def dispersion_waterfall(nchan,nsamp,noise,tsamp,bwchan,fch1,dm,amp,tau1,alpha,width,dmerr,offset,show):
    base = np.zeros((nchan, nsamp))
    time=np.arange(nsamp)*tsamp
    vif,chan_idx=freq_splitter_idx(nchan,0,nchan,bwchan,fch1)
    toas=np.array(delaypos(vif,bwchan,fch1,dm))
    bin=10
    matrix=np.ones((nsamp,bin))*np.linspace(-0.5,0.5,bin)*tsamp
    timematrix=(np.ones((nsamp,bin)).T*time).T
    global finergrid
    finergrid=(matrix+timematrix).flatten()
    ampx=amp
    for i in range(nchan):
        t0=nsamp//2*tsamp+toas[i]+offset
        #print (ampx)
        #scat_pulse(t,t0,tau1,dm,dmerr,sigma,alpha,a,vi)
        if tau1 !=0:
            base[i]+=np.mean(scat_pulse_smear(finergrid,t0,tau1,dm,0,width,alpha,10,vif[i]).reshape(nsamp,-1),axis=1)
        else: #single_pulse_smear(t,t0,dm,dmerr,sigma,a,vi)
            base[i]+=np.mean(single_pulse_smear(finergrid,t0,dm,0,width,10,vif[i]).reshape(nsamp,-1),axis=1)
    base_sn=quick_snr(base)
    base=base/base_sn*ampx + np.random.randn(nchan, nsamp)*noise
    if show:
        plt.imshow(base,aspect='auto')
        plt.yticks([0,335],[vif[0],vif[335]])
        plt.ylabel('Frequency (MHz)')
        plt.xlabel("Time Samples")
        plt.tight_layout()
        plt.show()

    # np.save(arr=base,file=output)
    if dmerr !=float(0):
        base2 = np.random.randn(nchan, nsamp)*noise
        time=np.arange(nsamp)*tsamp
        toas=np.array(delaypos(vif,bwchan,fch1,dm-dmerr))
        for i in range(nchan):
            t0=nsamp//2*tsamp+toas[i]
            #print (ampx)
            #scat_pulse(t,t0,tau1,dm,dmerr,sigma,alpha,a,vi)
            if tau1 !=0:
                base2[i]+=np.mean(scat_pulse_smear(finergrid,t0,tau1,dm,0,width,alpha,1,vif[i]).reshape(nsamp,-1),axis=1)
            else: #single_pulse_smear(t,t0,dm,dmerr,sigma,a,vi)
                base2[i]+=np.mean(single_pulse_smear(finergrid,t0,dm,0,width,1,vif[i]).reshape(nsamp,-1),axis=1)
        base_sn=quick_snr(base2)
        base2= base2/base_sn*ampx +np.random.randn(nchan, nsamp)*noise
        if show:
            plt.figure()
            plt.imshow(base2,aspect='auto')
            plt.yticks([0,335],[vif[0],vif[335]])
            plt.ylabel('Frequency (MHz)')
            plt.xlabel("Time Samples")
            plt.tight_layout()
            plt.show()
        return base2
        # np.save(arr=base2,file=output+"_shifted")

    else:
        return base

def dmdelays(dm,vi,ftop): ## dispersion time delay offset
    beta=2
    #ftop=1464/1000 ## MHz--->GHz
    #fbot=1128/1000 ## MHz--->GHz
    ### ftop is in GHz
    ti=4.15*dm*(vi**(-beta)-ftop**(-beta)) ### ms
    return ti ### ms
def delaypos(f,bwchan,fch1,dm):
    ftop=fch1/1000
    dw=float(bwchan)/1000
    #print(bwchan)
    times=[]
    for i in f:
        fchan=i/1000
        #print(ftop,i,dw)
        #print(fchan,dmdelays(dm,fchan,ftop)*(-1*dw/np.abs(dw)))
        times.append(dmdelays(dm,fchan,ftop))
    #print ("times generated, max delay is",times[-1])
    return times
def freq_splitter_idx(n,skip,end,bwchan,fch1):
    dw=(end-skip)/n
    #print(dw,bwchan)
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
### functions

def gaus_func(sigi,t0,t,ti):
    #ti=0### gaussian function
    sit=1/np.sqrt(np.pi*2*(sigi**2))*np.exp(-(t-t0-ti)**2/2*(sigi**2)) ### model 0 in ravi 2018
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
def scat(t,t0,tau1,alpha,v,ti):
    ###tau=tau1/1000 ## ms
    t_0=t0+ti
    flux=np.zeros(len(t))
    #flux[t>=t0]=np.exp(-(t[t>=t0]-t0)/(tau1*(v/fcentre)**(-alpha)))
    flux[t>=t_0]=np.exp(-(t[t>=t_0]-t_0)/(tau1*(v/fcentre)**(-alpha)))

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
    sf=pulse
    sn_norm=quick_snr(sf)
    flux=sf/sn_norm### normalise
    return a*flux
def single_pulse_smear(t,t0,dm,dmerr,sigma,a,vi):
    ### vi is MHz
    #dmerr=dmerr ### smaller
    ti=tidm(dmerr,vi) ##ms
    smear=delta_t(dm+dmerr,vi) ##ms
    width=np.sqrt(sigma**2+smear**2)
    pulse=gaus_func(width,t0,t,ti) ## create pulse
    sf=pulse
    sn_norm=quick_snr(sf)
    flux=sf/sn_norm### normalise
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
    sf=convolve(pulse,scat_corr,'same')
    sn_norm=quick_snr(sf)
    flux=sf/sn_norm### normalise
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
    # flux=convolve(scat_corr,pulse,'same')
    sf=convolve(pulse,scat_corr,'same')
    sn_norm=quick_snr(sf)
    flux=sf/sn_norm### normalise
    return a*flux

def quick_snr(sf):
    return np.sum(sf[sf>0]**2)**0.5


if __name__ == '__main__':
    _main()
