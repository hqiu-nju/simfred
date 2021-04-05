
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
    parser = ArgumentParser(description='Makes a single filterbank with a number of pulses of same property', formatter_class=ArgumentDefaultsHelpFormatter)
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
    parser.add_argument('-A', '--snfac',type=float, default=20,help='Define simulated S/N')
    parser.add_argument('-t', '--tau',type=float, default=0,help='millisecond tsamp')
    parser.add_argument('-a', '--alpha',type=float, default=0,help='millisecond tsamp')
    parser.add_argument('-I', '--spectralindex',type=float, default=0)
    parser.add_argument('-x','--offset',type=float,default=0.5, help='Offset within sample')
    parser.add_argument("--noise",type=int,default=1,help='noise level adjustment')
    parser.add_argument('--inject',action='store_true',help='make ASKAP filterbanks, only works on py27 at the moment')
    parser.add_argument('-N','--nfrb',type=int,default=1, help='how many FRBs to inject')
    parser.add_argument("--fbstd",type=int,default=18,help='filterbank units')
    parser.add_argument("--fbbase",type=int,default=128,help='filterbank baseline')
    parser.add_argument('-L','--nsamp',type=int,default=5000, help='data length')
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
 # GHz
    bandwidth=values.bw
    show=values.show
    output=values.output

    bwchan=bandwidth/nchan
    ftop=fch1/1000
    fcentre=fch1+bandwidth/2 ## MHz
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
    global burst
    print("creating waterfall")
    burst,dedisp_burst=dispersion_waterfall(nchan,nsamp,tsamp,bwchan,fch1,dm,amp,tau1,alpha,width,dmerr,offset,show) ### this is a noise free burst
    print(quick_snr(burst),quick_snr(dedisp_burst))
    if values.inject :
        mockheader=makeheader(fch1,bwchan,nchan,nsamp,dmerr)
        inject(mockheader,output,nsamp,nchan,fbstd,noise,base,N,burst,dedisp_burst,amp)
    else:
        filbank=burst+np.random.randn(nchan, nsamp)
        init_sn=mask_check_sn(filbank)
        print("initial sn",init_sn, amp/init_sn)
        finalfil=burst*amp/init_sn+np.random.randn(nchan, nsamp)
        print("final sn",mask_check_sn(finalfil))
        print("match filter sn",quick_snr(burst/init_sn*amp))
        np.save(arr=(finalfil*fbstd+base).astype(np.uint8),file=output)

def inject(mockheader,output,nsamp,nchan,fbstd,noise,base,nfrb,burst,dedisp_b,amp):
    filterbank=fbio.makefilterbank(output+".fil",header=mockheader)
    # filterbank=sgp.SigprocFile(output+'.fil','w',mockheader)
    # print filterbank.header
    filterbank.writenoise(5000,fbstd*noise,base)
    # noise=(np.random.randn(nchan, nsamp)*fbstd + fbbase).astype(np.uint8)
    # noise.T.tofile(filterbank.fin)
    # burst=dispersion_waterfall(nchan,nsamp,tsamp,bwchan,fch1,dm,amp,tau1,alpha,width,dmerr,offset,show=False)
    for i in range(nfrb):
        np.random.seed(i)
        filbank=burst+np.random.randn(nchan, nsamp)
        dedispfil=dedisp_b+np.random.randn(nchan, nsamp)
        init_sn=check_your_snr(dedispfil)
        print("dedisp vs disp sn",init_sn,mask_check_sn(filbank))

        finalfil=burst*amp/init_sn+np.random.randn(nchan, nsamp)

        print("final sn",mask_check_sn(finalfil),mask_check_sn(dedisp_b*amp/init_sn+np.random.randn(nchan, nsamp)))
        newburst=(finalfil*fbstd+base).astype(np.uint8)
        filterbank.writeblock(newburst)
        # burst.T.tofile(filterbank.fin)
    filterbank.writenoise(5000,fbstd*noise,base)
    filterbank.closefile()

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

######## snr functions here ########
def quick_snr(sf):
    return np.sum(sf[sf>0]**2)**0.5

def check_your_snr(sf):
    time_series=sf.sum(0)
#     maxpos=np.argmax(time_series)
    # print(snr_your(burst.sum(0),width))
#     print(mask.shape)
    basemean=np.mean(time_series[:len(time_series//3)])
#     print(basemean)
    burstcut=time_series-basemean
    std=np.std(burstcut[:len(time_series//3)])
#     print(burstcut.max()/std)
    return burstcut.max()/std

def mask_check_sn(sf):
    time_series=sf.sum(0)
    # print(snr_your(burst.sum(0),width))
#     print(mask.shape)
    basemean=np.mean(time_series[:len(time_series//3)])
    maxpulse=np.sum(np.max(sf,axis=1))-basemean
#     print(basemean)
    burstcut=time_series-basemean
    std=np.std(burstcut[:len(time_series//3)])
#     print(burstcut.max()/std)
    return maxpulse/std

############## waterfall ##############
def dispersion_waterfall(nchan,nsamp,tsamp,bwchan,fch1,dm,amp,tau1,alpha,width,dmerr,offset,show):
    time=np.arange(nsamp)*tsamp
    vif,chan_idx=freq_splitter_idx(nchan,0,nchan,bwchan,fch1)
    bin=10
    matrix=np.ones((nsamp,bin))*np.linspace(-0.5,0.5,bin)*tsamp
    timematrix=(np.ones((nsamp,bin)).T*time).T
    global finergrid
    finergrid=(matrix+timematrix).flatten()
    ampx=amp
    A=100
    base = np.zeros((nchan, nsamp))
    base2 = np.random.randn(nchan, nsamp)*0
    t02=nsamp//3*tsamp
    if dmerr == float(0):
        toas=np.array(delaypos(vif,bwchan,fch1,dm))
    else:
        toas=np.array(delaypos(vif,bwchan,fch1,dm+dmerr))
        # plt.plot(gaus_func(0.3,2000+toas[0],finergrid,0))
        # plt.show()
    for i in range(nchan):
        t0=nsamp//3+toas[i]+offset
        #print (ampx)
        #scat_pulse(t,t0,tau1,dm,dmerr,sigma,alpha,a,vi)
        # print("channel",i)
        if tau1 !=0:
            base[i]+=np.mean(scat_pulse_smear(finergrid,t0,tau1,dm,0,width,alpha,A,vif[i],fch1).reshape(nsamp,-1),axis=1)
            base2[i]+=np.mean(scat_pulse_smear(finergrid,t02,tau1,dm,0,width,alpha,A,vif[i],fch1).reshape(nsamp,-1),axis=1)
        else: #single_pulse_smear(t,t0,dm,dmerr,sigma,a,vi)
            base[i]+=np.mean(single_pulse_smear(finergrid,t0,dm,0,width,A,vif[i],fch1).reshape(nsamp,-1),axis=1)
            base2[i]+=np.mean(single_pulse_smear(finergrid,t02,dm,0,width,A,vif[i],fch1).reshape(nsamp,-1),axis=1)

                # print(quick_snr(base[i]))
    base_sn=quick_snr(base)
    base_sn2=quick_snr(base2)
    print(base_sn,base_sn2)
    base=base/base_sn2*ampx
    base2=base2/base_sn2*ampx

    if show:
        plt.imshow(base,aspect='auto')
        plt.yticks([0,335],[vif[0],vif[335]])
        plt.ylabel('Frequency (MHz)')
        plt.xlabel("Time Samples")
        plt.tight_layout()
        plt.show()
    return base,base2
        # np.save(arr=base2,file=output+"_shifted")

##########
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
    chan_idx=chan_idx.astype(np.int64)
    return vi,chan_idx


### basic functions here ###
### Basic functions should be general purpose use for all data formats,
### you should be able to copy them to other versions directly, please do not change the units in these functions
### functions

def gaus_func(sigi,t0,t,ti):
    #ti=0### gaussian function
    sit=1/np.sqrt(np.pi*2*(sigi**2))*np.exp(-(t-t0-ti)**2/2/(sigi**2)) ### model 0 in ravi 2018
    return sit

### adjust dm
def tidm(dmerr,vi,fch1): ## dispersion time delay offset
    vi=vi/1000 ## MHz---> GHz
    dmerr=dmerr
    beta=2
    ftop=fch1/1000
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


def single_pulse_smear(t,t0,dm,dmerr,sigma,a,vi,fch1):
    ### vi is MHz
    #dmerr=dmerr ### smaller
    ti=tidm(dmerr,vi,fch1) ##ms
#     print(ti)
    smear=delta_t(dm+dmerr,vi) ##ms
#     print(smear)
    width=np.sqrt(sigma**2+smear**2)
    # print(width)
    pulse=gaus_func(width,t0,t,ti) ## create pulse
    # plt.plot(t,pulse)
    # plt.show()
    sf=pulse
    sn_norm=quick_snr(sf)
    # print(sn_norm,a)
    flux=sf/sn_norm### normalise
    return a*flux

def scat_pulse_smear(t,t0,tau1,dm,dmerr,sigma,alpha,a,vi,fch1):
    ### vi is MHz
    vi=vi
    #dmerr=dmerr ## smaller
    ti=tidm(dmerr,vi,fch1)##ms
    smear=delta_t(dm+dmerr,vi) ##ms
    width=np.sqrt(sigma**2+smear**2)
    gt0=np.mean(t)
    pulse=gaus_func(width,t0,t,ti) ## create pulse
    scat_corr=scat(t,gt0,tau1,alpha,vi,ti) ## create scatter kernel
    # flux=convolve(scat_corr,pulse,'same')
    sf=convolve(pulse,scat_corr,'same')
    sn_norm=quick_snr(sf)
    flux=sf/sn_norm### normalise
    return a*flux

##########

if __name__ == '__main__':
    _main()
