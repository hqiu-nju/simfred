import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from numpy import convolve
from matplotlib.gridspec import GridSpec
import fbio


class spectra:
    def __init__(self,fch1=1100,nchan=336,bwchan=1,tsamp=1):
        """initiate function for creating a mock dynamic spectrum. This sets up the header.
        Parameters
        ----------
        fch1 : float
            First channel frequency (MHz)
        nchan : int
            Total number of channels
        bwchan : float
            channel bandwidth (MHz)
        tsamp : float
            time resolution (ms)
        """
        self.ftop=fch1-bwchan*0.5
        self.fch1=fch1
        self.bwchan=bwchan
        self.nchan=np.int(nchan)
        self.tsamp=tsamp
        self.header={'az_start': 0.0,
        'barycentric': None,
        'data_type': 1,
        'fch1':self.fch1,
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
        'tsamp': tsamp/1000,
        'tstart': 57946.52703893818,
        'za_start': 0.0}

    def create_filterbank(self,file_name,std=18,base=127):
        """Create a mock dynamic spectrum filterbank file.
        Parameters
        ----------
        file_name : string
            filename of filterbank
        std : float
            standard deviation of white noise
        base : float
            base level of array
        """
        self.filterbank=fbio.makefilterbank(file_name+".fil",header=self.header)
        self.fil_std=std
        self.fil_base=base

    def closefile(self):
        """Close writing filterbank

        """
        self.filterbank.closefile()

    def writenoise(self,nsamp=5000):
        """Write a block of white noise into the filterbank.
        Parameters
        ----------
        nsamp : int
            length of noise in units of tsamp
        """
        self.filterbank.writenoise(nsamp,self.fil_std,self.fil_base)

    def imprint(self,array):
        """Create a mock dynamic spectrum filterbank file.
        Parameters
        ----------
        nsamp : int
            length of noise in units of tsamp
        """
        bkg=np.random.randn(array.shape[0],array.shape[1])*self.fil_std+self.fil_base
        imprint=(bkg+array).astype(np.uint8)
        self.filterbank.writeblock(imprint)
        return imprint

    def simulate(self,array):
        bkg=np.random.randn(array.shape[0],array.shape[1])*self.fil_std+self.fil_base
        imprint=(bkg+array).astype(np.uint8)
        return imprint

    def burst(self,dm=200,width=1,A=20,t0=2500,nsamp=5000,mode="boxcar",show=False,tau=0.1,alpha=4,offset=0.5,fstart=0,fend=336):
        """Create a pulse.
        Parameters
        ----------
        mode : string
            enter pulse shape used for injection: boxcar,scat,single
            boxcar: dynspec.boxcar
            scat: dynspec.spectra.scat_pulse_smear
            single: dynspec.spectra.single_pulse_smear


        """
        self.dm=dm
        self.width=width
        self.amplitude=A
        self.t0=t0
        self.pulse_nsamp=nsamp
        tsamp=self.tsamp
        time=np.arange(nsamp)*tsamp
        vif,chan_idx=freq_splitter_idx(self.nchan,fstart,fend,self.bwchan,self.fch1)
        self.vif=vif
        bins=10
        matrix=np.ones((nsamp,bins))*np.linspace(-0.5,0.5,bins)*tsamp
        timematrix=(np.ones((nsamp,bins)).T*time).T
        finergrid=(matrix+timematrix).flatten()
        self.grid=finergrid
        self.x_time=time
        base = np.zeros((self.nchan, nsamp))
        base2 = np.zeros((self.nchan, nsamp))
        toas=np.array(delaypos(vif,self.bwchan,self.fch1,dm))
        for i in range(self.nchan):
            ti=t0+toas[i]+offset
            # print(t0)
            # print (ampx)
            # print("channel",i)
            if mode=='boxcar':
                smear=delta_t(dm,vif[i],self.bwchan)
                # print(vif[i],smear)
                box=np.sqrt(width**2+smear**2)
                # print(box)
                base[i]+=boxcar(time,ti,A,box)
                base2[i]+=boxcar(time,t0,A,box)

            elif mode=='scat':
                dm_p=np.mean(self.scat_pulse_smear(finergrid,ti,tau,dm,width,alpha,A,vif[i]).reshape(nsamp,-1),axis=1)
                dedisp_p=np.mean(self.scat_pulse_smear(finergrid,t0,tau,dm,width,alpha,A,vif[i]).reshape(nsamp,-1),axis=1)
                base[i]+=dm_p/np.max(dm_p)*A
                base2[i]+=dedisp_p/np.max(dedisp_p)*A
            elif mode=="single": #single_pulse_smear(t,t0,dm,dmerr,sigma,a,vi)
                dm_p=np.mean(self.single_pulse_smear(finergrid,ti,dm,width,A,vif[i]).reshape(nsamp,-1),axis=1)
                dedisp_p=np.mean(self.single_pulse_smear(finergrid,t0,dm,width,A,vif[i]).reshape(nsamp,-1),axis=1)
                base[i]+=dm_p/np.max(dm_p)*A
                base2[i]+=dedisp_p/np.max(dedisp_p)*A
        if show:
            plt.imshow(base,aspect='auto')
            plt.yticks([0,self.nchan-1],[vif[0],vif[self.nchan-1]])
            plt.ylabel('Frequency (MHz)')
            plt.xlabel("Time Samples")
            plt.tight_layout()
            plt.show()

        # snfactor=np.max(np.sum(base2,axis=0))
        # print(snfactor,snfactor2)
        self.burst_original=base#/snfactor*A
        self.burst_dedispersed=base2#/snfactor*A
        return self.burst_original,self.burst_dedispersed


    def single_pulse_smear(self,t,t0,dm,sigma,a,vi):
        ### vi is MHz
        #dmerr=dmerr ### smaller
        # fch1=self.fch1
        bwchan=self.bwchan
    #     print(ti)
        smear=delta_t(dm,vi,bwchan) ##msdelta_t(dm,v,bwchan)
    #     print(smear)
        width=np.sqrt(sigma**2+smear**2)
        # print(vi,width)
        pulse=gaus_func(t,t0,width/2) ## create pulse
        # plt.plot(t,pulse)
        # plt.show()
        sf=pulse
        # print("CHECK")
        # sn_norm=quick_snr(sf)
        sn_norm=np.max(sf)
        flux=sf/sn_norm*a### normalise
        # print(quick_snr(flux))
        return flux

    def scat_pulse_smear(self,t,t0,tau1,dm,sigma,alpha,a,vi):
        # fch1=self.fch1
        bwchan=self.bwchan
        ### vi is MHz
        vi=vi
        smear=delta_t(dm,vi,bwchan) ##ms
        # print (vi)
        width=np.sqrt(sigma**2+smear**2)
        gt0=np.mean(t)
        pulse=gaus_func(t,t0,width/2) ## create pulse
        scat_corr=scattering(t,gt0,tau1,alpha,vi) ## create scatter kernel
        # flux=convolve(scat_corr,pulse,'same')
        sf=convolve(pulse,scat_corr,'same')
        # sn_norm=quick_snr(sf)
        sn_norm=np.max(sf)

        flux=sf/sn_norm### normalise
        return a*flux
    # def model_pulse(self)

def simulate(array,std=18,base=127):
    bkg=np.random.randn(array.shape[0],array.shape[1])*std+base
    imprint=(bkg+array).astype(np.uint8)
    return imprint

def quick_snr(sf):
    # print("snr_check")
    return np.sum(sf[sf>0]**2)**0.5

def quad_sum(sf):
    # print("snr_check")
    return np.sum(sf**2)**0.5

def old_snr(sf):
    return np.sum(sf)/np.sqrt((sf.flatten()).shape[0])

def boxcar(t,t0,a,width):
    y=np.zeros(t.shape[0])
    mask=(t<t0+width)*(t>t0)
    y[mask]=a
    return y


def gaus_func(t,t0,sigi):
    #ti=0### gaussian function
    sit=1/np.sqrt(np.pi*2*(sigi**2))*np.exp(-(t-t0)**2/2/(sigi**2)) ### model 0 in ravi 2018
    return sit

def scattering(t,t_0,tau1,alpha,v):
    ###tau=tau1/1000 ## ms
    flux=np.zeros(len(t))
    flux[t>=t_0]=np.exp(-(t[t>=t_0]-t_0)/(tau1*(v/1000)**(-alpha)))
    return flux


def tidm(dmerr,vi,fch1): ## dispersion time delay offset
    vi=vi/1000 ## MHz---> GHz
    beta=2
    ftop=fch1/1000
    ### ftop is in GHz
    ti=4.15*dmerr*(vi**(-beta)-ftop**(-beta)) ### ms
    return ti ### ms

def delta_t(dm,v,bwchan): ### calculate dm smearing
    v=v/1000 ###MHz ---> GHz
    B=bwchan ###1MHz channels in fly's eye change if needed
    inverse_v=1/v #### 1/GHz
    #print(v)
    dt=8.3*dm*(inverse_v**3)*B/2 #### unit:us, 2 sigma---> 1 sigma
    return dt/1000 ###us ---> ms


#


def delaypos(f,bwchan,fch1,dm):
    ftop=fch1
    #print(bwchan)
    times=[]
    for i in f:
        #print(ftop,i,dw)
        #print(fchan,dmdelays(dm,fchan,ftop)*(-1*dw/np.abs(dw)))
        times.append(tidm(dm,i,ftop))
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

def fscrunch():
    np.sum((simdata.astype(np.float64)-base),axis=0)
