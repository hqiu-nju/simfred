import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
# from numpy import convolve
from scipy.signal import convolve
from matplotlib.gridspec import GridSpec
import fbio
import math as m

class TimeSeries:
    def __init__(self,tsamp=1,nsamp=1000,bins=10):
        """initiate function for creating a mock time series. This sets up the frequency.
        Parameters
        ----------
        tsamp : float
            time resolution (ms)
        nsamp : int
            This sets the length of the array. Must be long enough for scattering tail and dispersion track
        bins : int
            grid resolution of the array,
        """
        # self.fch=fch
        # self.bwchan=bwchan
        self.tsamp=tsamp
        self.nsamp=nsamp
        # self.alpha=scatindex
        tsamp=self.tsamp
        time=np.arange(nsamp)*tsamp
        matrix=np.ones((nsamp,bins))*np.linspace(-0.5,0.5,bins)*tsamp
        timematrix=(np.ones((nsamp,bins)).T*time).T
        finergrid=(matrix+timematrix).flatten()
        self.grid=finergrid
        self.x_time=time
    def boxcar(self,t0,width,a):
        tims=boxcar_func(self.x_time,t0,A,width)
        self.spectra=tims/np.max(tims)*a
        return self.spectra
    def pulse(self,t0,width,a):
        tims=np.mean(single_pulse(self.grid,t0,width,100).reshape(self.nsamp,-1),axis=1)
        self.spectra=tims/np.max(tims)*a
        return self.spectra
    def scatp(self,t0,width,a,tau):
        tims=np.mean(scat_pulse(self.grid,t0,tau,width,0,100,1000).reshape(self.nsamp,-1),axis=1)
        self.spectra=tims/np.max(tims)*a
        return self.spectra
    def inverse_scatp(self,t0,width,a,tau):
        tims=np.mean(invscat_pulse(self.grid,t0,tau,width,0,100,1000).reshape(self.nsamp,-1),axis=1)
        self.spectra=tims/np.max(tims)*a
        return self.spectra


class spectra:
    def __init__(self,fch1=1100,nchan=336,bwchan=1,tsamp=1,nbits=8,fbin=10,tbin=10):
        """initiate function for creating a mock dynamic spectrum data. This sets up the header.
        Parameters
        ----------
        fch1 : float
            First channel centre frequency (MHz)
        nchan : int
            Total number of channels
        bwchan : float
            channel bandwidth (MHz)
        tsamp : float
            time resolution (ms)
        """
        # self.ftop=fch1+bwchan*0.5
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
        'nbits': nbits,
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
        vif,chan_idx=freq_splitter_idx(self.nchan,0,self.nchan,self.bwchan,self.fch1)
         ### finer frequency grid use
        # ff,ff_idx=freq_splitter_idx(self.nchan*fbin,0,self.nchan*fbin,self.bwchan/fbin,self.fch1-self.bwchan*0.5)
        self.vif=vif
        self.fbin=fbin
        self.tbin=tbin
        # self.ff_vif=ff
        self.chan_idx=chan_idx
        # self.ff_idx=ff_idx

    def create_filterbank(self,file_name,std=np.sqrt(336),base=127):
        """Create a mock dynamic spectrum filterbank file.
        Parameters
        ----------
        file_name : string
            filename of filterbank
        std : float
            standard deviation of white noise, for normalised noise after fscrunching, set to sqrt(nchan)
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

    def inject(self,array):
        """Create a mock dynamic spectrum filterbank file.
        Parameters
        ----------
        array : numpy array object
            the burst array data to be injected into the filterbank object
        """
        scaledarray=array*self.fil_std/np.sqrt(array.shape[0])
        bkg=np.random.randn(array.shape[0],array.shape[1])*self.fil_std+self.fil_base
        imprint=(bkg+scaledarray).astype(np.uint8)
        self.filterbank.writeblock(imprint)
        self.injected_array=imprint


    def burst(self,t0=100,dm=200,width=1,A=20,nsamp=5000,mode="boxcar",show=False,tau=0.1,alpha=4,offset=0.,dmoff=0,bandfrac=None):
        """Create a dispersed pulse in noiseless data. Outputs both the dedispered and dedispered pulse
        Parameters
        ----------
        mode : string
            Enter pulse shape used for injection: boxcar,scat,single
            boxcar: dynspec.boxcar_func
            scat: dynspec.spectra.scat_pulse_smear
            single: dynspec.spectra.single_pulse_smear
        width : float
            This is the 1-sigma of the gaussian, in units of ms.
            Note: for the boxcar it is the full width of the boxcar
        nsamp : int
            This sets the length of the array. Must be long enough for the dispersion track.
        A : float
            This is now the channel amplitude of the pulse with no frequency variation applied. This is not the same for boxcar mode, this parameter decides the injected value of the boxcar.
        bandfrac : float
            array for channel spectra, input an array the length of channels to scale emission/mimic scintillation.

        """
        # t0=nsamp//3*self.tsamp
        # if self.bwchan > 0:
        #     t0=nsamp//3*2*self.tsamp
        self.dm=dm
        self.width=width
        self.amplitude=A
        self.t0=t0
        self.pulse_nsamp=nsamp
        # tsamp=self.tsamp
        time=np.arange(nsamp)*self.tsamp
        # bins=self.tbin
        # fbin=self.fbin
        scrunchbw=self.bwchan/self.fbin
        matrix=np.ones((nsamp,self.tbin))*np.linspace(-0.5,0.5,self.tbin)*self.tsamp
        timematrix=(np.ones((nsamp,self.tbin)).T*time).T
        finergrid=(matrix+timematrix).flatten()
        self.grid=finergrid
        self.x_time=time
        base = np.zeros((self.nchan*self.fbin, nsamp))
        # base2 = np.zeros((self.nchan*self.fbin, nsamp))
        # print(f"initialising grid, base shape is {np.shape(base)}")
        matrix=np.ones((self.nchan,self.fbin))*np.linspace(-0.5,0.5,self.fbin)*self.bwchan
        ffgrid=(np.ones((self.nchan,self.fbin)).T*self.vif).T
        vif=(ffgrid+matrix).flatten()
        self.ff_vif=vif
        smear=delta_t(dm,vif,scrunchbw) ## add the smear factor here to shift the pulse when considering smearing
        smeared=np.sqrt(smear**2+width**2)
        toas=np.array(delaypos(self.ff_vif,scrunchbw,self.fch1,dm+dmoff))
        self.ff_toas=toas
        # toas_withsmear=toas+np.sqrt(smear**2+width**2)
        self.toas=np.array(delaypos(self.vif,self.bwchan,self.fch1,dm+dmoff))
        effbw=self.bwchan
        if bandfrac is None:
            bandfrac=np.ones(self.nchan)
        for i in range(self.nchan*self.fbin):
            # excess=self.toas[i//self.fbin]-toas[i]
            shift=toas[i]
            # t_dedisp=t0+offset*self.tsamp+excess
            ti=t0+offset*self.tsamp+shift
            # print(t0)
            # print (ampx)
            # print("channel",i)
            if mode=='boxcar':
                # print(vif[i],smear)
                box=width
                # print(box)
                base[i]+=boxcar_func(time,ti,A,box)
                # base2[i]+=boxcar_func(time,t_dedisp,A,box)
            elif mode=='scat':

                dm_p=np.mean(scat_pulse(finergrid,ti,tau,smeared[i],alpha,10,vif[i]).reshape(nsamp,-1),axis=1)
                # dedisp_p=np.mean(scat_pulse(finergrid,t_dedisp,tau,smeared[i],alpha,10,vif[i]).reshape(nsamp,-1),axis=1)
                # base[i]+=dm_p/np.max(dm_p)*A
                # base2[i]+=dedisp_p/np.max(dedisp_p)*A
                base[i]+=dm_p/np.max(dm_p)*A

            elif mode=="single":
                dm_p=np.mean(single_pulse(finergrid,ti,smeared[i],A).reshape(nsamp,-1),axis=1)
                # print("dedisp file")
                # dedisp_p=np.mean(single_pulse(finergrid,t_dedisp,smeared[i],A).reshape(nsamp,-1),axis=1)
                # base[i]+=dm_p/np.max(dm_p)*A
                # base2[i]+=dedisp_p/np.max(dedisp_p)*A
                base[i]+=dm_p/np.max(dm_p)*A

            elif mode=="nosmear": #single_pulse_smear(t,t0,dm,dmerr,sigma,a,vi)
                dm_p=np.mean(single_pulse(finergrid,ti,width,A).reshape(nsamp,-1),axis=1)
                # print("dedisp file")
                # dedisp_p=np.mean(single_pulse(finergrid,t_dedisp,width,A).reshape(nsamp,-1),axis=1)
                base[i]+=dm_p/np.max(dm_p)*A
            # base2[i]+=dedisp_p/np.max(dedisp_p)*A

        if show:
            plt.imshow(base,aspect='auto')
            plt.yticks([0,self.nchan*self.fbin-1],[vif[0],vif[self.nchan-1]])
            plt.ylabel('Frequency (MHz)')
            plt.xlabel("Time Samples")
            plt.tight_layout()
            plt.show()

        # snfactor=np.max(np.sum(base2,axis=0))
        # print(snfactor,snfactor2)
        # dedispersed=np.empty(np.shape(base1))

        # self.base=base
        # p0s=np.argmax(self.base,axis=1)
        # for i in range(self.base.shape[0]):
        #     base2[i]=np.roll(self.base[i],t0-p0s[i])
        # self.base2=base2
        self.burst_original= np.transpose(base.reshape(self.nchan,self.fbin,nsamp).mean(1).T *bandfrac)#/snfactor*A
        # p0s=np.argmax(self.burst_original,axis=1)
        dedispersed=np.empty(np.shape(self.burst_original))
        for i in range(self.burst_original.shape[0]):
            dedispersed[i]=np.roll(self.burst_original[i],-self.toas[i].astype(np.int)) #-self.toas[i].astype(np.int)
        self.burst_dedispersed=dedispersed#/snfactor*A
        return self.burst_original,self.burst_dedispersed

    # def model_pulse(self)
    def write_snr(self):
        ### Harry's fscrunch and L2 snr script
        base2=self.burst_dedispersed
        quadsn=L2_clean(base2)
        fwhm=(m.sqrt(8.0*m.log(2.0)))*self.width

        # print(quadsn)
        # sf=base2

        return f"{self.dm};{self.width};{fwhm};{quadsn}\n",quadsn
    def write_flux(self):
        base2=self.burst_dedispersed
        flux=L2_flux(base2)
        return flux



class fgrid:
    def __init__(self,fch=1000,bwchan=1,nchan=336,tsamp=1,nsamp=1000,tbin=10,fbin=10):
        """Simulate a burst in a higher resolution grid. tgrid is the higher resolution while fgrid is the final dynamic higher resolution 
        Parameters
        ----------
        fch1 : float
            Upper frequency of array, not first channel frequency (MHz)
        bwchan: float
            Final product Channel bandwidth (MHz)
        nchan : int
            Number of channels of final product
        nsamp : int
            This sets the length of the array. Make it long enough for the dispersion/scattering track, use dynspec.tidm to estimate.
        tsamp : float
            time resolution (ms)
        nsamp : int
            This sets the length of the array. Must be long enough for scattering tail and dispersion track
        tbin : int
            grid time resolution of the tgrid array,
        fbin : int
            grid frequency resolution of the fgrid array,
        """
        self.fch=fch
        self.bwchan=bwchan
        self.nchan=nchan
        self.tsamp=tsamp
        self.nsamp=nsamp
        tsamp=self.tsamp
        time=np.arange(nsamp)*tsamp
        tims=TimeSeries(tsamp=tsamp,nsamp=nsamp,bins=tbin)
        self.tims=tims
        self.tgrid=tims.grid
        self.x_time=tims.x_time
        vif,chan_idx=freq_splitter_idx(self.nchan,0,self.nchan,self.bwchan,self.fch)
         ### finer frequency grid use
        # ff,ff_idx=freq_splitter_idx(self.nchan*fbin,0,self.nchan*fbin,self.bwchan/fbin,self.fch1-self.bwchan*0.5)
        self.vif=vif   ### frequency index
        self.fbin=fbin
        self.tbin=tbin
        # self.ff_vif=ff
        self.chan_idx=chan_idx  ### channel index
        vif,chan_idx=freq_splitter_idx(self.nchan*fbin,0,self.nchan*fbin,bwchan/fbin,self.fch-self.bwchan*0.5)
        self.fgrid=vif
        base = np.zeros((self.nchan*self.fbin, nsamp)) ## this is the grid
        self.array=base


    def pulse(self,t0,width,A,tau=10,alpha=4,dm=0,mode='gaussian'):
        """Simulates pulse in datagrid
        Parameters
        ----------
        mode : string
            Enter pulse shape used for injection: gaussian,scat,scat_r
            gaussian: TimeSeries.pulse, Gaussian pulse
            scat: TimeSeries.scatp
            scat_r: TimeSeries.inverse_scatp
        width : float
            This is the 1-sigma of the gaussian, in units of ms.
            Note: for the boxcar it is the full width of the boxcar
        A : float
            This is now the amplitude of the pulse with no frequency variation applied. This is not the same for boxcar mode, this parameter decides the injected value of the boxcar.
        t0 : float
            Pulse position. This is in the units of time not time sample.
        """
        for i in range(self.nchan*self.fbin):
            tstart=t0+tidm(dm,self.fgrid[i],self.fch)
            if mode == 'gaussian':
                self.array[i]=self.tims.pulse(tstart,width,A)
            elif mode =='scat':
                tscat=tau*(self.fgrid[i]/1000)**(-alpha)  ### tau is scaled to 1 GHz
                self.array[i]=self.tims.scatp(tstart,width,A,tscat)
        # self.model_original=base#/snfactor*A
        self.model_burst=np.mean(self.array.reshape(self.nchan,self.fbin,self.nsamp),axis=1)
        return self.model_burst

def scat_pulse_smear(t,t0,tau1,dm,sigma,alpha,a,vi,bwchan):
    # fch1=self.fch1
    # bwchan=self.bwchan
    ### vi is MHz
    vi=vi
    smear=delta_t(dm,vi,bwchan) ##ms smear is 1 sigma
    # print (vi)
    width=np.sqrt(sigma**2+smear**2)
    gt0=np.mean(t)
    pulse=gaus_func(t,t0,width) ## create pulse
    scat_corr=scattering(t,gt0,tau1,alpha,vi) ## create scatter kernel
    # flux=convolve(scat_corr,pulse,'same')
    sf=convolve(pulse,scat_corr,'same',method='fft')
    # sn_norm=quick_snr(sf)
    sn_norm=np.max(sf)

    flux=sf/sn_norm### normalise
    return a*flux

def single_pulse_smear(t,t0,dm,sigma,a,vi,bwchan):
    ### vi is MHz
    #dmerr=dmerr ### smaller
    # fch1=self.fch1
    # bwchan=self.bwchan
    # print (t0)
#     print(ti)
    smear=delta_t(dm,vi,bwchan) ##msdelta_t(dm,v,bwchan)
#     print(smear)
    width=np.sqrt(sigma**2+smear**2)
    # print(vi,width)
    pulse=gaus_func(t,t0,width) ## create pulse
    # plt.plot(t,pulse)
    # plt.show()
    sf=pulse
    # print("CHECK")
    # sn_norm=quick_snr(sf)
    sn_norm=np.sum(sf)
    if sn_norm==0.0:
        print(f"warning check these parameters channel{vi}, time={t0},width={width}, gives max intensity {sn_norm}, pulse may be located out of array limits")
        # plt.plot(t,pulse)
        # plt.show()

    # print(sn_norm)
    flux=sf/sn_norm*a### normalise
    # print(quick_snr(flux))
    return flux

def single_pulse(t,t0,sigma,a):
    ### vi is MHz
    #dmerr=dmerr ### smaller
    # fch1=self.fch1
    # print (t0)
#     print(ti)
#     print(smear)
    width=sigma
    # print(vi,width)
    pulse=gaus_func(t,t0,width) ## create pulse
    # plt.plot(t,pulse)
    # plt.show()
    sf=pulse
    # print("CHECK")
    # sn_norm=quick_snr(sf)
    sn_norm=np.max(sf)
    if sn_norm==0.0:
        print(f"warning check these parameters, time={t0},width={width}, gives max intensity {sn_norm}, pulse may be located out of array limits")
        # print(f"warning check these parameters{t0},{width}, gives max intensity {sn_norm}, pulse may be located out of array limits")
        # plt.plot(t,pulse)
        # plt.show()

    # print(sn_norm)
    flux=sf/sn_norm*a### normalise
    # print(quick_snr(flux))
    return flux

def scat_pulse(t,t0,tau1,sigma,alpha,a,vi):

    # print (vi)
    kernel=t
    gt0=0
    pulse=gaus_func(t,t0,sigma) ## create pulse
    scat_corr=scattering(kernel,gt0,tau1,alpha,vi) ## create scatter kernel
    # flux=convolve(scat_corr,pulse,'same')
    sf=convolve(pulse,scat_corr,'full')[:len(pulse)]
    # sn_norm=quick_snr(sf)
    sn_norm=np.max(sf)

    flux=sf/sn_norm### normalise
    return a*flux

def invscat_pulse(t,t0,tau1,sigma,alpha,a,vi):

    # print (vi)
    kernel=t
    gt0=0
    pulse=gaus_func(t,t0,sigma) ## create pulse
    scat_corr=inverse_scattering(kernel,gt0,tau1,alpha,vi) ## create scatter kernel
    # flux=convolve(scat_corr,pulse,'same')
    sf=convolve(pulse,scat_corr,'full')[:len(pulse)]
    # sn_norm=quick_snr(sf)
    sn_norm=np.max(sf)

    flux=sf/sn_norm### normalise
    return a*flux


def simulate(array,std=18,base=127,outtype=np.uint8):
    bkg=np.random.randn(array.shape[0],array.shape[1])*std+base
    imprint=(bkg+array).astype(outtype)
    return imprint

def quick_snr(sf):
    # print("snr_check")
    ## single band snr method
    return np.sum(sf[sf>0]**2)**0.5

def quad_sum(sf):
    # print("snr_check")
    ## repetition of quick_snr should be deprecated?
    return np.sum(sf**2)**0.5

def old_snr(sf):
    ## arbitary normalised snr method
    return np.sum(sf)/np.sqrt((sf.flatten()).shape[0])

def boxcar_func(t,t0,a,width):
    y=np.zeros(t.shape[0])
    samp_diff=np.diff(t)[0]
    hw=width/2
    p1=np.argmin(np.abs(t-t0+hw))
    y[p1:p1+np.int(width)]=a

    return y


def gaus_func(t,t0,sigi):
    #ti=0### gaussian function returns a function with amplitude of 10
    A=1/np.sqrt(np.pi*2*(sigi**2))
    sit=A*np.exp(-1/2*((t-t0)**2)/(sigi**2)) ### model 0 in ravi 2018

    ### normalisation factor is 1/np.sqrt(np.pi*2*(sigi**2)) replace A with this term for a total of 1 pdf
    return sit

def scattering(t,t_0,tau1,alpha=4,v=1000):
    ###tau=tau1/1000 ## ms
    flux=np.zeros(len(t))
    flux[t>=t_0]=np.exp(-(t[t>=t_0]-t_0)/(tau1*(v/1000)**(-alpha)))
    return flux
def inverse_scattering(t,t_0,tau1,alpha=4,v=1000):
    ###tau=tau1/1000 ## ms
    flux=np.zeros(len(t))
    flux[t<=t_0]=np.exp((t[t<=t_0]-t_0)/(tau1*(v/1000)**(-alpha)))
    return flux
def tidm(dmerr,vi,fch1): ## dispersion time delay offset
    vi=vi/1000 ## MHz---> GHz
    beta=2
    ftop=fch1/1000
    ### ftop is in GHz
    ti=4.1488*dmerr*(vi**(-beta)-ftop**(-beta)) ### ms
    return ti ### ms

def delta_t(dm,v,bwchan): ### calculate dm smearing
    """Calculate dm smearing
    Parameters
    ----------
    dm : float
        Dispersion measure value
    v : float
        frequency of channel MHz
    bwchan : float
        channel bandwidth MHz
    ---------
    Return smearing width in ms

    """
    v=v/1000 ###MHz ---> GHz
    B=bwchan ###1MHz channels in fly's eye change if needed
    inverse_v=1/v #### 1/GHz
    #print(v)
    dt=8.3*dm*(inverse_v**3)*B/2.355 #### unit:us, fwhm ---> 1 sigma
    return dt/1000 ###us ---> ms


#


def get_fwhm(i):
    fwhm=(m.sqrt(8.0*m.log(2.0)))*i
    return fwhm
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
    ### generates the frequency of channels and then group them into subbands, also returns an array that records the channel numbers of each subband
    ###
    dw=(end-skip)/n
    #print(dw,bwchan)
    vi=(np.arange(n)+0.5)*dw*bwchan
    base=fch1+skip*bwchan
    vi=base+vi
    chan_idx=np.arange(n)*dw+skip
    chan_idx=np.append(chan_idx,end)
    chan_idx=chan_idx.astype(np.int64)
    return vi,chan_idx

def fscrunch(array,prepost=2000):

    return fscrunched


def L2_snr(base2):
    ### Harry's fscrunch and L2 snr script
    simdata=simulate(base2,outtype=np.float64) #base2 is the clean burst array
    fscrunched=np.sum((simdata.astype(np.float64)),axis=0)
    fscrun_mean=np.mean(fscrunched)
    fscrun_median=np.median(fscrunched)
    # fscrun_rms=np.std(fscrunched)
    fscrun_mad=np.median(np.abs(fscrunched-fscrun_mean)) ##use MAD
    # print (fscrun_mad)
    mask=np.sum(base2,axis=0)/fscrun_mad>1 # no noise find where the pulse is after fscrunch
    # fwhm=(m.sqrt(8.0*m.log(2.0)))*self.width
    # print("rms",fscrun_rms)
    sf=((fscrunched-fscrun_median)/fscrun_mad)[mask]
    ### real snr here
    quadsn=(np.sum(sf**2)**0.5)
    # print(quadsn)
    # sf=base2

    return quadsn

def L2_clean(base2):
    ### Harry's fscrunch and L2 snr script with no noise, assume rms/std is 1
    ydata=base2 #base2 is the clean burst array
    fscrunched=np.mean((ydata.astype(np.float64)),axis=0)
    mask=np.mean(base2,axis=0)>0 # no noise find where the pulse is after fscrunch
    # fwhm=(m.sqrt(8.0*m.log(2.0)))*self.width
    # print("rms",fscrun_rms)
    sf=fscrunched[mask]
    ### real snr here
    quadsn=(np.sum(sf**2)**0.5)
    # print(quadsn)
    # sf=base2

    return quadsn

def triangle_snr(base2):
    simdata=simulate(base2,outtype=np.float64)
    fscrun_mean=np.mean(simdata)
    fscrun_median=np.median(simdata)
    fscrun_mad=np.median(np.abs(simdata-fscrun_mean))
    mask2d=(base2<=0)
    # mask2d_b=(base2>0)
    # medianresponse=np.sqrt(mask2d_b.sum(0)/336)
    # simdata2=np.copy(simdata)
    ydata=simdata-fscrun_median
    ydata[mask2d]=0
    sf=np.sum((ydata)/fscrun_mad,0)/np.sqrt(simdata.shape[0])
    ### real snr here
    quadsn=(np.sum(sf**2)**0.5)
    return quadsn

def triangle_clean(base2):
    simdata=base2
    effsig=np.sum(base2>0,0)
    mask=effsig>0

    sf=np.sum(simdata.astype(np.float64),0)[mask]
    ### real snr here
    quadsn=(np.sum(sf**2)**0.5)/np.sqrt(np.sum(effsig[mask]**2))
    return quadsn



def rollingbox(base2):
    ### rolling boxcar filter
    simdata=simulate(base2,outtype=np.float64) #base2 is the clean burst array
    fscrunched=np.sum((simdata.astype(np.float64)),axis=0)
    fscrun_mean=np.mean(fscrunched[:2000])
    fscrun_rms=np.std(fscrunched[:2000])
    mask=np.sum(base2,axis=0)/fscrun_rms>1 # no noise find where the pulse is after fscrunch
    p0=np.argmax(np.sum(base2,axis=0))
    snr=0
    for i in range(np.sum(mask)*2):
        width=i+1
    #     print(width)
    #     print(width//2,width-width//2)
        # print(f"{p0-width//2}:{p0+(width-width//2)}")
        box=fscrunched[p0-width//2:p0+(width-width//2)]
        box_snr=np.sum(box)/np.sqrt(width)
        if box_snr > snr:
            snr=box_snr

def L2_flux(base2):
    ### Harry's fscrunch and L2 snr script with no noise, assume rms/std is 1
    ydata=base2 #base2 is the clean burst array
    fscrunched=np.mean((ydata.astype(np.float64)),axis=0)
    mask=np.mean(base2,axis=0)>0 # no noise find where the pulse is after fscrunch
    # fwhm=(m.sqrt(8.0*m.log(2.0)))*self.width
    # print("rms",fscrun_rms)
    sf=fscrunched[mask]
    ### real snr here
    flux=np.sum(sf)
    # print(quadsn)
    # sf=base2

    return flux