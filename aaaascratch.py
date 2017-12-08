import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import sigproc as sgp
import time

def injector(frb,x,frbconvolvemap,normmap,toffset,nchan,tsamp,foff,froof,dm,amplitude,flu,w2):
    for c in xrange(nchan):
        ceil = froof + (c)*foff
        floor = froof + (c+1)*foff ###frequency of current channel
        roof_dispersion_delay_ms = 4.15*dm*(froof**-2 - ceil**-2)
        bott_dispersion_delay_ms = 4.15*dm*(froof**-2 - floor**-2)
        roof_dispersion_delay_samp = abs(roof_dispersion_delay_ms/tsamp)+toffset
        bott_dispersion_delay_samp = abs(bott_dispersion_delay_ms/tsamp)+toffset # offset by a bit so the first few samples aren't off the end of the thing
        length_scale = abs(int(roof_dispersion_delay_samp)-int(bott_dispersion_delay_samp))
        smearing_width_samp = abs((roof_dispersion_delay_samp)-(bott_dispersion_delay_samp))
        noise_start=np.round(roof_dispersion_delay_samp)
        noise_end=np.round(bott_dispersion_delay_samp)
        total_width2 = smearing_width_samp**2 + w2
        ac = flu/np.sqrt(total_width2)
        #snr+=ac/np.sqrt(noise_end-noise_start)*np.sqrt(nchan)
        #print ac,snr
        #x = np.exp(-(times-4)**2/total_width2/2.)*gauss_amp
        #x = amplitude*np.exp(-(times-roof_dispersion_delay_samp-smearing_width_samp)**2/total_width2/2.)
        #print smearing_width_samp
        if length_scale >= 1 :
            start_block=int(roof_dispersion_delay_samp)
            start_cut = roof_dispersion_delay_samp-int(roof_dispersion_delay_samp)
            frb[c,start_block] += abs(1-start_cut)*amplitude
            end_block=int(bott_dispersion_delay_samp)
            ###hence the full signal part is start_block+1 to end_block-1
            frb[c,start_block+1:end_block] += 1*amplitude
            end_cut = abs(bott_dispersion_delay_samp-int(bott_dispersion_delay_samp))
            frb[c,end_block] += (end_cut)*amplitude
        else:
            assert int(roof_dispersion_delay_samp)==int(bott_dispersion_delay_samp)
            frb[c,int(roof_dispersion_delay_samp)]=abs(roof_dispersion_delay_samp-bott_dispersion_delay_samp)
        if w2:
            convolved=np.convolve(x,frb[c])
        else:
            convolved=np.append(frb[c],np.zeros(99))
        normfac=np.sum(convolved)
        normmap[c]+=convolved/normfac
        frbconvolvemap[c]+=normmap[c]*flu


nchan = 336 # channels
fch1 = 1.448 # GHz
foff =  -1/1e3 # GHz
froof = 1.464 # GHz need change to upper limit
tsamp = 1.26646875 # milliseconds
amplitude=1
gauss_amp=amplitude # this has to be the same to achieve norm maximum?
t = np.random.rand()*2*0 + 5 # random seconds between 5 and 10 set at 5, currently 5 seconds
#dm = np.random.rand()*400
nsamp = int(t*1000/tsamp)
dataset = (np.random.randn(nchan, nsamp+99)*18 + 128).astype(np.uint8)
frb = np.zeros((nchan, nsamp))
times = np.arange(100)
frbconvolvemap=np.zeros((nchan,nsamp+99))
dm = 1000
fluence = 10  ###jy ms
fluence_fac = 18./2.  ### jansky  to units of sigma? which is 18
flu = fluence*fluence_fac   ###fluence in 8bit digital
###beamno=-1 ##########fake beam identity signature
frb = np.zeros((nchan, nsamp))
frbconvolvemap=np.zeros((nchan,nsamp+99))
normmap=np.zeros((nchan,nsamp+99))
#dm = dm+100
widthms = np.random.rand()*5 + 0.064
widthms= 1.0
print widthms
widthsamp = widthms/tsamp
#fluence= np.random.rand()*40 +10
#flu=fluence*fluence_fac
#toffset=50
#tinsert=(200/tsamp)
toffset=500/tsamp  ###500 time sample
###real sample time=(1+i)*nsamp+toffset
#d = file.readthismeanyseconds(t)
#print t,toffset/tsamp,widthms,dm
width2=widthsamp**2
sig=0
noise=0




x = gauss_amp*np.exp(-(times-50)**2/(2*width2))
idt = abs(4.15*dm*(froof**-2 - (froof+336*foff)**-2)/tsamp)
injector(frb,x,frbconvolvemap,normmap,toffset,nchan,tsamp,foff,froof,dm,amplitude,flu,width2)
for i in range(336):
    pulse=frbconvolvemap[i] >0
    a=(dataset+frbconvolvemap)[i]
    npos=a[pulse] > 137
    sig+=sum(frbconvolvemap[i][pulse])
    noise+=sum((dataset[i][pulse]*npos))

print sig/np.sqrt(noise)
