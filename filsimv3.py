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
import sigproc as sgp
import time

__author__ = "Harry Qiu"


def injector(frb,x,frbconvolvemap,normmap,toffset,nchan,tsamp,foff,froof,dm,amplitude,flu,w2,nsamp):
    nodm=np.zeros((nchan, nsamp))
    nodmcon=np.zeros((nchan,nsamp+99))
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
            nodm[c,int(toffset)] += abs(1-start_cut)*amplitude
            end_block=int(bott_dispersion_delay_samp)
            ###hence the full signal part is start_block+1 to end_block-1
            frb[c,start_block+1:end_block] += 1*amplitude
            nodm[c,int(toffset+1):int(toffset+length_scale)]+= 1*amplitude
            end_cut = abs(bott_dispersion_delay_samp-int(bott_dispersion_delay_samp))
            frb[c,end_block] += (end_cut)*amplitude
            nodm[c,int(toffset+length_scale)] += (end_cut)*amplitude
        else:
            assert int(roof_dispersion_delay_samp)==int(bott_dispersion_delay_samp)
            frb[c,int(roof_dispersion_delay_samp)]=abs(roof_dispersion_delay_samp-bott_dispersion_delay_samp)
            nodm[c,int(toffset)]+= abs(roof_dispersion_delay_samp-bott_dispersion_delay_samp)
        if w2:
            convolved=np.convolve(x,frb[c])
            ndcon=np.convolve(x,nodm[c])
        else:
            convolved=np.append(frb[c],np.zeros(99))
            ndcon=np.append(nodm[c],np.zeros(99))
        normfac=np.sum(convolved)
        normmap[c]+=convolved/normfac
        frbconvolvemap[c]+=normmap[c]*flu
        nodmcon[c]+=ndcon/normfac * flu
    return frbconvolvemap,normmap,nodmcon


date=time.strftime('%Y_%m_%d',time.localtime())
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
parser.add_argument('-d', '--dm',type=float, default=1000)
parser.add_argument('-w', '--width',type=float, default=0.1)
parser.add_argument('-p','--printer', action='store_true', help='Show progress')
#parser.add_argument('-m', '--menu',type=float, default=0, help='0 for test, 1 for list, 2 for random scramble')
#parser.add_argument('-l', '--list',type=str, default='')
parser.add_argument('-f', '--fluence',type=float, default=10)
parser.add_argument('-b', '--base',type=str, default='cut.fil')
parser.add_argument('-o', '--output',type=str, default='set_'+date+'.fil')
parser.add_argument('-s','--show', action='store_true', help='Show')
#parser.add_argument(dest='files', nargs='+')
parser.set_defaults(verbose=False)
values = parser.parse_args()
if values.verbose:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)
#date=time.strftime('%Y_%m_%d_%H_%M',time.localtime())
#date='2000'
inputname=values.base
outputname=values.output
f=open(outputname[:-4]+'.candlist','w')
f.write("##new file of truth of truth!!!!!"+outputname+"\n")
f.write("# S/N, sampno, secs from file start, boxcar, idt, dm, beamno\n")
fprint=values.printer
readin=sgp.SigprocFile(inputname)
readin.header['nchans']=336
mkout=sgp.SigprocFile(outputname,'w',readin.header)
mkout.seek_data()
readin.seek_data()
#datastring=np.fromfile(readin.fin,dtype=np.uint8)
#dataset=datastring.reshape(len(datastring)/336,1,336)
mkout.seek_data()
readin.seek_data()
#### these are system/constant parameters
nchan = readin.header['nchans']# channels
fch1 = 1.464 # GHz
foff =  -1/1e3 # GHz
froof = fch1 # GHz need change to upper limit
tsamp = float(readin.header['tsamp'])*1000
print (tsamp)# milliseconds
amplitude=1.
gauss_amp=amplitude # this has to be the same to achieve norm maximum?
t = np.random.rand()*2*0 + 5 # random seconds between 5 and 10 set at 5, currently 5 seconds
#dm = np.random.rand()*400
nsamp = int(t*1000/tsamp)
dataset = (np.random.randn(nchan, nsamp+99)*18 + 128).astype(np.uint8)
frb = np.zeros((nchan, nsamp))
times = np.arange(100)
frbconvolvemap=np.zeros((nchan,nsamp+99))
tblock=nsamp+99

##### these are model specific parameters

#dm = 1000  ### start dm
dm = values.dm
fluence = values.fluence  ###jy ms
fluence_fac = 18./2.  ### jansky  to units of sigma? which is 18 digits(1sigma) per 2 jansky
flu = fluence*fluence_fac/18.   ###fluence  units of jansky
beamno=-1 ##########fake beam identity signature
#dm = np.random.uniform()*2000 + 100
#widthms = np.random.rand()*1 + 0.064 +1
#widthms=1.0
widthms=values.width
widthsamp = widthms/tsamp
#fluence= np.random.rand()*40 +10
#flu=fluence*fluence_fac
#toffset=50
#tinsert=(200/tsamp)
toffset=100/tsamp  ###500 time sample
###real sample time=(1+i)*nsamp+toffset
#d = file.readthismeanyseconds(t)
#print t,toffset/tsamp,widthms,dm
width2=widthsamp**2
#print times
boxcar=int(widthsamp)+1
'''
widthms = np.random.rand()*5 + 0.064
widthms=20
widthsamp = widthms/tsamp
toffset=50
tinsert=50
#d = file.readthismeanyseconds(t)

print t,widthms,dm

width2=widthsamp**2
#print times
#x=np.array([0,50,100,50,0])
x = gauss_amp*np.exp(-(times-50)**2/(2*width2))
#print x
#pylab.figure()
#pylab.plot(x)
#pylab.show()
'''
print 'start'
print nsamp,flu,tblock
dataset.T.tofile(mkout.fin)  #### print first timeset of noise
#snr=50.0

#if values.menu == 0:
for i in xrange(6):
    dataset = (np.random.randn(nchan, nsamp+99) + 0)  #reset noise
    snr=0.
    snr_sig=0.
    snr_sig=336.*flu/np.sqrt(nchan)
    snr_bkg=0.
    sampn=0
    frb = np.zeros((nchan, nsamp))
    frbconvolvemap=np.zeros((nchan,nsamp+99))
    normmap=np.zeros((nchan,nsamp+99))
    #x=np.array([0,50,100,50,0])
    if width2:
        x = gauss_amp*np.exp(-(times-50)**2/(2*width2))
    else:
        x =1
    idt = abs(4.15*dm*(froof**-2 - (froof+336*foff)**-2)/tsamp)
    frbconvolvemap, normmap,dedisp = injector(frb,x,frbconvolvemap,normmap,toffset,nchan,tsamp,foff,froof,dm,amplitude,flu,width2,nsamp)
    #print i,t,toffset*tsamp,widthms,dm,flu
     # positions mask of burst on normmap
    datasetsum=(dataset.astype(float)+frbconvolvemap)*18+128
    #print dataset
    sat_mask=(datasetsum >= 255)
    datasetsum[sat_mask]=255
    dataset1= datasetsum.astype(np.uint8)
    dataset1.T.tofile(mkout.fin)
    d=dedisp
    #print('nosquare')
    #### snr calculation
    mask=(d>0.5)


    a=datasetsum
    pulse=normmap>0

    npos=a*pulse > 137
    #print(np.mean(abs(v1)),np.std(v1))

    snr = d[d > 0.5].sum()/np.sqrt((d > 0.5).sum())
    a=np.where(pulse[-1])
    tend=tblock*(i+1)+a[0][len(a[0])/2]
    print snr,tend
    #print snr_cal(dataset,normmap,frbconvolvemap,336)
    #print sn,np.sqrt(nchan)

    #print np.max(frbconvolvemap[0])
    # S/N, sampno, secs from file start, boxcar, idt, dm, beamno
    #f.write("%f %d %f %d %d %f %d\n"%(snr,toffset+(i+1)*(tblock)+idt,(toffset+(i+1)*nsamp+idt)*tsamp/1000,boxcar,idt,dm,beamno))
    #print("%d fluence=%f snr=%f samp=%d time=%f idt=%d dm=%d widthsamp=%f \n"%(i,fluence,snr,toffset+(i+1)*tblock+idt,(toffset+(i+1)*nsamp+idt)*tsamp/1000,idt,dm,widthsamp))
    f.write("%f %d %f %d %d %f %d\n"%(snr,tend,tend*tsamp/1000,boxcar,idt,dm,beamno))
    if fprint:
        print("%d fluence=%f snr=%f samp=%d time=%f idt=%d dm=%d widthsamp=%f \n"%(i,fluence,snr,tend,tend*tsamp/1000,idt,dm,widthsamp))

'''
if values.menu ==1 :
    for i in range(1000):
        dataset = (np.random.randn(nchan, nsamp+99)*18 + 128).astype(np.uint8)  #reset noise
        snr=0
        snr_sig=0
        snr_sig=336*flu/np.sqrt(nchan)
        snr_bkg=0
        frb = np.zeros((nchan, nsamp))
        frbconvolvemap=np.zeros((nchan,nsamp+99))
        normmap=np.zeros((nchan,nsamp+99))
        #dm = np.random.uniform()*2000 + 100
        widthms = np.random.rand()*1 + 0.064 +1
        #widthms=1.0
        widthms=tsamp
        widthsamp = widthms/tsamp
        #fluence= np.random.rand()*40 +10
        #flu=fluence*fluence_fac
        #toffset=50
        #tinsert=(200/tsamp)
        toffset=100/tsamp  ###500 time sample
        ###real sample time=(1+i)*nsamp+toffset
        #d = file.readthismeanyseconds(t)
        #print t,toffset/tsamp,widthms,dm
        width2=widthsamp**2
        #print times
        boxcar=int(widthsamp)+1

        #x=np.array([0,50,100,50,0])
        x = gauss_amp*np.exp(-(times-50)**2/(2*width2))
        idt = abs(4.15*dm*(froof**-2 - (froof+336*foff)**-2)/tsamp)
        injector(frb,x,frbconvolvemap,normmap,toffset,nchan,tsamp,foff,froof,dm,amplitude,flu,width2)
        #print i,t,toffset*tsamp,widthms,dm,flu
         # positions mask of burst on normmap
        datasetsum=dataset.astype(float)+frbconvolvemap
        sat_mask=(datasetsum >= 255)
        datasetsum[sat_mask]=255
        dataset1= datasetsum.astype(np.uint8)
        dataset1.T.tofile(mkout.fin)

        #### snr calculation
        for j in range(336):
            pulse=frbconvolvemap[j] >0
            a=(dataset+frbconvolvemap)[j]
            npos=a[pulse] > 137
            snr_sig+=sum(frbconvolvemap[j][pulse])
            snr_bkg+=sum((dataset[j][pulse]*npos))


        snr = snr_sig/np.sqrt(snr_bkg)
        print snr
        #print snr_cal(dataset,normmap,frbconvolvemap,336)
        #print sn,np.sqrt(nchan)

        #print np.max(frbconvolvemap[0])
        # S/N, sampno, secs from file start, boxcar, idt, dm, beamno
        f.write("%f %d %f %d %d %f %d\n"%(snr,toffset+(i+1)*nsamp+idt,(toffset+(i+1)*nsamp+idt)*tsamp/1000,boxcar,idt,dm,beamno))
        print("fluence=%f snr=%f samp=%d time=%f idt=%d dm=%d widthsamp=%f \n"%(fluence,snr,toffset+(i+1)*nsamp+idt,(toffset+(i+1)*nsamp+idt)*tsamp/1000,idt,dm,widthsamp))
'''

f.close()
mkout.fin.flush()
mkout.fin.close()
if values.show:
    pylab.figure()
    pylab.imshow(frbconvolvemap)
    pylab.figure()
    pylab.imshow(datasetsum)
    pylab.show()
    pylab.close()
