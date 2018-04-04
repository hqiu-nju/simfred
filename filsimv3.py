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

autoheader={'az_start': 0.0,
 'barycentric': None,
 'data_type': 1,
 'fch1': 1464.0,
 'fchannel': None,
 'foff': -1.0,
 'machine_id': 0,
 'nbits': 8,
 'nchans': 336,
 'nifs': 1,
 'nsamples': None,
 'period': None,
 'pulsarcentric': None,
 'rawdatafile': None,
 'refdm': None,
 'source_name': 'None',
 #'src_dej': -32441.084833752,
 #'src_raj': 215344.63079648,
 'src_raj':174540.1662,
 'src_dej':-290029.896,
 'telescope_id': 7,
 'tsamp': 0.00126646875,
 'tstart': 57946.52703893818,
 'za_start': 0.0}


def injector(frb,x,frbconvolvemap,normmap,tstart,nchan,tsamp,foff,froof,dm,amplitude,flu,w2,offset,diffsamp,index,ind_fac):
    toffset=int(tstart)+offset
    for c in xrange(nchan):
        ceil = froof + (c)*foff
        floor = froof + (c+1)*foff ###frequency of current channel
        roof_dispersion_delay_ms = 4.15*dm*(froof**-2 - ceil**-2)
        bott_dispersion_delay_ms = 4.15*dm*(froof**-2 - floor**-2)
        roof_dispersion_delay_samp = abs(roof_dispersion_delay_ms/tsamp)+toffset
        bott_dispersion_delay_samp = abs(bott_dispersion_delay_ms/tsamp)+toffset # offset by a bit so the first few samples aren't off the end of the thing
        length_scale = abs(int(roof_dispersion_delay_samp)-int(bott_dispersion_delay_samp))
        smearing_width_samp = abs((roof_dispersion_delay_samp)-(bott_dispersion_delay_samp))
        #noise_start=np.round(roof_dispersion_delay_samp)
        #noise_end=np.round(bott_dispersion_delay_samp)
        total_width2 = smearing_width_samp**2 + w2
        #ac = flu/np.sqrt(total_width2)
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
            convolved=np.append(frb[c],np.zeros(diffsamp))
        normfac=np.sum(convolved)
        normmap[c]+=convolved/normfac
        if index:
            frbconvolvemap[c]+=normmap[c]*flu*((ceil/froof)**ind_fac)
            normmap[c]=normmap[c]*((ceil/froof)**ind_fac)
        else:
            frbconvolvemap[c]+=normmap[c]*flu
        boxcar=np.sqrt(total_width2)+1
    return frbconvolvemap,normmap,boxcar


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
parser.add_argument('-B','--basefile', action='store_true', help='base filterbank')
parser.add_argument('-o', '--output',type=str, default='set_'+date)
parser.add_argument('-s','--show', action='store_true', help='Show')
parser.add_argument('--nchan',type=int,default=336)
parser.add_argument('-a', '--snr',action='store_true', help='vary snr')
parser.add_argument('-i', '--index',action='store_true', help='spectral index option')
parser.add_argument('-z', '--zero',action='store_true', help='zero noise')
parser.add_argument('-A', '--snfac',type=float, default=10)
parser.add_argument('-I', '--indexfac',type=float, default=0)
parser.add_argument('-x','--offset',type=float,default=0.5, help='Offset within sample')
parser.add_argument('-N','--number',type=int,default=10, help='how many FRBs')
#parser.add_argument(dest='files', nargs='+')
parser.set_defaults(verbose=False)
values = parser.parse_args()
if values.snr:
    mxsnr=values.snfac
if values.verbose:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)
ch_index=values.index
ch_fac=values.indexfac
#date=time.strftime('%Y_%m_%d_%H_%M',time.localtime())
#date='2000'
inputname=values.base
outputname=values.output+".fil"
f=open(outputname[:-4]+'.candlist','w')
f.write("##new file of truth of truth!!!!!"+outputname+"\n")
f.write("# S/N, sampno, secs from file start, boxcar, idt, dm, beamno, fluence, width, offset, spectral index\n")
fprint=values.printer
if values.basefile:
    readin=sgp.SigprocFile(inputname)
    readin.header['nchans']=values.nchan
    hdr_dic=readin.header
    readin.seek_data()
else:
    autoheader['nchans']=values.nchan
    hdr_dic=autoheader

mkout=sgp.SigprocFile(outputname,'w',hdr_dic)
mkout.seek_data()

#mkout.seek_data()
#readin.seek_data()
#### these are system/constant parameters

nchan = hdr_dic['nchans']# channels
fch1 = 1.464 # GHz
foff =  -1/1e3 # GHz
froof = fch1 # GHz need change to upper limit
tsamp = float(hdr_dic['tsamp'])*1000
#print (tsamp)# milliseconds
amplitude=1.
gauss_amp=amplitude # this has to be the same to achieve norm maximum?
t = 5 # random seconds between 5 and 10 set at 5, currently 5 seconds
#dm = np.random.rand()*400
nsamp = int(t*1000/tsamp)
if values.dm >2000:
    t=7
    nsamp = int(t*1000/tsamp)
if values.dm > 4000:
    tes=4.15*values.dm*((1.464-0.336)**-2 - 1.464**-2) #ms
    t=tes+2000
    nsamp = int(t/tsamp)

dataset = (np.random.randn(nchan, nsamp)*18 + 128).astype(np.uint8)
if values.zero:
    dataset= (np.random.randn(nchan, nsamp)*0 + 128).astype(np.uint8)
#frb = np.zeros((nchan, nsamp))
#frbconvolvemap=np.zeros((nchan,nsamp))
tblock=nsamp

##### these are model specific parameters

#dm = 1000  ### start dm
dm = values.dm
fluence = values.fluence  ###jy ms
fluence_fac = 18./2.  ### jansky  to units of sigma? which is 18 digits(1sigma) per 2 jansky
flu = fluence*fluence_fac/18.   ###fluence  units of sigma
beamno=-1 ##########fake beam identity signature
#dm = np.random.uniform()*2000 + 100
#widthms = np.random.rand()*1 + 0.064 +1
#widthms=1.0
widthms=values.width
widthsamp = widthms/tsamp
times = np.arange(int(widthsamp)*6+1)
#fluence= np.random.rand()*40 +10
#flu=fluence*fluence_fac
#toffset=50
#tinsert=(200/tsamp)
toffset=100/tsamp  ###100 ms by timesample units
###real sample time=(1+i)*nsamp+toffset
#d = file.readthismeanyseconds(t)
#print t,toffset/tsamp,widthms,dm
width2=widthsamp**2
#print times

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
if width2:
    x = gauss_amp*np.exp(-(times-int(widthsamp)*3)**2/(2*width2))
else:
    x =1
#print 'start'
#print nsamp,flu,tblock
dataset.T.tofile(mkout.fin)  #### print first timeset of noise
#snr=50.0
xoff=values.offset
realsamp=nsamp+int(widthsamp)*6
diffsamp=int(widthsamp)*6
#if values.menu == 0:
for i in xrange(values.number):
    dataset = (np.random.randn(nchan, realsamp) + 0)  #reset noise
    snr=0.
    snr_sig=0.
    snr_sig=336.*flu/np.sqrt(nchan)
    snr_bkg=0.
    sampn=0
    frb = np.zeros((nchan, nsamp))
    frbconvolvemap=np.zeros((nchan,realsamp))
    normmap=np.zeros((nchan,realsamp))
    #x=np.array([0,50,100,50,0])
    idt = abs(4.15*dm*(froof**-2 - (froof+336*foff)**-2)/tsamp)
    if ch_index:
        wrt_index=ch_fac
    else:
        wrt_index='0'
    frbconvolvemap, normmap,boxcar = injector(frb,x,frbconvolvemap,normmap,toffset,nchan,tsamp,foff,froof,dm,amplitude,flu,width2,xoff,diffsamp,ch_index,ch_fac)
    #print i,t,toffset*tsamp,widthms,dm,flu
    if values.snr:
        d = normmap.flatten()
        #print('nosquare')
        #### snr calculation
        pulse=normmap>0
        snr1 = d.sum()/np.sqrt((d > 0.0).sum())
        mfactor=mxsnr/snr1
        frbconvolvemap=normmap*mfactor
        d2 = frbconvolvemap.flatten()
        #print('nosquare')
        #### snr calculation
        snr = d2.sum()/np.sqrt((d2 > 0.0).sum())
        fluence=frbconvolvemap[1].sum()
    else:
        d = frbconvolvemap.flatten()
        #print('nosquare')
        #### snr calculation
        pulse=normmap>0
        snr = d.sum()/np.sqrt((d > 0.0).sum())

    #print(np.mean(abs(v1)),np.std(v1))

    # positions mask of burst on normmap
    if values.zero:
        datasetsum=(frbconvolvemap)*18+128
    else:
        datasetsum=(dataset.astype(float)+frbconvolvemap)*18+128
    #print dataset
    sat_mask=(datasetsum >= 255)
    datasetsum[sat_mask]=255
    dataset1= datasetsum.astype(np.uint8)
    dataset1.T.tofile(mkout.fin)


    #print (d>0.0).sum()
    a=np.where(pulse[-1])
    tend=nsamp+realsamp*i+a[0][len(a[0])/2]
    #print snr,tend
    #print snr_cal(dataset,normmap,frbconvolvemap,336)
    #print sn,np.sqrt(nchan)

    #print np.max(frbconvolvemap[0])
    # S/N, sampno, secs from file start, boxcar, idt, dm, beamno
    #f.write("%f %d %f %d %d %f %d\n"%(snr,toffset+(i+1)*(tblock)+idt,(toffset+(i+1)*nsamp+idt)*tsamp/1000,boxcar,idt,dm,beamno))
    #print("%d fluence=%f snr=%f samp=%d time=%f idt=%d dm=%d widthsamp=%f \n"%(i,fluence,snr,toffset+(i+1)*tblock+idt,(toffset+(i+1)*nsamp+idt)*tsamp/1000,idt,dm,widthsamp))
    f.write("%f %d %f %d %d %f %d %f %f %f %s\n"%(snr,tend,tend*tsamp/1000,boxcar+1,idt,dm,beamno,fluence,widthms,xoff,wrt_index))
    if fprint:
        print("%d fluence=%f snr=%f samp=%d time=%f idt=%d dm=%d widthsamp=%f boxcar=%f \n"%(i,fluence,snr,tend,tend*tsamp/1000,idt,dm,widthsamp,boxcar))
f.close()
#dataset = (np.random.randn(nchan, nsamp)*18 + 128).astype(np.uint8)
#dataset.T.tofile(mkout.fin)
mkout.fin.flush()
mkout.fin.close()
if values.show:
    pylab.figure(1)
    pylab.xlabel("Time Sample")
    pylab.ylabel("Channel")
    pylab.imshow(frbconvolvemap,aspect='auto')
    pylab.tight_layout()
    pylab.savefig("convolvemap.eps")
    pylab.figure(2)
    pylab.xlabel("Time Sample")
    pylab.ylabel("Channel")
    pylab.imshow(datasetsum,aspect='auto')
    pylab.tight_layout()
    pylab.savefig("data_sum.eps")
    pylab.show()
    pylab.close()
