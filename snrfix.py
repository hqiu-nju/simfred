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
    return frbconvolvemap,normmap

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

inputname=values.base
outputname=values.output
f=open(outputname[:-4]+'.candlist','w')
f.write("##new file of truth of truth!!!!!"+outputname+"\n")
f.write("# S/N, sampno, secs from file start, boxcar, idt, dm, beamno\n")
fprint=values.printer
readin=sgp.SigprocFile(inputname)
mkout=sgp.SigprocFile(outputname,'w',readin.header)
