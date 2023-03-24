import dynspec
import matplotlib.pyplot as plt
import numpy as np
import math as m






def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-l', '--label',type=str, default="test",help='Output File Name, no suffix')
    parser.add_argument('-m', '--mode',type=str, default="single",help='Injection modes: single, scat, boxcar')
    parser.add_argument('--snmode',type=str, default="fluence",help='Injection : fluence or snr')
    parser.add_argument('-A','--amplitude',type=float, default=50,help='Injection : fluence or snr')
    parser.add_argument('-t','--tbin',type=int, default=10,help='time samples per bin during simulation')
    parser.add_argument('-f','--fbin',type=int, default=10,help='freq channels per bin during simulation')
    parser.add_argument('-s','--samples',type=int, default=20000,help='injection block sample length')
    parser.add_argument('--nchan',type=int, default=336,help='number of channels')
    parser.add_argument('--tsamp',type=int, default=1,help='time resolution (ms)')
    parser.add_argument('--fch1',type=int, default=1100,help='first channel center freq')
    parser.add_argument('--bwchan',type=int, default=1,help='injection block sample length')
    parser.add_argument('-N','--npulse',type=int, default=100,help='number of pulses to inject')
    parser.add_argument('--dm_start',type=float, default=0,help='min dm (pc cm-3)')
    parser.add_argument('--step',type=float, default=50,help='step dm (pc cm-3)')
    parser.add_argument('--dm',type=float, default=2000,help='max dm (pc cm-3)')
    parser.add_argument('--sig_start',type=float, default=0.5,help='starting pulse width sigma (ms)')
    parser.add_argument('--sig_step',type=float, default=0.5,help='starting pulse width sigma (ms)')
    parser.add_argument('--sig',type=float, default=0.5,help='max pulse width sigma (ms)')
    values = parser.parse_args()
    sigmarange=np.arange(values.sig_start,values.sig+0.5*values.sig_step,values.sig_step)
    dmrange=np.arange(values.dm_start,values.dm+0.5*values.step,values.step)
    tbin=values.tbin
    fbin.values.fbin
    fch1=values.fch1
    bwchan=values.bwchan
    nchan=values.nchan
    tsamp=values.tsamp
    nsamp=values.samples
    mode=values.mode
    label=values.label
    npulse=values.npulse
    snr=values.amplitude
    if values.snmode == 'fluence':
        fluencebatch(fch1,bwchan,nchan,tsamp,mode,label,nsamp,npulse,sigmarange,dmrange,tbin,fbin)
    elif values.snmode == 'snr':
        snrbatch(fch1,bwchan,nchan,tsamp,mode,label,nsamp,npulse,sigmarange,dmrange,tbin,fbin)

def fluencebatch(fch1,bwchan,nchan,tsamp,mode,label,nsamp,npulse,sigmarange,dmrange,tbin,fbin):
    ## this script generates 1 pulse for each parameter
    model=dynspec.spectra(fch1=fch1,nchan=nchan,bwchan=bwchan,tsamp=tsamp,tbin=tbin,fbin=fbin)
    testname=f"{label}_{mode}"
    w=open(f"{testname}.txt",'w')
    ### create file
    printloop=0
    print("starting injection\n")
    for i in sigmarange:  ### intrinsic standard deviation sigma
        for j in dmrange:  ### DM
            model.create_filterbank(f"{testname}_dm{np.round(j,0)}_width{np.round(i,1)}",std=18,base=127)
            print(f"created file {testname}_dm{np.round(j,0)}_width{np.round(i,1)}", end = "\r")
            # w=open(f"{testname}_dm{np.round(j,0)}_width{np.round(i,1).txt",'w')
            # print (f"make DM{i} width{j}\n")
            xset=np.random.rand()-0.5
            model.writenoise(nsamp=nsamp)
            model.writenoise(nsamp=nsamp)
            base1,base2=model.burst(dm=j,A=50,width=i,mode=mode,nsamp=nsamp,offset=xset)
            # print(model.L2_snr())
            # print(i)
            # print(model.L2_snr()[0][:-2]+";"+str(dynspec.L2_snr(base2/model.L2_snr()[1]*50))+"\n")
            for printloop in range(npulse):  ### how many pulses in the data
                model.writenoise(nsamp=nsamp)
                model.inject(base1/model.write_flux()*50)
                w.write(model.write_snr()[0][:-2]+";"+str(dynspec.L2_snr(base2/model.write_snr()[1]*50))+f";{xset}"+"\n")
            # print(f" {np.round(i,1)}/11.1 completed", end = "\r")
                model.writenoise(nsamp=nsamp)
            model.writenoise(nsamp=nsamp) ## write noise
            model.closefile()
    w.close()
    # print(f" {i}/2000 of this file completed")
    print("\nFinished ")

def snrbatch(fch1,bwchan,nchan,tsamp,mode,label,nsamp,npulse,sigmarange,dmrange):
    ## this script generates 1 pulse for each parameter
    model=dynspec.spectra(fch1=fch1,nchan=nchan,bwchan=bwchan,tsamp=tsamp)
    testname=f"{label}_{mode}"
    w=open(f"{testname}.txt",'w')
    ### create file
    printloop=0
    print("starting injection\n")
    for i in sigmarange:  ### intrinsic standard deviation sigma
        for j in dmrange:  ### DM
            model.create_filterbank(f"{testname}_dm{np.round(j,0)}_width{np.round(i,1)}",std=18,base=127)
            print(f"created file {testname}_dm{np.round(j,0)}_width{np.round(i,1)}", end = "\r")
            # w=open(f"{testname}_dm{np.round(j,0)}_width{np.round(i,1).txt",'w')
            # print (f"make DM{i} width{j}\n")
            xset=np.random.rand()-0.5
            model.writenoise(nsamp=20000)
            model.writenoise(nsamp=20000)
            base1,base2=model.burst(dm=j,A=50,width=i,mode=mode,nsamp=20000,offset=xset)
            # print(model.L2_snr())
            # print(i)
            # print(model.L2_snr()[0][:-2]+";"+str(dynspec.L2_snr(base2/model.L2_snr()[1]*50))+"\n")
            for printloop in range(50):  ### how many pulses in the data
                model.writenoise(nsamp=20000)
                model.inject(base1/model.write_snr()*50)
                w.write(model.write_snr()[0][:-2]+";"+str(dynspec.L2_snr(base2/model.write_snr()[1]*50))+f";{xset}"+"\n")
            # print(f" {np.round(i,1)}/11.1 completed", end = "\r")
                model.writenoise(nsamp=20000)
            model.writenoise(nsamp=20000) ## write noise
            model.closefile()
    w.close()
    # print(f" {i}/2000 of this file completed")
    print("\nFinished ")
