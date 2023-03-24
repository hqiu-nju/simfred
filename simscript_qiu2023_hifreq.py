import dynspec
import matplotlib.pyplot as plt
import numpy as np
import math as m

## this script generates 1 pulse for each parameter
model=dynspec.spectra(fch1=1100+336,nchan=336,bwchan=-1,tsamp=1)
mode='single'
testname=f"hifreq_{mode}"
w=open(f"{testname}.txt",'w')
### create file

printloop=0
print("starting injection\n")
for i in np.arange(0.5,11.1,0.5):  ### intrinsic standard deviation sigma
    for j in np.arange(0,3001,50):  ### DM
        model.create_filterbank(f"{testname}_dm{np.round(j,0)}_width{np.round(i,1)}",std=18,base=127)
        print(f"created file {testname}_dm{np.round(j,0)}_width{np.round(i,1)}", end = "\r")
        # w=open(f"{testname}_dm{np.round(j,0)}_width{np.round(i,1).txt",'w')
        # print (f"make DM{i} width{j}\n")
        xset=np.random.rand()-0.5
        model.writenoise(nsamp=20000)
        model.writenoise(nsamp=20000)
        base1,base2=model.burst(dm=j,A=50,width=i,mode=mode,nsamp=10000,offset=xset)
        # print(model.L2_snr())
        # print(i)
        # print(model.L2_snr()[0][:-2]+";"+str(dynspec.L2_snr(base2/model.L2_snr()[1]*50))+"\n")
        for printloop in range(50):  ### how many pulses in the data
            model.writenoise(nsamp=20000)
            model.inject(base1/model.write_flux()*50)
            w.write(model.write_snr()[0][:-2]+";"+str(dynspec.L2_snr(base2/model.write_snr()[1]*50))+f";{xset}"+"\n")
        # print(f" {np.round(i,1)}/11.1 completed", end = "\r")
            model.writenoise(nsamp=20000)
        model.writenoise(nsamp=20000) ## write noise
        model.closefile()
w.close()
# print(f" {i}/2000 of this file completed")
print("\nFinished ")
