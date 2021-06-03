import dynspec
import matplotlib.pyplot as plt
import numpy as np
import math as m

### script for injecting a hundred pulses for each parameter, use for final simulation results
model=dynspec.spectra(fch1=1100,nchan=336,bwchan=1,tsamp=1)


w=open(f"test_DM_width.txt",'w')
for i in np.arange(10,2001,10):
    ### create file
    model.create_filterbank(f"test_DM{i}_width1",std=18,base=127)
    print (f"test_DM{i}_width1")

    model.writenoise() ## write noise
    base1,base2=model.burst(dm=i,A=50,width=1,mode='single',nsamp=10000)
    w.write(model.L2_snr())
    print("starting injection")
    for printloop in range(100):
        print(f"{printloop}% of this file completed", end = "\r")
        imprint=model.imprint(base1)
    print("\nsaving final file")
    model.writenoise()
    model.closefile()
w.close()
