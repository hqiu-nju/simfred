import dynspec
import matplotlib.pyplot as plt
import numpy as np
import math as m


model=dynspec.spectra(fch1=1100,nchan=336,bwchan=1,tsamp=1)


w=open(f"test_DM_width.txt",'w')
for i in np.arange(100,200,50):
    ### create file
    model.create_filterbank(f"test_DM{i}_width",std=18,base=127)

    model.writenoise() ## write noise
    base1,base2=model.burst(dm=i,A=50,width=1,mode='single',nsamp=10000)
    w.write(model.L2_snr())

    for printloop in range(5):
        imprint=model.imprint(base1)
    model.writenoise()
    model.closefile()
w.close()
