import dynspec
import matplotlib.pyplot as plt
import numpy as np
import math as m


model=dynspec.spectra(fch1=1100,nchan=336,bwchan=1,tsamp=1)



for i in np.arange(100,2000,50):
    model.create_filterbank(f"test_DM{i}_width1",std=18,base=127)
    model.writenoise()
    base1,base2=model.burst(dm=i,A=50,width=j,mode='single',nsamp=10000)
    for printloop in range(100):
        imprint=model.imprint(base1)
    model.writenoise()
    model.closefile()
