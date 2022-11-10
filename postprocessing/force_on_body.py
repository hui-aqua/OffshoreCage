import numpy as np
import matplotlib.pyplot as plt

num_seg = 110  # can be 55,110,220
run_time=250  #s [5,25,50,75,100,150,200,250,300]       # Total simulation run time, unit [s]   
dt = 1.0e-3     # Time steps, unit [s]

f1=np.load("ami2/New_T"+str(run_time)+"s_dt"+str(dt) +"s_DL"+str(int(1100/num_seg))+ "m_all_force.npy",allow_pickle=True)
fv1=np.array(f1.item().get('FV1'))
fv1_mag=np.linalg.norm(fv1,axis=1)/1000.0

plt.figure()

plt.plot(fv1_mag)
plt.yscale("log")
plt.show()