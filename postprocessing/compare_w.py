import numpy as np
import matplotlib.pyplot as plt

num_seg = 110  # can be 55,110,220
run_time=60.0  #s [5,25,50,75,100,150,200,250,300]       # Total simulation run time, unit [s]   
dt = 1.0e-2     # Time steps, unit [s]

i=0
time_list=[]
while round(dt*i,2)<run_time:
    time_list.append(round(dt*i,2))
    i+=40

    


plt.figure()

for w in [2,10,100,200]:
    data=np.loadtxt(f"temp_w{w}_T{run_time}s_dt{dt}s_DL5_results.txt")
    y_m=np.array(data)[1:,1]
    y_e=np.array(data)[1:,3]/2
    plt.plot(time_list[1:],y_m,label='Num_sub_step='+str(w))
    plt.fill_between(time_list[1:],y_m-y_e,y_m+y_e,alpha=0.2)
    
plt.legend()
plt.xlabel("Time (s)")
plt.ylabel("Mean nodal velocity (m/s)")
plt.yscale('log')
plt.savefig('Effect of num_sub_iteration.jpg',dpi=1000)
plt.show()