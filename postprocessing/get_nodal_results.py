import numpy as np
import json
import matplotlib.pyplot as plt


run_time=60.0  #s [5,25,50,75,100,150,200,250,300]       # Total simulation run time, unit [s]   
dt = 1.0e-2     # Time steps, unit [s]

i=0
time_list=[]
while round(dt*i,2)<run_time:
    time_list.append(round(dt*i,2))
    i+=40


for w in [2,10,100,200]:
    f=open(f"results/2single_w{w}_T{run_time}s_dt{dt}s_DL5_results.json")
    resu=json.load(f)
    f.close()
    print("read finish")
    v_max_mean_min_std=[]

    for item in time_list:
        v_node=np.linalg.norm(np.array(resu["velocity"+str(item)]),axis=1)
        v_max_mean_min_std.append([np.max(v_node),np.mean(v_node),np.min(v_node),np.std(v_node)])
    np.savetxt(f"temp_w{w}_T{run_time}s_dt{dt}s_DL5_results.txt",np.array(v_max_mean_min_std))
    print("save finish")

plt.figure()
plt.plot(time_list,np.array(v_max_mean_min_std))
plt.yscale('log')
plt.show()