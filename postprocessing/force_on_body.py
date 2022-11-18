import numpy as np
import json
import matplotlib.pyplot as plt

num_seg = 110  # can be 55,110,220
run_time=30.0  #s [5,25,50,75,100,150,200,250,300]       # Total simulation run time, unit [s]   
dt = 1.0e-3     # Time steps, unit [s]


BP1 = num_seg                       # Body attached point of fairlead 1 location
ML1a = BP1 - 1                      # 1st point on mooring line 1
ML1b = BP1 + num_seg                # 1st point on mooring line 2
BP2 = num_seg + (1*num_seg*2) + 1   # Body attached point of fairlead 2 location
ML2a = BP2 - 1                      # 1st point on mooring line 3
ML2b = BP2 + num_seg                # 1st point on mooring line 4
BP3 = num_seg + (2*num_seg*2) + 2   # Body attached point of fairlead 3 location
ML3a = BP3 - 1                      # 1st point on mooring line 5
ML3b = BP3 + num_seg                # 1st point on mooring line 6
BP4 = num_seg + (3*num_seg*2) + 3   # Body attached point of fairlead 4 location
ML4a = BP4 - 1                      # 1st point on mooring line 7
ML4b = BP4 + num_seg                # 1st point on mooring line 8

BP_14=[BP1,BP2,BP3,BP4]
ML_18=[ML1a,ML1b,ML2a,ML2b,ML3a,ML3b,ML4a,ML4b]  

# force= k_fiber*(np.linalg.norm(position[BP_14[j//2]]-position[ML_18[j]])-seg_length)    # Calculate force magnitude by mooring line 1
# forceVector = np.array(position[BP_14[j//2]])-np.array(position[ML_18[j]])              # Calculate line vector of BP1 to ML1
# force_unit_vector=forceVector/np.linalg.norm(forceVector)
# all_forces['FV'+str(j)].append(force*force_unit_vector)

# resu=np.load("New_w"+str(num_sub_step)+"_T"+str(run_time)+"s_dt"+str(dt) +"s_DL"+str(int(1100/num_seg))+ "results.npy",allow_pickle=True)
# f=open("New_w"+str(num_sub_step)+"_T"+str(run_time)+"s_dt"+str(dt) +"s_DL"+str(int(1100/num_seg))+ "results.json")
i=0
time_list=[]
while round(dt*i,2)<run_time:
    time_list.append(round(dt*i,2))
    i+=10


for w in [1,2,5,10,15,20,25,50]:
    f=open("results/New_w"+str(w)+"_T30.0s_dt0.001s_DL5results.json")
    resu=json.load(f)
    f.close()
    print("read finish")
    v_max_mean_min_std=[]

    for item in time_list:
        v_node=np.linalg.norm(np.array(resu["velocity"+str(item)]),axis=1)
        v_max_mean_min_std.append([np.max(v_node),np.mean(v_node),np.min(v_node),np.std(v_node)])
    np.savetxt("temp_w"+str(w)+".txt",np.array(v_max_mean_min_std))
    print("save finish")

plt.figure()
plt.plot(time_list,np.array(v_max_mean_min_std))
plt.yscale('log')
plt.show()