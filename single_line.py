# this case is to validiat with Orcaflex and sima.
import sys
import numpy as np
import json
import src.geoMaker.oceanFarm_mooring as geo
import src.visualization.saveVtk as sv
import src.element.line as l


seg_length =5.0  #total length of mooring is 1100 m 
num_seg = int(1100/seg_length)  # can be 220,
run_time=60.0  #s [5,25,50,75,100,150,200,250,300]       # Total simulation run time, unit [s]   
dt = 1.0e-2     # Time steps, unit [s]
num_sub_step=sys.argv[1]

nodes=np.linspace([-1070,0,-150],[0,0,0],num_seg+1)
chain_line_con=[]
for i in range(num_seg):
    chain_line_con.append([i,i+1])
mass= np.array([[147*1100/len(nodes)]]*len(nodes))
fixed_point=[0,num_seg]
# define structural properties


EA_chain= 7e8
k_chain = EA_chain/seg_length

chain_line=l.lines(chain_line_con,k_chain,0.088)   # Axial stiffness[MN] 680.81 (Chain) 235.44 (Fiber)

chain_line.assign_length(seg_length)

## setting
gravity=np.array([0,0,-9.81])

## initialization 
position=np.array(nodes)
velocity=np.zeros_like(position)
results={}

i=0
while round(dt*i,2)<run_time:
# for i in range(int(run_time/dt)):       
    if i % round(float(0.4)/float(dt)) == 0:                # Write vtk result per 0.04s, 25fps
        print("save results at {:.2f} s".format(round(dt*i,2)))
        sv.write_line_vtk("ami2/3single_w"+str(num_sub_step)+"_T"+str(run_time)+"s_dt"+str(dt) +"s_DL"+str(int(1100/num_seg))+ "m"+str(i),
                          point=position.tolist(),line=chain_line_con)
        results["position"+str(round(dt*i,2))]=position.tolist()
        results["velocity"+str(round(dt*i,2))]=velocity.tolist()
        print(np.sum(chain_line.line_length))

    ### Forward Euler (Explicit)
    ## External loads
    pre_position=position.copy()
    
    # Gravity force
    velocity += dt*gravity *min(1,i/1e4)
    
    ## boundary condition
    velocity[fixed_point] *= 0.0  # velocity restriction
    velocity[position[:,2]<-150]*=np.array([1,1,0])# ground
    position += dt*velocity
    position[fixed_point]=np.array(nodes)[fixed_point]
    
    for w in range(int(num_sub_step)):
    ### constraint function 
        position+=chain_line.pbd_edge_constraint(position,mass,dt)
        
    
    ### velocity correction
    velocity=(position-pre_position)/dt  
    i+=1
    
    
json=json.dumps(results)    
# open file for writing, "w" 
f = open("results/single_w"+str(num_sub_step)+"_T"+str(run_time)+"s_dt"+str(dt) +"s_DL"+str(int(1100/num_seg))+ "results.json","w")
# write json object to file
f.write(json)
# close file
f.close()
