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

# the initial position has huge influence on the convergence
# nodes=np.linspace([-1070,0,-150],[0,0,0],num_seg+1)
z=np.sqrt(1100*1100-1070*1070)
nodes=np.linspace([-1070,0,z],[0,0,0],num_seg+1)

chain_line_con=[]
for i in range(num_seg):
    chain_line_con.append([i,i+1])
mass= np.array([[147*1100/len(nodes)]]*len(nodes))
fixed_point=[0,num_seg]
# define structural properties


EA_chain= 7e7
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
        sv.write_line_vtk(f"ami2/k1single_w{num_sub_step}_T{run_time}s_dt{dt}s_DL{int(1100/num_seg)}m{i}",
                          point=position.tolist(),line=chain_line_con)
        results[f"position{round(dt*i,2)}"]=position.tolist()
        results[f"velocity{round(dt*i,2)}"]=velocity.tolist()
        print(np.sum(chain_line.line_length))

    ### Forward Euler (Explicit)
    ## External loads
    pre_position=position.copy()
    
    # Gravity force
    velocity += dt*gravity *min(1,i/1e4)
    
    ## boundary condition
    # velocity[fixed_point] *= 0.0  # velocity restriction
    velocity[position[:,2]<-150]*=np.array([1,1,0])# ground
    position += dt*velocity
    position[-1]=np.array(nodes)[-1]
    position[0,:1]=np.array(nodes)[0,:1]
    
    for w in range(int(num_sub_step)):
    ### constraint function 
        position+=chain_line.pbd_edge_constraint(position,mass,dt)
        
    
    ### velocity correction
    velocity=(position-pre_position)/dt  
    i+=1
    
    
json=json.dumps(results)    
# open file for writing, "w" 
f = open(f"results/2single_w{num_sub_step}_T{run_time}s_dt{dt}s_DL{int(1100/num_seg)}_results.json","w")
# write json object to file
f.write(json)
# close file
f.close()
