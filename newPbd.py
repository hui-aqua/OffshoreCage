import sys
import numpy as np
import json
import src.geoMaker.oceanFarm_mooring as geo
import src.visualization.saveVtk as sv
import src.element.line as l

num_seg = 110  # can be 55,110,220
run_time=30.0  #s [5,25,50,75,100,150,200,250,300]       # Total simulation run time, unit [s]   
dt = 1.0e-1     # Time steps, unit [s]
num_sub_step=sys.argv[1]

nodes,fiber_line,chain_line,mass,fixed_point=geo.generate_model(num_seg)
line=fiber_line+chain_line
# define structural properties
seg_length = 1100.0/num_seg #total length of mooring is 1100 m 

EA_fiber= 235.44e6  #  Bore and Fossan's thesis (2015)
k_fiber = EA_fiber/seg_length

EA_chain= 680.81e6
k_chain = EA_chain/seg_length

chain_line=l.lines(chain_line,k_chain,0.088)   # Axial stiffness[MN] 680.81 (Chain) 235.44 (Fiber)
fiber_line=l.lines(fiber_line,k_fiber,0.160)

chain_line.assign_length(seg_length)
fiber_line.assign_length(seg_length)

## setting
gravity=np.array([0,0,-9.81])
mass_matrix = np.array(mass).reshape(len(mass),1)

all_forces={'FV1':[],'FV2':[],'FV3':[],'FV4':[],'FV5':[],'FV6':[],'FV7':[],'FV0':[]}

## initialization 
position=np.array(nodes)
velocity=np.zeros_like(position)
results={}

i=0
while round(dt*i,2)<run_time:
# for i in range(int(run_time/dt)):       
    if i % round(float(0.4)/float(dt)) == 0:                # Write vtk result per 0.04s, 25fps
        print("save results at {:.2f} s".format(round(dt*i,2)))
        sv.write_line_vtk("ami2/New_w"+str(num_sub_step)+"_T"+str(run_time)+"s_dt"+str(dt) +"s_DL"+str(int(1100/num_seg))+ "m"+str(i),
                          point=position.tolist(),line=line)
    if i % round(float(0.1)/float(dt)) == 0:                # Write force result per 0.01s
        results["position"+str(round(dt*i,2))]=position.tolist()
        results["velocity"+str(round(dt*i,2))]=velocity.tolist()

    ### Forward Euler (Explicit)
    ## External loads
    pre_position=position.copy()
    
    # Gravity force
    velocity += dt*gravity
    
    ## boundary condition
    velocity[fixed_point] *= 0.0  # velocity restriction
    velocity[position[:,2]<-150]*=np.array([1,1,0])# ground
    position += dt*velocity
    position[fixed_point]=np.array(nodes)[fixed_point]
    
    for w in range(int(num_sub_step)):
    ### constraint function 
        position+=chain_line.pbd_edge_constraint(position,mass_matrix,dt)
        position+=fiber_line.pbd_edge_constraint(position,mass_matrix,dt)
    
    ### velocity correction
    velocity=(position-pre_position)/dt  
    i+=1
    
    
json=json.dumps(results)    
# open file for writing, "w" 
f = open("results/New_w"+str(num_sub_step)+"_T"+str(run_time)+"s_dt"+str(dt) +"s_DL"+str(int(1100/num_seg))+ "results.json","w")
# write json object to file
f.write(json)
# close file
f.close()
