import numpy as np
import src.geoMaker.oceanFarm_mooring as geo
import src.visualization.saveVtk as sv
import src.element.line as l

num_seg = 220  # can be 55,110,220
run_time=250  #s [5,25,50,75,100,150,200,250,300]       # Total simulation run time, unit [s]   
dt = 1.0e-3     # Time steps, unit [s]

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

for i in range(int(run_time/dt)):       
    if i % round(float(0.04)/float(dt)) == 0:                # Write vtk result per 0.04s, 25fps
        sv.write_line_vtk("ami2/New_T"+str(run_time)+"s_dt"+str(dt) +"s_DL"+str(int(1100/num_seg))+ "m_mooring_line"+str(i),point=position.tolist(),line=line)
    if i % round(float(0.01)/float(dt)) == 0:                # Write force result per 0.01s
        for j in range(8):
            force= k_fiber*(np.linalg.norm(position[BP_14[j//2]]-position[ML_18[j]])-seg_length)    # Calculate force magnitude by mooring line 1
            forceVector = np.array(position[BP_14[j//2]])-np.array(position[ML_18[j]])              # Calculate line vector of BP1 to ML1
            force_unit_vector=forceVector/np.linalg.norm(forceVector)
            all_forces['FV'+str(j)].append((force*force_unit_vector).tolist())

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
    
    ### constraint function 
    position+=chain_line.pbd_edge_constraint(position,mass_matrix,dt)
    position+=fiber_line.pbd_edge_constraint(position,mass_matrix,dt)
    
    ### velocity correction
    velocity=(position-pre_position)/dt  

np.save("ami2/New_T"+str(run_time)+"s_dt"+str(dt) +"s_DL"+str(int(1100/num_seg))+ "m_all_force.npy",all_forces)
np.savetxt("ami2/New_T"+str(run_time)+"s_dt"+str(dt) +"s_DL"+str(int(1100/num_seg))+ "m_final_position.txt",position)
