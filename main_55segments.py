import numpy as np
import src.geoMaker.oceanFarm1_55 as geo
import src.visualization.saveVtk as sv
import src.element.line as l
import src.element.quad as q
import math as m
import matplotlib.pyplot as plt

nodes   = geo.mooring_point_new
line    = geo.mooring_line_new
face    = geo.netFace
num_seg = geo.num_seg

sv.write_vtk("initial_mooring_line",point=nodes,line=line)
sv.write_vtk("initial_cage",point=geo.nodes,line=geo.all_line,face=geo.netFace)

# define structural properties
seg_length = 1100/num_seg #total length of mooring is 1100 m 

E_fiber = 11.7e9 
A_fiber = 0.020103051
EA_fiber=235e6
k_fiber = EA_fiber/seg_length

E_chain = 60e9
A_chain = 0.012170404
EA_chain=687e6
k_chain = EA_chain/seg_length

chain_line=l.lines(geo.mooring_line_chain,k_chain,0.088)   # Axial stiffness[MN] 680.81 (Chain) 235.44 (Fiber)
fiber_line=l.lines(geo.mooring_line_fiber,k_fiber,0.160)

chain_line.assign_length(20.0)
fiber_line.assign_length(20.0)

## setting
fixed_point =  geo.fixed_point          # anchor point
fixed_point += geo.body_attached_point  # fish cage body

gravity=np.array([0,0,-9.81])
mass_matrix = np.array(geo.mass_mooring_line_new).reshape(len(geo.mass_mooring_line_new),1)
total_time = [5,25,50,75,100,150,200,250,300]       # Total simulation runn time, unit [s]

all_forces={'FV1':[],'FV2':[],'FV3':[],'FV4':[],'FV5':[],'FV6':[],'FV7':[],'FV0':[]}

## initialization 
position=np.array(nodes)
velocity=np.zeros_like(position)

BP1 = num_seg                       # Body attached point of fairlead 1 location
ML1a = BP1 - 1                      # 1st point on mooring line 1
ML1b = BP1 + num_seg                # 1st point on mooring line 2
BP2 = num_seg + (1*num_seg*2) + 1   # Body attached point of fairlead 2 location
ML2a = BP1 - 1                      # 1st point on mooring line 3
ML2b = BP1 + num_seg                # 1st point on mooring line 4
BP3 = num_seg + (2*num_seg*2) + 2   # Body attached point of fairlead 3 location
ML3a = BP1 - 1                      # 1st point on mooring line 5
ML3b = BP1 + num_seg                # 1st point on mooring line 6
BP4 = num_seg + (3*num_seg*2) + 3   # Body attached point of fairlead 4 location
ML4a = BP1 - 1                      # 1st point on mooring line 7
ML4b = BP1 + num_seg                # 1st point on mooring line 8

BP_14=[BP1,BP2,BP3,BP4]
ML_18=[ML1a,ML1b,ML2a,ML2b,ML3a,ML3b,ML4a,ML4b]
    
run_time=total_time[6]
dt = 1.0e-4                         # Time steps, unit [s]

for i in range(int(run_time/dt)):       
    if i % 400 == 0:                # Write vtk result per 0.04s, 25fps
        sv.write_line_vtk("ami2/"+str(run_time)+"W"+str(dt) + "mooring_line"+str(i),point=position.tolist(),line=line)
    if i % 100 == 0:                # Write force result per 0.01s
        for j in range(8):
            force= k_fiber*(np.linalg.norm(position[BP_14[j//2]]-position[ML_18[j]])-seg_length)    # Calculate force magnitude by mooring line 1
            forceVector = np.array(position[BP_14[j//2]])-np.array(position[ML_18[j]])              # Calculate line vector of BP1 to ML1
            force_unit_vector=forceVector/np.linalg.norm(forceVector)
            all_forces['FV'+str(j)].append(force*force_unit_vector)

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
    

np.save('all_force.npy',all_forces)
np.savetxt('position.txt',position)
