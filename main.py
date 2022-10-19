import numpy as np
import pyvista as pv
import src.geoMaker.oceanFarm1 as geo
import src.visualization.saveVtk as sv
import src.element.line as l
import src.element.quad as q


nodes   = geo.mooring_point_new
line    = geo.mooring_line_new
face    = geo.netFace

sv.write_vtk("initial_mooring_line",point=nodes,line=line)
sv.write_vtk("initial_cage",point=geo.nodes,line=geo.all_line,face=geo.netFace)

# define structural properties
num_mooring_line_seg       = geo.num_seg #50
l1=l.lines(geo.mooring_line,680.81e6,0.088) # Axial stiffness[MN] 680.81 (Chain) 235.44 (Fiber)
l1.assign_length(10.0)

## setting
fixed_point=[num_mooring_line_seg*i for i in range(8)]
fixed_point+=[219,879,659,439]

gravity=np.array([0,0,-9.81])
mass_matrix = np.array(geo.mass_mooring_line).reshape(len(geo.mass_mooring_line),1)
run_time = 10  # unit [s]
dt = 2e-2    # unit [s]


## initialization 
position=np.array(nodes)
velocity=np.zeros_like(position)






for i in range(int(run_time/dt)):       
    if i % 5 == 0:

        sv.write_line_vtk("ami2/"+"mooring_line"+str(i),point=position.tolist(),line=line)
               
    ### forward Euler (explicit)
    ## external loads
    # gravity force
    pre_position=position.copy()
    velocity += gravity*dt
    
    ## boundary condition
    velocity[fixed_point] *= 0.0  # velocity restriction
    velocity[position[:,2]<-150]*=np.array([1,1,0])
    position += velocity*dt
    ### constraint function 
    position+=l1.pbd_edge_constraint(position,mass_matrix,dt)
    
    ### velocity correction
    velocity=(position-pre_position)/dt
