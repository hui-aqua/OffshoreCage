import numpy as np
# import pyvista as pv
import src.geoMaker.oceanFarm1 as geo
import src.visualization.saveVtk as sv
import src.element.line as l
import src.element.quad as q

# for circular cage
# import src.case.circularCage as case
# case.main()
nodes=geo.nodes
line=geo.all_line
face=geo.netFace

# sv.write_vtk('initial',point=nodes,line=line,face=face)
sv.write_vtk("initial",point=nodes,face=geo.netFace)
sv.write_line_vtk("initial_main_frame",point=nodes,line=geo.main_frame)
sv.write_line_vtk("initial_mooring_line",point=nodes,line=geo.mooring_line)

# define structural properties
num_bodyPoint=77 #31
num_seg=110 #50

m_body=5600000 #[kg] 4980760 
l1=l.lines(geo.mooring_line,680e6,0.088) # Axial stiffness[MN] 680.81 (Chain) 235.44 (Fiber)
l1.assign_length("10")

# Main frame properties
l2=l.lines(geo.top_cross_beam,10e9,2.05)
l3=l.lines(geo.top_hor_beam,10e9,2.29)
l4=l.lines(geo.mid_hor_beam,10e9,1)
l5=l.lines(geo.bottom_hor_beam,10e9,2.05)
l6=l.lines(geo.dia_beam,10e9,1)
l7=l.lines(geo.side_column_wo_pontoon,10e9,2.8)
l8=l.lines(geo.bottom_rad_beam,10e9,1.75)
l9=l.lines(geo.side_column_w_pontoon,10e9,3.56)
l10=l.lines(geo.pontoon_cylinder,10e9,12)
l11=l.lines(geo.center_column,10e9,3.56)
l12=l.lines(geo.center_pontoon_cylinder,10e9,17)


# 
fixed_point=[num_bodyPoint+num_seg*i for i in range(8)]

xyz=np.array(nodes)
dxyz=np.zeros_like(xyz)
gravity=np.array([0,0,-9.81])


run_time = 10  # unit [s]
dt = 2e-2    # unit [s]

# forward Euler (explicit)
for i in range(int(run_time/dt)):       
       
    # gravity force
    dxyz += dt*gravity
    
    # boundary condition
    dxyz[fixed_point] *= 0.0  # velocity restriction
    
    dxyz[xyz[:,2]<-150.0]*=np.array([1.0,1.0,0.0]) 
    # velocity[fixed_point] *= np.array([1.0,1.0,0.0])  # fixed on xy plane
    xyz += dt*dxyz
    # print(dxyz)
    # print(position0)
    nodes = xyz.tolist()
    # print(nodes)
    if i % 5 == 0:
        # sv.write_vtk('initial',point=nodes,line=line,face=face)
        sv.write_vtk("ami2/"+"fa"+str(i),point=nodes,face=face)
        sv.write_line_vtk("ami2/"+"main_frame"+str(i),point=nodes,line=geo.main_frame)
        sv.write_line_vtk("ami2/"+"pontoon_cylinder"+str(i),point=nodes,line=geo.pontoon_cylinder)
        sv.write_line_vtk("ami2/"+"pontoon_cone"+str(i),point=nodes,line=geo.pontoon_cone)
        #sv.write_line_vtk("ami2/"+"l5"+str(i),point=nodes,line=geo.l5)
        sv.write_line_vtk("ami2/"+"mooring_line"+str(i),point=nodes,line=geo.mooring_line)
        
'''
# for circular cage
# import src.case.circularCage as case
# case.main()
nodes=geo.nodes
line=geo.l_all
face=geo.netFace
# sv.write_vtk('initial',point=nodes,line=line,face=face)
sv.write_vtk("initial",point=nodes,face=geo.netFace)
sv.write_line_vtk("initial_l1",point=nodes,line=geo.l1)
sv.write_line_vtk("initial_l2",point=nodes,line=geo.l2)
sv.write_line_vtk("initial_l3",point=nodes,line=geo.l3)
sv.write_line_vtk("initial_l5",point=nodes,line=geo.l5)
sv.write_line_vtk("initial_lm",point=nodes,line=geo.lm)

# define structural properties
num_bodyPoint=31
num_seg=50

m_body=4980760 #[kg]
# l1=l.lines(geo.l1,10e9,0.5)
# l2=l.lines(geo.l2,10e9,1.0)
# l3=l.lines(geo.l3,10e9,1.5)
# l5=l.lines(geo.l5,10e9,2.5)
# lm=l.lines(geo.lm,10e9,0.01)




# 
fixed_point=[num_bodyPoint+num_seg*i for i in range(8)]

xyz=np.array(nodes)
dxyz=np.zeros_like(xyz)
gravity=np.array([0,0,-9.81])


run_time = 10  # unit [s]
dt = 2e-2    # unit [s]

# forward Euler (explicit)
for i in range(int(run_time/dt)):       
       
    # gravity force
    dxyz += dt*gravity
    
    # boundary condition
    dxyz[fixed_point] *= 0.0  # velocity restriction
    dxyz[xyz[:,2]<-22.0]*=np.array([1.0,1.0,0.0]) 
    # velocity[fixed_point] *= np.array([1.0,1.0,0.0])  # fixed on xy plane
    xyz += dt*dxyz
    # print(dxyz)
    # print(position0)
    nodes = xyz.tolist()
    # print(nodes)
    if i % 5 == 0:
        # sv.write_vtk('initial',point=nodes,line=line,face=face)
        sv.write_vtk("ami2/"+"fa"+str(i),point=nodes,face=face)
        sv.write_line_vtk("ami2/"+"l1"+str(i),point=nodes,line=geo.l1)
        sv.write_line_vtk("ami2/"+"l2"+str(i),point=nodes,line=geo.l2)
        sv.write_line_vtk("ami2/"+"l3"+str(i),point=nodes,line=geo.l3)
        sv.write_line_vtk("ami2/"+"l5"+str(i),point=nodes,line=geo.l5)
        sv.write_line_vtk("ami2/"+"lm"+str(i),point=nodes,line=geo.lm)       
'''