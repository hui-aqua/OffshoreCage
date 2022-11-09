import numpy as np
import src.geoMaker.oceanFarm1 as geo
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
k_fiber = E_fiber*A_fiber/seg_length

E_chain = 60e9
A_chain = 0.012170404
k_chain = E_chain*A_chain/seg_length

num_mooring_line_seg = geo.num_seg 
l1=l.lines(geo.mooring_line_new,k_chain,0.088) # Axial stiffness[MN] 680.81 (Chain) 235.44 (Fiber)
l1f=l.lines(geo.mooring_line_fiber,k_fiber,0.160)
l1.assign_length(10.0)

## setting
fixed_point =  geo.fixed_point          # anchor point
fixed_point += geo.body_attached_point  # fish cage body

gravity=np.array([0,0,-9.81])
current=np.array([[1.0,0,0]]*len(nodes))
mass_matrix = np.array(geo.mass_mooring_line_new).reshape(len(geo.mass_mooring_line_new),1)

total_time = [5,25,50,75,100,150,200,250,300]       # Total simulation runn time, unit [s]
dt = 0.0001                         # Time steps, unit [s]
force_on_cage=[]
force_on_BP1 = []
force_on_cage = []
force_on_BP1  = []
## initialization 
position=np.array(nodes)
velocity=np.zeros_like(position)

# print(l1f.k)
for run_time in total_time:
# for dt in time_step:

    for i in range(int(run_time/dt)):       
        if i % (round(1/dt/30,-2)) == 0:        # write result for 30 frames per second  

        sv.write_line_vtk("ami2/"+str(run_time)+"_with_"+str(dt) + "mooring_line"+str(i),point=position.tolist(),line=line)
               
    ### forward Euler (explicit)
    ## external loads
    # gravity force
    pre_position=position.copy()
    velocity += gravity*dt
    
    ## boundary condition
    velocity[fixed_point] *= 0.0  # velocity restriction
    velocity[position[:,2]<-150]*=np.array([1,1,0])# ground
    position += velocity*dt
    position[fixed_point]=np.array(nodes)[fixed_point]
    
    ### constraint function 
    position+=l1.pbd_edge_constraint(position,mass_matrix,dt)
    
    ### velocity correction
    velocity=(position-pre_position)/dt
    
    BP1 = 110
    ML1a = BP1 - 1
    ML1b = BP1 + 110

    BP2 = 331
    ML2a = BP1 - 1
    ML2b = BP1 + 110

    BP3 = 552
    ML3a = BP1 - 1
    ML3b = BP1 + 110

    BP4 = 773
    ML4a = BP1 - 1
    ML4b = BP1 + 110

    force1 = 235.44e6*(np.linalg.norm(position[BP1]-position[ML1a]) -10)
    forceVector1 = np.array(position[BP1])-np.array(position[ML1a])
    VectorMagnitude1 = np.linalg.norm(forceVector1)
    FV1x = force1*forceVector1[0]/VectorMagnitude1
    FV1y = force1*forceVector1[1]/VectorMagnitude1
    FV1z = force1*forceVector1[2]/VectorMagnitude1
    FV1 = [FV1x, FV1y, FV1z]
    
    force2 = 235.44e6*(np.linalg.norm(position[BP1]-position[ML1b]) -10)
    forceVector2 = np.array(position[BP1])-np.array(position[ML1b])
    VectorMagnitude2 = np.linalg.norm(forceVector2)
    FV2x = force2*forceVector2[0]/VectorMagnitude2
    FV2y = force2*forceVector2[1]/VectorMagnitude2
    FV2z = force2*forceVector2[2]/VectorMagnitude2
    FV2 = [FV2x, FV2y, FV2z]

    force3 = 235.44e6*(np.linalg.norm(position[BP2]-position[ML2a]) -10)
    forceVector3 = np.array(position[BP2])-np.array(position[ML2a])
    VectorMagnitude3 = np.linalg.norm(forceVector3)
    FV3x = force3*forceVector3[0]/VectorMagnitude3
    FV3y = force3*forceVector3[1]/VectorMagnitude3
    FV3z = force3*forceVector3[2]/VectorMagnitude3
    FV3 = [FV3x, FV3y, FV3z]

    force4 = 235.44e6*(np.linalg.norm(position[BP2]-position[ML2b]) -10)
    forceVector4 = np.array(position[BP2])-np.array(position[ML2b])
    VectorMagnitude4 = np.linalg.norm(forceVector4)
    FV4x = force4*forceVector3[0]/VectorMagnitude4
    FV4y = force4*forceVector3[1]/VectorMagnitude4
    FV4z = force4*forceVector3[2]/VectorMagnitude4
    FV4 = [FV4x, FV4y, FV4z]

    force5 = 235.44e6*(np.linalg.norm(position[BP3]-position[ML3a]) -10)
    forceVector5 = np.array(position[BP3])-np.array(position[ML3a])
    VectorMagnitude5 = np.linalg.norm(forceVector5)
    FV5x = force5*forceVector5[0]/VectorMagnitude5
    FV5y = force5*forceVector5[1]/VectorMagnitude5
    FV5z = force5*forceVector5[2]/VectorMagnitude5
    FV5 = [FV5x, FV5y, FV5z]
    
    force6 = 235.44e6*(np.linalg.norm(position[BP3]-position[ML3b]) -10)
    forceVector6 = np.array(position[BP3])-np.array(position[ML3b])
    VectorMagnitude6 = np.linalg.norm(forceVector6)
    FV6x = force6*forceVector6[0]/VectorMagnitude6
    FV6y = force6*forceVector6[1]/VectorMagnitude6
    FV6z = force6*forceVector6[2]/VectorMagnitude6
    FV6 = [FV6x, FV6y, FV6z]

    force7 = 235.44e6*(np.linalg.norm(position[BP4]-position[ML4a]) -10)
    forceVector7 = np.array(position[BP4])-np.array(position[ML4a])
    VectorMagnitude7 = np.linalg.norm(forceVector7)
    FV7x = force7*forceVector7[0]/VectorMagnitude7
    FV7y = force7*forceVector7[1]/VectorMagnitude7
    FV7z = force3*forceVector7[2]/VectorMagnitude7
    FV7 = [FV7x, FV7y, FV7z]

    force8 = 235.44e6*(np.linalg.norm(position[BP4]-position[ML4b]) -10)
    forceVector8 = np.array(position[BP4])-np.array(position[ML4b])
    VectorMagnitude8 = np.linalg.norm(forceVector8)
    FV8x = force8*forceVector8[0]/VectorMagnitude8
    FV8y = force8*forceVector8[1]/VectorMagnitude8
    FV8z = force4*forceVector8[2]/VectorMagnitude8
    FV8 = [FV8x, FV8y, FV8z]

    forceBP1 = [FV1x+FV2x, FV1y+FV2y, FV1z+FV2z]
    forceBP2 = [FV3x+FV4x, FV3y+FV4y, FV3z+FV4z]
    forceBP3 = [FV5x+FV6x, FV5y+FV6y, FV5z+FV6z]
    forceBP4 = [FV7x+FV8x, FV8y+FV8y, FV7z+FV8z]

    force_on_cage.append([forceBP1,forceBP2,forceBP3,forceBP4])
    force_on_BP1.append([FV1x+FV2x, FV1y+FV2y, FV1z+FV2z])


np.savetxt( str(run_time) + "s_with_dt-" + str(dt) + 's_ForceOnBP1.csv', force_on_BP1)

        force3 = k_fiber*(np.linalg.norm(position[BP2]-position[ML2a]) -10)
        forceVector3 = np.array(position[BP2])-np.array(position[ML2a])
        VectorMagnitude3 = np.linalg.norm(forceVector3)
        FV3x = force3*forceVector3[0]/VectorMagnitude3
        FV3y = force3*forceVector3[1]/VectorMagnitude3
        FV3z = force3*forceVector3[2]/VectorMagnitude3
        FV3 = [FV3x, FV3y, FV3z]

        force4 = k_fiber*(np.linalg.norm(position[BP2]-position[ML2b]) -10)
        forceVector4 = np.array(position[BP2])-np.array(position[ML2b])
        VectorMagnitude4 = np.linalg.norm(forceVector4)
        FV4x = force4*forceVector3[0]/VectorMagnitude4
        FV4y = force4*forceVector3[1]/VectorMagnitude4
        FV4z = force4*forceVector3[2]/VectorMagnitude4
        FV4 = [FV4x, FV4y, FV4z]

        force5 = k_fiber*(np.linalg.norm(position[BP3]-position[ML3a]) -10)
        forceVector5 = np.array(position[BP3])-np.array(position[ML3a])
        VectorMagnitude5 = np.linalg.norm(forceVector5)
        FV5x = force5*forceVector5[0]/VectorMagnitude5
        FV5y = force5*forceVector5[1]/VectorMagnitude5
        FV5z = force5*forceVector5[2]/VectorMagnitude5
        FV5 = [FV5x, FV5y, FV5z]
        
        force6 = k_fiber*(np.linalg.norm(position[BP3]-position[ML3b]) -10)
        forceVector6 = np.array(position[BP3])-np.array(position[ML3b])
        VectorMagnitude6 = np.linalg.norm(forceVector6)
        FV6x = force6*forceVector6[0]/VectorMagnitude6
        FV6y = force6*forceVector6[1]/VectorMagnitude6
        FV6z = force6*forceVector6[2]/VectorMagnitude6
        FV6 = [FV6x, FV6y, FV6z]

        force7 = k_fiber*(np.linalg.norm(position[BP4]-position[ML4a]) -10)
        forceVector7 = np.array(position[BP4])-np.array(position[ML4a])
        VectorMagnitude7 = np.linalg.norm(forceVector7)
        FV7x = force7*forceVector7[0]/VectorMagnitude7
        FV7y = force7*forceVector7[1]/VectorMagnitude7
        FV7z = force3*forceVector7[2]/VectorMagnitude7
        FV7 = [FV7x, FV7y, FV7z]

        force8 = k_fiber*(np.linalg.norm(position[BP4]-position[ML4b]) -10)
        forceVector8 = np.array(position[BP4])-np.array(position[ML4b])
        VectorMagnitude8 = np.linalg.norm(forceVector8)
        FV8x = force8*forceVector8[0]/VectorMagnitude8
        FV8y = force8*forceVector8[1]/VectorMagnitude8
        FV8z = force4*forceVector8[2]/VectorMagnitude8
        FV8 = [FV8x, FV8y, FV8z]

        forceBP1 = [FV1x+FV2x, FV1y+FV2y, FV1z+FV2z]
        forceBP2 = [FV3x+FV4x, FV3y+FV4y, FV3z+FV4z]
        forceBP3 = [FV5x+FV6x, FV5y+FV6y, FV5z+FV6z]
        forceBP4 = [FV7x+FV8x, FV8y+FV8y, FV7z+FV8z]

        force_on_cage.append([forceBP1,forceBP2,forceBP3,forceBP4])
        
        if i % (run_time / dt / 10000) == 0: # write force histroy for force every 1000 iteration
            force_on_BP1.append([i/10000, m.sqrt((FV1x+FV2x)**2 + (FV1y+FV2y)**2 + (FV1z+FV2z)**2), FV1x+FV2x, FV1y+FV2y, FV1z+FV2z])
        
    np.savetxt( str(run_time) + "s_with_dt-" + str(dt) + 's_ForceOnBP1.csv', force_on_BP1, delimiter=",") # Save file

  
# # x axis values
# x = force_on_BP1[:,0]
# # corresponding y axis values
# y = force_on_BP1[:,1]

# print(x)

# # plotting the points 
# plt.plot(x, y)
  
# # naming the x axis
# plt.xlabel('Run time (s)')
# # naming the y axis
# plt.ylabel('Force magnitude (N)')
# plt.yscale("log")

# # giving a title to my graph
# plt.title('Force time history')
  
# # function to show the plot
# plt.show()