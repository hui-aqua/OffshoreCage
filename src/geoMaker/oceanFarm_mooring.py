#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np


anchor_point = [[429.356, 1017.798, -150.000], # A Node 38
                [997.667, 449.488, -150.000], # B Node 38
                [997.667, -449.488, -150.000], # C Node 42
                [429.356, -1017.798, -150.000], # D Node 42
                [-429.356, -1017.798, -150.000], # E Node 44
                [-997.667, -449.488, -150.000], # F Node 44
                [-997.667, 449.488, -150.000], # G Node 48
                [-429.356, 1017.798, -150.000]] # H Node 48

attached_point = [[27.5,  47.631, -27.0], # a and b
                  [27.5, -47.631, -27.0], # c and d
                  [-27.5,-47.631, -27.0], # e and f
                  [-27.5, 47.631, -27.0]] # g and h


mooring_length = 1100

def generate_model(num_seg=110):
    # automatic generate the mooring lines according to number of segment per each line
    # num_seg should be an integer and at least be 11. 
    # can also be 55,110,220 for this study.
    # the position of nodes, element and fixed position are return in this function
    
    mooring_line_point=[]
    for i in range(8):
        ml=np.linspace(anchor_point[i],attached_point[i//2],num_seg,endpoint=False)
        mooring_line_point+=ml.tolist()
        if i%2==0:
            mooring_line_point+=[attached_point[i//2]]

    num_p=num_seg*2+1
    fixed_point= [0, 
                  (1*num_seg)+1, 
                  (2*num_seg)+1, 
                  (3*num_seg)+2, 
                  (4*num_seg)+2, 
                  (5*num_seg)+3, 
                  (6*num_seg)+3, 
                  (7*num_seg)+4]  # anchor point
    body_attached_point=[num_seg+num_p*i for i in range(4)]
    fixed_point+=body_attached_point
    mooring_line_fiber =[]
    for j in range(4):
        for i in range(int(num_seg/11)):
            mooring_line_fiber.append([num_seg-i+num_p*j,num_seg-1-i+num_p*j])
        mooring_line_fiber.append([num_seg+num_p*j,2*num_seg+num_p*j])
        for i in range(int(num_seg/11)-1):
            mooring_line_fiber.append([2*num_seg-i+num_p*j,2*num_seg-1-i+num_p*j])

    mooring_line_chain=[]
    for j in range(4):
        for i in range(int(num_seg/11*10)):
            mooring_line_chain.append([i+num_p*j,i+1+num_p*j])
        for i in range(int(num_seg/11*10)):
            mooring_line_chain.append([num_seg+1+i+num_p*j,num_seg+1+i+1+num_p*j])

    mooring_line_mass =[147*mooring_length/num_seg/2]
    mooring_line_mass+=[147*mooring_length/num_seg]*int(num_seg/11*10-1)
    mooring_line_mass+=[147*mooring_length/num_seg/2 + 4*mooring_length/num_seg/2]
    mooring_line_mass+=[4*mooring_length/num_seg]*int(num_seg/11)

    mooring_line_mass+=[147*mooring_length/num_seg/2]
    mooring_line_mass+=[147*mooring_length/num_seg]*int(num_seg/11*10-1)
    mooring_line_mass+=[147*mooring_length/num_seg/2 + 4*mooring_length/num_seg/2]
    mooring_line_mass+=[4*mooring_length/num_seg]*int(num_seg/11-1)

    mooring_line_mass+=mooring_line_mass
    mooring_line_mass+=mooring_line_mass

    return mooring_line_point,mooring_line_fiber,mooring_line_chain,mooring_line_mass, fixed_point

if __name__ == '__main__':
    pass