#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
# import pyvista as pyv
import src.visualization.saveVtk as sv
import csv

# node on main structure
nodes = [[0.000,0.000,4.800],
[0.000,55.000,4.800],
[27.500,47.631,4.800],
[47.631,27.500,4.800],
[55.000,0.000,4.800],
[47.631,-27.500,4.800],
[27.500,-47.631,4.800],
[0.000,-55.000,4.800],
[-27.500,-47.631,4.800],
[-47.631,-27.500,4.800],
[-55.000,0.000,4.800],
[-47.631,27.500,4.800],
[-27.500,47.631,4.800],
[0.000,55.000,2.855],
[27.500,47.631,2.855],
[47.631,27.500,2.855],
[55.000,0.000,2.855],
[47.631,-27.500,2.855],
[27.500,-47.631,2.855],
[0.000,-55.000,2.855],
[-27.500,-47.631,2.855],
[-47.631,-27.500,2.855],
[-55.000,0.000,2.855],
[-47.631,27.500,2.855],
[-27.500,47.631,2.855],
[0.000,55.000,-13.573],
[27.500,47.631,-13.573],
[47.631,27.500,-13.573],
[55.000,0.000,-13.572],
[47.631,-27.500,-13.573],
[27.500,-47.631,-13.573],
[0.000,-55.000,-13.572],
[-27.500,-47.631,-13.573],
[-47.631,-27.500,-13.573],
[-55.000,0.000,-13.573],
[-47.631,27.500,-13.573],
[-27.500,47.631,-13.573],
[0.000,55.000,-30.000],
[27.500,47.631,-30.000],
[47.631,27.500,-30.000],
[55.000,0.000,-30.000],
[47.631,-27.500,-30.000],
[27.500,-47.631,-30.000],
[0.000,-55.000,-30.000],
[-27.500,-47.631,-30.000],
[-47.631,-27.500,-30.000],
[-55.000,0.000,-30.000],
[-47.631,27.500,-30.000],
[-27.500,47.631,-30.000],
[0.000,55.000,-33.000],
[27.500,47.631,-33.000],
[47.631,27.500,-33.000],
[55.000,0.000,-33.000],
[47.631,-27.500,-33.000],
[27.500,-47.631,-33.000],
[0.000,-55.000,-33.000],
[-27.500,-47.631,-33.000],
[-47.631,-27.500,-33.000],
[-55.000,0.000,-33.000],
[-47.631,27.500,-33.000],
[-27.500,47.631,-33.000],
[0.000,0.000,-35.000],
[0.000,0.000,-39.000],
[27.500,47.631,-39.000],
[55.000,0.000,-39.000],
[27.500,-47.631,-39.000],
[-27.500,-47.631,-39.000],
[-55.000,0.000,-39.000],
[-27.500,47.631,-39.000],
[0.000,0.000,-46.000],
[27.500,47.631,-46.000],
[55.000,0.000,-46.000],
[27.500,-47.631,-46.000],
[-27.500,-47.631,-46.000],
[-55.000,0.000,-46.000],
[-27.500,47.631,-46.000],
[0.000,0.000,-38.500]]
num_body_point=len(nodes)
anchor_point = [[429.356,1017.798,-150.000], # A Node 38
                [997.667,449.488,-150.000], # B Node 38
                [997.667,-449.488,-150.000], # C Node 42
                [429.356,-1017.798,-150.000], # D Node 42
                [-429.356,-1017.798,-150.000], # D Node 44
                [-997.667,-449.488,-150.000], # E Node 44
                [-997.667,449.488,-150.000], # F Node 48
                [-429.356,1017.798,-150.000]] # G Node 48
        

attached_point = [38,38,42,42,44,44,48,48]

num_seg=110
mooring_point=[]
for i in range(8):
    ml=np.linspace(anchor_point[i],nodes[attached_point[i]],num_seg,endpoint=False)
    mooring_point+=ml.tolist()

nodes += mooring_point

top_cross_beam = [[15,21],
                  [18,24],
                  [13,19]]

top_hor_beam = [[13,14],
               [14,15],
               [15,16],
               [16,17],
               [17,18],
               [18,19],
               [19,20],
               [20,21],
               [21,22],
               [22,23],
               [23,24],
               [24,13]]

mid_hor_beam = [[25,26],
               [26,27],
               [27,28],
               [28,29],
               [29,30],
               [30,31],
               [31,32],
               [32,33],
               [33,34],
               [34,35],
               [35,36],
               [36,25]]

bottom_hor_beam = [[37,38],
                   [38,39],
                   [39,40],
                   [40,41],
                   [41,42],
                   [42,43],
                   [43,44],
                   [44,45],
                   [45,46],
                   [46,47],
                   [47,48],
                   [48,37]]

dia_beam = [[13,38],
            [38,15],
            [15,40],
            [40,17],
            [17,42],
            [42,19],
            [19,44],
            [44,21],
            [21,46],
            [46,23],
            [23,48],
            [48,13]]

side_column_wo_pontoon = [[1,49],
                          [3,51],
                          [5,53],
                          [7,55],
                          [9,57],
                          [11,59]]

bottom_rad_beam = [[76,37],
                   [76,38],
                   [76,39],
                   [76,40],
                   [76,41],
                   [76,42],
                   [76,43],
                   [76,44],
                   [76,45],
                   [76,46],
                   [76,47],
                   [76,48]]

side_column_w_pontoon = [[2,50],
                         [4,52],
                         [6,54],
                         [8,56],
                         [10,58],
                         [12,60]]

pontoon_cylinder = [[63,70],
                    [64,71],
                    [65,72],
                    [66,73],
                    [67,74],
                    [68,75]]

pontoon_cone = [[50,63],
                [52,64],
                [54,65],
                [56,66],
                [58,67],
                [60,68]]

center_column = [[0,61]]
center_pontoon_cylinder = [[62,69]]
center_pontoon_cone = [[61,62]]

mooring_line=[]
for i in range(8):
    for j in range(num_seg-1):
        mooring_line.append([i*num_seg+num_body_point+j,i*num_seg+num_body_point+j+1])
    mooring_line.append([i*num_seg+num_body_point+j+1,attached_point[i]])

sideColumnWPontoon   = side_column_w_pontoon + pontoon_cylinder + pontoon_cone
centerColumnWPontoon = center_column + center_pontoon_cylinder + center_pontoon_cone
main_frame           = sideColumnWPontoon + centerColumnWPontoon + top_cross_beam + top_hor_beam + mid_hor_beam + bottom_hor_beam + dia_beam + side_column_wo_pontoon + bottom_rad_beam
all_line             = top_cross_beam + top_hor_beam + mid_hor_beam + bottom_hor_beam + dia_beam + side_column_w_pontoon + side_column_wo_pontoon + pontoon_cylinder + pontoon_cone + mooring_line + bottom_rad_beam

netFace = [[13,25,14,26],
        [14,26,15,27],
        [15,27,16,28],
        [16,28,17,29],
        [17,29,18,30],
        [18,30,19,31],
        [19,31,20,32],
        [20,32,21,33],
        [21,33,22,34],
        [22,34,23,35],
        [23,35,24,36],
        [24,36,13,25],
        [25,37,26,38],
        [26,38,27,39],
        [27,39,28,40],
        [28,40,29,41],
        [29,41,30,42],
        [30,42,31,43],
        [31,43,32,44],
        [32,44,33,45],
        [33,45,34,46],
        [34,46,35,47],
        [35,47,36,48],
        [36,48,25,37],
        [37,38,76],
        [38,39,76],
        [39,40,76],
        [40,41,76],
        [41,42,76],
        [42,43,76],
        [43,44,76],
        [44,45,76],
        [45,46,76],
        [46,47,76],
        [47,48,76],
        [48,37,76]]

if __name__ == '__main__':
    pass

'''
# private function
def __gen_points():
    with open('point.csv', 'r') as pf:
        next(pf)
        point_data = csv.reader(pf)
        point_one_cage = []
        for line in point_data:
            point_one_cage.append([float(x) for x in line])
    return point_one_cage


def __gen_lines():    
    with open('lines.csv', 'r') as lf:
        next(lf)
        lines_data = csv.reader(lf)
        lines = []
        for line in lines_data:
            lines.append([int(x) for x in line])
    return lines


def __gen_surfs():
    with open('surfs.csv', 'r') as sf:
        next(sf)
        surfs_data = csv.reader(sf)
        surf_element = []
        for line in surfs_data:
            surf_element.append([int(x) for x in line])
    return surf_element

# public function
def gen_lines():
    lines=__gen_lines()
    return lines

def gen_cage():
    points = __gen_points()
    surfs = __gen_surfs()
    return points, surfs
'''