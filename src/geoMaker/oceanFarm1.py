#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pyvista as pyv
import src.visualization.saveVtk as sv
import csv

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