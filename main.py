import numpy as np
import src.visualization.vtkVisualization as pv
import src.geoMaker.oceanFarm1 as geo
import src.visualization.saveVtk as sv
import src.element.line as l
import src.visualization.showMatrix as sm
import src.element.quad as q
import src.waterWave.regularWaves as ww

point=geo.__gen_points()

print(point)
sv.write_vtkPoint(point,"debug")