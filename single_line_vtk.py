import numpy as np
import src.visualization.saveVtk as sv
import src.element.line as le
from geometry import Orcaflex_dl5_w147 as op

# print(op.dl2)
dL5=[]
for i in range(int(1100/5)):
    dL5.append([i,i+1])
dL2=[]
for i in range(int(1100/2)):
    dL2.append([i,i+1])

sv.write_line_vtk("ami1/orcaflex_dL5",point=op.dL5,line=dL5)    
sv.write_line_vtk("ami1/orcaflex_dL2",point=op.dL2,line=dL2)    


line5=le.lines(dL5,1e9,0.1)
line5.assign_length(np.array(op.dL5))
print(np.sum(line5.line_length))


line2=le.lines(dL2,1e9,0.1)
line2.assign_length(np.array(op.dL2))
print(np.sum(line2.line_length))