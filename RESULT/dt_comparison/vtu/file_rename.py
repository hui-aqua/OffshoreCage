import os
import glob

path = os.getcwd()
filenames = os.listdir(path)

filenames= glob.glob('*.vtu')

   
for file in filenames:
    old_name = file
    new_name1 = file.replace("_with_","w")
    # new_name2 = new_name1.replace("_","-")
    os.rename(os.path.join(path, old_name), os.path.join(path, new_name1))