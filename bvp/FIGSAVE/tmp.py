import os
import shutil
for fname in os.listdir('.'):
    i = fname.find('singlegrid');
    if i != 0:
        continue;
    ffname = fname.replace('singlegrid', 'sinpiya2e4', 1)
    cmd = 'mv '+ fname + ' '+ffname
    print(cmd)
    os.system(cmd)
    
    
