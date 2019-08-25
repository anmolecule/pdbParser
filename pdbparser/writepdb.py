#See COPYING for license 

import numpy as np

def writeca(div,file):
    with open(file) as fout:
        for d in div:
            print(d)
            atnr,atnam,resname,resnr,x,y,z,segname=d
            fout.write("{:6s}{:5d} {:^4s}{:4s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:4s}\n".format('ATOM',atnr,atnam,resname,resnr,x,y,z,0.0,0.0,segname))
