#See COPYING for license 

import numpy as np
import logging

def getch(coord,chlist):
    flatchlist=[]
    for mol in chlist:
        for ch in mol:
            flatchlist.append(ch)
    delch=coord[np.in1d(coord['ch'],flatchlist)]
    return delch

def getca(coord,altloc='A',chlist):
    subcoord=getch(coord,chlist)
    delalter=subcoord[(subcoord['altloc'] == altloc) | (subcoord['altloc'] == '')]
    logging.warning('Cleaning alternative locations if present')
    logging.warning('Default alternative location is A')
    ca=delalter[np.in1d(delalter['atname'],'CA')]
    return ca
