#See COPYING for license 
import numpy as np
import logging

def missinginfo(r465,chlist,ca):
    flatchlist=[]
    for mol in chlist:
        for ch in mol:
            flatchlist.append(ch)

    missinginfo={}
    if len(r465) != 0:
        logging.info('Identified missing residues.')
        logging.info('Checking the integrity of the structure')
        for ch in flatchlist:
            rmis=r465[r465['ch']==ch]['rid'].tolist()
            max=ca[ca['ch']== ch ]['resnr'].tolist()[-1]
            min=ca[ca['ch']== ch ]['resnr'].tolist()[0]
            for i in rmis:
                if i < max and i > min:
                    logging.warning('If there are multiple biological assemblies, I will try to pick the one that is not broken. Else, I will quit')
                    #we set the chain as not complete if the missing residues are at the termini only.
                    #0 broken, 1 complete
                    missinginfo[ch]=0
                    continue
    else:
        logging.info('No missing residues detected.')
        for ch in flatchlist:
            missinginfo[ch]=1
    #This is to set the remaining chains as broken. 
    # Perhaps modifying the loop would help, but when the below commands are within the first *if* loop chain info is overwritten. 
    for ch in flatchlist:
        if ch not in missinginfo.keys():
            missinginfo[ch]=1
    if all(ch == 0 for ch in missinginfo.values()) is True:
        return None
    else:
        return missinginfo
