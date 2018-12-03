#See COPYING for license 

from urllib2 import urlopen
import readpdb, clean_pdb, divide_mer
import missing, writepdb
import logging

def pdbTitle(pdb):
    for line in pdb:
        if 'TITLE' in line:
            title=line.split()
            break
    for i in title:
        if 'CHIMERA' in title or 'FUSED' in title:
            return True 
        else:
            return False

def pdbParser(pdb,pdbid,mer,altloc,chlist):
    logging.info('Retriving CA coordinates')
    logging.info('Checking for missing residues')
    atomlines=readpdb.readatom(pdb)
    coords=readpdb.coord(atomlines)
    ca=clean_pdb.getca(coords,altloc,chlist)
    return ca

    #Above loop is the part where I check the missing residues. I have a feeling we should do this before. But if we want a sequence alignment between two structure and retrive a region automatically for eBDIMS then doing it after is better so that we can also use the remarks and seqres to chop and align the sequences.
