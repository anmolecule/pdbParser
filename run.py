#!/usr/bin/env python

import argparse,logging
from os import getcwd
from pdbParser import pdbParser as pP
from pdbParser import readpdb as rp
from pdbParser import writepdb as wp
#from pdbparser.pdbParser import pdbParser
#from pdbparser.readpdb import getpdb
#from pdbparser.readpdb import checkmulti
#from pdbparser.writepdb import writeca
from pdbParser import alignment as a
from pdbParser import prepENS as pE
#from align.alignment import getaligned, multialigned
#import prepENS.prepENS

parser = argparse.ArgumentParser(description='Identification of Missing residues')
parser.add_argument('--start', metavar='PDB File', nargs=1 , help='Starting structure')
parser.add_argument('--target', metavar='PDB File', nargs=1 , help='Target structure')
parser.set_defaults(local=False)
parser.add_argument('--schains',dest='schains',help='Supply the chain ids of the starting structure: A-E or A,B,C')
parser.set_defaults(schains='A')
parser.add_argument('--echains',dest='echains',help='Supply the chain ids of the target structure: A-E or A,B,C')
parser.set_defaults(echains='A')
parser.add_argument('--multimeric', dest='mer', help='If the protein is multimeric provide the number of chains. Default is 1', type=int)
parser.set_defaults(mer=1)
parser.add_argument('--dir', dest='cwd', help='The directory to save the output files. Default is current work directory.')
parser.set_defaults(cwd=getcwd())
parser.add_argument('--altloc', dest='altloc', help='Alternative location to be extracted, default is A')
parser.set_defaults(altloc='A')
parser.add_argument('--prepENS',dest='prepENS',metavar='Uniport ID',nargs=1 ,help='This argument accepts a Uniprot ID, and returns the PDB Files that are complete. PDB IDs with missing residues are written with the flat broken_ . When this argument is given, all the other arguments except for multimeric is not used.')
parser.add_argument('--exclude',dest='exclude',metavar='PDB IDs to skip',help='Provide a list of PDBs to skip: in the form of ID1A,ID1B,...')
parser.add_argument('--clustal',dest='clustal',metavar='Path to clustal program')
parser.set_defaults(clustal=None)
parser.set_defaults(prepENS=None)
parser.set_defaults(exclude=None)

args=parser.parse_args()
try:
    sid=args.start[0]
    eid=args.target[0]
except TypeError:
    if args.prepENS is None:
        parser.print_help()
        exit()
    else:
        pass
logging.info('Output directory is %s.\n\t If the input files are not given with full path, current working directory is used to search.' %args.cwd)

#Below does the prepENS and then exits. So the rest of the commands are ignored.
if args.prepENS is not None and args.mer is not None:
    if args.exclude is not None:
        exclude=args.exclude.split()
    else:
        exclude=None
    info=pE.PDBInfo(args.prepENS[0],args.mer,exclude)
    if info.result:
        pE.downloadPDB(info,args.cwd)
        clustal=args.clustal
        pE.msa(info,args.cwd,clustal)
        pE.getcore(info,args.cwd)
    else:
        logging.error('No structure is available')
    exit()
else:
    parser.print_help()
    exit()

def chlistsplit(chlist):
    finlist=[]
    tmplist=chlist.split(",")
    for sub in tmplist:
        if len(sub)==3 and sub[1]=="-":
            sublist=sub.split("-")
            firstch=sublist[0]
            lastch=sublist[1]
            finsublist=[chr(i) for i in range(ord(firstch),ord(lastch)+1)]
            for ch in finsublist:finlist.append(ch)
        elif len(sub)==1 and isinstance(sub, str):
            finlist.append(sub)
        else:
            logging.critical("Please supply the chain ids in these formats: A-E or A-C,E or A,B")
            exit()
    if args.mer == len(finlist):
        return[finlist]
    if args.mer > len(finlist):
        if len(finlist) % args.mer == 0:
            moldivision=[]
            for i in xrange(0,len(finliist),args.mer):
                moldivision.append(finlist[i:i+args.mer])
            return moldivision
        else:
            logging.critical("Chain labels are not equal or multiple of total number of chains given.")
            exit()

schains=chlistsplit(args.schains)
echains=chlistsplit(args.echains)

#toAlign=True

outf=args.cwd+'/error.dat'
logging.basicConfig(level=logging.INFO,filename=outf)
if args.mer == 1:
    for ch1,ch2 in zip(schains,echains):
        start=rp.getpdb(sid)
        rp.checkmulti(start)
        logging.info('Processing PDB files')
        sca=pP.pdbParser(start,sid,args.mer,args.altloc,ch1)
        end=rp.getpdb(eid)
        rp.checkmulti(end)
        logging.info('Processing PDB files')
        eca=pP.pdbParser(end,eid,args.mer,args.altloc,ch2)
        logging.info('Extracting the core region.')
        score,ecore,correct=a.getaligned(sca,eca)
        if correct is True:
            wp.writeca(score,args.cwd+'/start.pdb')
            wp.writeca(ecore,args.cwd+'/target.pdb')
            break
        else:
            continue

elif args.mer !=1:
    start=rp.getpdb(sid)
    rp.checkmulti(start)
    logging.info('Processing PDB files')
    sca=pP.pdbParser(start,sid,args.mer,args.altloc,schains)
    end=rp.getpdb(eid)
    rp.checkmulti(end)
    logging.info('Processing PDB files')
    eca=pP.pdbParser(end,eid,args.mer,args.altloc,echains)
    logging.info('Extracting the core region.')
    score,ecore,correct=a.multialigned(sca,eca,args.mer)
    if correct is True:
        wp.writeca(score,args.cwd+'/start.pdb')
        wp.writeca(ecore,args.cwd+'/target.pdb')

