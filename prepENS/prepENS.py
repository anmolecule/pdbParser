
import requests
import logging
import urllib
import sys,os
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath('prepENS.py'))))
from pdbparser.pdbParser import pdbParser, pdbTitle
from pdbparser.writepdb import writeca
from pdbparser.clean_pdb import getca
from align.alignment import msa_clustal, getseq

class PDBInfo():
    def __init__(self,query,mer):
        self.query=query
        self.mer=mer
        self.result,self.refseq=self.get_pdbinfo()

    def get_pdbinfo(self):
        URLbase = ('http://www.uniprot.org/uniprot/')

        idparam = {
            'query': 'ID:{}'.format(self.query),
            'format': 'tab',
            'columns': 'database(PDB),sequence'
        }


        result1 = requests.get(URLbase,params=idparam).text
        if len(result1) > 0:
            pdbids,refseq=result1.split('\n')[1].split('\t')
            pdbids=['{}'.format(i) for i in pdbids.split(';') if len(i)>1]
            chainids={z.split(';')[1].strip():z.split(';')[4].split('=')[0].strip() for z in [i for i in urllib.urlopen(URLbase+self.query+'.txt').readlines() if i.startswith('DR   PDB;')] if z.split()[3]=='X-ray;'}
        else:
            logging.critical('Cannot retrive the information for query number %s' %(self.query))
            return(None,None)

        returninfo={}
        for pdb in pdbids:
            try:
                count=0
                try:
                    chains=chainids[pdb].split('/')
                except AttributeError:
                    continue
                nchain=len(chains)
                if nchain == self.mer:
                    returninfo[pdb]=[count+1,chains]
                elif nchain > self.mer:
                    if nchain % self.mer == 0:
                        newchains=[]
                        for chnr in xrange(0,nchain,self.mer):
                            newchains.append(chains[chnr:chnr+self.mer])
                            count=count+1
                        returninfo[pdb]=[count,newchains]
                    else:
                        logging.warning('Cannot process PDB id %s. It does not contain complete set' %pdb)
                        chainids.pop(pdb)
                        pdbids.remove(pdb)
                else:
                    logging.warning('Cannot process PDB id %s. It does not contain complete set' %pdb)
                    chainids.pop(pdb)
                    pdbids.remove(pdb)
            except KeyError:
                logging.warning('PDB ID %s has no chain information' %pdb)
                pdbids.remove(pdb)
        return(returninfo,refseq)


def downloadPDB(pdblist,mer,altloc,refseq,query,cwd):
    outseq=open(cwd+'/'+query+'_seq.txt','w')
    outresmap=open(cwd+'/'+query+'_residmap.txt','w')
    outseq.write('>refseq'+'\n'+refseq+'\n')
    for pdb in pdblist.keys():
        urllib.urlretrieve('http://files.rcsb.org/download/%s.pdb' %pdb, cwd+'/'+pdb+'.pdb')
        #time.sleep(-1)
        for mol in range(0,pdblist[pdb][0]):
            pdblines=open(cwd+'/'+pdb+'.pdb').readlines()
            if pdbTitle(pdblines) is True:
                pdblist.pop(pdb)
                logging.critical('PDB ID %s is a chimera, skipping this file' %pdb)
                continue
            coord=pdbParser(pdblines,pdb,mer,altloc,[pdblist[pdb][1][mol]])
            coord=getca(coord,altloc,[pdblist[pdb][1][mol]])
            writeca(coord,cwd+'/'+pdb+'_'+str(mol+1)+'.pdb')
            for ch in pdblist[pdb][1][mol]:
                ca=getca(coord,altloc,ch)
                seq,map=getseq(ca)
                outseq.write('>'+pdb+'_'+str(mol+1)+'.pdb'+'|'+ch+'|'+'\n'+seq+'\n')
                code,name,nr=zip(*map)
                outresmap.write('>'+pdb+'_'+str(mol+1)+'.pdb'+'|'+ch+'|'+'\n'+'-'.join([str(i) for i in nr])+'\n')
        os.remove(cwd+'/'+pdb+'.pdb')
    outseq.close()
    outresmap.close()

def msa(seqfile,resmap,query,cwd,clustalopath,merinfo,totmer):
    outfile=cwd+'/'+query+'_msa.aln.fasta'
    seqfile=cwd+'/'+seqfile
    complete,resids=msa_clustal(seqfile,resmap,outfile,clustalopath,cwd,merinfo,query,totmer)
    altloc='A'
    print complete
    for pdb in complete:
        chains=complete[pdb]
        pdblines=open(pdb,'r').readlines()
        ca=pdbParser(pdblines,pdb,totmer,altloc,[chains])
        newca=None
        for ch in chains:
            nter,cter=resids[pdb+'|'+ch+'|']
            if newca is None:
                newca=ca[(ca['ch']==ch)&(ca['resnr']>=nter) & (ca['resnr']<=cter)]
            else:
                newca=np.concatenate([newca,ca[(ca['ch']==ch)&(ca['resnr']>=nter) & (ca['resnr']<=cter)]])
        writeca(ca,cwd+'/'+'correct_'+pdb)
        os.remove(cwd+'/'+pdb)





