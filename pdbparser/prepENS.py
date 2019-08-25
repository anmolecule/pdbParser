import requests
import logging
import urllib
import sys,os
import numpy as np
from pdbParser import pdbParser as pP
from pdbParser import writepdb as wp
from pdbParser import clean_pdb as cp
from pdbParser import alignment as a

class PDBInfo():
    def __init__(self,query,mer,exclude=None):
        self.exclude=exclude
        self.query=query
        self.mer=mer
        self.result,self.refseq=self.get_pdbinfo()
        self.broken=None
        self.seqfilename=query+'_seq.txt'
        self.residmapfilename=query+'_residmap.txt'
        self.alnfasta=query+'_init.aln.txt'
        self.coremer=None
        self.coreresids=None

    def get_pdbinfo(self):
        URLbase = ('http://www.uniprot.org/uniprot/')

        idparam = {
            'query': 'ID:{}'.format(self.query),
            'format': 'tab',
            'columns': 'database(PDB),sequence'
        }
        idparam=urllib.parse.urlencode(idparam)
        result1 = urllib.request.Request(URLbase+'?'+idparam)
        try:
            result1=urllib.request.urlopen(result1).readlines()[1]
        except:
            logging.error('Cannot retrive the information for query number %s' %(self.query))
            return(None,None)

        if len(result1) > 0:
            pdbids,refseq=result1.split(b'\t')
            refseq=str(refseq)
            pdbids=['{}'.format(i.decode('utf-8')) for i in pdbids.split(b';') if len(i)>1]

            if self.exclude is not None:
                try:
                    for ex in self.exclude:
                        pdbids.remove(ex)
                        logging.info('Removing PDB ID %s' %ex)
                except KeyError:
                    logging.warning('Did not find the PDB ID %s' %ex)

            chainids={}
            for i in urllib.request.urlopen(URLbase+self.query+'.txt').readlines():
                if i.startswith (b'DR   PDB;') and i.split(b';')[3] not in ['NMR;','model;']:
                    chainids[i.split(b';')[1].strip().decode('utf-8')]=i.split(b';')[4].split(b'=')[0].strip().decode('utf-8')
        print(pdbids)
        print(chainids)
        returninfo={}
        tmpids=[*pdbids]
        for pdb in tmpids:
            count=0
            try:
                chains=chainids[pdb]
            except KeyError:
                logging.warning('PDB ID %s is either an NMR structure or a model. Skipping' %pdb)
                pdbids.remove(pdb)
                continue
            try:
                chains=chains.split('/')
            except AttributeError:
                continue
            else:
                nchain=len(chains)
                if nchain == self.mer:
                    returninfo[pdb]=[count+1,[chains]]
                elif nchain > self.mer:
                    if nchain % self.mer == 0:
                        newchains=[]
                        for chnr in range(0,nchain,self.mer):
                            newchains.append(chains[chnr:chnr+self.mer])
                            count=count+1
                        returninfo[pdb]=[count,newchains]
                    else:
                        logging.critical('Cannot process PDB id %s. It does not contain complete set' %pdb)
                        chainids.pop(pdb)
                        pdbids.remove(pdb)
                else:
                    logging.critical('Cannot process PDB id %s. It does not contain complete set' %pdb)
                    chainids.pop(pdb)
                    pdbids.remove(pdb)
        if len(returninfo) == 0:
            return(None,None)
        return(returninfo,refseq)


def downloadPDB(info,cwd):
    query=info.query
    pdblist=info.result
    mer=info.mer
    refseq=info.refseq
    altloc="A"
    outseq=open(cwd+'/'+info.seqfilename,'w')
    outresmap=open(cwd+'/'+info.residmapfilename,'w')
    outseq.write('>refseq'+'\n'+refseq+'\n')
    for pdb in pdblist.keys():
        urllib.request.urlretrieve('http://files.rcsb.org/download/%s.pdb' %pdb, cwd+'/'+pdb+'.pdb')
        #time.sleep(-1)
        for mol in range(0,pdblist[pdb][0]):
            pdblines=open(cwd+'/'+pdb+'.pdb').readlines()
            if pP.pdbTitle(pdblines) is True:
                try:
                    pdblist.pop(pdb)
                    logging.critical('PDB ID %s is a chimera, skipping this file' %pdb)
                    continue
                except KeyError:
                    logging.info('This structure was already removed. I am an example of bad programming. Nothing to worry about')
                    continue
            coord=pP.pdbParser(pdblines,pdb,mer,altloc,[pdblist[pdb][1][mol]])
            coord=cp.getca(coord,altloc,[pdblist[pdb][1][mol]])
            wp.writeca(coord,cwd+'/'+pdb+'_'+str(mol+1)+'.pdb')
            for ch in pdblist[pdb][1][mol]:
                ca=cp.getca(coord,altloc,ch)
                seq,map=a.getseq(ca)
                outseq.write('>'+pdb+'_'+str(mol+1)+'.pdb'+'|'+ch+'|'+'\n'+seq+'\n')
                code,name,nr=zip(*map)
                outresmap.write('>'+pdb+'_'+str(mol+1)+'.pdb'+'|'+ch+'|'+'\n'+'-'.join([str(i) for i in nr])+'\n')
        os.remove(cwd+'/'+pdb+'.pdb')
    outseq.close()
    outresmap.close()

def msa(info,cwd,clustalopath,alnf=None):
    query=info.query
    seqfile=info.seqfilename
    resmap=info.residmapfilename
    merinfo=info.result
    totmer=info.mer
    outfile=cwd+'/'+info.alnfasta
    seqfile=cwd+'/'+seqfile
    complete,resids,broken=a.msa_clustal(seqfile,resmap,outfile,clustalopath,cwd,merinfo,query,totmer,alnf)
    info.broken=broken
    info.coremer=complete
    info.coreresids=resids
    
def getcore(info,cwd):
    altloc='A'
    complete=info.coremer
    resids=info.coreresids
    broken=info.broken
    totmer=info.mer
    for pdb in complete:
        if pdb in broken:
            try:
                os.rename(cwd+'/'+pdb,cwd+'/'+"broken_"+pdb)
                continue
            except (OSError,IOError):
                continue
        chains=complete[pdb]
        try:
            pdblines=open(cwd+'/'+pdb,'r').readlines()
        except (OSError,IOError):
            logging.warning('File does not exist: '+pdb+' skipped')
            continue
        ca=pP.pdbParser(pdblines,pdb,totmer,altloc,[chains])
        newca=None
        for ch in chains:
            nter,cter=[int(i) for i in resids[pdb+'|'+ch+'|']]
            if newca is None:
                newca=ca[(ca['ch']==ch)&(ca['resnr']>=nter) & (ca['resnr']<=cter)]
            else:
                newca=np.concatenate([newca,ca[(ca['ch']==ch)&(ca['resnr']>=nter) & (ca['resnr']<=cter)]])
        wp.writeca(newca,cwd+'/'+'correct_'+pdb)
        os.remove(cwd+'/'+pdb)
