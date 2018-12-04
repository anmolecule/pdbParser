import requests
import logging
import urllib
import sys,os
import numpy as np
import pdbParser as pP
import writepdb as wp
import clean_pdb as cp
import alignment as a

class PDBInfo():
    def __init__(self,query,mer):
        self.query=query
        self.mer=mer
        self.result,self.refseq=self.get_pdbinfo()
        self.broken=None
        self.seqfilename=query+'_seq.txt'
        self.residmapfilename=query+'_seq.txt'
        self.alnfasta=query+'_init.aln.txt'

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
                    returninfo[pdb]=[count+1,[chains]]
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


def downloadPDB(info,cwd):
    query=info.query
    pdblist=info.result
    mer=info.mer
    refseq=info.refseq
    altloc="A"
    outseq=open(cwd+'/'+info.seqfilename,'w')
    outresmap=open(cwd+'/'+info.seqfilename,'w')
    outseq.write('>refseq'+'\n'+refseq+'\n')
    for pdb in pdblist.keys():
        urllib.urlretrieve('http://files.rcsb.org/download/%s.pdb' %pdb, cwd+'/'+pdb+'.pdb')
        #time.sleep(-1)
        for mol in range(0,pdblist[pdb][0]):
            pdblines=open(cwd+'/'+pdb+'.pdb').readlines()
            if pP.pdbTitle(pdblines) is True:
                pdblist.pop(pdb)
                logging.critical('PDB ID %s is a chimera, skipping this file' %pdb)
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

def msa(seqfile,resmap,query,cwd,clustalopath,merinfo,totmer):
    outfile=cwd+'/'+info.alnfasta
    #wd+'/'+query+'_msa.aln.fasta'
    seqfile=cwd+'/'+seqfile
    complete,resids=a.msa_clustal(seqfile,resmap,outfile,clustalopath,cwd,merinfo,query,totmer)
    altloc='A'
    print complete
    for pdb in complete:
        chains=complete[pdb]
        try:
            pdblines=open(pdb,'r').readlines()
        except IOError:
            print pdb+' skipped'
            continue
        ca=pP.pdbParser(pdblines,pdb,totmer,altloc,[chains])
        newca=None
        for ch in chains:
            nter,cter=resids[pdb+'|'+ch+'|']
            if newca is None:
                newca=ca[(ca['ch']==ch)&(ca['resnr']>=nter) & (ca['resnr']<=cter)]
            else:
                newca=np.concatenate([newca,ca[(ca['ch']==ch)&(ca['resnr']>=nter) & (ca['resnr']<=cter)]])
        wp.writeca(ca,cwd+'/'+'correct_'+pdb)
        os.remove(cwd+'/'+pdb)





