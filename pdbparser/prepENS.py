import requests
import logging
import urllib as ub

logger = logging.getLogger()
logger.setLevel(logging.INFO)

query='Q7NDN8'
totchain=5

class PDBInfo():
    def __init__(self,query,totchain):
        self.query=query
        self.totchain=totchain
        self.result=self.get_pdbinfo()
        
    def get_pdbinfo(self):
        URLbase = ('http://www.uniprot.org/uniprot/')

        idparam = {
            'query': 'ID:{}'.format(self.query),
            'format': 'tab',
            'columns': 'database(PDB),sequence'
        }


        result1 = requests.get(URLbase,params=idparam).text
        if len(result1) > 0:
            pdbids,fasta=result1.split('\n')[1].split('\t')
            pdbids=['{}'.format(i) for i in pdbids.split(';') if len(i)>1]
            chainids={z.split(';')[1].strip():z.split(';')[4].split('=')[0].strip() for z in [i for i in ub.urlopen(URLbase+query+'.txt').readlines() if i.startswith('DR   PDB;')]}
        else:
            print('Cannot retrive the information for query number %s' %(query))
            return('None')

        returninfo={}
        for pdb in pdbids:
            try:
                count=0
                try:
                    chains=chainids[pdb].split('/')
                except AttributeError:
                    continue
                nchain=len(chains)
                if nchain == self.totchain:
                    returninfo[pdb]=[count,chains]
                elif nchain > self.totchain:
                    if nchain % self.totchain == 0:
                        newchains=[]
                        for chnr in xrange(0,nchain,self.totchain):
                            newchains.append(chains[chnr:chnr+self.totchain])
                            print chnr
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

        return(returninfo)
    
    def downloadPDB(self,pdblist):
        for pdb in pdblist:
            ub.urlretrieve('http://files.rcsb.org/download/%s.pdb' %pdb, pdb+'.pdb')
            time.sleep(5)

info=PDBInfo(query,totchain)
for pdb in info.result:
    if info.result[pdb][0]>1:
        print info.result[pdb]

info.downloadPDB(info.result.keys())





