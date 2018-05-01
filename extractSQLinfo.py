import logging, os,sys,sqlalchemy as sql
from collections import defaultdict

log = logging.getLogger('')
log.setLevel(logging.DEBUG)
format = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(format)
log.addHandler(ch)

def getUniqueCoreAccIds(queryDB,refDB,queryOrgs,refOrgs,cvg=0.5,eval=1e-4,RNA=False):
    ## first get the set of candidate genes from the query list

    engine = sql.create_engine('sqlite:///'+queryDB)
    csr = engine.connect()
    queryCts = dict()
    query = csr.execute('Select orgname,hmmhit,seqid from HMMhits WHERE orgname in ("{}") '
                        'AND hmmcov >= {:.4f} '
                        'AND evalue <= {:.4f}'.format('", "'.join(queryOrgs),cvg,eval))
    for orgname, hmmhit,seqid in query:
        hitDict = queryCts.setdefault(hmmhit,defaultdict(list))
        hitDict[orgname].append(seqid)
    csr.close()

    candidateGenes = {k:v for k,v in queryCts.items() if all(len(v[x]) == 1 for x in queryOrgs)}
    engine = sql.create_engine('sqlite:///'+refDB)
    csr = engine.connect()

    refCts = dict()
    query = csr.execute('Select orgname,hmmhit,seqid from HMMhits WHERE orgname in ("{}") '
                        'AND hmmhit in ("{}") '
                        'AND hmmcov >= {:.4f} '
                        'AND evalue <= {:.4f}'.format('", "'.join(refOrgs),'" , "'.join(candidateGenes.keys()), cvg, eval))
    for orgname,hmmhit,seqid in query:
        hitDict = refCts.setdefault(hmmhit,defaultdict(list))
        hitDict[orgname].append(seqid)
    csr.close()

    mlstGenes = {k:v for k, v in refCts.items() if all(len(v[x]) == 1 for x in refOrgs)}
    log.info('Found {} candidate mlst genes'.format(len(mlstGenes)))

    idDict = {}
    for hmmhit in mlstGenes.keys():
        idDict[hmmhit] = defaultdict(list)
        for queryOrg in queryOrgs:
            idDict[hmmhit]['query'].extend(candidateGenes[hmmhit][queryOrg])
        for refOrg in refOrgs:
            idDict[hmmhit]['ref'].extend(mlstGenes[hmmhit][refOrg])
    return idDict

def cleanName(name):
    name = name.split('(')[0].strip()
    name = name.replace(' ','_')
    return(name)

def writeFastas(idDict,queryDB,refDB,refOrgs,outputFolder='.'):
    ## First Construct the ref org dict
    orgDict = dict()
    engine = sql.create_engine('sqlite:///' + refDB)
    csr = engine.connect()
    query = csr.execute('Select asmID,speciesName,type_strain FROM organisms '
                        'WHERE asmID in ("{}")'.format('" , "'.join(refOrgs)))
    for asmID,name,ts in query:
        orgDict[asmID] = (cleanName(name),ts)
    csr.close()

    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)
    for hmmhit in idDict:
        log.info('Writing {} Seqs'.format(hmmhit))
        seqs = []
        refOrgs = set()
        ## get the query sequences from the
        engine = sql.create_engine('sqlite:///' + queryDB)
        csr = engine.connect()
        query = csr.execute('Select orgname,aaseq FROM Seqs WHERE seqid in ({})'.format(
            ' , '.join(str(x) for x in idDict[hmmhit]['query'])))
        for orgname,aaseq in query:
            seqs.append(('QS--'+orgname,aaseq))
        csr.close()

        engine = sql.create_engine('sqlite:///' + refDB)
        csr = engine.connect()
        query = csr.execute('Select orgname,aaseq FROM Seqs WHERE seqid in ({})'.format(
            ' , '.join(str(x) for x in idDict[hmmhit]['ref'])))
        for orgname,aaseq in query:
            name,ts = orgDict.get(orgname,[orgname,False])
            if ts:
                refName = 'TS--{}:{}'.format(orgname,name)
            else:
                refName = '{}:{}'.format(orgname, name)
            seqs.append((refName,aaseq))
        csr.close()

        with open(os.path.join(outputFolder,'{}.fasta'.format(hmmhit)),'w') as outfile:
            for name,seq in seqs:
                outfile.write('>{}\n{}\n'.format(name,seq))
    return