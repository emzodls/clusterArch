import logging, os,sys,sqlalchemy as sql
from itertools import chain
from collections import defaultdict

log = logging.getLogger('')
log.setLevel(logging.DEBUG)
format = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(format)
log.addHandler(ch)

def getUniqueCoreAccIds(queryDB,refDB,queryOrgs,refOrgs,cvg=0.5,eval=1e-4):
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

    log.info('Found {} candidate mlst genes in query organism'.format(len(candidateGenes)))
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

def writeFastas(idDict,queryDB,refDB,reforgs,ogorgs=[],addAA=False,RNA=False,outputFolder='.'):
    ## First Construct the ref org dict
    orgDict = dict()
    engine = sql.create_engine('sqlite:///' + refDB)
    csr = engine.connect()
    query = csr.execute('Select asmID,speciesName,type_strain FROM organisms '
                        'WHERE asmID in ("{}")'.format('" , "'.join(chain(reforgs,ogorgs))))
    for asmID,name,ts in query:
        orgDict[asmID] = (cleanName(name),ts)
    csr.close()
    queryOrgs = set()
    refOrgs = dict()
    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)
    for hmmhit in idDict:
        log.info('Writing {} Seqs'.format(hmmhit))
        seqs = []
        ## get the query sequences from the
        engine = sql.create_engine('sqlite:///' + queryDB)
        csr = engine.connect()
        query = csr.execute('Select orgname,naseq,aaseq FROM Seqs WHERE seqid in ({})'.format(
            ' , '.join(str(x) for x in idDict[hmmhit]['query'])))
        for orgname,naseq,aaseq in query:
            seqs.append(('QS--'+orgname,naseq,aaseq))
            queryOrgs.add('QS--'+orgname)
        csr.close()

        engine = sql.create_engine('sqlite:///' + refDB)
        csr = engine.connect()
        query = csr.execute('Select orgname,naseq,aaseq FROM Seqs WHERE seqid in ({})'.format(
            ' , '.join(str(x) for x in idDict[hmmhit]['ref'])))
        for asmID,naseq,aaseq in query:
            name,ts = orgDict.get(asmID,[asmID,False])
            if asmID in ogorgs:
                if ts:
                    refName = 'OG-TS--{}-{}'.format(asmID,name)
                else:
                    refName = 'OG--{}-{}'.format(asmID, name)
            elif ts:
                refName = 'TS--{}-{}'.format(asmID,name)
            else:
                refName = '{}-{}'.format(asmID, name)
            seqs.append((refName,naseq,aaseq))
            refOrgs[asmID] = refName
        csr.close()

        if addAA:
            with open(os.path.join(outputFolder,'{}.fna'.format(hmmhit)),'w') as ntfasta, \
                open(os.path.join(outputFolder,'{}.faa'.format(hmmhit)),'w') as aafasta:
                for name,naseq,aaseq in seqs:
                    ntfasta.write('>{}\n{}\n'.format(name,naseq))
                    aafasta.write('>{}\n{}\n'.format(name,aaseq))
        else:
            with open(os.path.join(outputFolder, '{}.fna'.format(hmmhit)), 'w') as ntfasta:
                for name,naseq,aaseq in seqs:
                    ntfasta.write('>{}\n{}\n'.format(name,naseq))
    if RNA:
        engine = sql.create_engine('sqlite:///' + queryDB)
        csr = engine.connect()
        query = csr.execute('Select orgname,seq FROM RNAHits '
                            'WHERE rnahit="16S_rRNA" AND '
                            'orgname in ("{}")'.format('" , "'.join(queryOrgs)))
        queryRNAs = [x for x in query]
        csr.close()
        ## Only do main DB query if there is 16s information for all of the query sequences

        engine = sql.create_engine('sqlite:///' + refDB)
        csr = engine.connect()
        query = csr.execute('Select orgname,seq FROM RNAHits '
                                'WHERE rnahit="16S_rRNA" AND '
                                'orgname in ("{}")'.format('" , "'.join(refOrgs.keys())))
        refRNAs = [x for x in query]

        if not len(refRNAs) == len(refOrgs):
            log.warning('Did not find 16s information for all of the reference organisms')
        if not len(queryRNAs) == len(queryOrgs):
            log.warning('Did not find 16s information for all of the query organisms')
        with open(os.path.join(outputFolder, '16s.fna'), 'w') as rnafasta:
            for orgname, seq in queryRNAs:
                rnafasta.write('>{}\n{}\n'.format(orgname, seq))
            for orgname, seq in refRNAs:
                rnafasta.write('>{}\n{}\n'.format(refOrgs[orgname], seq))
    else:
        log.info('Skipping 16s genes...')
    return