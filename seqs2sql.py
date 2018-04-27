import re, sys, os, time,gzip, logging,sqlalchemy as sql
from sqlalchemy import MetaData
from Bio import Seq,SeqIO
from Bio.Seq import Seq
from bx.intervals.intersection import IntervalTree, Interval
from sortedcontainers import SortedDict
from operator import itemgetter
from Bio.Alphabet import generic_protein,generic_dna
from sqlalchemy.orm import sessionmaker

from logging import handlers
from logging.handlers import RotatingFileHandler

log = logging.getLogger('')
log.setLevel(logging.DEBUG)
format = logging.Formatter("%(asctime)s - %(name)s - %(levelnamee)s - %(message)s")
ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(format)
log.addHandler(ch)

def addrRnaHits(db,barrnapFastas,pythonSort=False):
    '''
    :param db: path to sql database
    :param barrnapFastas: list of fastafiles created from
    :return:
    '''

    records = []

    objHandler = MetaData()
    engine = sql.create_engine('sqlite:///'+db)
    csr = engine.connect()
    objHandler.reflect(bind=engine)
    rRnaTable = objHandler.tables['RNAHits']
    rnaHits = []
    if not pythonSort:
        for fastaFile in barrnapFastas:
            fastaHandle = SeqIO.parse(fastaFile,'fasta')
            for rnaHit in fastaHandle:
                orgname, rnahit, eval, location, strand,partial = rnaHit.id.strip().split(':')
                eval = float(eval)
                loc_start,loc_end = location.split('-')
                loc_start = int(loc_start)
                loc_end  = int(loc_end)
                dir_dict = {'1':'+','-1':'-'}
                loc_strand = dir_dict[strand]
                seq = str(rnaHit.seq)
                if partial == 'partial':
                    partialFlag = 1
                else:
                    partialFlag = 0
                rnaHits.append(
                    {'rnahit': rnahit, 'orgname': orgname, 'eval': eval, 'loc_start': loc_start, 'loc_end': loc_end,
                     'loc_strand': loc_strand, 'seq': seq,'partialFlag':partialFlag})
        csr.execute(rRnaTable.insert(), rnaHits)
        ## get the best evalues for each of the RNA subunits, prioritize non-partial hits
        ## NOTE: this does not guarantee unique sequences but this is a barrnap limitation as there are dissimilar sequences
        ## that have the same score
        csr.execute('DELETE FROM "RNAHits" WHERE (seqid,orgname,rnahit,eval,partialFlag) '
                    'NOT IN (SELECT t.seqid,t.orgname,t.rnahit,MIN(t.eval),MIN(t.partialFlag) '
                    'FROM "RNAhits" t GROUP BY rnahit, orgname);')
        csr.close()
    else:
        for fastaFile in barrnapFastas:
            rnaRecs = dict()
            fastaHandle = SeqIO.parse(fastaFile,'fasta')
            for rnaHit in fastaHandle:
                orgname, rnahit, eval, location, strand, partial = rnaHit.id.strip().split(':')
                candSSUs = rnaRecs.setdefault(rnahit,[])
                eval = float(eval)
                loc_start, loc_end = location.split('-')
                loc_start = int(loc_start)
                loc_end = int(loc_end)
                dir_dict = {'1': '+', '-1': '-'}
                loc_strand = dir_dict[strand]
                seq = str(rnaHit.seq)
                if partial == 'partial':
                    partialFlag = 1
                else:
                    partialFlag = 0
                candSSUs.append(
                    {'rnahit': rnahit, 'orgname': orgname, 'eval': eval, 'loc_start': loc_start, 'loc_end': loc_end,
                     'loc_strand': loc_strand, 'seq': seq,'partialFlag':partialFlag})
            for candSSUs in rnaRecs.values():
                candSSUs.sort(key=lambda x: x['eval'])
                candSSUs.sort(key=lambda x: sum(1 for x in x['seq'] if x == 'N'))
                candSSUs.sort(key=lambda x: x['partialFlag'])
                rnaHits.append(candSSUs[0])
        csr.execute(rRnaTable.insert(), rnaHits)
    ## remove duplicates pick latest entry
    csr.execute('DELETE FROM "RNAHits" WHERE (seqid) '
                    'NOT IN (SELECT MAX(t.seqid)'
                    'FROM "RNAhits" t GROUP BY rnahit, orgname,loc_start);')
    csr.close()
    return

def addCoreHmmHits(db,hmmHits):
    objHandler = MetaData()
    engine = sql.create_engine('sqlite:///'+db)
    csr = engine.connect()
    objHandler.reflect(bind=engine)
    hmmTable = objHandler.tables['HMMhits']
    csr.execute(hmmTable.insert(), hmmHits)
    ### remove duplicates pick latest entry
    csr.execute(
        'DELETE FROM "HMMhits" WHERE hmmhitID NOT IN (SELECT max(t.hmmhitID) FROM "HMMhits" t '
        'GROUP BY hmmhit,seqid,hmmstart,hmmend)')
    csr.close()


def calculate_window(coordinates):
    return coordinates[-1]-coordinates[0]+1.

def resolve_conflicts(pfam_hit_dict,minDomSize = 9,verbose=False):
    '''
    :param pfam_hit_dict: dictionary of hits for the gene in the following format
    hit start,hit end : int
    hit id : str
    score, model coverage percent : float
    {(hit start,hit end):('hit id',score,model coverage percent)}
    :param minDomSize: int, the minimum window size that will be considered a domain
    :return:
    a sorted dictionary with the position of the hit as the keys and ('hit id',score,model coverage percent)
    '''
    # initialize output
    gene_hits = SortedDict()
    redoFlag = True
    while redoFlag:
        if verbose: print("Sorting through intervals", pfam_hit_dict)
        redoFlag = False
        intervals_scores = [(key,value['score']) for key,value in pfam_hit_dict.items()]
        # sort intervals from pfam hits by score and place the highest score first
        intervals_scores.sort(key=itemgetter(1),reverse=True)
        # initialize intersect tree for quick overlap search
        intersectTree = IntervalTree()
        #add the intervals with the highest scores first
        for (interval,score) in intervals_scores:
            intervalStart = interval[0]
            intervalEnd = interval[1]
            intervalLength = intervalEnd-intervalStart+1
            # if the interval is less than the minimum domain size don't bother
            if intervalLength > minDomSize:
                intersectingIntervals = [(x.start,x.end) for x in intersectTree.find(intervalStart,intervalEnd)]
                overLapFlag = False
                # for every interval that you're adding resolve the overlapping intervals
                while len(intersectingIntervals) > 0 and intervalLength > 1:

                    start,end = intersectingIntervals[0]

                    # interval completely covers existing coverage, break up into two intervals and redo the process
                    if (intervalStart < start and intervalEnd > end):
                        if verbose: print("Split Interval", interval,intersectingIntervals, pfam_hit_dict[interval])
                        left_scale = calculate_window((intervalStart,start-1))/intervalLength
                        right_scale = calculate_window((end+1,intervalEnd))/intervalLength
                        pfam_hit_dict[(intervalStart,start-1)] = (pfam_hit_dict[interval][0],
                                                                  pfam_hit_dict[interval][1],
                                                                  pfam_hit_dict[interval][2] * left_scale)
                        pfam_hit_dict[(end+1,intervalEnd)] = (pfam_hit_dict[interval][0],
                                                              pfam_hit_dict[interval][1],
                                                              pfam_hit_dict[interval][2] * right_scale)
                        # delete original hit and iterate
                        del pfam_hit_dict[interval]
                        redoFlag = True
                        break
                    else:
                        #completely in the interval
                        if (intervalStart >= start and intervalEnd <= end):
                            #if completely overlapping then ignore since we already sorted by score
                            overLapFlag = True
                            break
                        #intersection covers the left hand side of the interval
                        elif intervalStart >= start:
                            intervalStart = end + 1
                        #intersection covers the right hand side of the interval
                        elif intervalEnd <= end:
                            intervalEnd = start - 1
                            # recalculate the interval length and see if there are still intersecting intervals
                        intervalLength = intervalEnd-intervalStart+1
                        intersectingIntervals = [(x.start,x.end) for x in intersectTree.find(intervalStart,intervalEnd)]

                if redoFlag:
                    if verbose: print("Exiting For Loop to Reinitialize",pfam_hit_dict)
                    break
                # if loop did not break because of an overlap add the annotation after resolving overlap,
                # check for minimum length after you merge intervals
                elif not overLapFlag and intervalLength > minDomSize:
                    if verbose: print("Adding Hit",(intervalStart,intervalEnd),pfam_hit_dict[interval][0])
                    # scale the hitCoverage based on the reduction this works since interval is a tuple and isn't mutated
                    hitCoverage = pfam_hit_dict[interval][2]*(intervalLength/(interval[1]-interval[0]+1.))
                    gene_hits[(intervalStart,intervalEnd)] = (pfam_hit_dict[interval][0],
                                                              pfam_hit_dict[interval][1],
                                                              hitCoverage)
                    intersectTree.add_interval(Interval(float(intervalStart),intervalEnd))
    if verbose: print("Merging Hits")
    # Merge Windows Right Next to one another that have the same pFam ID,
    # redoFlag: need to restart the process after a successful merge
    redoFlag = True
    while redoFlag:
        for idx in range(len(gene_hits)-1):
            left_hit = gene_hits.keys()[idx]
            right_hit = gene_hits.keys()[idx+1]
            left_window_size = calculate_window(left_hit)
            right_window_size = calculate_window(right_hit)
            merged_window_size = calculate_window((left_hit[0],right_hit[1]))
            new_coverage = (gene_hits[left_hit][2] + gene_hits[right_hit][2])*\
                           (left_window_size+ right_window_size)/merged_window_size
            # Will merge a hit under the following conditions:
            # 1. Gap between the two hits is less than the minimum domain
            # 2. Cumulative coverage of the two hits is less than 1 (this avoids merging repeats together)
            if right_hit[0]-left_hit[1] < minDomSize and gene_hits[left_hit][0] == gene_hits[right_hit][0] \
                    and new_coverage < 1:
                gene_hits[(left_hit[0],right_hit[1])] = (gene_hits[left_hit][0],
                                                         left_window_size/merged_window_size * gene_hits[left_hit][1] +
                                                         right_window_size/merged_window_size * gene_hits[right_hit][1],
                                                         new_coverage)
                redoFlag = True
                del gene_hits[left_hit]
                del gene_hits[right_hit]
                if verbose: print("Merged", left_hit,right_hit)
                break
        else:
            redoFlag = False
    if verbose: print("Deleting Domains Under Minimum Domain Size")
    # Finally check if any of the domains are less than the minimum domain size
    keysToDelete = [coordinates for coordinates in gene_hits.keys() if calculate_window(coordinates) < minDomSize]
    for key in keysToDelete:
        del gene_hits[key]
        if verbose: print("Deleting",key)
    if verbose: print("Final Annotation", gene_hits)
    return gene_hits

def parseHmmSearchDomtbl(domtblPath,minDomSize,hmmLib,eval_cutoff=1e-5,cutoff_score=25,verbose = False,speciesFilter = None,
                         proteinFilter = None):
    """
    :param path: path to pfam hits that you want parsed

    :return: {protein:gene hit dictionary (output of resolve_conflicts)}
    hit dictionary of pfam hits that can be linked to proteins for their pfam domain annotation
    """
    hit_dict = dict()
    with open(domtblPath) as hmmfile:
        for line in hmmfile:
            if 'Query file' in line:
                hmmLib = os.path.split(line.strip().split()[1])[1].split('.hmm')[0]
            elif not line.startswith('#'):
                hit = line.strip().split()
                score = float(hit[13])
                orgname,seqid = hit[0].split('|')
                evalue = float(hit[12])

                if score >= cutoff_score  and evalue <= eval_cutoff:
                    #check if you're done parsing the hits of the old gene and annotate the gene

                    # check if any hits pass the filter
                    hmmhit = hit[3]
                    hmmstart = int(hit[15])
                    hmmend = int(hit[16])
                    hmmlen = int(hit[5])
                    hmmcov = (hmmend-hmmstart+1.)/hmmlen
                    geneStart = int(hit[17])
                    geneEnd = int(hit[18])
                    genelen = int(hit[2])
                    genecov = (geneEnd-geneStart+1.)/genelen
                    if (orgname,seqid,geneStart,geneEnd,hmmhit) in hit_dict and score > hit_dict[(orgname,seqid,geneStart,geneEnd,hmmhit)]['score']:
                        hit_dict[(orgname,seqid,geneStart,geneEnd,hmmhit)] = {'hmmhit':hmmhit,'score':score,'hmmstart':hmmstart,
                                                         'hmmend':hmmend,'hmmlen':hmmlen,'hmmcov':hmmcov,'evalue':evalue,
                                                         'geneStart':geneStart,'geneEnd':geneEnd,'genecov':genecov,
                                                         'hmmLib':hmmLib,'orgname':orgname,'seqid':seqid,'genelen':genelen}
                    elif (orgname,seqid,geneStart,geneEnd,hmmhit) not in hit_dict:
                        hit_dict[(orgname,seqid,geneStart,geneEnd,hmmhit)] = {'hmmhit':hmmhit,'score':score,'hmmstart':hmmstart,
                                                         'hmmend':hmmend,'hmmlen':hmmlen,'hmmcov':hmmcov,'evalue':evalue,
                                                         'geneStart':geneStart,'geneEnd':geneEnd,'genecov':genecov,
                                                         'hmmLib':hmmLib,'orgname':orgname,'seqid':seqid,'genelen':genelen}
    return list(hit_dict.values())


def writeFasta(db,fastaHandle,nuc=False,latest=False):
    objHandler = MetaData()
    engine = sql.create_engine('sqlite:///' + db)
    csr = engine.connect()
    objHandler.reflect(bind=engine)
    seqsTable = objHandler.tables['Seqs']
    if latest:
        maxval = sql.select([sql.func.max(objHandler.tables['Seqs'].c.lastscan)])
        result = csr.execute(maxval)
        maxval = result.fetchone()[0]
    else:
        maxval = 0

    if nuc:
        seqsRequest = sql.select([seqsTable.c.orgname,seqsTable.c.seqid,seqsTable.c.naseq]).where(seqsTable.c.lastscan >= maxval)
    else:
        seqsRequest = sql.select([seqsTable.c.orgname, seqsTable.c.seqid, seqsTable.c.aaseq]).where(
            seqsTable.c.lastscan >= maxval)

    seqs = csr.execute(seqsRequest)
    with open(fastaHandle,'a') as outfile:
        for result in seqs:
            asmID,seqid,sequence = result
            outfile.write('>{}|{}\n{}\n'.format(asmID,seqid,sequence))

def importProdigalGeneFile(geneFilePath,database,orgname=False):
    '''
    imports prodigal CDS fasta into the database,
    :param geneFilePath: fasta output
    :param database: mlst database
    :param orgname: defaults to filename, better to specify
    :return:
    '''
    objHandler = MetaData()
    engine = sql.create_engine('sqlite:///'+database)
    csr = engine.connect()
    objHandler.reflect(bind=engine)

    orgTable = objHandler.tables['organisms']
    seqsTable = objHandler.tables['Seqs']

    timeStamp = int(time.time())
    recs = []
    source = os.path.split(geneFilePath)[1]
    sequence = ''
    if not orgname:
        orgname = os.path.splitext(os.path.split(geneFilePath)[1])[0]
    for line in open(geneFilePath):
        if line.startswith('>'):
            if sequence:
                #log.debug(sequence)
                dnaSeq = Seq(sequence,generic_dna)
                proteinSeq = dnaSeq.translate()
                recDict = dict(orgname=orgname,gene=internal_id,description=descr,source=source,loc_start=gene_start,
                               loc_end=gene_end,loc_strand=direction_id,acc=acc,lastscan=timeStamp,
                               naseq = str(dnaSeq),aaseq = str(proteinSeq))
                recs.append(recDict)
            sequence = ''
            lineParse = line[1:].split('#')
            rawCDSid = lineParse[0]
            cdsInfoParse = rawCDSid.split('_')
            acc = '_'.join(cdsInfoParse[:-1])
            cds_ctr = int(cdsInfoParse[-1])
            gene_start = int(lineParse[1])
            gene_end = int(lineParse[2])
            if lineParse[3] == '1':
                direction_id = '+'
            else:
                direction_id = '-'

            internal_id = "{}_CDS_{:05d}".format(acc, cds_ctr)
            descr = None
        else:
            sequence += line.strip()

    dnaSeq = Seq(sequence, generic_dna)
    proteinSeq = dnaSeq.translate()
    recDict = dict(orgname=orgname, gene=internal_id, description=descr, source=source, loc_start=gene_start,
                   loc_end=gene_end, loc_strand=direction_id, acc=acc, lastscan=timeStamp,
                   naseq=str(dnaSeq), aaseq=str(proteinSeq))
    recs.append(recDict)

    if recs:
    ## First check if organism is already in Organisms if not add
        if not csr.execute(sql.select([orgTable]).where(orgTable.c.asmID == orgname)).fetchone():
            csr.execute(orgTable.insert().values(asmID=orgname,asmVersion='1'))
        csr.execute(seqsTable.insert(),recs)
        log.info('Added {:d} Records'.format(len(recs)))
        csr.execute(
         'DELETE FROM "Seqs" WHERE seqid NOT IN (SELECT min(t.seqid) FROM "Seqs" t GROUP BY orgname,gene,description,loc_start,loc_end,loc_strand)')
        csr.close()

    return recs

def importAsmCDSfile(asmCDSpth,database,orgname=False):
    '''
    will populate the SQL database with the sequences from cds_from_genomic fasta
    :param asmCDSpth: path to file
    :param database: SQL Database Path
    :param orgname: Override organism name
    :return:
    '''
    objHandler = MetaData()
    engine = sql.create_engine('sqlite:///'+database)
    csr = engine.connect()
    objHandler.reflect(bind=engine)

    orgTable = objHandler.tables['organisms']
    seqsTable = objHandler.tables['Seqs']

    timeStamp = int(time.time())
    recs = []
    source = os.path.split(asmCDSpth)[1]
    sequence = ''
    if not orgname:
        # grabs the asmID base assuming you don't change the file name from the assembly download
        orgname = source.split('.')[0]

    for line in gzip.open(asmCDSpth):
        line = line.decode()
        if line.startswith('>'):
            if sequence:
                #log.debug(sequence)
                dnaSeq = Seq(sequence,generic_dna)
                proteinSeq = dnaSeq.translate()
                recDict = dict(orgname=orgname,gene=internal_id,description=descr,source=source,loc_start=gene_start,
                               loc_end=gene_end,loc_strand=direction_id,acc=acc,lastscan=timeStamp,
                               naseq = str(dnaSeq),aaseq = str(proteinSeq))
                recs.append(recDict)
            sequence = ''
            lineParse = line.split()
            rawCDSid = lineParse[0].split('|')[1]
            cdsInfoParse = rawCDSid.split('_')
            acc = rawCDSid.split('_cds_')[0]
            cds_ctr = int(cdsInfoParse[-1])
            descriptors = re.findall('\[{1}\w+\={1}[^=]*\]', line)
            descriptors_dict = dict()
            for match in descriptors:
                try:
                    key, value = match[1:-1].split('=')
                    descriptors_dict[key] = value
                except ValueError:
                    log.info("Wasn't able to parse descriptor, {}".format(match[1:-1]))
                    pass
            location = re.findall('\d+', descriptors_dict['location'])
            gene_start = int(location[0])
            gene_end = int(location[-1])
            if 'complement' in descriptors_dict['location']:
                direction_id = '-'
            else:
                direction_id = '+'
            internal_id = "{}_CDS_{:05d}".format(acc, cds_ctr)
            if 'protein_id' in descriptors_dict.keys():
                descr = descriptors_dict['protein_id']
            elif 'gene' in descriptors_dict.keys():
                descr = descriptors_dict['gene']
            elif 'locus_tag' in descriptors_dict.keys():
                descr = descriptors_dict['locus_tag']
            else:
                descr = None
        else:
            sequence += line.strip()
    ## Flush final buffer
    dnaSeq = Seq(sequence, generic_dna)
    proteinSeq = dnaSeq.translate()
    recDict = dict(orgname=orgname, gene=internal_id, description=descr, source=source, loc_start=gene_start,
                   loc_end=gene_end, loc_strand=direction_id, acc=acc, lastscan=timeStamp,
                   naseq=str(dnaSeq), aaseq=str(proteinSeq))
    recs.append(recDict)
    log.info("Finished parsing file found {} new coding sequences".format(len(recs)))
    log.info("Inserting into Database".format(len(recs)))
    if recs:
    ## First check if organism is already in Organisms if not add
        if not csr.execute(sql.select([orgTable]).where(orgTable.c.asmID == orgname)).fetchone():
            csr.execute(orgTable.insert().values(asmID=orgname))
        csr.execute(seqsTable.insert(),recs)
        log.info('Added {:d} Records'.format(len(recs)))
        csr.execute(
         'DELETE FROM "Seqs" WHERE seqid NOT IN (SELECT min(t.seqid) FROM "Seqs" t GROUP BY orgname,gene,description,loc_start,loc_end,loc_strand)')
        csr.close()

    return