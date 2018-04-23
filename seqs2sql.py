import re, sys, os, time,gzip, logging,sqlalchemy as sql
from sqlalchemy import MetaData
from Bio import Seq
from Bio.Seq import Seq
from initialize_db import Organism,CDS
from Bio.Alphabet import generic_protein,generic_dna
from sqlalchemy.orm import sessionmaker

from logging import handlers
from logging.handlers import RotatingFileHandler

log = logging.getLogger('')
log.setLevel(logging.DEBUG)
format = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(format)
log.addHandler(ch)


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

    log.info("Finished parsing file found {} new coding sequences".format(len(recs)))
    log.info("Inserting into Database".format(len(recs)))
    if recs:
    ## First check if organism is already in Organisms if not add
        if not csr.execute(sql.select([orgTable]).where(orgTable.c.name == orgname)).fetchone():
            csr.execute(orgTable.insert().values(name=orgname))
        csr.execute(seqsTable.insert(),recs)
        log.info('Added {:d} Records'.format(len(recs)))
        csr.execute(
         'DELETE FROM "Seqs" WHERE seqid NOT IN (SELECT min(t.seqid) FROM "Seqs" t GROUP BY orgname,gene,description,loc_start,loc_end,loc_strand)')
        csr.close()

    return