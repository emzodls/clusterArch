import re, sys, os, time,gzip, setlog, sqlalchemy as sql
from Bio import Seq
from Bio.Alphabet import generic_protein,generic_dna
from sqlalchemy.orm import sessionmaker

global log
log = setlog.init(toconsole=True)

def importAsmCDSfile(asmCDSpth,database,orgname=False):
    '''
    will populate the SQL database with the sequences from cds_from_genomic fasta
    :param asmCDSpth: path to file
    :param database: SQL Database Path
    :param orgname: Override organism name
    :return:
    '''

    engine = sql.create_engine('sqlite:///'+database)
    csr = engine.connect()
    maxid = csr.execute('Select Max(seqid) from "Seqs"').fetchone()[0]
    timeStamp = int(time.time())
    recs = []
    source = os.path.split(asmCDSpth)[1]
    if not orgname:
        # grabs the asmID base assuming you don't change the file name from the assembly download
        orgname = source.split('.')[0]
    for line in gzip.open(asmCDSpth):
        if '>' in line:
            startFlag = False
            id = maxid + 1
            if startFlag and sequence:
                dnaSeq = Seq(sequence,generic_dna)
                proteinSeq = dnaSeq.translate()
                recs.append(id,orgname,internal_id,descr,source,gene_start,gene_end,direction_id,
                            acc,timeStamp,dnaSeq.seq,proteinSeq.seq)
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
                    print(match[1:-1])
                    pass
            location = re.findall('\d+', descriptors_dict['location'])
            gene_start = int(location[0])
            gene_end = int(location[-1])
            if 'complement' in descriptors_dict['location']:
                direction_id = '-'
            else:
                direction_id = '+'
            internal_id = "%s_CDS_%.5i" % (species_id, cds_ctr)
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

    ## First check if organism is already in Organisms if not add
    if not csr.execute('Select name,id from organisms where organisms.name = "{}"'.format(orgname)):
        maxOrgID = csr.execute('Select Max(seqid) from "organisms"').fetchone()[0]
        csr.execute('Insert into "organisms" (name,id) VALUES ("{}",{})'.format(orgname,maxOrgID+1))
