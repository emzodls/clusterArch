from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys,subprocess

def execute(commands, input=None):
    "Execute commands in a system-independent manner"

    if input is not None:
        stdin_redir = subprocess.PIPE
    else:
        stdin_redir = None

    try:
        proc = subprocess.Popen(commands, stdin=stdin_redir,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate(input=input)
        retcode = proc.returncode
        return out, err, retcode
    except OSError as e:
        print("%r %r returned %r" % (commands, input[:40] if input is not None else None, e))
        raise

def parseSeqFile(SeqFilePath,geneDict):
    extension = SeqFilePath.split('.')[-1]
    # Bring through different parsing workflows based on extension
    genesToAdd = []
    # genbank
    if extension == 'gbk':
        genbank_entries = SeqIO.parse(open(SeqFilePath), "genbank")
        cds_ctr = 0
        for genbank_entry in genbank_entries:
            CDS_list = (feature for feature in genbank_entry.features if feature.type == 'CDS')
            species_id = genbank_entry.id
            for CDS in CDS_list:
                cds_ctr += 1
                direction = CDS.location.strand
                # Ensure that you don't get negative values, Biopython parser will not ignore slices that are greater
                # than the entry so you don't need to worry about the other direction
                internal_id = "%s_CDS_%.5i" % (species_id, cds_ctr)
                protein_id = internal_id
                genbank_seq = CDS.location.extract(genbank_entry)

                # Try to find a common name for the promoter, otherwise just use the internal ID
                if 'protein_id' in CDS.qualifiers.keys():
                    protein_id = CDS.qualifiers['protein_id'][0]
                else:
                    for feature in genbank_seq.features:
                        if 'locus_tag' in feature.qualifiers:
                            protein_id = feature.qualifiers['locus_tag'][0]
                existingGenes = geneDict.keys()
                if protein_id in existingGenes:
                    protein_id = protein_id + '_' + internal_id
                genesToAdd.append(protein_id)
                geneDict[protein_id] = str(genbank_seq)

    # fasta
    elif extension == 'fasta' or extension == 'fa':
        genes = SeqIO.parse(open(SeqFilePath), "fasta")
        for gene in genes:
            existingGenes = geneDict.keys()
            if gene.id in existingGenes:
                geneName = gene.id + '_' + SeqFilePath.split(extension)[-1].split('/')[-1]
            else:
                geneName = gene.id
            genesToAdd.append(geneName)
            geneDict[geneName] = str(gene.seq)
    return genesToAdd,geneDict

def parseHMMfile(HMMfilePath,HMMdict):
    hmmsToAdd = [x.strip().split()[-1] for x in open(HMMfilePath) if 'NAME' in x]
    for hmm in hmmsToAdd:
        HMMdict[hmm] = HMMfilePath
    return hmmsToAdd,HMMdict

def MakeBlastDB(makeblastdbExec,dbPath):
    command = [makeblastdbExec, "-in", dbPath, "-dbtype", "prot"]
    out, err, retcode = execute(command)
    if retcode != 0:
        print('makeblastDB failed with retcode %d: %r' % (retcode, err))
    return out,err,retcode

def runBLAST(blastExec,inputFastas,outFolder,dbPath,eValue='1E-05'):
    command = [blastExec, "-db", dbPath, "-query", inputFastas, "-outfmt", "6", "-max_target_seqs", "10000", "-evalue",
               eValue, "-out", outFolder + "_blast.out"]
    out, err, retcode = execute(command)
    if retcode != 0:
        print('BLAST failed with retcode %d: %r' % (retcode, err))
    return out,err,retcode

def runHmmsearch(hmmSearchExec,hmmDBase,outFolder,dbPath,eValue='1E-05'):
    command = [hmmSearchExec,'--domtblout', outFolder+'_hmmSearch.out', '--noali',
               '-E', eValue, hmmDBase, dbPath]
    out, err, retcode = execute(command)
    if retcode != 0:
        print('hmmsearch failed with retcode %d: %r' % (retcode, err))
    return out,err,retcode

def generateInputFasta(geneList,geneDict,outFolder):
    with open('%s/gene_queries.fa' % outFolder,'w') as outfile:
        for gene in geneList:
            prot_entry = SeqRecord(geneDict[gene], id=gene,
                               description='%s' % (gene))
            SeqIO.write(prot_entry,outfile,'fasta')

def generateHMMdb(hmmFetchExec,hmmDict,hmmSet,outFolder):
    errFlag = False
    failedToFetch = set()
    with open('%s/hmmDB.hmm' % outFolder,'wb') as outfile:
        for hmm in hmmSet:
            out, err, retcode = execute([hmmFetchExec, hmmDict[hmm], hmmDict[hmm]])
            if retcode == 0:
                outfile.write(out)
            else:
                print('hmmfetch failed with retcode %d: %r' % (retcode, err))
                errFlag = True
                failedToFetch.add(hmm)
    return errFlag,failedToFetch
