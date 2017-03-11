# Copyright (C) 2017 Emmanuel LC. de los Santos
# University of Warwick
# Warwick Integrative Synthetic Biology Centre
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
'''
    This file is part of clusterTools.

    clusterTools is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    clusterTools is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with clusterTools.  If not, see <http://www.gnu.org/licenses/>.
'''
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein,generic_dna
from clusterTools import clusterAnalysis
from random import random
import sys,subprocess,os,platform
import gzip,re,urllib

### To fix the file paths in windows ###
# from http://stackoverflow.com/questions/23598289/how-to-get-windows-short-file-name-in-python
if platform.system() == 'Windows':
    import ctypes
    from ctypes import wintypes
    _GetShortPathNameW = ctypes.windll.kernel32.GetShortPathNameW
    _GetShortPathNameW.argtypes = [wintypes.LPCWSTR, wintypes.LPWSTR, wintypes.DWORD]
    _GetShortPathNameW.restype = wintypes.DWORD

    def get_short_path_name(long_name):
        """
        Gets the short path name of a given long path.
        http://stackoverflow.com/a/23598461/200291
        """
        output_buf_size = 0
        while True:
            output_buf = ctypes.create_unicode_buffer(output_buf_size)
            needed = _GetShortPathNameW(long_name, output_buf, output_buf_size)
            if output_buf_size >= needed:
                return output_buf.value
            else:
                output_buf_size = needed
###########################################

def execute(commands, input=None,shell=False):
    "Execute commands in a system-independent manner"

    if input is not None:
        stdin_redir = subprocess.PIPE
    else:
        stdin_redir = None

    try:
        proc = subprocess.Popen(commands, stdin=stdin_redir,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,shell=False)
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
                if protein_id not in existingGenes:
                    genesToAdd.append(protein_id)
                    geneDict[protein_id] = str(genbank_seq)

    # fasta
    elif extension == 'fasta' or extension == 'fa':
        genes = SeqIO.parse(open(SeqFilePath), "fasta")
        for gene in genes:
            existingGenes = geneDict.keys()
            if gene.id in existingGenes:
                geneName = gene.id + '.' + SeqFilePath.split(extension)[0].split('/')[-1]
            else:
                geneName = gene.id
            if geneName not in existingGenes:
                genesToAdd.append(geneName)
                geneDict[geneName] = str(gene.seq)
    return genesToAdd,geneDict

def parseHMMfile(HMMfilePath,HMMdict):
    hmmsToAdd = [x.strip().split()[-1] for x in open(HMMfilePath) if 'NAME' in x]
    for hmm in hmmsToAdd:
        if hmm not in HMMdict.keys():
            HMMdict[hmm] = HMMfilePath

    return hmmsToAdd,HMMdict

def MakeBlastDB(makeblastdbExec,dbPath,outputDir,outDBName):
    if platform.system() == 'Windows':
        dbPath = get_short_path_name(dbPath)
        outputDir = get_short_path_name(outputDir)
    command = [makeblastdbExec, '-in', dbPath, '-dbtype', "prot",'-out',os.path.join(outputDir,outDBName)]
    out, err, retcode = execute(command)
    if retcode != 0:
        print('makeblastDB failed with retcode %d: %r' % (retcode, err))
    return out,err,retcode

def runBLAST(blastExec,inputFastas,outputDir,searchName,dbPath,eValue='1E-05'):
    if platform.system() == 'Windows':
        path, outputDBname = os.path.split(dbPath)
        dbPath = os.path.join(get_short_path_name(path),outputDBname)
        inputFastas = get_short_path_name(inputFastas)
        outputDir = get_short_path_name(outputDir)
    command = [blastExec, "-db", dbPath, "-query", inputFastas, "-outfmt", "6", "-max_target_seqs", "10000", "-evalue",
               eValue, "-out", os.path.join(outputDir,"{}_blast_results.out".format(searchName))]
    out, err, retcode = execute(command)
    if retcode != 0:
        print('BLAST failed with retcode %d: %r' % (retcode, err))
    return out,err,retcode

def runHmmCheck(hmmSearchExec,hmmDBase):
    if platform.system() == 'Windows':
        hmmDBase = get_short_path_name(hmmDBase)
    command = [hmmSearchExec,hmmDBase,'testHMM.fasta']
    out,err,retcode = execute(command)
    if retcode != 0:
        return False
    else:
        return True

def runHmmBuild(hmmBuildExec,inFile,outFile):
    path,hmmName = os.path.split(outFile)
    hmmName = hmmName.split('.hmm')[0]
    if platform.system() == 'Windows':
        inFile = get_short_path_name(inFile)
        outFile = get_short_path_name(outFile)
    command = [hmmBuildExec,'-n',hmmName,outFile,inFile]
    if platform.system() == 'Windows':
        print(platform.system())
        out,err,retcode = execute(command,shell=True)
    else:
        out, err, retcode = execute(command)
    if retcode != 0:
        print(err,retcode)
        return False

    else:
        return True

def runHmmsearch(hmmSearchExec,hmmDBase,outputDir,searchName,dbPath,eValue='1E-05'):
    if platform.system() == 'Windows':
        hmmDBase = get_short_path_name(hmmDBase)
        outputDir = get_short_path_name(outputDir)
        dbPath = get_short_path_name(dbPath)
    command = [hmmSearchExec,'--domtblout', os.path.join(outputDir,'{}_hmmSearch.out'.format(searchName)), '--noali',
               '-E', eValue, hmmDBase, dbPath]
    out, err, retcode = execute(command)
    if retcode != 0:
        print('hmmsearch failed with retcode %d: %r' % (retcode, err))
    return out,err,retcode

def generateInputFasta(forBLAST,outputDir,searchName):
    with open(os.path.join('%s' % outputDir,'{}_gene_queries.fa'.format(searchName)),'w') as outfile:
        for gene in forBLAST.keys():
            prot_entry = SeqRecord(Seq(forBLAST[gene],generic_protein), id=gene,
                               description='%s' % (gene))
            SeqIO.write(prot_entry,outfile,'fasta')

def generateHMMdb(hmmFetchExec,hmmDict,hmmSet,outputDir,searchName):
    errFlag = False
    failedToFetch = set()
    if platform.system() == 'Windows':
        outputDir = get_short_path_name(outputDir)
    with open(os.path.join('%s' % outputDir,'{}_hmmDB.hmm'.format(searchName)),'wb') as outfile:
        for hmm in hmmSet:
            if platform.system() == 'Windows':
                hmmSource = get_short_path_name(hmmDict[hmm])
            else:
                hmmSource = hmmDict[hmm]
            out, err, retcode = execute([hmmFetchExec, hmmSource, hmm])
            if retcode == 0:
                outfile.write(out)
            else:
                print('hmmfetch failed with retcode %d: %r' % (retcode, err))
                errFlag = True
                failedToFetch.add(hmm)
    return errFlag,failedToFetch

def processSearchListOptionalHits(requiredBlastList,requiredHmmList,blastOutFile,blastEval,hmmOutFile,hmmScore, hmmDomLen,windowSize,
                                  totalHitsRequired,additionalBlastList=[],additionalHmmList=[]):
    # Gather all of the proteins, might be a memory issue...code memory friendly version with sequential filters (?)
    prots = dict()
    if requiredBlastList or additionalBlastList:
        prots = clusterAnalysis.parseBLAST(blastOutFile,prots,swapQuery=True,evalCutoff=blastEval)
    if requiredHmmList or additionalHmmList:
        prots = clusterAnalysis.parse_hmmsearch_domtbl_anot(hmmOutFile,hmmDomLen,'hmm',prots,cutoff_score=hmmScore)

    requiredBlastHitDict = dict()
    requiredHmmHitDict = dict()

    additionalBlastHitDict = dict()
    additionalHmmHitDict = dict()

    if requiredBlastList:
        requiredBlastHitDict = {hitName:set(protein for protein in prots.values() if hitName in protein.hit_dict['blast'].hits)
                  for hitName in requiredBlastList}
    if requiredHmmList:
        requiredHmmHitDict = {hmms: set(protein for protein in prots.values() if len(set(hmms) & protein.getAnnotations('hmm')) == len(hmms))
               for hmms in requiredHmmList}
    if additionalBlastList:
        additionalBlastHitDict = {hitName:set(protein for protein in prots.values() if hitName in protein.hit_dict['blast'].hits)
                  for hitName in additionalBlastList}
    if additionalHmmList:
        additionalHmmHitDict = {hmms: set(protein for protein in prots.values() if len(set(hmms) & protein.getAnnotations('hmm')) == len(hmms))
               for hmms in additionalHmmList}
    requiredHitDict = {**requiredBlastHitDict,**requiredHmmHitDict}
    additionalHitDict = {**additionalBlastHitDict, **additionalHmmHitDict}

    #need this for repeat domains

    requiredHitList = requiredBlastList+requiredHmmList
    additionalHitList = additionalBlastList+additionalHmmList
    numReqHits = len(requiredHitList)

    hitDict = {**requiredHitDict, **additionalHitDict}

    hitProteins = set()
    hitProteins.update(*requiredHitDict.values())
    hitProteins.update(*additionalHitDict.values())

    putativeClusters = clusterAnalysis.cluster_proteins(hitProteins,windowSize)

    assert totalHitsRequired >= numReqHits

    numExtraHitsNeeded = totalHitsRequired - numReqHits

    filteredClusters = dict()
    for species,clusters in putativeClusters.items():
        for cluster in clusters:
            clusterProts = set(protein for protein in cluster)
            # first term checks required hits, second term checks extra hits 3rd term checks that there are enough genes for a unique values
            if (sum(1 for hitID in requiredHitList if len(clusterProts & requiredHitDict[hitID]) >= 1) == numReqHits) \
                    and (sum(1 for hitID in additionalHitList if len(clusterProts & additionalHitDict[hitID]) >= 1) >= numExtraHitsNeeded)\
                    and (len(clusterProts) >= totalHitsRequired):
                filteredClusters[(species,cluster.location[0],cluster.location[1])] = \
                {hitQuery:[protein.name for protein in (hitSet & clusterProts)] for hitQuery,hitSet in hitDict.items()}
    return filteredClusters

def processGbkDivFile(gbkDivFile,database,guiSignal=None):
    ## unzip gbkDivFile
    try:
        genbankHandle = SeqIO.parse(gzip.open(gbkDivFile,mode='rt'),'genbank')
        entryIDlist = set()
        entryCtr = 1
        for genbankEntry in genbankHandle:
            genesToWrite = []
            species_id = genbankEntry.name
            # Make sure there are no collisions in dictionary
            if species_id in entryIDlist:
                species_id = '{}.clusterTools{}'.format(species_id,entryCtr)
                entryCtr += 1
            entryIDlist.add(species_id)
            cds_ctr = 0
            CDS_list = (feature for feature in genbankEntry.features if feature.type == 'CDS')
            if guiSignal:
                guiSignal.emit(species_id)
            for CDS in CDS_list:
                cds_ctr += 1
                direction = CDS.location.strand
                internal_id = "%s_CDS_%.5i".format(species_id, cds_ctr)
                protein_id = internal_id
                # Ensure that you don't get negative values, Biopython parser will not ignore slices that are greater
                # than the entry so you don't need to worry about the other direction
                internal_id = "%s_CDS_%.5i" % (species_id, cds_ctr)
                protein_id = internal_id

                gene_start = max(0, CDS.location.nofuzzy_start)
                gene_end = max(0, CDS.location.nofuzzy_end)
                # Try to find a common name for the promoter, otherwise just use the internal ID
                if 'protein_id' in CDS.qualifiers.keys():
                    protein_id = CDS.qualifiers['protein_id'][0]
                elif 'locus_tag' in CDS.qualifiers.keys():
                    protein_id = CDS.qualifiers['locus_tag'][0]

                if 'translation' in CDS.qualifiers.keys():
                    prot_seq = Seq(CDS.qualifiers['translation'][0])
                    if direction == 1:
                        direction_id = '+'
                    else:
                        direction_id = '-'
                else:
                    genbank_seq = CDS.location.extract(genbankEntry)
                    if protein_id == internal_id:
                        for feature in genbank_seq.features:
                            if 'locus_tag' in feature.qualifiers:
                                protein_id = feature.qualifiers['locus_tag'][0]
                    nt_seq = genbank_seq.seq
                    if direction == 1:
                        direction_id = '+'
                        # for protein sequence if it is at the start of the entry assume that end of sequence is in frame
                        # if it is at the end of the genbank entry assume that the start of the sequence is in frame
                        if gene_start == 0:
                            if len(nt_seq) % 3 == 0:
                                prot_seq = nt_seq.translate()
                            elif len(nt_seq) % 3 == 1:
                                prot_seq = nt_seq[1:].translate()
                            else:
                                prot_seq = nt_seq[2:].translate()
                        else:
                            prot_seq = nt_seq.translate()
                    if direction == -1:
                        direction_id = '-'
                        nt_seq = genbank_seq.seq
                        if gene_start == 0:
                            prot_seq = nt_seq.translate()
                        else:
                            if len(nt_seq) % 3 == 0:
                                prot_seq = nt_seq.translate()
                            elif len(nt_seq) % 3 == 1:
                                prot_seq = nt_seq[:-1].translate()
                            else:
                                prot_seq = nt_seq[:-2].reverse_complement().translate()
                # Write protein file
                if len(prot_seq) > 0 and '*' not in prot_seq[:-1]:
                    prot_entry = SeqRecord(prot_seq, id='%s|%i-%i|%s|%s|%s' % (species_id, gene_start + 1,
                                                                               gene_end, direction_id,
                                                                               internal_id, protein_id),
                                           description='%s in %s' % (protein_id, species_id))
                    genesToWrite.append(prot_entry)
            with open(database,'a') as outfileHandle:
                SeqIO.write(genesToWrite,outfileHandle,'fasta')
        return entryIDlist
    except Exception as e:
        print(e)
        if guiSignal:
            guiSignal.emit('Failed')
        raise Exception

def proccessGbks(taskList,outputDir,guiSignal=None):
    # make sure species list is unique
    speciesList = set()
    failedToProcess = []
    for gbkFile in taskList:
        try:
            genbank_entries = SeqIO.parse(open(gbkFile), "genbank")
            path,fileName = os.path.split(gbkFile)
            species_base,ext = os.path.splitext(fileName)
            if guiSignal:
                guiSignal.emit(fileName)
            CDS_prot_outfile_name = outputDir
            cds_ctr = 0
            entry_ctr = 1
            # See if user wants a different name
            for genbank_entry in genbank_entries:
                clusterNumber = None
                prot_seqs = []
                species_id = genbank_entry.id
                if species_id in speciesList:
                    species_id = species_base + '.entry%.3i' % entry_ctr
                # check for uniqueness, if there is already an entry on the list insert random number
                ## Check if it is an antismash file
                clusters = [cluster for cluster in genbank_entry.features if cluster.type == 'cluster']
                ## if it is specifically only has 1 antismash cluster, tag it as such with the species ID and clustertype
                if len(clusters) == 1 and entry_ctr == 1:
                    cluster = clusters[0]
                    try:
                        clusterNumber = cluster.qualifiers['note'][0].split(':')[1].strip()
                    except:
                        clusterNumber = None
                    if clusterNumber:
                        species_id += '.antismashCluster{}'.format(clusterNumber)
                    if 'product' in cluster.qualifiers.keys():
                        productID = cluster.qualifiers['product'][0].split()
                        productID = ''.join(productID)
                        # only get top 4
                        productID = productID.split('-')
                        endIdx = min(len(productID),4)
                        productID = '-'.join(productID[:endIdx])
                        species_id += '.{}'.format(productID).strip()
                species_id = species_id.strip()
                if species_id in speciesList:
                    splitSpecies = species_id.split('.entry')[0]
                    species_id = splitSpecies + '%.5i.entry%.3i'.format((random()*10000),entry_ctr)
                    speciesList.add(species_id)
                else:
                    speciesList.add(species_id)

                CDS_list = (feature for feature in genbank_entry.features if feature.type == 'CDS')
                for CDS in CDS_list:
                    cds_ctr += 1
                    direction = CDS.location.strand
                    # Ensure that you don't get negative values, Biopython parser will not ignore slices that are greater
                    # than the entry so you don't need to worry about the other direction
                    internal_id = "%s_CDS_%.5i" % (species_id, cds_ctr)
                    protein_id = internal_id

                    gene_start = max(0, CDS.location.nofuzzy_start)
                    gene_end = max(0, CDS.location.nofuzzy_end)

                    # Try to find a common name for the promoter, otherwise just use the internal ID
                    if 'protein_id' in CDS.qualifiers.keys():
                        protein_id = CDS.qualifiers['protein_id'][0]
                    elif 'locus_tag' in CDS.qualifiers.keys():
                        protein_id = CDS.qualifiers['locus_tag'][0]
                    if 'translation' in CDS.qualifiers.keys():
                        prot_seq = Seq(CDS.qualifiers['translation'][0])
                        if direction == 1:
                            direction_id = '+'
                        else:
                            direction_id = '-'
                    else:
                        genbank_seq = CDS.location.extract(genbank_entry)
                        nt_seq = genbank_seq.seq
                        if direction == 1:
                            direction_id = '+'
                            # for protein sequence if it is at the start of the entry assume that end of sequence is in frame
                            # if it is at the end of the genbank entry assume that the start of the sequence is in frame
                            if gene_start == 0:
                                if len(nt_seq) % 3 == 0:
                                    prot_seq = nt_seq.translate()
                                elif len(nt_seq) % 3 == 1:
                                    prot_seq = nt_seq[1:].translate()
                                else:
                                    prot_seq = nt_seq[2:].translate()
                            else:
                                prot_seq = nt_seq.translate()
                        if direction == -1:
                            direction_id = '-'

                            nt_seq = genbank_seq.seq

                            if gene_start == 0:
                                prot_seq = nt_seq.translate()
                            else:
                                if len(nt_seq) % 3 == 0:
                                    prot_seq = nt_seq.translate()
                                elif len(nt_seq) % 3 == 1:
                                    prot_seq = nt_seq[:-1].translate()
                                else:
                                    prot_seq = nt_seq[:-2].reverse_complement().translate()

                    # Write protein file
                    if len(prot_seq) > 0:
                        prot_entry = SeqRecord(prot_seq, id='%s|%i-%i|%s|%s|%s' % (species_id, gene_start + 1,
                                                                                   gene_end, direction_id,
                                                                                   internal_id, protein_id))
                        prot_seqs.append(prot_entry)
                with open(CDS_prot_outfile_name, 'a') as outfile_handle:
                    SeqIO.write(prot_seqs, outfile_handle, 'fasta')
                entry_ctr += 1
        except Exception as e:
            if guiSignal:
                guiSignal.emit('Error Reading {}'.format(fileName))
            failedToProcess.append(fileName)
            print(gbkFile,e)
            pass
    return failedToProcess

def ncbiGenomeFastaParser(fastaHandle):
    # returns a fasta dictionary with entry as title and sequence as output
    sequence = ''
    id = ''
    fastaDict = dict()
    for line in fastaHandle:
        if '>' in line:
            if id and sequence:
                dnaSeq = Seq(sequence,generic_dna)
                proteinSeq = dnaSeq.translate()
                if '*' in proteinSeq[:-1]:
                    print(id,proteinSeq)
                fastaDict[id] = dnaSeq.translate()
            sequence = ''
            lineParse = line.split()
            rawCDSid = lineParse[0].split('|')[1]
            cdsInfoParse = rawCDSid.split('_')
            species_id = rawCDSid.split('_cds_')[0]
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
                protein_id = descriptors_dict['protein_id']
            elif 'gene' in descriptors_dict.keys():
                protein_id = descriptors_dict['gene']
            elif 'locus_tag' in descriptors_dict.keys():
                protein_id = descriptors_dict['locus_tag']
            else:
                protein_id = internal_id
            id = '%s|%i-%i|%s|%s|%s' % (species_id, gene_start + 1,
                                    gene_end, direction_id,
                                    internal_id, protein_id)
        else:
            sequence += line.strip()
    return fastaDict

def fastaDictToSeqRecs(fastaDict):
    return [SeqRecord(seq,id=id) for id,seq in fastaDict.items()]
def writeSeqRecs(handle,SeqRecs):
    try:
        SeqIO.write(SeqRecs,handle,'fasta')
        return True
    except:
        return False

### stack exchange http://stackoverflow.com/questions/12523586/python-format-size-application-converting-b-to-kb-mb-gb-tb
def humanbytes(B):
   'Return the given bytes as a human friendly KB, MB, GB, or TB string'
   B = float(B)
   KB = float(1024)
   MB = float(KB ** 2) # 1,048,576
   GB = float(KB ** 3) # 1,073,741,824
   TB = float(KB ** 4) # 1,099,511,627,776

   if B < KB:
      return '{0} {1}'.format(B,'Bytes' if 0 == B > 1 else 'Byte')
   elif KB <= B < MB:
      return '{0:.2f} KB'.format(B/KB)
   elif MB <= B < GB:
      return '{0:.2f} MB'.format(B/MB)
   elif GB <= B < TB:
      return '{0:.2f} GB'.format(B/GB)
   elif TB <= B:
      return '{0:.2f} TB'.format(B/TB)

if __name__ == "__main__":

    print(processSearchListOptionalHits([],[('AT_DOMAIN', 'KS_DOMAIN'),
                                      ('AT_DOMAIN', 'KS_DOMAIN'),('AT_DOMAIN', 'KS_DOMAIN'),('AT_DOMAIN', 'KS_DOMAIN')],'',
                                  1e5,'/Users/emzodls/Downloads/testClusterTools/hmmSearch.out',30,15,50000,4,additionalHmmList=[('AT_DOMAIN', 'KS_DOMAIN')]))