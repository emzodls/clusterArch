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
from Bio.Alphabet import generic_protein
from clusterTools import clusterAnalysis
import sys,subprocess,os

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

def MakeBlastDB(makeblastdbExec,dbPath,outputDir,outDBName):
    command = [makeblastdbExec, "-in", dbPath, "-dbtype", "prot","-out",os.path.join(outputDir,outDBName)]
    out, err, retcode = execute(command)
    if retcode != 0:
        print('makeblastDB failed with retcode %d: %r' % (retcode, err))
    return out,err,retcode

def runBLAST(blastExec,inputFastas,outputDir,dbPath,eValue='1E-05'):
    command = [blastExec, "-db", dbPath, "-query", inputFastas, "-outfmt", "6", "-max_target_seqs", "10000", "-evalue",
               eValue, "-out", outputDir + "/blast_results.out"]
    out, err, retcode = execute(command)
    if retcode != 0:
        print('BLAST failed with retcode %d: %r' % (retcode, err))
    return out,err,retcode

def runHmmsearch(hmmSearchExec,hmmDBase,outputDir,dbPath,eValue='1E-05'):
    command = [hmmSearchExec,'--domtblout', outputDir+'/hmmSearch.out', '--noali',
               '-E', eValue, hmmDBase, dbPath]
    out, err, retcode = execute(command)
    if retcode != 0:
        print('hmmsearch failed with retcode %d: %r' % (retcode, err))
    return out,err,retcode

def generateInputFasta(forBLAST,outputDir):
    with open('%s/gene_queries.fa' % outputDir,'w') as outfile:
        for gene in forBLAST.keys():
            prot_entry = SeqRecord(Seq(forBLAST[gene],generic_protein), id=gene,
                               description='%s' % (gene))
            SeqIO.write(prot_entry,outfile,'fasta')

def generateHMMdb(hmmFetchExec,hmmDict,hmmSet,outputDir):
    errFlag = False
    failedToFetch = set()
    with open('%s/hmmDB.hmm' % outputDir,'wb') as outfile:
        for hmm in hmmSet:
            out, err, retcode = execute([hmmFetchExec, hmmDict[hmm], hmm])
            if retcode == 0:
                outfile.write(out)
            else:
                print('hmmfetch failed with retcode %d: %r' % (retcode, err))
                errFlag = True
                failedToFetch.add(hmm)
    return errFlag,failedToFetch

def processSearchList(blastList,hmmList,blastOutFile,hmmOutFile,windowSize = 50000):
    # Gather all of the proteins, might be a memory issue...code memory friendly version with sequential filters (?)
    prots = dict()
    prots = clusterAnalysis.parseBLAST(blastOutFile,prots,swapQuery=True)
    prots = clusterAnalysis.parse_hmmsearch_domtbl_anot(hmmOutFile,15,'hmm',prots)

    blastDict = {hitName:set(protein for protein in prots.values() if hitName in protein.hit_dict['blast'].hits)
                  for hitName in blastList}
    hmmDict = {hmms: set(protein for protein in prots.values() if len(set(hmms) & protein.getAnnotations('hmm')) == len(hmms))
               for hmms in hmmList}
    hitDict = {**blastDict,**hmmDict}

    putativeClusters = clusterAnalysis.cluster_proteins(prots.values(),windowSize)

    filteredClusters = dict()
    for species,clusters in putativeClusters.items():
        for cluster in clusters:
            clusterProts = set(protein for protein in cluster)
            if sum(1 for hitSet in hitDict.values() if len(clusterProts & hitSet) >= 1) == len(hitDict.keys()):
                filteredClusters[(species,cluster.location[0],cluster.location[1])] = \
                {hitQuery:[protein.name for protein in (hitSet & clusterProts)] for hitQuery,hitSet in hitDict.items()}
    return filteredClusters

