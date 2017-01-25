from Bio import SeqIO

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
