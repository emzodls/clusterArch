

with open('/Volumes/Data/lola_AS3/antimycin/antimycin_hits_filtered.csv','w') as outfile:
    for line in open('/Volumes/Data/lola_AS3/antimycin/antimycin_hits.csv'):
        line = line.strip()
        lineparse = line.split(',')
        txtE = set(x.strip() for x in lineparse[4].split(';'))
        orf07489 = set(x.strip() for x in lineparse[5].split(';'))
        orf07487 = set(x.strip() for x in lineparse[6].split(';'))
        C_Dom = set(x.strip() for x in lineparse[7].split(';'))
        if len(orf07489|orf07487|C_Dom|txtE) >=4 and len(orf07489|orf07487) >= 2 and \
                        len(orf07489|C_Dom) >= 2 and len(orf07487|C_Dom) >=2 and len(txtE|orf07487) >= 2 and len(txtE|orf07489) >= 2\
                and len(C_Dom|txtE) >= 2:
            outfile.write(line + '\n')
