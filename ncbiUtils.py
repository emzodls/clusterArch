import urllib.request
import xml.etree.ElementTree as etree
import os,timeit

def ncbiQuery(keyword,organism,accession,minLength=0,maxLength=10000000000,retmax=1000):
    if not keyword and not organism and not accession:
        print("Nothing")
    else:
        # Accession will override any search
        organism = organism.strip().replace(" ","%20")
        keyword = keyword.strip().replace(" ", "%20")
        if accession:
            esearchURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term={}[accn]%20AND%20{}:{}[SLEN]" \
                         "&retmode=text&retmax={}&tool=clusterTools&email=e.de-los-santos@warwick.ac.uk".format(accession,minLength,maxLength,retmax)
        # otherwise parse the keywords and organisms entries
        else:
            baseURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term="
            if organism and keyword:
                baseURL += "%22{}%22[orgn]%20AND%20{}%20AND%20{}:{}[SLEN]".format(organism,keyword,minLength,maxLength)
            elif organism:
                baseURL += "%22{}%22[orgn]%20AND%20{}:{}[SLEN]".format(organism,minLength,maxLength)
            elif keyword:
                baseURL += "%22{}%22%20AND%20{}:{}[SLEN]".format(keyword,minLength,maxLength)

            esearchURL = baseURL + "&retmode=text&retmax={}&tool=clusterTools&email=e.de-los-santos@warwick.ac.uk".format(retmax)

        try:
            ncbiRequest = urllib.request.urlopen(esearchURL)
            parseRequest = etree.parse(ncbiRequest)
            numHits = int(parseRequest.find('Count').text)
            ncbiIDsList = [idElem.text for idElem in parseRequest.iter('Id')]
        except:
            print("Error Correcting to NCBI Server")
            raise
    return numHits,ncbiIDsList

def ncbiSummary(idList,chunkSize = 250):
    ncbiDict = dict()
    idChunks = [idList[x:x+chunkSize] for x in range(0,len(idList),chunkSize)]
    for dataChunk in idChunks:
        esummaryURL ="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&" \
                     "id={}&retmode=xml&retmax=251&tool=clusterTools&email=e.de-los-santos@warwick.ac.uk".format(','.join(dataChunk))
        try:
            ncbiRequest = urllib.request.urlopen(esummaryURL)
            parseRequest = etree.parse(ncbiRequest)
        except:
            print("Error Correcting to NCBI Server")
            pass
        for entry in parseRequest.findall('DocSum'):
            try:
                gi = int(entry.find("*[@Name='Gi']").text)
                acc = entry.find("*[@Name='Caption']").text
                descr = entry.find("*[@Name='Title']").text
                ncbiDict[gi] = (gi,acc,descr)
            except Exception as e:
                print(gi,acc,descr,str(e))
                pass
    return ncbiDict

def ncbiFetch(idList,outputFolder,chunkSize = 100,batchMode=True,ncbiDict={}):
    if batchMode:
        idChunks = [idList[x:x+chunkSize] for x in range(0,len(idList),chunkSize)]
        for dataChunk in idChunks:
            fetchURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}" \
                       "&retmode=text&rettype=gbwithparts&retmax=251&tool=clusterTools&email=e.de-los-santos@warwick.ac.uk".format(','.join(dataChunk))

            ## from Eli Korvigo in biostars (https://www.biostars.org/p/66921/) modified because i'm using urllib
            def extract_records(records_handle):
                buffer = []
                for line in records_handle:
                    line = line.decode()
                    if line.startswith("LOCUS") and buffer:
                        # yield accession number and record
                        yield buffer[0].split()[1], "".join(buffer)
                        buffer = [line]
                    else:
                        buffer.append(line)
                yield buffer[0].split()[1], "".join(buffer)

            fetchReq = urllib.request.urlopen(fetchURL)

            for accession,record in extract_records(fetchReq):
                with open(os.path.join(outputFolder,'{}.gbk'.format(accession)),'w') as output:
                    output.write(record)
    else:
        for gi in idList:
            fetchURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}" \
                       "&retmode=text&rettype=gbwithparts&retmax=251&tool=clusterTools&email=e.de-los-santos@warwick.ac.uk".format(gi)
            fetchReq = urllib.request.urlopen(fetchURL)
            with open(os.path.join(outputFolder,'{}.gbk'.format(ncbiDict[int(gi)][1])),'wb') as output:
                for line in fetchReq:
                    output.write(line)

def fetchBatch():
    l, idList = ncbiQuery('biosynthesis', 'Streptomyces', '')
    idList = idList[:100]
    ncbiDict = ncbiSummary(idList, chunkSize=250)
    ncbiFetch(idList,'testFetch/')

def fetchSingle():
    l, idList = ncbiQuery('biosynthesis', 'Streptomyces', '')
    idList = idList[:100]
    ncbiDict = ncbiSummary(idList, chunkSize=250)
    ncbiFetch(idList,'testFetch/',batchMode=False,ncbiDict=ncbiDict)

if __name__ == '__main__':
    #l, idList = ncbiQuery('biosynthesis', 'Streptomyces', '')
    #print(idList)
    print(timeit.timeit("fetchSingle()", setup="from __main__ import fetchSingle", number=2))
    print(timeit.timeit("fetchBatch()",setup="from __main__ import fetchBatch",number=2))

