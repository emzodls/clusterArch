import urllib.request
import xml.etree.ElementTree as etree
import os,timeit
import wget
from ftplib import FTP

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
    acc2gi = dict()
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
                acc = entry.find("*[@Name='Caption']").text.strip()
                descr = entry.find("*[@Name='Title']").text.strip()
                ncbiDict[gi] = (gi,acc,descr)
                acc2gi[acc] = gi
            except Exception as e:
                print(gi,acc,descr,str(e))
                pass
    return ncbiDict,acc2gi

def ncbiFetch(idList,outputFolder,chunkSize = 100,batchMode=True,ncbiDict={},guiSignal=None):
    if batchMode:
        idChunks = [idList[x:x+chunkSize] for x in range(0,len(idList),chunkSize)]
        for dataChunk in idChunks:
            print(dataChunk)
            fetchURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}" \
                       "&retmode=text&rettype=gbwithparts&retmax=251&tool=clusterTools&email=e.de-los-santos@warwick.ac.uk".format(','.join(str(gi) for gi in dataChunk))
            print(fetchURL)
            ## from Eli Korvigo in biostars (https://www.biostars.org/p/66921/) modified because i'm using urllib
            def extract_records(records_handle):
                buffer = []
                for line in records_handle:
                    line = line.decode()
                    if line.startswith("LOCUS") and buffer:
                        # yield accession number and record
                        currentAcc = buffer[0].split()[1]
                        if guiSignal:
                            guiSignal.emit(currentAcc)
                        yield currentAcc, "".join(buffer)
                        buffer = [line]
                    else:
                        buffer.append(line)
                currentAcc = buffer[0].split()[1]
                if guiSignal:
                    guiSignal.emit(currentAcc)
                yield currentAcc, "".join(buffer)

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
### Functions from Marnix Modified
def ftp_connect(guiSignal=None):
  #Connect to NCBI FTP site and change to genbank directory
  try:
    ftp = FTP('ftp.ncbi.nlm.nih.gov')   # connect to host, default port
    ftp.login()               # user anonymous, passwd anonymous@
    ftp.cwd("genbank")
    if guiSignal:
        guiSignal.emit(True)
    return ftp
  except:
    try:
      ftp = FTP('bio-mirror.net')   # connect to host, default port
      ftp.login()               # user anonymous, passwd anonymous@
      ftp.cwd("biomirror")
      ftp.cwd("genbank")
      if guiSignal:
          guiSignal.emit(True)
      return ftp
    except:
      guiSignal.emit(False)

def getGbkDlList(keyword):

    def filterList(line,fileList,keyword):
        fileheader = 'gb' + keyword.lower()
        fileheader.strip()
        lineParse = line.split()
        fileName  = lineParse[-1].strip()
        if fileName.endswith('.gz') and fileName.startswith(fileheader):
            estSize = float(lineParse[4])
            fileList.append((fileName,estSize))
    fileList = []
    ftpHandle = ftp_connect()
    ftpHandle.retrlines('LIST',callback=lambda x: filterList(x,fileList,keyword))

    return fileList

def parseAssemblyReportFile(database,division,keyword=None):
    try:
        assemblyHandle = urllib.request.urlopen('https://ftp.ncbi.nlm.nih.gov/genomes/{}/{}/assembly_summary.txt'.format(database,division.lower()))
        genomeDict = dict()
        for line in assemblyHandle:
            line = line.decode().strip()
            if '#' in line:
                pass
            else:
                lineParse = line.split('\t')
                try:
                    if lineParse[10] == 'latest':
                        species = lineParse[7]
                        if not keyword or keyword.lower() in species.lower():
                            acc = lineParse[0]
                            url = lineParse[19]
                            genomeDict[url] = (acc,species)
                except IndexError:
                    pass
        return genomeDict
    except Exception as e:
        print(e)
        return



if __name__ == '__main__':
    test = parseAssemblyReportFile('bacteria',keyword='burkholderia')
    print(len(test),test)
    # names,sizes = zip(*getGbkDlList('ENV'))
    # print(names,sum(sizes))
    #print(sum(x[1] for x in test))
    # dl = wget.download('ftp://ftp.ncbi.nlm.nih.gov/genbank/{}'.format(test[2][0]),test[0][0])
   # print(dl)