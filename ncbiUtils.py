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
import urllib.request
import xml.etree.ElementTree as etree
import os
from ftplib import FTP
from math import floor,ceil

def ncbiQuery(keyword,organism,accession,minLength=0,maxLength=10000000000,retmax=1000,db='nuccore'):
    if not keyword and not organism and not accession:
        print("Nothing")
    else:
        # Accession will override any search
        organism = organism.strip().replace(" ","%20")
        keyword = keyword.strip().replace(" ", "%20")
        if accession:
            esearchURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={}&term={}[accn]%20AND%20{}:{}[SLEN]" \
                         "&retmode=text&retmax={}&tool=clusterTools&email=e.de-los-santos@warwick.ac.uk".format(db,accession,minLength,maxLength,retmax)
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
            print("Error Connecting to NCBI Server")
            raise
    return numHits,ncbiIDsList

def ncbiSummary(idList,db='nuccore',chunkSize = 250):
    ncbiDict = dict()
    acc2gi = dict()
    idChunks = [idList[x:x+chunkSize] for x in range(0,len(idList),chunkSize)]
    for dataChunk in idChunks:
        esummaryURL ="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db={}&" \
                     "id={}&retmode=xml&retmax=251&tool=clusterTools&email=e.de-los-santos@warwick.ac.uk".format(db,','.join(dataChunk))
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

def ncbiFetch(idList,outputFolder,chunkSize = 100,batchMode=True,ncbiDict={},guiSignal=None,db='nuccore'):
    if batchMode:
        idChunks = [idList[x:x+chunkSize] for x in range(0,len(idList),chunkSize)]
        for dataChunk in idChunks:
            fetchURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={}&id={}" \
                       "&retmode=text&rettype=gbwithparts&retmax=251&tool=clusterTools&email=e.de-los-santos@warwick.ac.uk".format(db,','.join(str(gi) for gi in dataChunk))
            ## from Eli Korvigo in biostars (https://www.biostars.org/p/66921/) modified because i'm using urllib
            def extract_records(records_handle):
                buffer = []
                for line in records_handle:
                    line = line.decode("utf-8", "ignore")
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

def accToFasta(accList,db,outputfolder):
    acc2gi = dict()
    for acc in accList:
        esearchURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={}&term={}[accn]" \
                     "&retmode=text&retmax=10&tool=clusterTools&email=e.de-los-santos@warwick.ac.uk".format(db,acc)
        try:
            ncbiRequest = urllib.request.urlopen(esearchURL)
            parseRequest = etree.parse(ncbiRequest)
            ncbiID = [idElem.text for idElem in parseRequest.iter('Id')][0]
            acc2gi[acc] = ncbiID
        except:
            print("Error Connecting to NCBI Server for acc {}".format(acc))
            pass
    idList = list(acc2gi.values())
    idChunks = [idList[x:x + 100] for x in range(0, len(acc2gi), 100)]
    for dataChunk in idChunks:
        fetchURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={}&id={}" \
                   "&retmode=text&rettype=fasta&retmax=251&tool=clusterTools&email=e.de-los-santos@warwick.ac.uk".format(
            db, ','.join(str(gi) for gi in dataChunk))
        ## from Eli Korvigo in biostars (https://www.biostars.org/p/66921/) modified because i'm using urllib
        def extract_records(records_handle):
            buffer = []
            for line in records_handle:
                line = line.decode("utf-8", "ignore")
                if line.startswith(">") and buffer:
                    # yield accession number and record
                    currentAcc = buffer[0].split()[0][1:].split('.')[0]
                    yield currentAcc, "".join(buffer)
                    buffer = [line]
                else:
                    buffer.append(line)
            currentAcc = buffer[0].split()[0][1:].split('.')[0]
            yield currentAcc, "".join(buffer)

        fetchReq = urllib.request.urlopen(fetchURL)

        for accession, record in extract_records(fetchReq):
            with open(os.path.join(outputfolder, '{}.fasta'.format(accession)), 'w') as output:
                output.write(record)

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

def parseAssemblyReportFile(database,division,keyword=None,categoryFilters=None):
    try:
        assemblyHandle = urllib.request.urlopen('https://ftp.ncbi.nlm.nih.gov/genomes/{}/{}/assembly_summary.txt'.format(database,division.lower()))
        genomeDict = dict()
        for line in assemblyHandle:
            line = line.decode("utf-8", "ignore").strip()
            if '#' in line:
                pass
            else:
                lineParse = line.split('\t')
                try:
                    category = lineParse[4]
                    if lineParse[10] == 'latest':
                        if not categoryFilters or any(catFilter in category for catFilter in categoryFilters):
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

def fetchGbksWithAcc(clusterList,window,outputFolder,guiSignal=None):
    # cluster list format (acc,(start,end))
    dlList = set()
    acc2ncbi = dict()
    for acc,coordinates in clusterList:
        esearchURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={}&term={}[accn]" \
                     "&retmode=text&retmax=10&tool=clusterTools&email=e.de-los-santos@warwick.ac.uk".format('nuccore',acc)
        try:
            if acc not in acc2ncbi.keys():
                ncbiRequest = urllib.request.urlopen(esearchURL)
                parseRequest = etree.parse(ncbiRequest)
                ncbiID = [idElem.text for idElem in parseRequest.iter('Id')][0]
                acc2ncbi[acc] = ncbiID
            else:
                ncbiID = acc2ncbi[acc]
            clusterMidpoint = sum(coordinates)/2
            window_start = floor(max(1,clusterMidpoint- window/2))
            window_end = ceil(clusterMidpoint + window/2)
            fetchURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}" \
                   "&retmode=text&rettype=gbwithparts&strand=1&seq_start={}&seq_stop={}" \
                   "&retmax=251&tool=clusterTools&email=e.de-los-santos@warwick.ac.uk".format(ncbiID,window_start,window_end)
            fetchReq = urllib.request.urlopen(fetchURL)
            buffer = []
            for line in fetchReq:
                line = line.decode("utf-8", "ignore")
                buffer.append(line)
                if line.startswith('DEFINITION'):
                    ## try to parse the species name from the definition
                    summary = line.split('DEFINITION')[1].strip().split(',')[0].strip()

            with open(os.path.join(outputFolder,'{}-{}-{}.gbk'.format(acc,window_start,window_end)),'w') as gbkFile:
                gbkFile.write(''.join(buffer))
            if guiSignal:
                guiSignal.emit((acc,ncbiID,summary))
            dlList.add((acc,ncbiID,summary))
        except:
            if guiSignal:
                guiSignal.emit((acc,None,None))
            dlList.add((acc,None,None))
    dlSuccess = set(entry for entry in dlList if entry[2])
    dlFail = set(entry for entry in dlList if not entry[2])
    with open(os.path.join(outputFolder, 'downloadSummary.txt'), 'w') as summaryFile:
        summaryFile.write('Failed to Resolve:\n')
        for entry,blank,blank2 in dlFail:
            summaryFile.write('{}\n'.format(entry))
        summaryFile.write('Downloaded:\n')
        for fileNameAcc,gi,summary in dlSuccess:
            summaryFile.write('{},{},{}\n'.format(fileNameAcc, gi,summary))
    return dlList

def protToFasta(proteinList,outfile):
    protID2Acc = dict()
    acc2ProtID = dict()
    acc2gi = dict()
    failedToFetch = list()
    success = list()
    for protein in proteinList:
        acc = protein.name.split('.')[0]
        proteinID = protein.fastaID
        esearchURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term={}[accn]" \
                     "&retmode=text&retmax=10&tool=clusterTools&email=e.de-los-santos@warwick.ac.uk".format(acc)
        try:
            ncbiRequest = urllib.request.urlopen(esearchURL)
            parseRequest = etree.parse(ncbiRequest)
            ncbiID = [idElem.text for idElem in parseRequest.iter('Id')][0]
            acc2gi[acc] = ncbiID
            protID2Acc[proteinID] = acc
            acc2ProtID[acc] = proteinID
            success.append(acc)
        except:
            print("Error Connecting to NCBI Server for acc {}".format(acc))
            failedToFetch.append(acc)
            pass
    idList = list(acc2gi.values())
    idChunks = [idList[x:x + 100] for x in range(0, len(acc2gi), 100)]
    if os.path.isfile(outfile):
        os.remove(outfile)

    for dataChunk in idChunks:
        fetchURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={}" \
                   "&retmode=text&rettype=fasta&retmax=251&tool=clusterTools&email=e.de-los-santos@warwick.ac.uk".format(','.join(str(gi) for gi in dataChunk))
        def extract_records(records_handle):
            buffer = []
            for line in records_handle:
                line = line.decode("utf-8", "ignore")
                if line:
                    if line.startswith(">"):
                        accID = line.split()[0][1:].split('.')[0]

                        buffer.append('>{}\n'.format(acc2ProtID[accID]))
                    else:
                        buffer.append(line)
            return "".join(buffer)
        fetchReq = urllib.request.urlopen(fetchURL)
        with open(outfile, 'a') as output:
            output.write(extract_records(fetchReq))

if __name__ == '__main__':
    clusterList = []
    #
    # for line in open('/Volumes/Data/lola_AS3/antimycin_ozmN/hits_genomesGB.csv'):
    #     if '##' in line:
    #         pass
    #     else:
    #         lineParse = line.split(',')
    #         clusterList.append((lineParse[0],(int(lineParse[1]),int(lineParse[2]))))
    # print(clusterList)
    clusterList = [('EF552687', (20505, 87657)), ('CP006871', (6578021, 6611567)), ('CP006871', (1405296, 1502258)), ('CP006871', (9179933, 9267062)), ('CP007574', (6209940, 6264994)), ('CP000249', (2320878, 2324816))]

    test = fetchGbksWithAcc(clusterList[:10], 100000, '/Volumes/Data/lola_AS3/antimycin_ozmN/testDL', guiSignal=None)
    print(test)
    # names,sizes = zip(*getGbkDlList('ENV'))
    # print(names,sum(sizes))
    #print(sum(x[1] for x in test))
    # dl = wget.download('ftp://ftp.ncbi.nlm.nih.gov/genbank/{}'.format(test[2][0]),test[0][0])
   # print(dl)