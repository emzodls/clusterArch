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

import os,csv
from PyQt5.QtCore import (QRegularExpression,
        Qt,QCoreApplication,QThread,pyqtSignal,pyqtSlot,QObject)
from PyQt5 import QtGui
from PyQt5.QtWidgets import (QApplication, QFileDialog, QMainWindow,
        QMessageBox,QCheckBox,QTableWidgetItem,QWidget,QListWidgetItem)
from utils import parseSeqFile,parseHMMfile,generateInputFasta,generateHMMdb,MakeBlastDB,runBLAST,\
    runHmmsearch,processSearchList,proccessGbks,processSearchListOptionalHits,humanbytes,processGbkDivFile
from collections import Counter
from glob import iglob,glob
from itertools import chain
from ncbiUtils import ncbiQuery,ncbiSummary,ncbiFetch,getGbkDlList

import mainGuiNCBI,resultsWindow,status,parameters,sys,createDbStatus,addGene,downloadNcbiFilesWin,wget

class ncbiQueryWorker(QObject):
    start = pyqtSignal()
    connectionErr = pyqtSignal()
    result = pyqtSignal(int,dict,dict)
    def __init__(self,keyword,organism,accession,minLength,maxLength,retmax):
        self.keyword = keyword
        self.organism = organism
        self.accession = accession
        self.minLength = minLength
        self.maxLength = maxLength
        self.retmax = retmax
        super(ncbiQueryWorker, self).__init__()

    @pyqtSlot()
    @pyqtSlot(int,dict,dict)
    def run(self):
        try:
            numHits, ncbiIDsList = ncbiQuery(self.keyword,self.organism,self.accession,
                                             self.minLength,self.maxLength,self.retmax)
            ncbiDict,acc2gi = ncbiSummary(ncbiIDsList, chunkSize=250)
            self.result.emit(numHits, ncbiDict, acc2gi)
        except:
            self.connectionErr.emit()

class downloadNcbiFilesWorker(QObject):
    start = pyqtSignal()
    currentFile = pyqtSignal(str)
    finished = pyqtSignal()
    def __init__(self,idList,outputFolder):
        self.idList = idList
        self.outputFolder = outputFolder
        super(downloadNcbiFilesWorker,self).__init__()
    @pyqtSlot()
    @pyqtSlot(str)
    def run(self):
        try:
            ncbiFetch(self.idList, self.outputFolder,guiSignal=self.currentFile)
            self.currentFile.emit("Successfully Downloaded All Files.")
            self.finished.emit()
        except Exception as e:
            self.currentFile.emit("There was a Problem with the Download. \n Error: %s" % str(e))
            self.finished.emit()

class downloadGbDivFilesWorker(QObject):
    start = pyqtSignal()
    currentFile = pyqtSignal(str)
    finished = pyqtSignal(list)
    def __init__(self,filelist,outputFolder):
        self.filelist = filelist
        self.outputFolder = outputFolder
        self.failedToFetch = []
        super(downloadGbDivFilesWorker,self).__init__()
    @pyqtSlot()
    @pyqtSlot(str)
    @pyqtSlot(list)
    def run(self):
        for gbkDivFile in self.filelist:
            try:
                self.currentFile.emit(gbkDivFile)
                wget.download('ftp://ftp.ncbi.nlm.nih.gov/genbank/{}'.format(gbkDivFile),
                              os.path.join(self.outputFolder,'{}.clusterToolsDB'.format(gbkDivFile)), False)
            # if it fails try mirror
            except:
                try:
                    self.currentFile.emit(gbkDivFile)
                    wget.download('ftp://bio-mirror.net/biomirror/genbank/{}'.format(gbkDivFile),
                                  os.path.join(self.outputFolder, '{}.clusterToolsDB'.format(gbkDivFile)), False)
                except:
                    self.failedToFetch.append(gbkDivFile)
                    pass
        self.finished.emit(self.failedToFetch)

class processGbDivFilesWorker(QObject):
    start = pyqtSignal()
    currentFile = pyqtSignal(str)
    currentSpecies = pyqtSignal(str)
    finished = pyqtSignal(list)

    def __init__(self,filelist,targetfolder):
        self.filelist = filelist
        self.targetfolder = targetfolder
        self.failedToProcess = []
        super(processGbDivFilesWorker,self).__init__()
    @pyqtSlot()
    @pyqtSlot(str)
    @pyqtSlot(list)
    def run(self):
        for gbkDivFile in self.filelist:
            try:
                print(gbkDivFile)
                path,gbkFileName = os.path.split(gbkDivFile)
                self.currentFile.emit(gbkFileName)
                database = os.path.join(self.targetfolder,'{}.clusterTools.fasta'.format(gbkDivFile))
                processGbkDivFile(gbkDivFile,database,self.currentSpecies)
                os.remove(gbkDivFile)
            except Exception as e:
                print(e)
                self.failedToProcess.append(gbkFileName)
                pass
        clusterToolsFastaFiles = glob(os.path.join(self.targetfolder,'*.clusterTools.fasta'))
        with open(os.path.join(self.targetfolder,'gbDivDB.clusterTools.fasta'), 'a') as outfile:
            for clusterToolsFastaFile in clusterToolsFastaFiles:
                with open(clusterToolsFastaFile) as infile:
                    for line in infile:
                        outfile.write(line)
        self.finished.emit(self.failedToProcess)

class createDbWorker(QObject):
    start = pyqtSignal()
    currentGbk = pyqtSignal(str)
    finished = pyqtSignal()
    def __init__(self,taskList,dbPath):
        self.taskList = taskList
        self.dbPath = dbPath
        super(createDbWorker,self).__init__()
    @pyqtSlot()
    def run(self):
        try:
            proccessGbks(self.taskList,self.dbPath,self.currentGbk)
            self.currentGbk.emit("Database Creation Successful!")
            self.finished.emit()
        except:
            self.finished.emit()

class runChecksWorker(QObject):
    checkError = pyqtSignal(str)
    start = pyqtSignal(str)
    finished = pyqtSignal(bool)
    def __init__(self,outputDir,pathToDatabase):
        self.outputDir = outputDir
        self.pathToDatabase = pathToDatabase
        super(runChecksWorker,self).__init__()

    @pyqtSlot(str)
    def run(self):
        # First check that you can read and write to the required directories
        if not self.outputDir:
            self.checkError.emit('noOutputDir')
            self.finished.emit(False)
        else:
            if not os.access(self.outputDir, os.W_OK):
                self.checkError.emit('noWriteOutputDir')
                self.finished.emit(False)
            else:
                if not self.pathToDatabase:
                    self.checkError.emit('noDBpath')
                    self.finished.emit(False)
                else:
                    if not os.access(self.pathToDatabase, os.R_OK):
                        self.checkError.emit('noReadDB')
                        self.finished.emit(False)
                    else:
                    # Now check if the specified database is in the required format
                        with open(self.pathToDatabase) as dbHandle:
                            hasValidEntriesFlag = False
                            someInvalidEntriesFlag = False
                            for line in dbHandle:
                                if '>' in line:
                                    try:
                                        line_parse = line.split('|')
                                        coordinates = [int(x) for x in line_parse[1].split('-')]
                                        direction = line_parse[2]
                                        queryIntID = line_parse[3]
                                        queryProtID = line_parse[4]
                                        queryIntIdx = int(queryIntID.split('_')[-1])
                                        hasValidEntriesFlag = True
                                    except (IndexError,ValueError):
                                        someInvalidEntriesFlag = True
                                        pass
                                else:
                                    pass
                        if not hasValidEntriesFlag:
                            self.checkError.emit('invalidDB')
                            self.finished.emit(False)
                        else:
                            if someInvalidEntriesFlag:
                                self.checkError.emit('entriesWarn')
                                self.finished.emit(True)
                            else:
                                self.checkError.emit('good')
                                self.finished.emit(True)

class runBlastWorker(QObject):
    inputGenerated = pyqtSignal(bool)
    dbCreated = pyqtSignal(bool)
    makeDB = pyqtSignal(bool)
    doneBLAST = pyqtSignal(bool)
    start = pyqtSignal()

    def __init__(self,forBLAST,outputDir,pathToDatabase,makeblastdbExec,blastpExec,blastEval):
        self.forBLAST = forBLAST
        self.outputDir = outputDir
        self.pathToDatabase = pathToDatabase
        self.makeblastdbExec = makeblastdbExec
        self.blastpExec = blastpExec

        self.blastEval = blastEval
        super(runBlastWorker,self).__init__()
        # QThread.__init__(self)
    @pyqtSlot()
    @pyqtSlot(bool)
    def run(self):
        ### Generate Input Fasta
        generateInputFasta(self.forBLAST,self.outputDir)
        self.inputGenerated.emit(True)
        baseDir, ext = os.path.splitext(self.pathToDatabase)
        phrCheck = iglob(baseDir + '*phr')
        pinCheck = iglob(baseDir + '*pin')
        psqCheck = iglob(baseDir + '*psq')

        outdirPhrCheck = iglob(os.path.join(self.outputDir,'*phr'))
        outdirPinCheck = iglob(os.path.join(self.outputDir, '*pin'))
        outdirPsqCheck = iglob(os.path.join(self.outputDir, '*psq'))

        ### Check for Blastp formatted Database
        path, outputDBname = os.path.split(self.pathToDatabase)
        ## BLAST formatted database is in database path
        if any(True for _ in phrCheck) and any(True for _ in pinCheck) and any(True for _ in psqCheck):
            self.makeDB.emit(False)
            dbInOutputFolder = False
        ## BLAST formatted database is in output directory
        elif any(True for _ in outdirPhrCheck) and any(True for _ in outdirPinCheck) and any(True for _ in outdirPsqCheck):
            self.makeDB.emit(False)
            dbInOutputFolder = True
        else:
            self.makeDB.emit(True)
            dbInOutputFolder = True

            out, err, retcode = MakeBlastDB(self.makeblastdbExec, self.pathToDatabase, self.outputDir,outputDBname)
            if retcode != 0:
                self.dbCreated.emit(False)
                self.terminate()
                return
            else:
                self.dbCreated.emit(True)
        ### Run Blastp
        inputFastas = os.path.join(self.outputDir, 'gene_queries.fa')

        if dbInOutputFolder:
            out, err, retcode = runBLAST(self.blastpExec, inputFastas, self.outputDir,
                                         os.path.join(self.outputDir, outputDBname),
                                         eValue=str(self.blastEval))
        else:
            out, err, retcode = runBLAST(self.blastpExec, inputFastas, self.outputDir,
                                         self.pathToDatabase, eValue=str(self.blastEval))
        if retcode != 0:
            self.doneBLAST.emit(False)
            return
        else:
            self.doneBLAST.emit(True)
            return

class runHmmerWorker(QObject):
    start = pyqtSignal()
    hmmDBbuilt = pyqtSignal(bool,list)
    hmmSearchComplete = pyqtSignal(bool)

    def __init__(self,hmmFetchExec,hmmSearchExec,hmmDict,hmmSet,outputDir,pathToDatabase,hmmEval):
        self.hmmFetchExec = hmmFetchExec
        self.hmmSearchExec = hmmSearchExec
        self.hmmDict = hmmDict
        self.hmmSet = hmmSet
        self.outputDir = outputDir
        self.pathToDatabase = pathToDatabase
        self.hmmEval = hmmEval
        super(runHmmerWorker, self).__init__()
    @pyqtSlot()
    @pyqtSlot(bool)
    @pyqtSlot(bool,list)
    def run(self):
        errFlag, failedToFetch = generateHMMdb(self.hmmFetchExec, self.hmmDict, self.hmmSet, self.outputDir)
        if errFlag:
            self.hmmDBbuilt.emit(False,failedToFetch)
            return
        else:
            self.hmmDBbuilt.emit(True,[])
            hmmDBase = os.path.join(self.outputDir, 'hmmDB.hmm')
            out, err, retcode = runHmmsearch(self.hmmSearchExec, hmmDBase, self.outputDir,
                                             self.pathToDatabase, eValue=str(self.hmmEval))
            if retcode != 0:
                self.hmmSearchComplete.emit(False)
                return
            else:
                self.hmmSearchComplete.emit(True)
                return

class runProcessSearchList(QObject):
    start = pyqtSignal()
    result = pyqtSignal(dict)
    def __init__(self,blastList,hmmList,blastOutFile,hmmOutFile,hmmScore,hmmDomLen,windowSize):
        self.blastList = blastList
        self.hmmList = hmmList
        self.blastOutFile = blastOutFile
        self.hmmOutFile = hmmOutFile

        self.hmmScore = hmmScore
        self.hmmDomLen = hmmDomLen
        self.windowSize = windowSize

        super(runProcessSearchList, self).__init__()
    @pyqtSlot()
    @pyqtSlot(dict)
    def run(self):
        filteredClusters = processSearchList(self.blastList, self.hmmList, self.blastOutFile, self.hmmOutFile,
                                             self.hmmScore,self.hmmDomLen,self.windowSize)
        self.result.emit(filteredClusters)

class runProcessSearchListOptionalHits(QObject):
    start = pyqtSignal()
    result = pyqtSignal(dict)
    def __init__(self,requiredBlastList,requiredHmmList,blastOutFile,hmmOutFile,hmmScore,hmmDomLen,windowSize,
                 totalHitsRequired, additionalBlastList=[], additionalHmmList=[]):
        self.requiredBlastList = requiredBlastList
        self.requiredHmmList = requiredHmmList

        self.additionalBlastList = additionalBlastList
        self.additionalHmmList = additionalHmmList

        self.blastOutFile = blastOutFile
        self.hmmOutFile = hmmOutFile

        self.hmmScore = hmmScore
        self.hmmDomLen = hmmDomLen
        self.windowSize = windowSize

        self.totalHitsRequired = totalHitsRequired


        super(runProcessSearchListOptionalHits, self).__init__()
    @pyqtSlot()
    @pyqtSlot(dict)
    def run(self):
        filteredClusters = processSearchListOptionalHits(self.requiredBlastList, self.requiredHmmList,
                                                         self.blastOutFile, self.hmmOutFile,
                                                         self.hmmScore,self.hmmDomLen,self.windowSize,
                                                         self.totalHitsRequired,self.additionalBlastList,self.additionalHmmList)
        self.result.emit(filteredClusters)

class downloadNcbiFilesWindow(QWidget,downloadNcbiFilesWin.Ui_downloadNcbiWin):
    def __init__(self):
        super(self.__class__,self).__init__()
        self.setupUi(self)

class addGeneWindow(QWidget,addGene.Ui_addGeneWindow):
    def __init__(self):
        super(self.__class__,self).__init__()
        self.setupUi(self)
        regex = QRegularExpression('^[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]+$')
        self.AminoAcidValidator = QtGui.QRegularExpressionValidator(regex)

class statusWindow(QWidget,status.Ui_runSearch):
    def __init__(self):
        super(self.__class__,self).__init__()
        self.setupUi(self)

class createDbStatusWindow(QWidget,createDbStatus.Ui_createDbWin):
    def __init__(self):
        super(self.__class__,self).__init__()
        self.setupUi(self)

class paramsWindow(QWidget,parameters.Ui_parametersWin):
    def __init__(self):
        super(self.__class__,self).__init__()
        self.setupUi(self)
        regex = QRegularExpression('(\d*\.?\d+(E-|E+|E|e-|e\+|e|\d+)\d+|\d*\.?\d+)')
        SciNoteValidator = QtGui.QRegularExpressionValidator(regex)
        posIntValidator = QtGui.QIntValidator()
        posIntValidator.setBottom(1)

        self.blastEval.setValidator(SciNoteValidator)
        self.hmmEval.setValidator(SciNoteValidator)

        self.hmmScore.setValidator(posIntValidator)
        self.hmmDomLen.setValidator(posIntValidator)
        self.windowSize.setValidator(posIntValidator)

class resultsWindow(QWidget,resultsWindow.Ui_Results):
    def __init__(self):
        super(self.__class__,self).__init__()
        self.setupUi(self)

class mainApp(QMainWindow, mainGuiNCBI.Ui_clusterArch):
    def __init__(self,makeblastdbExec,blastExec,hmmFetchExec,hmmSearchExec,verbose=False):

        super(self.__class__, self).__init__()
        self.setupUi(self)

        ##### Set up NCBI File Interface
        self.ncbiDict = dict()
        self.acc2gi = dict()
        self.ncbiDLdict = dict()

        posIntValidator = QtGui.QIntValidator(0,1000000000)

        maxEntriesValidator = QtGui.QIntValidator(1,1000)

        self.minLength.setValidator(posIntValidator)
        self.maxLength.setValidator(posIntValidator)
        self.numEntries.setValidator(maxEntriesValidator)

        self.minLength.setText('0')
        self.maxLength.setText('1000000000')
        self.numEntries.setText('1000')

        self.searchNcbiBtn.clicked.connect(self.ncbiQuery)
        self.addNcbiFilesBtn.clicked.connect(self.addNcbiFiles)
        self.removeNcbiFilesBtn.clicked.connect(self.removeNcbiFiles)
        self.selectNcbiFileDirBtn.clicked.connect(self.setNcbiFileDlDirectory)
        self.dlNcbiFilesBtn.clicked.connect(self.downloadNcbiFiles)

        self.searchFilter.editingFinished.connect(self.selectFilteredItems)
        self.clearSearchSelectionBtn.clicked.connect(self.ncbiFileSearchResults.clearSelection)
        ####  Set up Genbank DB Interface ####
        self.gbDivFilesToDlList = []
        self.checkBoxList = [self.bctChk,self.envChk,self.gssChk,self.htgChk,self.invChk,self.mamChk,self.patChk,self.plnChk,
                             self.priChk,self.rodChk,self.stsChk,self.synChk,self.tsaChk,self.vrlChk,self.vrtChk]
        self.gbDivList = ['BCT','ENV','GSS','HTG','INV','MAM','PAT','PLN','PRI','ROD','STS','SYN','TSA','VRL','VRT']

        self.checkRequestBtn.clicked.connect(self.checkGbDivDlRequest)
        self.selectGbDivDirBtn.clicked.connect(self.selectGbDivDir)
        self.downloadGbDivFileBtn.clicked.connect(self.downloadGenbankDivFiles)
        # when state of checkbox changes, force user to reestimate the download
        for checkbox in self.checkBoxList:
            checkbox.stateChanged.connect(self.redoGbDlEstimate)
        ### Set Up Parameters and Shared Data Structures
        self.nameToParamDict = {'blastEval': 1e-5, 'hmmEval': 1e-5, 'hmmScore': 30,
                                'hmmDomLen': 15, 'windowSize': 50000}

        self.geneDict = dict()
        self.hmmDict = dict()
        self.forBLAST = dict()
        self.forHmmer = dict()
        self.pathToDatabase = ''
        self.outputDir  = ''
        self.runnerThread = None

        self.makeblastdbExec = makeblastdbExec
        self.blastExec = blastExec
        self.hmmFetchExec = hmmFetchExec
        self.hmmSearchExec = hmmSearchExec
        self.verbose = verbose

        self.currentGeneSelected = ''

        ### Set Up Workflow Variables

        self.checkSuccessful = False
        self.blastSuccessful = False
        self.hmmerSuccessful = False

        ### Set Up Button Connections
        self.genbankDirBtn.clicked.connect(self.loadGbkDir)
        self.addGenbankBtn.clicked.connect(self.addGbkFiles)
        self.deleteGenbankBtn.clicked.connect(self.removeGbkFiles)
        self.chooseDBoutputDirBtn.clicked.connect(self.chooseDBoutDir)
        self.createDbBtn.clicked.connect(self.createDB)

        self.addGenefilePathBtn.clicked.connect(self.openGene)
        self.addHMMfilePathBtn.clicked.connect(self.openHmm)
        self.addGeneSeqBtn.clicked.connect(self.showAddGeneWin)
        self.addGenefileBtn.clicked.connect(self.loadGeneFile)
        self.addHMMfilebtn.clicked.connect(self.loadHMMFile)

        self.addGeneBtn.clicked.connect(self.loadGene)
        self.addHMMbtn.clicked.connect(self.loadHMM)

        self.removeSearchTermBtn.clicked.connect(self.removeSearchTerm)

        self.clearGeneBtn.clicked.connect(self.clearGene)
        self.clearHMMbtn.clicked.connect(self.clearHMM)

        self.runSearch.clicked.connect(self.runChecks)

        self.closeBtn.clicked.connect(QCoreApplication.instance().quit)

        self.dataBaseSelector.currentIndexChanged.connect(self.databaseFunction)
        self.dataBaseSelector.activated.connect(self.databaseFunction)
        self.outputDirectorySelector.currentIndexChanged.connect(self.outdirFunction)
        self.outputDirectorySelector.activated.connect(self.outdirFunction)

        #### Rename Genes in Gene List on Double Click
        self.geneList.currentRowChanged.connect(self.updateSelectedGene)
        self.geneList.itemDelegate().commitData.connect(self.updateGeneName)
        ### For Required Genes
        self.searchList.itemChanged.connect(self.updateSpinBox)
        self.totalHitsRequired = 0
        self.hitsNeededSpinBox.valueChanged.connect(self.updateTotalHitsRequired)
        ### Initialize settings for parameters window
        self.paramsWin = paramsWindow()

        self.paramsWin.blastEval.textChanged.connect(self.activateApply)
        self.paramsWin.hmmEval.textChanged.connect(self.activateApply)
        self.paramsWin.hmmDomLen.textChanged.connect(self.activateApply)
        self.paramsWin.hmmScore.textChanged.connect(self.activateApply)
        self.paramsWin.windowSize.textChanged.connect(self.activateApply)
        self.paramsWin.applyBtn.clicked.connect(self.updateParams)
        self.setParamsBtn.clicked.connect(self.showParamWindow)
    ##### NCBI File Query Functions #########
    def ncbiQuery(self):
        keyword = self.ncbiKeyword.text()
        organism = self.ncbiOrganism.text()
        accession = self.ncbiAccession.text()
        minLength = int(self.minLength.text())
        maxLength = int(self.maxLength.text())
        retmax = self.numEntries.text()

        if not keyword and not organism and not accession:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Please Specify at Least One Search Parameter")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
        else:
            if minLength >= maxLength:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setText("Min Length > Max Length, Using defaults.")
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()

            self.runnerThread = QThread()
            self.runnerThread.start()
            self.ncbiQueryWorker = ncbiQueryWorker(keyword,organism,accession,minLength,maxLength,retmax)
            self.ncbiQueryWorker.connectionErr.connect(self.connectionError)
            self.ncbiQueryWorker.moveToThread(self.runnerThread)
            self.ncbiQueryWorker.start.connect(self.ncbiQueryWorker.run)
            self.ncbiQueryWorker.start.emit()
            ### Reset Search Results
            self.ncbiFileSearchResults.clear()
            self.ncbiDict = dict()
            self.acc2gi = {k:v for k,v in self.acc2gi.items() if v in self.ncbiDLdict.keys()}
            self.ncbiFileSearchResults.addItem('Querying NCBI Database...')
            self.ncbiQueryWorker.result.connect(self.updateSearchResults)

    def connectionError(self):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Problems Connecting to NCBI Database")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec()

    def updateSearchResults(self,numHits,ncbiDict,acc2gi):
        self.ncbiFileSearchResults.takeItem(0)
        if numHits == 0:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("No Results Found")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
        elif numHits >= 1000:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("More than 1000 Hits, showing only Top %s results" % self.numEntries.text())
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()

        self.totalSearchResults.setText(str(numHits))
        self.ncbiDict = ncbiDict
        self.acc2gi = {**self.acc2gi,**acc2gi}
        for gi,entry in ncbiDict.items():
            self.ncbiFileSearchResults.addItem("{}: {}".format(entry[1],entry[2]))
        if self.runnerThread:
            self.runnerThread.terminate()

    def selectFilteredItems(self):
        selectedItems = self.ncbiFileSearchResults.findItems(self.searchFilter.text(),Qt.MatchContains)
        for item in selectedItems:
            item.setSelected(True)

    def addNcbiFiles(self):
        for ncbiFile in self.ncbiFileSearchResults.selectedItems():
            acc = ncbiFile.text().split(':')[0].strip()
            gi = self.acc2gi[acc]
            if gi not in self.ncbiDLdict.keys():
                self.ncbiDLdict[gi] = self.ncbiDict[gi]
                self.ncbiFileDLlist.addItem(ncbiFile.text())
        self.ncbiDLtotalFilesCt.setText(str(len(self.ncbiDLdict)))

    def removeNcbiFiles(self):
        for ncbiFile in self.ncbiFileDLlist.selectedItems():
            acc = ncbiFile.text().split(':')[0].strip()
            gi  = self.acc2gi[acc]
            del self.ncbiDLdict[gi]
            self.ncbiFileDLlist.takeItem(self.ncbiFileDLlist.row(ncbiFile))
        self.ncbiDLtotalFilesCt.setText(str(len(self.ncbiDLdict)))

    def setNcbiFileDlDirectory(self):
        dirName = QFileDialog.getExistingDirectory()
        if dirName:
            if os.access(dirName, os.W_OK):
               self.ncbiFileDlDir.setText(dirName)
            else:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText("Can't Write to that Folder. Please Specify Another Directory")
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()

    def downloadNcbiFiles(self):
        ncbiFileDlDir = self.ncbiFileDlDir.text()
        taskList = set(self.ncbiDLdict.keys())
        if not ncbiFileDlDir:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("No Output Directory Specified")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
        elif len(taskList) == 0:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Nothing in Download List")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
        else:
            alreadyDownloadedAcc = set(os.path.split(x)[1].split('.gbk')[0] for x in glob(os.path.join(ncbiFileDlDir,'*.gbk')))
            alreadyDownloadedGI = set(self.acc2gi[acc] for acc in alreadyDownloadedAcc if acc in self.acc2gi.keys())
            alreadyDlGIDict = {gi:self.ncbiDLdict[gi][1] for gi in alreadyDownloadedGI}
            if self.verbose:
                print(taskList,alreadyDownloadedGI,taskList&alreadyDownloadedGI)
            if len(taskList&alreadyDownloadedGI) > 1:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Information)
                msg.setText("Already Genbanks for the following Accession IDs {}. "
                            "Skipping These.".format(', '.join(alreadyDlGIDict[gi] for gi in taskList&alreadyDownloadedGI)))
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()
            if not self.runnerThread:
                self.runnerThread = QThread()
            self.runnerThread.start()
            if self.verbose:
                print(taskList-alreadyDownloadedGI)
            self.downloadNcbiFilesWorker = downloadNcbiFilesWorker(list(taskList-alreadyDownloadedGI),ncbiFileDlDir)
            self.downloadNcbiFilesWin = downloadNcbiFilesWindow()
            self.downloadNcbiFilesWin.doneBtn.clicked.connect(self.downloadNcbiFilesWin.close)
            self.downloadNcbiFilesWin.cancelBtn.clicked.connect(lambda: self.abortSearch(self.downloadNcbiFilesWin,self.runnerThread))
            self.downloadNcbiFilesWin.progressBar.setMaximum(int(len(taskList-alreadyDownloadedGI)))
            self.downloadNcbiFilesWin.progressBar.setValue(0)
            self.downloadNcbiFilesWorker.moveToThread(self.runnerThread)
            self.downloadNcbiFilesWorker.start.connect(self.downloadNcbiFilesWorker.run)
            self.downloadNcbiFilesWorker.currentFile.connect(self.updateCurrentFile)
            self.downloadNcbiFilesWorker.finished.connect(self.ncbiFileDlFinished)
            self.downloadNcbiFilesWorker.start.emit()

            self.downloadNcbiFilesWin.show()

    def updateCurrentFile(self,currentFile):
        self.downloadNcbiFilesWin.currentTask.setText(currentFile)
        self.downloadNcbiFilesWin.progressBar.setValue(self.downloadNcbiFilesWin.progressBar.value()+1)
        self.downloadNcbiFilesWin.progressLabel.setText('{}/{}'.format(self.downloadNcbiFilesWin.progressBar.value(),
                                                                      self.downloadNcbiFilesWin.progressBar.maximum()))
    def ncbiFileDlFinished(self):
        self.ncbiFileDLlist.clear()
        self.ncbiDLdict = dict()
        self.downloadNcbiFilesWin.doneBtn.setEnabled(True)
        if self.runnerThread:
            self.runnerThread.terminate()
    ### Genbank Download Functions #####
    def redoGbDlEstimate(self):
        self.checkRequestBtn.setEnabled(True)
        self.totalFilesGbDiv.setText('TBD')
        self.estSize.setText('TBD')
        self.downloadGbDivFileBtn.setEnabled(False)
    def checkGbDivDlRequest(self):
        gbDivsToDl = [self.gbDivList[idx] for idx,checkbox in enumerate(self.checkBoxList) if checkbox.checkState() == 2]
        self.gbDivFilesToDlList = []
        totalSize = 0
        if self.verbose:
            print(gbDivsToDl)
        for keyword in gbDivsToDl:
            fileList,sizes = zip(*getGbkDlList(keyword))
            if self.verbose:
                print(fileList)
            self.gbDivFilesToDlList.extend(fileList)
            totalSize += sum(sizes)
        self.totalFilesGbDiv.setText(str(len(self.gbDivFilesToDlList)))
        self.estSize.setText(humanbytes(totalSize))
        self.checkRequestBtn.setEnabled(False)
        self.downloadGbDivFileBtn.setEnabled(True)
    def selectGbDivDir(self):
        dirName = QFileDialog.getExistingDirectory()
        if dirName:
            if os.access(dirName, os.W_OK):
               self.gbDivDLoutputDir.setText(dirName)
            else:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText("Can't Write to that Folder. Please Specify Another Directory")
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()
    # will ignore files of a similar name that were already downloaded
    def downloadGenbankDivFiles(self):
        gbDivFileDir = self.gbDivDLoutputDir.text()
        tasklist = set(x + '.clusterToolsDB' for x in self.gbDivFilesToDlList)
        if not gbDivFileDir:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("No Output Directory Specified")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
        elif len(tasklist) == 0:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Nothing in Download List")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
        else:
            alreadyDownloaded = set(os.listdir(gbDivFileDir))
            targetlist = [x.split('.clusterToolsDB')[0] for x in tasklist-alreadyDownloaded]
            if self.verbose:
                print(targetlist)
            if not self.runnerThread:
                self.runnerThread = QThread()
            self.runnerThread.start()
            self.genbankDivFilesStatusWin = downloadNcbiFilesWindow()

            self.downloadGbDivFilesWorker = downloadGbDivFilesWorker(targetlist,gbDivFileDir)
            self.genbankDivFilesStatusWin.doneBtn.clicked.connect(self.genbankDivFilesStatusWin.close)

            self.genbankDivFilesStatusWin.cancelBtn.clicked.connect(
                lambda: self.abortSearch(self.genbankDivFilesStatusWin, self.runnerThread))

            self.genbankDivFilesStatusWin.progressBar.setMaximum(int(len(tasklist - alreadyDownloaded)))
            self.genbankDivFilesStatusWin.progressBar.setValue(0)
            self.genbankDivFilesStatusWin.show()
            self.downloadGbDivFilesWorker.moveToThread(self.runnerThread)
            self.downloadGbDivFilesWorker.start.connect(self.downloadGbDivFilesWorker.run)
            self.downloadGbDivFilesWorker.currentFile.connect(self.updateCurrentGbDlFile)
            self.downloadGbDivFilesWorker.finished.connect(lambda x: self.gbDivFileDlFinished(x,gbDivFileDir))
            self.downloadGbDivFilesWorker.start.emit()

    def gbDivFileDlFinished(self,dlFailedList,targetDir):
        if dlFailedList:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("Failed to Download Some Files: {}\n"
                        "Proceeding with DB Creation With Downloaded Files".format('\n'.join(dlFailedList)))
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()

        filesToProcess = glob(os.path.join(targetDir,'*.gz.clusterToolsDB'))
        if self.verbose:
            print(filesToProcess)
        self.genbankDivFilesStatusWin.title.setText('Processing Genbank Division File:')
        self.genbankDivFilesStatusWin.currentTask.setText(' : ')
        self.genbankDivFilesStatusWin.progressBar.setMaximum(int(len(filesToProcess)))
        self.genbankDivFilesStatusWin.progressBar.setValue(0)

        if not self.runnerThread:
            self.runnerThread = QThread()
        self.runnerThread.start()
        self.processGbDivFilesWorker = processGbDivFilesWorker(filesToProcess,targetDir)
        self.processGbDivFilesWorker.moveToThread(self.runnerThread)
        self.processGbDivFilesWorker.start.connect(self.processGbDivFilesWorker.run)
        self.processGbDivFilesWorker.currentFile.connect(self.updateGbProcessFile)
        self.processGbDivFilesWorker.currentSpecies.connect(self.updateGbProcessSpecies)
        self.processGbDivFilesWorker.finished.connect(lambda x: self.gbProcessFileFinished(x,dlFailedList))
        self.processGbDivFilesWorker.start.emit()

    def updateCurrentGbDlFile(self,currentFile):
        self.genbankDivFilesStatusWin.currentTask.setText(currentFile)
        self.genbankDivFilesStatusWin.progressBar.setValue(self.genbankDivFilesStatusWin.progressBar.value() + 1)
        self.genbankDivFilesStatusWin.progressLabel.setText('{}/{}'.format(self.genbankDivFilesStatusWin.progressBar.value(),
                                                                       self.genbankDivFilesStatusWin.progressBar.maximum()))

    def updateGbProcessFile(self,currentFile):
        self.genbankDivFilesStatusWin.currentTask.setText('{} : {}'.format(currentFile,' '))
        self.genbankDivFilesStatusWin.progressBar.setValue(self.genbankDivFilesStatusWin.progressBar.value()+1)
        self.genbankDivFilesStatusWin.progressLabel.setText('{}/{}'.format(self.genbankDivFilesStatusWin.progressBar.value(),
                                                                      self.genbankDivFilesStatusWin.progressBar.maximum()))

    def updateGbProcessSpecies(self,currentSpecies):
        gbkFile, gbkSpecies = self.genbankDivFilesStatusWin.currentTask.text().split(' : ')
        self.genbankDivFilesStatusWin.currentTask.setText('{} : {}'.format(gbkFile,currentSpecies))

    def gbProcessFileFinished(self,processFailed,dlFailed):
        self.genbankDivFilesStatusWin.title.setText('Database Completed!')
        self.genbankDivFilesStatusWin.doneBtn.setEnabled(True)
        failedFiles = ''
        if dlFailed:
            failedFiles += 'Failed to Download: {}\n'.format('\n'.join(dlFailed))
        if processFailed:
            failedFiles += 'Failed to Process: {}'.format('\n'.join(processFailed))
        if failedFiles:
            self.genbankDivFilesStatusWin.currentTask.setText(failedFiles)
        if self.runnerThread:
            self.runnerThread.terminate()

    ###########################
    def showAddGeneWin(self):
        self.addGeneWin = addGeneWindow()
        self.addGeneWin.addGeneBtn.clicked.connect(self.addGeneSeq)
        self.addGeneWin.show()
    def addGeneSeq(self):
        geneName = self.addGeneWin.geneName.text()
        geneSeqToAdd = self.addGeneWin.geneSequence.toPlainText()
        geneSeqToAdd = geneSeqToAdd.replace('\n','')
        geneSeqToAdd = geneSeqToAdd.replace('\t', '')
        geneSeqToAdd = geneSeqToAdd.strip()
        geneSeqToAdd = geneSeqToAdd.upper()
        if geneName:
            if geneSeqToAdd:
                if geneName in self.geneDict.keys():
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Warning)
                    msg.setText('There is already a Gene Named %s. Please choose another name.' % geneName)
                    msg.setStandardButtons(QMessageBox.Ok)
                    msg.exec()
                    self.addGeneWin.geneName.setText('')
                else:
                    if self.verbose:
                        print(self.addGeneWin.AminoAcidValidator.validate(geneSeqToAdd,0))
                    if self.addGeneWin.AminoAcidValidator.validate(geneSeqToAdd,0)[0] == 2:
                        self.geneDict[geneName] = str(geneSeqToAdd)
                        geneItem = QListWidgetItem()
                        geneItem.setText(geneName)
                        geneItem.setFlags(Qt.ItemIsEditable|Qt.ItemIsSelectable | Qt.ItemIsEnabled)
                        self.geneList.addItem(geneItem)
                        self.addGeneWin.close()
                    else:
                        msg = QMessageBox()
                        msg.setIcon(QMessageBox.Warning)
                        msg.setText('Please Enter A Valid Amino Acid Sequence')
                        msg.setStandardButtons(QMessageBox.Ok)
                        msg.exec()
                        self.addGeneWin.geneSequence.setText('')
            else:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setText('No Amino Acid Sequence Detected')
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()
                self.addGeneWin.geneSequence.setText('')
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText('No Gene Name Specified')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
            self.addGeneWin.geneSequence.setText('')
    def updateTotalHitsRequired(self,value):
        self.totalHitsRequired = value
        if self.verbose:
            print("Total Hits Required: ",self.totalHitsRequired)
    def updateSelectedGene(self,rowNum):
        if self.geneList.currentItem():
            self.currentGeneSelected = self.geneList.currentItem().text()
        else:
            self.currentGeneSelected = ''
        if self.verbose:
            print("Updated Gene Selected", self.currentGeneSelected)
    def updateGeneName(self):
        if self.verbose:
            print("Modifying:",self.currentGeneSelected,self.geneDict[self.currentGeneSelected])
        newName = self.geneList.currentItem().text()
        ## Check that you actually edited the gene
        if newName != self.currentGeneSelected:
            ## Check if there is already a gene with that name
            if newName not in self.geneDict.keys():
                if self.verbose: print("Changing:", self.currentGeneSelected,"to ",newName)
                self.geneDict[newName] = self.geneDict[self.currentGeneSelected]
                searchListIdxs = []
                for idx in range(self.searchList.rowCount()):
                    if self.searchList.item(idx, 0).text() == 'GENE' and self.searchList.item(idx,1).text() == self.currentGeneSelected:
                        searchListIdxs.append(idx)
                if searchListIdxs:
                    self.forBLAST[newName] = self.forBLAST[self.currentGeneSelected]
                    del self.forBLAST[self.currentGeneSelected]
                    for idx in searchListIdxs:
                        self.searchList.item(idx, 1).setText(newName)
                self.updateSpinBox()

                del self.geneDict[self.currentGeneSelected]
                self.currentGeneSelected = newName
                if self.verbose: print("Dict Value Changed: ",
                                       self.currentGeneSelected,self.geneDict[self.currentGeneSelected])
            else:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setText('There is already a gene named %s. Please Choose Another Name' % newName)
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()
                self.geneList.currentItem().setText(self.currentGeneSelected)
    ### Utility Functions
    def updateParams(self):
        parameterObjs = [self.paramsWin.blastEval,self.paramsWin.hmmEval,self.paramsWin.hmmDomLen,
                         self.paramsWin.hmmScore,self.paramsWin.windowSize]
        print(self.nameToParamDict)
        for parameter in parameterObjs:
            parameterName = parameter.objectName()
            validator = parameter.validator()
            state = validator.validate(parameter.text(),0)[0]
            if state == QtGui.QValidator.Acceptable:
                if parameterName in {'blastEval','hmmEval'}:
                    self.nameToParamDict[parameterName] = float(parameter.text())
                else:
                    self.nameToParamDict[parameterName] = int(parameter.text())
            else:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setText('Please Input a Valid Value for %s' % parameterName)
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()
                parameter.setText(str(self.nameToParamDict[parameterName]))
        print(self.nameToParamDict)
        self.paramsWin.applyBtn.setEnabled(False)

    def showParamWindow(self):
        self.paramsWin.blastEval.setText(str(self.nameToParamDict[self.paramsWin.blastEval.objectName()]))
        self.paramsWin.hmmEval.setText(str(self.nameToParamDict[self.paramsWin.hmmEval.objectName()]))
        self.paramsWin.hmmDomLen.setText(str(self.nameToParamDict[self.paramsWin.hmmDomLen.objectName()]))
        self.paramsWin.hmmScore.setText(str(self.nameToParamDict[self.paramsWin.hmmScore.objectName()]))
        self.paramsWin.windowSize.setText(str(self.nameToParamDict[self.paramsWin.windowSize.objectName()]))
        self.paramsWin.show()

    def activateApply(self):
        self.paramsWin.applyBtn.setEnabled(True)

    def updateCheckSuccessful(self,flag):
        self.checkSuccessful = flag

    ### Genbank Panel Functions
    def loadGbkDir(self):
        if sys.platform == 'win32':
            dirName = QFileDialog.getExistingDirectory(options=QFileDialog.DontUseNativeDialog)
        else:
            dirName = QFileDialog.getExistingDirectory()
        if dirName and os.access(dirName, os.R_OK):
            self.genbankDirPath.setText(dirName)
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Can't Read Files In Folder")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
    def addGbkFiles(self):
        gbkPath = self.genbankDirPath.text()
        if gbkPath:
            gbkList = chain(iglob(gbkPath+'/*.gbk'),iglob(gbkPath+'/*.gb'))
            for gbkFile in gbkList:
                self.genbankList.addItem(gbkFile)
            self.totalFiles.setText(str(self.genbankList.count()))
            self.genbankDirPath.clear()
    def removeGbkFiles(self):
        gbksToRemove = self.genbankList.selectedItems()
        for gbkFile in gbksToRemove:
            self.genbankList.takeItem(self.genbankList.row(gbkFile))
        self.totalFiles.setText(str(self.genbankList.count()))
    def chooseDBoutDir(self):
        dbName , _ = QFileDialog.getSaveFileName(caption="Specify Database Name",filter='*.fasta')
        if dbName:
            dirName,fileName = os.path.split(dbName)
            if not fileName.endswith('.fasta'):
                fileName += '.fasta'
            if dirName and os.access(dirName, os.W_OK):
                self.dbOutputDir.setText(os.path.join(dirName,fileName))
            else:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText("Can't Write to that Folder")
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()
    def createDB(self):
        dbPath = self.dbOutputDir.text()
        taskList = [genbankFile.text() for genbankFile in (self.genbankList.item(x) for x in range(self.genbankList.count()))]

        if not dbPath:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("No Output Directory Specified")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
        elif len(taskList) == 0:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("No Genbank Files to Process")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
        else:
            self.runnerThread = QThread()
            self.runnerThread.start()
            if self.verbose:
                print(taskList)
            self.createDbWorker = createDbWorker(taskList,dbPath)
            self.createDbStatusWin = createDbStatusWindow()
            self.createDbStatusWin.doneBtn.clicked.connect(self.createDbStatusWin.close)
            self.createDbStatusWin.cancelBtn.clicked.connect(lambda: self.abortSearch(self.createDbStatusWin,
                                                                                      self.runnerThread))
            self.createDbStatusWin.progressBar.setMaximum(int(self.totalFiles.text()))
            self.createDbWorker.moveToThread(self.runnerThread)
            self.createDbWorker.start.connect(self.createDbWorker.run)
            self.createDbWorker.currentGbk.connect(self.updateCurrentGbk)
            self.createDbWorker.finished.connect(self.createDbSuccessful)
            self.createDbWorker.start.emit()
            self.createDbStatusWin.show()

    def updateCurrentGbk(self,currentGbk):
        self.createDbStatusWin.currentTask.setText(currentGbk)
        self.createDbStatusWin.progressBar.setValue(self.createDbStatusWin.progressBar.value()+1)
        self.createDbStatusWin.percentLabel.setText("{}/{}".format(self.createDbStatusWin.progressBar.value(),
                                                                   self.totalFiles.text()))
    def createDbSuccessful(self):
        self.genbankList.clear()
        self.createDbStatusWin.doneBtn.setEnabled(True)
        if self.runnerThread:
            self.runnerThread.terminate()
    ### Status Window Methods
    def updateStatusWinText(self,statusWin,currentTaskText,percentCmpLabelText,percentCmpBarValue):
        statusWin.currentTask.setText(currentTaskText)
        statusWin.percentCmpLabel.setText(percentCmpLabelText)
        statusWin.percentCmpBar.setValue(percentCmpBarValue)
    def abortSearch(self,statusWin,currentThread):
        statusWin.close()
        currentThread.terminate()
    ### Results Window Methods
    def showResultsWindow(self,blastList,hmmList,filteredClusters):
        self.statusWin.close()
        self.checkSuccessful = False
        self.resultsWin = resultsWindow()

        self.resultsWin.resultsList.setColumnCount(4 + len(blastList) + len(hmmList))
        self.resultsWin.resultsList.setHorizontalHeaderItem(0, QTableWidgetItem("Species/Accession ID"))
        self.resultsWin.resultsList.setHorizontalHeaderItem(1, QTableWidgetItem("Start"))
        self.resultsWin.resultsList.setHorizontalHeaderItem(2, QTableWidgetItem("End"))
        self.resultsWin.resultsList.setHorizontalHeaderItem(3, QTableWidgetItem("Size"))
        for idx, blastHit in enumerate(blastList):
            if self.verbose:
                print(idx, blastHit)
            self.resultsWin.resultsList.setHorizontalHeaderItem(idx + 4, QTableWidgetItem(blastHit))
        for idx, hmmHit in enumerate(hmmList):
            if self.verbose: print(idx,' and '.join(hmmHit))
            self.resultsWin.resultsList.setHorizontalHeaderItem(idx + len(blastList) + 4,
                                                          QTableWidgetItem(' and '.join(hmmHit)))
        for cluster, hitDict in filteredClusters.items():
            currentRowCount = self.resultsWin.resultsList.rowCount()
            self.resultsWin.resultsList.insertRow(currentRowCount)
            self.resultsWin.resultsList.setItem(currentRowCount, 0, QTableWidgetItem(cluster[0]))
            self.resultsWin.resultsList.setItem(currentRowCount, 1, QTableWidgetItem(str(cluster[1])))
            self.resultsWin.resultsList.setItem(currentRowCount, 2, QTableWidgetItem(str(cluster[2])))
            self.resultsWin.resultsList.setItem(currentRowCount, 3, QTableWidgetItem(str(cluster[2]-cluster[1]+1)))
            for idx, blastHit in enumerate(blastList):
                self.resultsWin.resultsList.setItem(currentRowCount, idx + 4,
                                              QTableWidgetItem('; '.join(filteredClusters[cluster].get(blastHit,''))))
            for idx, hmmHit in enumerate(hmmList):
                self.resultsWin.resultsList.setItem(currentRowCount, idx + len(blastList) + 4,
                                              QTableWidgetItem('; '.join(filteredClusters[cluster].get(hmmHit,''))))

        self.resultsWin.resultsList.resizeColumnsToContents()
        self.resultsWin.exportResults.clicked.connect(self.exportResults)
        self.resultsWin.show()
    def exportResults(self):
        totalFilePath, _ = QFileDialog.getSaveFileName(caption="Export Search Results",filter='*.csv')
        if totalFilePath:
            path,fileName = os.path.split(totalFilePath)
            if os.access(path, os.W_OK):
                with open(totalFilePath, 'w') as stream:
                    writer = csv.writer(stream)
                    rowData = [self.resultsWin.resultsList.horizontalHeaderItem(idx).text().strip()
                               for idx in range(self.resultsWin.resultsList.columnCount())]
                    writer.writerow(rowData)
                    for rowIdx in range(self.resultsWin.resultsList.rowCount()):
                        rowData = [self.resultsWin.resultsList.item(rowIdx, colIdx).text().strip()
                                   for colIdx in range(self.resultsWin.resultsList.columnCount())]
                        writer.writerow(rowData)
            else:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText("Can't Write to Directory")
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()

    ### Run Permissions Checks
    def runChecks(self):
        if self.searchList.rowCount() >= 1:
            if not self.checkSuccessful:
                self.runnerThread = QThread()
                self.runnerThread.start()

                self.checkSuccessful = False
                self.blastSuccessful = False
                self.hmmerSuccessful = False

                self.checksWorker = runChecksWorker(self.outputDir,self.pathToDatabase)
                self.checksWorker.checkError.connect(lambda x: self.doneCheck(x,self.statusWin))
                self.statusWin = statusWindow()
                self.statusWin.percentCmpBar.setMaximum(6)
                self.updateStatusWinText(self.statusWin, 'Running Database and Output Checks', '1/6', 1)
                self.statusWin.cancelBtn.clicked.connect(lambda: self.abortSearch(self.statusWin,self.runnerThread))
                self.statusWin.show()
                self.checksWorker.moveToThread(self.runnerThread)
                self.checksWorker.start.connect(self.checksWorker.run)
                self.checksWorker.start.emit('whee')
                self.checksWorker.finished.connect(self.updateCheckSuccessful)
            else:
                self.runBlast()
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('No Entries in Search List')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
    @pyqtSlot(str)
    def doneCheck(self,checkErr,statusWin):
        if checkErr == 'noOutputDir':
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('Please Specify An Output Directory')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
            statusWin.close()
            if self.runnerThread:
                self.runnerThread.terminate()
        elif checkErr == 'noDBpath':
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('Please Specify A Database to Query')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
            statusWin.close()
            if self.runnerThread:
                self.runnerThread.terminate()
        elif checkErr == 'noWriteOutputDir':
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('Cannot Write to Output Folder Specified')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
            statusWin.close()
            if self.runnerThread:
                self.runnerThread.terminate()
        elif checkErr == 'noReadDB':
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('Cannot Read Database')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
            statusWin.close()
            if self.runnerThread:
                self.runnerThread.terminate()
        elif checkErr == 'invalidDB':
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('No Valid Entries in Database')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
            statusWin.close()
            if self.runnerThread:
                self.runnerThread.terminate()
        elif checkErr == 'entriesWarn':
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText('Warning: Some Entries in Database are Incorrectly Formatted, these will be ignored.')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
            self.checkSuccessful = True
            self.runBlast()
            return
        elif checkErr == 'good':
            self.checkSuccessful = True
            self.runBlast()
            return
    ### Run BLAST Methods
    def runBlast(self):
        self.updateStatusWinText(self.statusWin, 'Running BLAST Search: Generating Input Fasta', '2/6',
                                 self.statusWin.percentCmpBar.value() + 1)
        if self.checkSuccessful:
            if not self.blastSuccessful:
                if self.forBLAST:
                    self.blastWorker = runBlastWorker(self.forBLAST, self.outputDir, self.pathToDatabase,
                                                      self.makeblastdbExec, self.blastExec,self.nameToParamDict['blastEval'])

                    self.blastWorker.inputGenerated.connect(lambda x: self.generateFasta(x, self.statusWin))
                    self.blastWorker.makeDB.connect(lambda x: self.checkDB(x, self.statusWin))
                    self.blastWorker.dbCreated.connect(lambda x: self.makeBlastDB(x, self.statusWin))
                    self.blastWorker.doneBLAST.connect(lambda x: self.doneBLAST(x, self.statusWin))
                    if not self.runnerThread:
                        self.runnerThread = QThread()
                    self.runnerThread.start()
                    self.blastWorker.moveToThread(self.runnerThread)
                    self.blastWorker.start.connect(self.blastWorker.run)
                    self.blastWorker.start.emit()
                else:
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Information)
                    msg.setText("No genes in search list, skipping BLAST step...")
                    msg.setStandardButtons(QMessageBox.Ok)
                    msg.exec()
                    self.blastSuccessful = True
                    self.runHmmer()
            else:
                self.runHmmer()
    @pyqtSlot(bool)
    def generateFasta(self,flag,statusWin):
        if flag:
            self.updateStatusWinText(statusWin, 'Running BLAST Search: Checking for BLAST Formatted DB', '2/6',
                                     statusWin.percentCmpBar.value())
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('Generating Input Fasta File Failed')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
            statusWin.close()
            if self.runnerThread:
                self.runnerThread.terminate()
    @pyqtSlot(bool)
    def checkDB(self,flag,statusWin):
        if flag:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("No Blastp formatted database found, creating one.")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
            self.updateStatusWinText(statusWin, 'Running BLAST Search: Creating Blastp database', '2/6',
                                     statusWin.percentCmpBar.value())
        else:
            self.updateStatusWinText(statusWin,
                                     'Running BLAST Search: Found Blastp database, \n'
                                     'Performing Blastp with Query Proteins \n (E-Value: %.1e)' % self.nameToParamDict['blastEval'],
                                     '2/6',
                                     statusWin.percentCmpBar.value())
    @pyqtSlot(bool)
    def makeBlastDB(self,flag,statusWin):
        if flag:
            self.updateStatusWinText(statusWin,
                                     'Running BLAST Search: '
                                    'Performing Blastp with Query Proteins \n(E-Value: %.1e)' % self.nameToParamDict['blastEval'],
                                     '2/6',
                                     statusWin.percentCmpBar.value())
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('Creating Blast Database Failed')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
            statusWin.close()
            if self.runnerThread:
                self.runnerThread.terminate()
    @pyqtSlot(bool)
    def doneBLAST(self,flag,statusWin):
        if flag:
            self.updateStatusWinText(statusWin,'Running HMMSearch: Extracting HMMs to make reference file',
                                     '3/6',statusWin.percentCmpBar.value()+1)
            self.blastSuccessful = True
            self.runHmmer()
            return
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('Blastp Failed')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
            statusWin.close()
            if self.runnerThread:
                self.runnerThread.terminate()
            return
    ### Run Hmmer Methods
    def runHmmer(self):
        if self.checkSuccessful:
            self.updateStatusWinText(self.statusWin, 'Running HMMSearch: Extracting HMMs to make reference file',
                                     '3/6', self.statusWin.percentCmpBar.value() + 1)
            if not self.hmmerSuccessful:
                hmmSet = set()
                hmmSet.update(*self.forHmmer.keys())
                if hmmSet:
                    self.hmmerWorker = runHmmerWorker(self.hmmFetchExec, self.hmmSearchExec,
                                                      self.hmmDict, hmmSet, self.outputDir,self.pathToDatabase,
                                                      self.nameToParamDict['hmmEval'])
                    self.hmmerWorker.hmmDBbuilt.connect(lambda x, y: self.hmmDBbuiltCheck(x, y, self.statusWin))
                    self.hmmerWorker.hmmSearchComplete.connect(lambda x: self.hmmSuccessCheck(x, self.statusWin))
                    if not self.runnerThread:
                        self.runnerThread = QThread()
                    self.runnerThread.start()
                    self.hmmerWorker.moveToThread(self.runnerThread)
                    self.hmmerWorker.start.connect(self.hmmerWorker.run)
                    self.hmmerWorker.start.emit()
                else:
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Information)
                    msg.setText("No HMMs to search, skipping HmmSearch Step")
                    msg.setStandardButtons(QMessageBox.Ok)
                    msg.exec()
                    self.hmmerSuccessful = True
                    self.processResultsOptionalHits()
            else:
                self.processResultsOptionalHits()
    @pyqtSlot(bool,list)
    def hmmDBbuiltCheck(self,flag,hmmsNotFetched,statusWin):
        if flag:
            self.updateStatusWinText(self.statusWin,
                                     'Running HMMSearch: Performing HmmSearch (E-Value: %.1e)' % self.nameToParamDict['hmmEval'],
                                     '3/6',
                                     self.statusWin.percentCmpBar.value())
            return
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('Failed Finding Following HMMs: %s' % ','.join(hmmsNotFetched))
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
            msg.exec()
            statusWin.close()
            if self.runnerThread:
                self.runnerThread.terminate()
            return
    @pyqtSlot(bool)
    def hmmSuccessCheck(self,flag,statusWin):
        if flag:
            self.hmmerSuccessful = True
            self.processResultsOptionalHits()
            return
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('Hmmsearch Failed')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
            statusWin.close()
            if self.runnerThread:
                self.runnerThread.terminate()
    ### Process Filtered Clusters
    def processResultsOptionalHits(self):
        self.updateStatusWinText(self.statusWin, 'Parsing Output Files: Creating Search Lists', '4/6', 4)
        if self.checkSuccessful and self.blastSuccessful and self.hmmerSuccessful:
            requiredBlastList = []
            requiredHmmList = []

            optionalBlastList = []
            optionalHmmList = []

            blastOutFile = self.outputDir + "/blast_results.out"
            hmmOutFile = self.outputDir + '/hmmSearch.out'

            for idx in range(self.searchList.rowCount()):
                if self.searchList.item(idx, 0).text() == 'HMMER Hits':
                    if self.searchList.item(idx,2).checkState()== 2:
                        requiredHmmList.append(tuple(sorted(self.searchList.item(idx, 1).text().split(' and '))))
                    else:
                        optionalHmmList.append(tuple(sorted(self.searchList.item(idx, 1).text().split(' and '))))
                elif self.searchList.item(idx, 0).text() == 'GENE':
                    if self.searchList.item(idx, 2).checkState() == 2:
                        requiredBlastList.append(self.searchList.item(idx, 1).text())
                    else:
                        optionalBlastList.append(self.searchList.item(idx, 1).text())

            passBLASTcheck = False
            passHMMcheck = False

            ### Final Checks before Passing on to Parsers
            if len(requiredBlastList) + len(requiredHmmList) + len(optionalBlastList) + len(optionalHmmList) == 0:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText('No HMM Hits or Genes Included in Search List')
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()
                self.statusWin.close()
                self.checkSuccessful = False
                if self.runnerThread:
                    self.runnerThread.terminate()
                return
            else:
                if requiredBlastList or optionalBlastList:
                    try:
                        assert os.path.isfile(blastOutFile), "BLAST Hits included in search but no " \
                                                             "BLAST output file found"
                        passBLASTcheck = True
                    except AssertionError:
                        msg = QMessageBox()
                        msg.setIcon(QMessageBox.Critical)
                        msg.setText('BLAST Hits included in search but no BLAST output file found!')
                        msg.setStandardButtons(QMessageBox.Ok)
                        msg.exec()
                        self.statusWin.close()
                        self.checkSuccessful = False
                        if self.runnerThread:
                            self.runnerThread.terminate()
                        return
                else:
                    passBLASTcheck = True
                if requiredHmmList or optionalHmmList:
                    try:
                        assert os.path.isfile(hmmOutFile), "HMMer Hits included in search but no " \
                                                                             "HMMer output file found"
                        passHMMcheck = True
                    except AssertionError:
                        msg = QMessageBox()
                        msg.setIcon(QMessageBox.Critical)
                        msg.setText('HMMer Hits requested in search but no HMMer output file found!')
                        msg.setStandardButtons(QMessageBox.Ok)
                        msg.exec()
                        self.statusWin.close()
                        self.checkSuccessful = False
                        if self.runnerThread:
                            self.runnerThread.terminate()
                        return

                else:
                    passHMMcheck = True
            if passBLASTcheck and passHMMcheck:
                self.updateStatusWinText(self.statusWin,
                                         'Output Files Parsed. Searching for Clusters.\n'
                                         '(Minimum Domain Score: %i, Minimum Domain Length: %i Window Size: %i)'
                    % (self.nameToParamDict['hmmScore'],self.nameToParamDict['hmmDomLen'],self.nameToParamDict['windowSize']),
                                         '5/6', 5)
                if not self.runnerThread:
                    self.runnerThread = QThread()
                self.runnerThread.start()
                if self.totalHitsRequired == len(requiredBlastList) + len(requiredHmmList) :
                    self.processSearchListWorker = runProcessSearchList(requiredBlastList, requiredHmmList, blastOutFile, hmmOutFile,
                                                                        self.nameToParamDict['hmmScore'],
                                                                        self.nameToParamDict['hmmDomLen'],
                                                                        self.nameToParamDict['windowSize'])
                    self.processSearchListWorker.result.connect(
                        lambda x: self.generateResults(x, self.statusWin, requiredBlastList, requiredHmmList))
                    self.processSearchListWorker.moveToThread(self.runnerThread)
                    self.processSearchListWorker.start.connect(self.processSearchListWorker.run)
                    self.processSearchListWorker.start.emit()
                else:
                    self.processSearchListWorker = runProcessSearchListOptionalHits(requiredBlastList,requiredHmmList,
                                                                                    blastOutFile,hmmOutFile,
                                                                                    self.nameToParamDict['hmmScore'],
                                                                                    self.nameToParamDict['hmmDomLen'],
                                                                                    self.nameToParamDict['windowSize'],
                                                                                    self.totalHitsRequired,
                                                                                    optionalBlastList,optionalHmmList)
                    self.processSearchListWorker.result.connect(
                        lambda x: self.generateResults(x, self.statusWin,requiredBlastList+optionalBlastList,
                                                       requiredHmmList+optionalHmmList))
                    self.processSearchListWorker.moveToThread(self.runnerThread)
                    self.processSearchListWorker.start.connect(self.processSearchListWorker.run)
                    self.processSearchListWorker.start.emit()

    def processResults(self):
        self.updateStatusWinText(self.statusWin, 'Parsing Output Files: Creating Search Lists', '4/6', 4)
        if self.checkSuccessful and self.blastSuccessful and self.hmmerSuccessful:
            blastList = []
            hmmList = []

            blastOutFile = self.outputDir + "/blast_results.out"
            hmmOutFile = self.outputDir + '/hmmSearch.out'

            for idx in range(self.searchList.rowCount()):
                if self.searchList.item(idx, 0).text() == 'HMMER Hits':
                    hmmList.append(tuple(sorted(self.searchList.item(idx, 1).text().split(' and '))))
                elif self.searchList.item(idx, 0).text() == 'GENE':
                    blastList.append(self.searchList.item(idx, 1).text())

            passBLASTcheck = False
            passHMMcheck = False

            ### Final Checks before Passing on to Parsers
            if len(blastList) + len(hmmList) == 0:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText('No HMM Hits or Genes Included in Search List')
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()
                self.statusWin.close()
                self.checkSuccessful = False
                if self.runnerThread:
                    self.runnerThread.terminate()
                return
            else:
                if blastList:
                    try:
                        assert os.path.isfile(blastOutFile), "BLAST Hits included in search but no " \
                                                             "BLAST output file found"
                        passBLASTcheck = True
                    except AssertionError:
                        msg = QMessageBox()
                        msg.setIcon(QMessageBox.Critical)
                        msg.setText('BLAST Hits included in search but no BLAST output file found!')
                        msg.setStandardButtons(QMessageBox.Ok)
                        msg.exec()
                        self.statusWin.close()
                        self.checkSuccessful = False
                        if self.runnerThread:
                            self.runnerThread.terminate()
                        return
                else:
                    passBLASTcheck = True
                if hmmList:
                    try:
                        assert os.path.isfile(hmmOutFile), "HMMer Hits included in search but no " \
                                                                             "HMMer output file found"
                        passHMMcheck = True
                    except AssertionError:
                        msg = QMessageBox()
                        msg.setIcon(QMessageBox.Critical)
                        msg.setText('HMMer Hits requested in search but no HMMer output file found!')
                        msg.setStandardButtons(QMessageBox.Ok)
                        msg.exec()
                        self.statusWin.close()
                        self.checkSuccessful = False
                        if self.runnerThread:
                            self.runnerThread.terminate()
                        return

                else:
                    passHMMcheck = True
            if passBLASTcheck and passHMMcheck:
                self.updateStatusWinText(self.statusWin,
                                         'Output Files Parsed. Searching for Clusters.\n'
                                         '(Minimum Domain Score: %i, Minimum Domain Length: %i Window Size: %i)'
                    % (self.nameToParamDict['hmmScore'],self.nameToParamDict['hmmDomLen'],self.nameToParamDict['windowSize']),
                                         '5/6', 5)
                if not self.runnerThread:
                    self.runnerThread = QThread()
                self.runnerThread.start()
                self.processSearchListWorker = runProcessSearchList(blastList, hmmList, blastOutFile, hmmOutFile,
                                                                    self.nameToParamDict['hmmScore'],
                                                                    self.nameToParamDict['hmmDomLen'],
                                                                    self.nameToParamDict['windowSize'])
                self.processSearchListWorker.result.connect(
                    lambda x: self.generateResults(x, self.statusWin, blastList, hmmList))
                self.processSearchListWorker.moveToThread(self.runnerThread)
                self.processSearchListWorker.start.connect(self.processSearchListWorker.run)
                self.processSearchListWorker.start.emit()
    @pyqtSlot(dict)
    def generateResults(self,filteredClusters,statusWin,blastList,hmmList):
        if filteredClusters:
            self.updateStatusWinText(statusWin, 'Search Complete',
                                     '6/6', 6)
            statusWin.viewResultsBtn.setEnabled(True)
            statusWin.viewResultsBtn.clicked.connect(lambda: self.showResultsWindow(blastList,hmmList,filteredClusters))
            if self.runnerThread:
                self.runnerThread.terminate()
        else:
            self.updateStatusWinText(statusWin, 'Search Complete: No Results Found',
                                     '6/6', 6)
            self.checkSuccessful = False
            if self.runnerThread:
                self.runnerThread.terminate()
            return
    #### GUI Functions #####
    def updateSpinBox(self):
        if self.searchList.rowCount() == 0:
            self.hitsNeededSpinBox.setEnabled(False)
        else:
            maxHits = self.searchList.rowCount()
            minHits = sum(1 for idx in range(self.searchList.rowCount())
                          if self.searchList.item(idx,2) and
                          self.searchList.item(idx,2).checkState()==2)
            if self.verbose:
                print("Setting range: %i-%i" % (minHits,maxHits))
            self.hitsNeededSpinBox.setRange(minHits,maxHits)
            self.hitsNeededSpinBox.setValue(minHits)
            self.hitsNeededSpinBox.setEnabled(True)
        self.totalHitsRequired = self.hitsNeededSpinBox.value()
        if self.verbose:
            print("Total Hits Required: ",self.totalHitsRequired)

    def removeSearchTerm(self):
        itemType =  self.searchList.item(self.searchList.currentRow(),0).text()
        termToRemove = self.searchList.item(self.searchList.currentRow(),1).text()
        if itemType == 'GENE':
            if self.verbose:
                print("Deleted %s" % termToRemove)
            del self.forBLAST[termToRemove]
            self.searchList.removeRow(self.searchList.currentRow())
        elif itemType == 'HMMER Hits':
            hmmKey = tuple(sorted(termToRemove.split(' and ')))
            hmmerSearchTerms = Counter()
            for idx in range(self.searchList.rowCount()):
                if self.searchList.item(idx,0).text() == 'HMMER Hits':
                    hmmerSearchTerms[tuple(sorted(self.searchList.item(idx,1).text().split(' and ')))] += 1
            if hmmerSearchTerms[hmmKey] <= 1:
                del self.forHmmer[hmmKey]
            if self.verbose: print(self.forHmmer)
            self.searchList.removeRow(self.searchList.currentRow())
        self.updateSpinBox()
    def outdirFunction(self):
        if 'Select Directory...' in self.outputDirectorySelector.currentText():
            dirName = QFileDialog.getExistingDirectory(self)
            if dirName:
                self.outputDirectorySelector.addItem(dirName)
                self.outputDirectorySelector.setCurrentIndex(self.outputDirectorySelector.count()-1)
        if 'Select Directory...' not in self.outputDirectorySelector.currentText():
            self.outputDir = self.outputDirectorySelector.currentText()
            if self.verbose: print(self.outputDir)
    def databaseFunction(self):
        if 'add database...' in self.dataBaseSelector.currentText():
            fileName, _ = QFileDialog.getOpenFileName(self)
            if fileName:
                self.dataBaseSelector.addItem(fileName)
                self.dataBaseSelector.setCurrentIndex(self.dataBaseSelector.count()-1)
        if 'add database...' not in self.dataBaseSelector.currentText():
            self.pathToDatabase = self.dataBaseSelector.currentText()
            if self.verbose: print(self.pathToDatabase)
    def openGene(self):
        fileName, _ = QFileDialog.getOpenFileName(self)
        if self.verbose: print(fileName)
        if fileName:
            self.GeneFilePath.setText(fileName)
    def loadGeneFile(self):
        if self.GeneFilePath.text():
            valid_extensions = set(['gbk','fa','fasta'])
            extension = self.GeneFilePath.text().split('.')[-1]
            if extension in valid_extensions:
                genesToAdd,self.geneDict = parseSeqFile(self.GeneFilePath.text(),self.geneDict)
                for gene in genesToAdd:
                    geneName = gene
                    geneItem = QListWidgetItem()
                    geneItem.setText(geneName)
                    geneItem.setFlags(Qt.ItemIsEditable|Qt.ItemIsSelectable | Qt.ItemIsEnabled)
                    self.geneList.addItem(geneItem)
                self.clearGeneFilePath()
            else:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setText('File Type Not Recognized')
                msg.setStandardButtons(QMessageBox.Ok)
                msg.buttonClicked.connect(self.clearGeneFilePath)
                msg.exec()
    def openHmm(self):
        fileName, _ = QFileDialog.getOpenFileName(self)
        if self.verbose: print(fileName)
        if fileName:
            self.hmmFilePath.setText(fileName)
    def loadHMMFile(self):
        if self.hmmFilePath.text():
            hmmsToAdd,self.hmmDict = parseHMMfile(self.hmmFilePath.text(),self.hmmDict)
            for hmm in hmmsToAdd:
                self.hmmList.addItem(hmm)
            self.clearHmmFilePath()
    def loadGene(self):
        genesToAdd = self.geneList.selectedItems()
        for gene in genesToAdd:
            if self.verbose: (gene.text())
            if gene.text() not in self.forBLAST.keys():
                currentRowCount = self.searchList.rowCount()
                checkbox = QTableWidgetItem()
                checkbox.setFlags(Qt.ItemIsUserCheckable|Qt.ItemIsEditable|Qt.ItemIsEnabled)
                checkbox.setCheckState(2)
                self.searchList.insertRow(currentRowCount)
                self.searchList.setItem(currentRowCount,0,QTableWidgetItem('GENE'))
                self.searchList.setItem(currentRowCount,1,QTableWidgetItem(gene.text()))
                self.searchList.setItem(currentRowCount,2,checkbox)
                self.searchList.resizeColumnsToContents()
                self.forBLAST[gene.text()] = self.geneDict[gene.text()]
                self.updateSpinBox()
    def loadHMM(self):
        # hmmDict -> {(hmmset):set(paths to hmms)}
        hmmsToAdd = tuple(sorted(hmmObject.text() for hmmObject in self.hmmList.selectedItems()))
        if hmmsToAdd:
            if self.searchList.rowCount() == 0:
                self.hitsNeededSpinBox.setRange(1, 1)
                self.hitsNeededSpinBox.setEnabled(True)
            self.forHmmer[hmmsToAdd] = set(self.hmmDict[hmm] for hmm in hmmsToAdd)
            currentRowCount = self.searchList.rowCount()
            checkbox = QTableWidgetItem()
            checkbox.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEditable | Qt.ItemIsEnabled)
            checkbox.setCheckState(2)
            self.searchList.insertRow(currentRowCount)
            self.searchList.setItem(currentRowCount, 0, QTableWidgetItem('HMMER Hits'))
            self.searchList.setItem(currentRowCount, 1, QTableWidgetItem(' and '.join(hmmsToAdd)))
            self.searchList.setItem(currentRowCount, 2, checkbox)
            self.searchList.resizeColumnsToContents()
            self.updateSpinBox()
            if self.verbose:
                print(self.forHmmer)
            self.hmmList.clearSelection()
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText('No HMMs Selected')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
    def clearGeneFilePath(self):
        self.GeneFilePath.setText('')
    def clearHmmFilePath(self):
        self.hmmFilePath.setText('')
    def clearGene(self):
        self.geneList.clear()
        self.geneDict = dict()
    def clearHMM(self):
        self.hmmList.clear()
        self.hmmDict = dict()

def main():

    hmmFetchExec = 'hmmfetch'
    makeblastdbExec = 'makeblastdb'
    blastExec = 'blastp'
    hmmSearchExec = 'hmmsearch'

    app = QApplication(sys.argv)
    form = mainApp(makeblastdbExec,blastExec,hmmFetchExec,hmmSearchExec,verbose=True)
    form.show()
    app.exec_()


if __name__ == '__main__':              # if we're running file directly and not importing it
    main()