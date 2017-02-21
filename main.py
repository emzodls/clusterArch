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
    runHmmsearch,processSearchList,proccessGbks,processSearchListOptionalHits
from collections import Counter
from glob import iglob
from itertools import chain
from ncbiUtils import ncbiQuery,ncbiSummary

import mainGuiNCBI,resultsWindow,status,parameters,sys,createDbStatus,addGene

class ncbiQueryWorker(QObject):
    start = pyqtSignal()
    connectionErr = pyqtSignal()
    result = pyqtSignal(int,dict)
    def __init__(self,keyword,organism,accession,minLength,maxLength,retmax):
        self.keyword = keyword
        self.organism = organism
        self.accession = accession
        self.minLength = minLength
        self.maxLength = maxLength
        self.retmax = retmax
        super(ncbiQueryWorker, self).__init__()

    @pyqtSlot()
    @pyqtSlot(int,dict)
    def run(self):
        try:
            numHits, ncbiIDsList = ncbiQuery(self.keyword,self.organism,self.accession,
                                             self.minLength,self.maxLength,self.retmax)
            ncbiDict = ncbiSummary(ncbiIDsList, chunkSize=250)
        except:
            self.connectionErr.emit()
        self.result.emit(numHits,ncbiDict)


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
        ### Check for Blastp formatted Database
        if any(True for _ in phrCheck) and any(True for _ in pinCheck) and any(True for _ in psqCheck):
            self.makeDB.emit(False)
            makeDBintl = False
        else:
            self.makeDB.emit(True)
            makeDBintl = True
            path, outputDBname = os.path.split(self.pathToDatabase)
            out, err, retcode = MakeBlastDB(self.makeblastdbExec, self.pathToDatabase, self.outputDir,outputDBname)
            if retcode != 0:
                self.dbCreated.emit(False)
                self.terminate()
                return
            else:
                self.dbCreated.emit(True)
        ### Run Blastp
        inputFastas = os.path.join(self.outputDir, 'gene_queries.fa')
        if makeDBintl:
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

        ##### Set up NCBI Interface
        self.ncbiDict = dict()

        posIntValidator = QtGui.QIntValidator(0,1000000000)

        maxEntriesValidator = QtGui.QIntValidator(1,1000)

        self.minLength.setValidator(posIntValidator)
        self.maxLength.setValidator(posIntValidator)
        self.numEntries.setValidator(maxEntriesValidator)

        self.minLength.setText('0')
        self.maxLength.setText('1000000000')
        self.numEntries.setText('1000')

        self.searchNcbiBtn.clicked.connect(self.ncbiQuery)

        ### Set Up Parameters and Shared Data Structures
        self.nameToParamDict = {'blastEval': 1e-5, 'hmmEval': 1e-5, 'hmmScore': 30,
                                'hmmDomLen': 15, 'windowSize': 50000}

        self.geneDict = dict()
        self.hmmDict = dict()
        self.forBLAST = dict()
        self.forHmmer = dict()
        self.pathToDatabase = ''
        self.outputDir  = ''

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
    ##### NCBI Query Functions #########
    def ncbiQuery(self):
        keyword = self.ncbiKeyword.text()
        organism = self.ncbiOrganism.text()
        accession = self.ncbiAccession.text()
        minLength = self.minLength.text()
        maxLength = self.maxLength.text()
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
            self.ncbiQueryWorker.result.connect(self.updateSearchResults)

    def connectionError(self):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Problems Connecting to NCBI Database")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec()

    def updateSearchResults(self,numHits,ncbiDict):
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

        self.ncbiDict = ncbiDict

        for gi,entry in ncbiDict.items():
            self.ncbiFileSearchResults.addItem("{}:{}".format(entry[1],entry[2]))


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
        dirName = QFileDialog.getExistingDirectory()
        if dirName and os.access(dirName, os.W_OK):
            self.dbOutputDir.setText(dirName)
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
            self.createDbStatusWin.cancelBtn.clicked.connect(lambda: self.abortSearch(self.createDbStatusWin,self.runnerThread))
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
        self.createDbStatusWin.percentLabel.setText(str("%.0f%%" % (100*self.createDbStatusWin.progressBar.value()
                                                                    /int(self.totalFiles.text()))))
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
                    rowData = [self.resultsWin.resultsList.horizontalHeaderItem(idx).text()
                               for idx in range(self.resultsWin.resultsList.columnCount())]
                    writer.writerow(rowData)
                    for rowIdx in range(self.resultsWin.resultsList.rowCount()):
                        rowData = [self.resultsWin.resultsList.item(rowIdx, colIdx).text()
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