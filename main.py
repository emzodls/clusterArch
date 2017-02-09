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

import os
from PyQt5.QtCore import (QFile, QFileInfo, QPoint, QRect, QSettings, QSize,QRegularExpression,
        Qt, QTextStream,QCoreApplication,QThread,pyqtSignal,pyqtSlot,QObject)
from PyQt5 import QtGui
from PyQt5.QtWidgets import (QAction, QApplication, QFileDialog, QMainWindow,
        QMessageBox, QTextEdit,QTableWidgetItem,QWidget,QSizePolicy)
from utils import parseSeqFile,parseHMMfile,generateInputFasta,generateHMMdb,MakeBlastDB,runBLAST,\
    runHmmsearch,processSearchList,proccessGbks
from collections import Counter
from glob import iglob
from itertools import chain

# Handle back to sequence of gene {(filepath,geneName):sequence}


import myGui_Beta,resultsWindow,status,parameters,sys # This file holds our MainWindow and all design related things
              # it also keeps events etc that we defined in Qt Designer

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
        proccessGbks(self.taskList,self.dbPath,self.currentGbk)
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

class statusWindow(QWidget,status.Ui_runSearch):
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

class mainApp(QMainWindow, myGui_Beta.Ui_clusterArch):
    def __init__(self):
        # Explaining super is out of the scope of this article
        # So please google it if you're not familar with it
        # Simple reason why we use it here is that it allows us to
        # access variables, methods etc in the design.py file
        super(self.__class__, self).__init__()
        self.setupUi(self)  # This is defined in design.py file automatically
                            # It sets up layout and widgets that are defined
        ### Process Genbank Panel
        self.genbankDirBtn.clicked.connect(self.loadGbkDir)
        self.addGenbankBtn.clicked.connect(self.addGbkFiles)
        self.deleteGenbankBtn.clicked.connect(self.removeGbkFiles)
        self.chooseDBoutputDirBtn.clicked.connect(self.chooseDBoutDir)
        self.createDbBtn.clicked.connect(self.createDB)

        self.addGenefilePathBtn.clicked.connect(self.openGene)
        self.addHMMfilePathBtn.clicked.connect(self.openHmm)

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
        self.checkSuccessful = False
        self.blastSuccessful = False
        self.hmmerSuccessful = False
        ### Initialize settings for parameters window
        self.paramsWin = paramsWindow()

        self.paramsWin.blastEval.textChanged.connect(self.activateApply)
        self.paramsWin.hmmEval.textChanged.connect(self.activateApply)
        self.paramsWin.hmmDomLen.textChanged.connect(self.activateApply)
        self.paramsWin.hmmScore.textChanged.connect(self.activateApply)
        self.paramsWin.windowSize.textChanged.connect(self.activateApply)
        self.paramsWin.applyBtn.clicked.connect(self.updateParams)
        self.setParamsBtn.clicked.connect(self.showParamWindow)

        self.nameToParamDict = {'blastEval':1e-5,'hmmEval':1e-5,'hmmScore':30,
                                'hmmDomLen':15,'windowSize':50000}
    ### Genbank Panel Functions
    def loadGbkDir(self):
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
            self.genbankDirPath.clear()

    def removeGbkFiles(self):
        gbksToRemove = self.genbankList.selectedItems()
        for gbkFile in gbksToRemove:
            self.genbankList.takeItem(self.genbankList.row(gbkFile))

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
        taskList = [genbankFile.text() for genbankFile in self.genbankList.items()]

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

            self.createDbWorker = createDbWorker(taskList,dbPath)
            self.createDbWorker.moveToThread(self.runnerThread)
            self.createDbWorker.start.connect(self.createDbWorker.run)
            self.createDbWorker.start.emit()

            self.createDbWorker.finished.connect(self.createDbSuccessful)

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
    ### Status Window Methods
    def updateStatusWinText(self,statusWin,currentTaskText,percentCmpLabelText,percentCmpBarValue):
        statusWin.currentTask.setText(currentTaskText)
        statusWin.percentCmpLabel.setText(percentCmpLabelText)
        statusWin.percentCmpBar.setValue(percentCmpBarValue)
    def abortSearch(self,statusWin,currentThread):
        statusWin.close()
        currentThread.terminate()
    ### Run Permissions Checks
    def runChecks(self):
        if self.searchList.rowCount() >= 1:
            if not self.checkSuccessful:
                self.runnerThread = QThread()
                self.runnerThread.start()

                self.checkSuccessful = False
                self.blastSuccessful = False
                self.hmmerSuccessful = False

                self.checksWorker = runChecksWorker(outputDir,pathToDatabase)
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
            self.runnerThread.terminate()
        elif checkErr == 'noDBpath':
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('Please Specify A Database to Query')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
            statusWin.close()
            self.runnerThread.terminate()
        elif checkErr == 'noWriteOutputDir':
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('Cannot Write to Output Folder Specified')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
            statusWin.close()
            self.runnerThread.terminate()
        elif checkErr == 'noReadDB':
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('Cannot Read Database')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
            statusWin.close()
            self.runnerThread.terminate()
        elif checkErr == 'invalidDB':
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('No Valid Entries in Database')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
            statusWin.close()
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
        global forBLAST, outputDir, pathToDatabase, makeblastdbExec, blastExec
        self.updateStatusWinText(self.statusWin, 'Running BLAST Search: Generating Input Fasta', '2/6',
                                 self.statusWin.percentCmpBar.value() + 1)
        if self.checkSuccessful:
            if not self.blastSuccessful:
                if forBLAST:
                    self.blastWorker = runBlastWorker(forBLAST, outputDir, pathToDatabase,
                                                      makeblastdbExec, blastExec,self.nameToParamDict['blastEval'])
                    # self.currentThread = self.blastThread
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
            self.runnerThread.terminate()
            return
    ### Run Hmmer Methods
    def runHmmer(self):
        global forHmmer, hmmFetchExec, hmmSearchExec, hmmDict, outputDir,pathToDatabase
        if self.checkSuccessful:
            self.updateStatusWinText(self.statusWin, 'Running HMMSearch: Extracting HMMs to make reference file',
                                     '3/6', self.statusWin.percentCmpBar.value() + 1)

            if not self.hmmerSuccessful:
                hmmSet = set()
                hmmSet.update(*forHmmer.keys())
                if hmmSet:
                    self.hmmerWorker = runHmmerWorker(hmmFetchExec, hmmSearchExec,
                                                      hmmDict, hmmSet, outputDir,pathToDatabase,self.nameToParamDict['hmmEval'])
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
                    self.processResults()
            else:
                self.processResults()
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
            self.runnerThread.terminate()
            return
    @pyqtSlot(bool)
    def hmmSuccessCheck(self,flag,statusWin):
        if flag:
            self.hmmerSuccessful = True
            self.processResults()
            return
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('Hmmsearch Failed')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
            statusWin.close()
            self.runnerThread.terminate()
    ### Process Filtered Clusters
    def processResults(self):
        global outputDir,verbose
        self.updateStatusWinText(self.statusWin, 'Parsing Output Files: Creating Search Lists', '4/6', 4)

        if self.checkSuccessful and self.blastSuccessful and self.hmmerSuccessful:
            blastList = []
            hmmList = []

            blastOutFile = outputDir + "/blast_results.out"
            hmmOutFile = outputDir + '/hmmSearch.out'

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
                        self.runnerThread.terminate()
                        return
                else:
                    passBLASTcheck = True
                if hmmList:
                    try:
                        assert os.path.isfile(outputDir + '/hmmSearch.out'), "HMMer Hits included in search but no " \
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
    @pyqtSlot()
    def showResultsWindow(self,blastList,hmmList,filteredClusters):
        self.statusWin.close()
        self.checkSuccessful = False
        self.resultsWin = QWidget()
        resultsUi = resultsWindow.Ui_Results()
        resultsUi.setupUi(self.resultsWin)
        resultsUi.resultsList.setColumnCount(3 + len(blastList) + len(hmmList))
        resultsUi.resultsList.setHorizontalHeaderItem(0, QTableWidgetItem("Species/Accession ID"))
        resultsUi.resultsList.setHorizontalHeaderItem(1, QTableWidgetItem("Start"))
        resultsUi.resultsList.setHorizontalHeaderItem(2, QTableWidgetItem("End"))
        for idx, blastHit in enumerate(blastList):
            resultsUi.resultsList.setHorizontalHeaderItem(idx + 3, QTableWidgetItem(blastHit))
        for idx, hmmHit in enumerate(hmmList):
            resultsUi.resultsList.setHorizontalHeaderItem(idx + len(blastList) + 3,
                                                          QTableWidgetItem(' and '.join(hmmHit)))
        for cluster, hitDict in filteredClusters.items():
            currentRowCount = resultsUi.resultsList.rowCount()
            resultsUi.resultsList.insertRow(currentRowCount)
            resultsUi.resultsList.setItem(currentRowCount, 0, QTableWidgetItem(cluster[0]))
            resultsUi.resultsList.setItem(currentRowCount, 1, QTableWidgetItem(str(cluster[1])))
            resultsUi.resultsList.setItem(currentRowCount, 2, QTableWidgetItem(str(cluster[2])))

            for idx, blastHit in enumerate(blastList):
                resultsUi.resultsList.setItem(currentRowCount, idx + 3,
                                              QTableWidgetItem('; '.join(filteredClusters[cluster][blastHit])))
            for idx, hmmHit in enumerate(hmmList):
                resultsUi.resultsList.setItem(currentRowCount, idx + len(blastList) + 3,
                                              QTableWidgetItem('; '.join(filteredClusters[cluster][hmmHit])))

        resultsUi.resultsList.resizeColumnsToContents()
        resultsUi.resultsList.update()

        self.resultsWin.updateGeometry()
        self.resultsWin.show()
    @pyqtSlot(dict)
    def generateResults(self,filteredClusters,statusWin,blastList,hmmList):
        if filteredClusters:
            self.updateStatusWinText(statusWin, 'Search Complete',
                                     '6/6', 6)
            statusWin.viewResultsBtn.setEnabled(True)
            statusWin.viewResultsBtn.clicked.connect(lambda: self.showResultsWindow(blastList,hmmList,filteredClusters))
        else:
            self.updateStatusWinText(statusWin, 'Search Complete: No Results Found',
                                     '6/6', 6)
            self.checkSuccessful = False
            return
    def removeSearchTerm(self):
        global forHmmer,forBLAST,hmmDict,geneDict,verbose

        itemType =  self.searchList.item(self.searchList.currentRow(),0).text()
        termToRemove = self.searchList.item(self.searchList.currentRow(),1).text()

        if itemType == 'GENE':
            if verbose: print("Deleted %s" % termToRemove)
            del forBLAST[termToRemove]
            self.searchList.removeRow(self.searchList.currentRow())
        elif itemType == 'HMMER Hits':
            hmmKey = tuple(sorted(termToRemove.split(' and ')))
            hmmerSearchTerms = Counter()
            for idx in range(self.searchList.rowCount()):
                if self.searchList.item(idx,0).text() == 'HMMER Hits':
                    hmmerSearchTerms[tuple(sorted(self.searchList.item(idx,1).text().split(' and ')))] += 1
            if hmmerSearchTerms[hmmKey] <= 1:
                del forHmmer[hmmKey]
            if verbose: print(forHmmer)
            self.searchList.removeRow(self.searchList.currentRow())
    def outdirFunction(self):
        global outputDir,verbose
        if 'Select Directory...' in self.outputDirectorySelector.currentText():
            dirName = QFileDialog.getExistingDirectory(self)
            if dirName:
                self.outputDirectorySelector.addItem(dirName)
                self.outputDirectorySelector.setCurrentIndex(self.outputDirectorySelector.count()-1)
        if 'Select Directory...' not in self.outputDirectorySelector.currentText():
            outputDir = self.outputDirectorySelector.currentText()
            if verbose: print(outputDir)
    def databaseFunction(self):
        global pathToDatabase,verbose
        if 'add database...' in self.dataBaseSelector.currentText():
            fileName, _ = QFileDialog.getOpenFileName(self)
            if fileName:
                self.dataBaseSelector.addItem(fileName)
                self.dataBaseSelector.setCurrentIndex(self.dataBaseSelector.count()-1)
        if 'add database...' not in self.dataBaseSelector.currentText():
            pathToDatabase = self.dataBaseSelector.currentText()
            if verbose: print(pathToDatabase)
    def openGene(self):
        global verbose
        fileName, _ = QFileDialog.getOpenFileName(self)
        if verbose: print(fileName)
        if fileName:
            self.GeneFilePath.setText(fileName)
    def loadGeneFile(self):
        global geneDict
        if self.GeneFilePath.text():
            valid_extensions = set(['gbk','fa','fasta'])
            extension = self.GeneFilePath.text().split('.')[-1]
            if extension in valid_extensions:
                genesToAdd,geneDict = parseSeqFile(self.GeneFilePath.text(),geneDict)
                for gene in genesToAdd:
                    self.geneList.addItem(gene)
                self.clearGeneFilePath()
            else:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setText('File Type Not Recognized')
                msg.setStandardButtons(QMessageBox.Ok)
                msg.buttonClicked.connect(self.clearGeneFilePath)
                msg.exec()
    def openHmm(self):
        global verbose
        fileName, _ = QFileDialog.getOpenFileName(self)
        if verbose: print(fileName)
        if fileName:
            self.hmmFilePath.setText(fileName)
    def loadHMMFile(self):
        global hmmDict
        if self.hmmFilePath.text():
            hmmsToAdd,hmmDict = parseHMMfile(self.hmmFilePath.text(),hmmDict)
            for hmm in hmmsToAdd:
                self.hmmList.addItem(hmm)
            self.clearHmmFilePath()
    def loadGene(self):
        global forBLAST,geneDict,verbose
        genesToAdd = self.geneList.selectedItems()
        for gene in genesToAdd:
            if verbose: (gene.text())
            if gene.text() not in forBLAST.keys():
                currentRowCount = self.searchList.rowCount()
                self.searchList.insertRow(currentRowCount)
                self.searchList.setItem(currentRowCount,0,QTableWidgetItem('GENE'))
                self.searchList.setItem(currentRowCount,1,QTableWidgetItem(gene.text()))
                self.searchList.resizeColumnsToContents()
                forBLAST[gene.text()] = geneDict[gene.text()]
    def loadHMM(self):
        # hmmDict -> {(hmmset):set(paths to hmms)}
        global hmmDict,forHmmer,verbose
        hmmsToAdd = tuple(sorted(hmmObject.text() for hmmObject in self.hmmList.selectedItems()))
        if hmmsToAdd:
            forHmmer[hmmsToAdd] = set(hmmDict[hmm] for hmm in hmmsToAdd)
            currentRowCount = self.searchList.rowCount()
            self.searchList.insertRow(currentRowCount)
            self.searchList.setItem(currentRowCount, 0, QTableWidgetItem('HMMER Hits'))
            self.searchList.setItem(currentRowCount, 1, QTableWidgetItem(' and '.join(hmmsToAdd)))
            self.searchList.resizeColumnsToContents()
            if verbose:
                print(forHmmer)
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
        global geneDict
        self.geneList.clear()
        geneDict = dict()
    def clearHMM(self):
        global hmmDict
        self.hmmList.clear()
        hmmDict = dict()

def main():

    global geneDict,hmmDict,forBLAST,forHmmer,pathToDatabase,blastExecutable,hmmerExecutable,verbose\
        ,outputDir,hmmFetchExec,makeblastdbExec,hmmSearchExec,blastExec

    verbose = True

    geneDict = dict()
    hmmDict = dict()
    forBLAST = dict()
    forHmmer = dict()

    hmmFetchExec = 'hmmfetch'
    makeblastdbExec = 'makeblastdb'
    blastExec = 'blastp'
    hmmSearchExec = 'hmmsearch'

    pathToDatabase = ''
    outputDir = ''

    app = QApplication(sys.argv)
    form = mainApp()
    form.show()
    app.exec_()


if __name__ == '__main__':              # if we're running file directly and not importing it
    main()