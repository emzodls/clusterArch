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
from PyQt5.QtCore import (QFile, QFileInfo, QPoint, QRect, QSettings, QSize,
        Qt, QTextStream,QCoreApplication)
from PyQt5 import QtGui
from PyQt5.QtWidgets import (QAction, QApplication, QFileDialog, QMainWindow,
        QMessageBox, QTextEdit,QTableWidgetItem,QWidget,QSizePolicy)
from utils import parseSeqFile,parseHMMfile,generateInputFasta,generateHMMdb,MakeBlastDB,runBLAST,\
    runHmmsearch,processSearchList
from collections import Counter
from glob import iglob

# Handle back to sequence of gene {(filepath,geneName):sequence}


import myGui_Beta,resultsWindow,sys # This file holds our MainWindow and all design related things
              # it also keeps events etc that we defined in Qt Designer

class mainApp(QMainWindow, myGui_Beta.Ui_clusterArch):
    def __init__(self):
        # Explaining super is out of the scope of this article
        # So please google it if you're not familar with it
        # Simple reason why we use it here is that it allows us to
        # access variables, methods etc in the design.py file
        super(self.__class__, self).__init__()
        self.setupUi(self)  # This is defined in design.py file automatically
                            # It sets up layout and widgets that are defined

        self.addGenefilePathBtn.clicked.connect(self.openGene)
        self.addHMMfilePathBtn.clicked.connect(self.openHmm)

        self.addGenefileBtn.clicked.connect(self.loadGeneFile)
        self.addHMMfilebtn.clicked.connect(self.loadHMMFile)

        self.addGeneBtn.clicked.connect(self.loadGene)
        self.addHMMbtn.clicked.connect(self.loadHMM)

        self.removeSearchTermBtn.clicked.connect(self.removeSearchTerm)

        self.clearGeneBtn.clicked.connect(self.clearGene)
        self.clearHMMbtn.clicked.connect(self.clearHMM)

        self.runSearch.clicked.connect(self.runSearchFunction)

        self.closeBtn.clicked.connect(QCoreApplication.instance().quit)

        self.dataBaseSelector.currentIndexChanged.connect(self.databaseFunction)
        self.outputDirectorySelector.currentIndexChanged.connect(self.outdirFunction)

    def runSearchFunction(self):
        global outputDir,hmmDict,pathToDatabase, forHmmer, forBLAST,hmmFetchExec,makeblastdbExec,blastExec,hmmSearchExec

        ## Run Checks

        if not outputDir:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('Please Specify An Output Directory')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
        elif not pathToDatabase:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('Please Specify A Database to Query')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
        elif not os.access(outputDir, os.W_OK):
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('Cannot Write to Output Folder Specified')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
        elif not os.access(pathToDatabase, os.R_OK):
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText('Cannot Read Database')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()

        ## Run BLAST
        if forBLAST:
            madeDBflag = False
            generateInputFasta(forBLAST,outputDir)
            print('Successfully Generated Input Fasta')
            baseDir, ext = os.path.splitext(pathToDatabase)
            phrCheck = iglob(baseDir + '*phr')
            pinCheck = iglob(baseDir + '*pin')
            psqCheck = iglob(baseDir + '*psq')

            if any(True for _ in phrCheck) and any(True for _ in pinCheck) and any(True for _ in psqCheck):
                pass
            else:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Information)
                msg.setText("BLAST Database not Found, creating BLAST Database in Output directory...")
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()
                madeDBflag = True
                path, outputDBname = os.path.split(pathToDatabase)
                out, err, retcode = MakeBlastDB(makeblastdbExec, pathToDatabase, outputDir,outputDBname)
                if retcode != 0:
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Critical)
                    msg.setText('No BLAST DB found, creating one failed')
                    msg.setStandardButtons(QMessageBox.Ok)
                    msg.exec()
                    return
            inputFastas = os.path.join(outputDir, 'gene_queries.fa')
            if madeDBflag:
                out, err, retcode = runBLAST(blastExec, inputFastas, outputDir, os.path.join(outputDir, outputDBname),
                                             eValue='1E-05')
            else:
                out, err, retcode = runBLAST(blastExec, inputFastas, outputDir, pathToDatabase, eValue='1E-05')
            if retcode != 0:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText('Blastp Failed')
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()
                return

        ## Run Hmmer
        hmmSet = set()
        hmmSet.update(*forHmmer.keys())
        if hmmSet:
            errFlag, failedToFetch = generateHMMdb(hmmFetchExec,hmmDict,hmmSet,outputDir)
            if errFlag:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText('Failed Finding Following HMMs: %s' % ','.join(failedToFetch))
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()
                return
            print('Successfully Generated Input HmmFile')
            hmmDBase = os.path.join(outputDir,'hmmDB.hmm')
            out, err, retcode = runHmmsearch(hmmSearchExec, hmmDBase, outputDir, pathToDatabase, eValue='1E-05')
            if retcode != 0:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText('Hmmsearch Failed')
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()
                return

        ## Create Task Search Lists from BLAST and HMMER requests and pass it off to helper function
        blastList = []
        hmmList =[]

        passBLASTcheck = False
        passHMMcheck = False

        blastOutFile = outputDir + "/blast_results.out"
        hmmOutFile = outputDir+'/hmmSearch.out'

        for idx in range(self.searchList.rowCount()):
            if self.searchList.item(idx, 0).text() == 'HMMER Hits':
                hmmList.append(tuple(sorted(self.searchList.item(idx,1).text().split(' and '))))
            elif self.searchList.item(idx, 0).text() == 'GENE':
                blastList.append(self.searchList.item(idx,1).text())
        # Check that outputfiles are successfully generated
        if blastList:
            try:
                assert os.path.isfile(blastOutFile) ,"BLAST Hits included in search but no " \
                                                                        "BLAST output file found"
                passBLASTcheck = True
            except AssertionError:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText('BLAST Hits included in search but no BLAST output file found!')
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()
                return
        else:
            passBLASTcheck = True
        if hmmList:
            try:
                assert os.path.isfile(outputDir+'/hmmSearch.out'), "HMMer Hits included in search but no " \
                                                                        "HMMer output file found"
                passHMMcheck = True
            except AssertionError:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText('HMMer Hits requested in search but no HMMer output file found!')
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()
                return
        else:
            passHMMcheck = True

        if passBLASTcheck and passHMMcheck:
            filteredClusters = processSearchList(blastList,hmmList,blastOutFile,hmmOutFile)

            if filteredClusters:
                print(filteredClusters)
                self.resultsWin = QWidget()

                # set up results to print

                resultsUi = resultsWindow.Ui_Results()
                resultsUi.setupUi(self.resultsWin)

                resultsUi.resultsList.setColumnCount(3+len(blastList)+len(hmmList))
                resultsUi.resultsList.setHorizontalHeaderItem(0,QTableWidgetItem("Species/Accession ID"))
                resultsUi.resultsList.setHorizontalHeaderItem(1, QTableWidgetItem("Start"))
                resultsUi.resultsList.setHorizontalHeaderItem(2, QTableWidgetItem("End"))

                for idx,blastHit in enumerate(blastList):
                    resultsUi.resultsList.setHorizontalHeaderItem(idx+ 3, QTableWidgetItem(blastHit))

                for idx,hmmHit in enumerate(hmmList):
                    resultsUi.resultsList.setHorizontalHeaderItem(idx+len(blastList)+3,
                                                                  QTableWidgetItem(' and '.join(hmmHit)))

                for cluster,hitDict in filteredClusters.items():
                    currentRowCount = resultsUi.resultsList.rowCount()
                    resultsUi.resultsList.insertRow(currentRowCount)
                    resultsUi.resultsList.setItem(currentRowCount, 0, QTableWidgetItem(cluster[0]))
                    resultsUi.resultsList.setItem(currentRowCount, 1, QTableWidgetItem(str(cluster[1])))
                    resultsUi.resultsList.setItem(currentRowCount, 2, QTableWidgetItem(str(cluster[2])))

                    for idx,blastHit in enumerate(blastList):
                        resultsUi.resultsList.setItem(currentRowCount, idx+3,
                                                      QTableWidgetItem('; '.join(filteredClusters[cluster][blastHit])))
                    for idx,hmmHit in enumerate(hmmList):
                        resultsUi.resultsList.setItem(currentRowCount, idx + len(blastList) + 3,
                                                      QTableWidgetItem('; '.join(filteredClusters[cluster][hmmHit])))

                resultsUi.resultsList.resizeColumnsToContents()
                resultsUi.resultsList.update()

                self.resultsWin.updateGeometry()
                self.resultsWin.show()
            else:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Information)
                msg.setText('No Hits Found.')
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()
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

    app = QApplication(sys.argv)  # A new instance of QApplication
    form = mainApp()                 # We set the form to be our ExampleApp (design)
    form.show()                         # Show the form
    app.exec_()                         # and execute the app


if __name__ == '__main__':              # if we're running file directly and not importing it
    main()