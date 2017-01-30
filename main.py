import os
from PyQt5.QtCore import (QFile, QFileInfo, QPoint, QRect, QSettings, QSize,
        Qt, QTextStream,QCoreApplication)
from PyQt5 import QtGui
from PyQt5.QtWidgets import (QAction, QApplication, QFileDialog, QMainWindow,
        QMessageBox, QTextEdit,QTableWidgetItem)
from utils import parseSeqFile,parseHMMfile
from collections import Counter

# Handle back to sequence of gene {(filepath,geneName):sequence}


import myGui_Beta,sys # This file holds our MainWindow and all design related things
              # it also keeps events etc that we defined in Qt Designer

class ExampleApp(QMainWindow, myGui_Beta.Ui_clusterArch):
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
        global outputDir,hmmDict,geneDict,hmmDict,pathToDatabase

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

        ## Get Tasklist from searchList

        ## See if BLAST DB exists, make one if not

        ## Generate Fasta query file

        ## Generate HMM Query FIle

        ## Blast

        ## Hmmer

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

    global geneDict,hmmDict,forBLAST,forHmmer,pathToDatabase,blastExecutable,hmmerExecutable,verbose,outputDir

    verbose = True

    geneDict = dict()
    hmmDict = dict()
    forBLAST = dict()
    forHmmer = dict()
    pathToDatabase = ''
    outputDir = ''

    app = QApplication(sys.argv)  # A new instance of QApplication
    form = ExampleApp()                 # We set the form to be our ExampleApp (design)
    form.show()                         # Show the form
    app.exec_()                         # and execute the app


if __name__ == '__main__':              # if we're running file directly and not importing it
    main()