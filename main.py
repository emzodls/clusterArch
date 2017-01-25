from PyQt5.QtCore import (QFile, QFileInfo, QPoint, QRect, QSettings, QSize,
        Qt, QTextStream,QCoreApplication)
from PyQt5 import QtGui
from PyQt5.QtWidgets import (QAction, QApplication, QFileDialog, QMainWindow,
        QMessageBox, QTextEdit,QTableWidgetItem)
from utils import parseSeqFile,parseHMMfile

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

        self.clearGeneBtn.clicked.connect(self.clearGene)
        self.clearHMMbtn.clicked.connect(self.clearHMM)

        self.closeBtn.clicked.connect(QCoreApplication.instance().quit)

    def openGene(self):
        fileName, _ = QFileDialog.getOpenFileName(self)
        print(fileName)
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
        fileName, _ = QFileDialog.getOpenFileName(self)
        print(fileName)
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
        global forBLAST,geneDict
        genesToAdd = self.geneList.selectedItems()
        for gene in genesToAdd:
            print(gene.text())
            if gene.text() not in forBLAST.keys():
                currentRowCount = self.searchList.rowCount()
                self.searchList.insertRow(currentRowCount)
                self.searchList.setItem(currentRowCount,0,QTableWidgetItem('GENE'))
                self.searchList.setItem(currentRowCount,1,QTableWidgetItem(gene.text()))
                forBLAST[gene.text()] = geneDict[gene.text()]

    def loadHMM(self):
        # hmmDict -> {hmmFile:hmmsToConsider}
        global hmmDict,forHmmer
        hmmsToAdd = self.hmmList.selectedItems()
        hmmList = []
        for hmmObject in hmmsToAdd:
            hmmName = hmmObject.text()
            hmmList.append(hmmName)
            hmmSourceFile = hmmDict[hmmName]
            hmmsToConsider = forHmmer.get(hmmSourceFile,set())
            hmmsToConsider.add(hmmName)
            forHmmer[hmmSourceFile] = hmmsToConsider
        currentRowCount = self.searchList.rowCount()
        self.searchList.insertRow(currentRowCount)
        self.searchList.setItem(currentRowCount, 0, QTableWidgetItem('HMMER Hits'))
        self.searchList.setItem(currentRowCount, 1, QTableWidgetItem(' and '.join(hmmList)))
        self.hmmList.clearSelection()


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
    global geneDict,hmmDict,forBLAST,forHmmer
    geneDict = dict()
    hmmDict = dict()
    forBLAST = dict()
    forHmmer = dict()
    app = QApplication(sys.argv)  # A new instance of QApplication
    form = ExampleApp()                 # We set the form to be our ExampleApp (design)
    form.show()                         # Show the form
    app.exec_()                         # and execute the app


if __name__ == '__main__':              # if we're running file directly and not importing it
    main()