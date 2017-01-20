from PyQt5.QtCore import (QFile, QFileInfo, QPoint, QRect, QSettings, QSize,
        Qt, QTextStream,QCoreApplication)
from PyQt5 import QtGui
from PyQt5.QtWidgets import (QAction, QApplication, QFileDialog, QMainWindow,
        QMessageBox, QTextEdit)

from Bio import SeqIO
# Handle back to sequence of gene {(filepath,geneName):sequence}
geneDict = dict()


import myGui_Beta,sys # This file holds our MainWindow and all design related things
              # it also keeps events etc that we defined in Qt Designer

class ExampleApp(QMainWindow, myGui_Beta.Ui_MainWindow):
    def __init__(self):
        # Explaining super is out of the scope of this article
        # So please google it if you're not familar with it
        # Simple reason why we use it here is that it allows us to
        # access variables, methods etc in the design.py file
        super(self.__class__, self).__init__()
        self.setupUi(self)  # This is defined in design.py file automatically
                            # It sets up layout and widgets that are defined

        self.addDNAfilePathBtn.clicked.connect(self.openGene)
        self.addDNAfileBtn.clicked.connect(self.loadGene)
        self.clearGeneBtn.clicked.connect(self.clearGene)



        self.closeBtn.clicked.connect(QCoreApplication.instance().quit)
    def openGene(self):
        fileName, _ = QFileDialog.getOpenFileName(self)
        print(fileName)
        if fileName:
            self.GeneFilePath.setText(fileName)
    def loadGene(self):
        if self.GeneFilePath.text():
            extension = self.GeneFilePath.text().split('.')[-1]
            # Bring through different parsing workflows based on extension

            # genbank
            if extension == 'gbk':
                genbank_entries = SeqIO.parse(open(self.GeneFilePath.text()), "genbank")
                genesToAdd = []
                cds_ctr = 0
                for genbank_entry in genbank_entries:
                    CDS_list = (feature for feature in genbank_entry.features if feature.type == 'CDS')
                    species_id = genbank_entry.id
                    for CDS in CDS_list:
                        cds_ctr += 1
                        direction = CDS.location.strand
                        # Ensure that you don't get negative values, Biopython parser will not ignore slices that are greater
                        # than the entry so you don't need to worry about the other direction
                        internal_id = "%s_CDS_%.5i" % (species_id, cds_ctr)
                        protein_id = internal_id
                        genbank_seq = CDS.location.extract(genbank_entry)

                        # Try to find a common name for the promoter, otherwise just use the internal ID
                        if 'protein_id' in CDS.qualifiers.keys():
                            protein_id = CDS.qualifiers['protein_id'][0]
                        else:
                            for feature in genbank_seq.features:
                                if 'locus_tag' in feature.qualifiers:
                                    protein_id = feature.qualifiers['locus_tag'][0]
                        genesToAdd.append(protein_id)
                        geneDict[(self.GeneFilePath.text(),protein_id)] = str(genbank_seq)
                for gene in genesToAdd:
                    self.geneList.addItem(gene)
            # fasta
            elif extension == 'fasta' or extension == 'fa':
                genes = SeqIO.parse(open(self.GeneFilePath.text()), "fasta")
                genesToAdd = []
                for gene in genes:
                    genesToAdd.append(gene.id)
                    geneDict[(self.GeneFilePath.text(),gene.id)] = str(gene.seq)
                for gene in genesToAdd:
                    self.geneList.addItem(gene)
            else:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setText('File Type Not Recognized')
                msg.setStandardButtons(QMessageBox.Ok)
                msg.buttonClicked.connect(self.clearFilePath)
                msg.exec()
    def clearFilePath(self):
        self.GeneFilePath.setText('')
    def clearGene(self):
        self.geneList.clear()
        geneDict = dict()


def main():
    global geneDict
    app = QApplication(sys.argv)  # A new instance of QApplication
    form = ExampleApp()                 # We set the form to be our ExampleApp (design)
    form.show()                         # Show the form
    app.exec_()                         # and execute the app


if __name__ == '__main__':              # if we're running file directly and not importing it
    main()