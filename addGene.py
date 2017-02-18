# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'addGene.ui'
#
# Created by: PyQt5 UI code generator 5.7.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_addGeneWindow(object):
    def setupUi(self, addGeneWindow):
        addGeneWindow.setObjectName("addGeneWindow")
        addGeneWindow.resize(429, 244)
        self.verticalLayout = QtWidgets.QVBoxLayout(addGeneWindow)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(addGeneWindow)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        self.geneName = QtWidgets.QLineEdit(addGeneWindow)
        self.geneName.setObjectName("geneName")
        self.horizontalLayout.addWidget(self.geneName)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.geneSequence = QtWidgets.QTextEdit(addGeneWindow)
        self.geneSequence.setAcceptRichText(False)
        self.geneSequence.setObjectName("geneSequence")
        self.horizontalLayout_2.addWidget(self.geneSequence)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.addGeneBtn = QtWidgets.QPushButton(addGeneWindow)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.addGeneBtn.sizePolicy().hasHeightForWidth())
        self.addGeneBtn.setSizePolicy(sizePolicy)
        self.addGeneBtn.setObjectName("addGeneBtn")
        self.horizontalLayout_3.addWidget(self.addGeneBtn)
        self.closeBtn = QtWidgets.QPushButton(addGeneWindow)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.closeBtn.sizePolicy().hasHeightForWidth())
        self.closeBtn.setSizePolicy(sizePolicy)
        self.closeBtn.setObjectName("closeBtn")
        self.horizontalLayout_3.addWidget(self.closeBtn)
        self.verticalLayout.addLayout(self.horizontalLayout_3)

        self.retranslateUi(addGeneWindow)
        self.closeBtn.clicked.connect(addGeneWindow.close)
        QtCore.QMetaObject.connectSlotsByName(addGeneWindow)

    def retranslateUi(self, addGeneWindow):
        _translate = QtCore.QCoreApplication.translate
        addGeneWindow.setWindowTitle(_translate("addGeneWindow", "Add Gene"))
        self.label.setText(_translate("addGeneWindow", "Gene Name"))
        self.addGeneBtn.setText(_translate("addGeneWindow", "Add Gene Sequence"))
        self.closeBtn.setText(_translate("addGeneWindow", "Close"))

