# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'buildHmmWin.ui'
#
# Created by: PyQt5 UI code generator 5.7.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_buildHmmWin(object):
    def setupUi(self, buildHmmWin):
        buildHmmWin.setObjectName("buildHmmWin")
        buildHmmWin.resize(433, 142)
        self.verticalLayout = QtWidgets.QVBoxLayout(buildHmmWin)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(buildHmmWin)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        self.inFilePath = QtWidgets.QLineEdit(buildHmmWin)
        self.inFilePath.setAcceptDrops(False)
        self.inFilePath.setReadOnly(True)
        self.inFilePath.setObjectName("inFilePath")
        self.horizontalLayout.addWidget(self.inFilePath)
        self.selectInFileBtn = QtWidgets.QPushButton(buildHmmWin)
        self.selectInFileBtn.setObjectName("selectInFileBtn")
        self.horizontalLayout.addWidget(self.selectInFileBtn)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_2 = QtWidgets.QLabel(buildHmmWin)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout_2.addWidget(self.label_2)
        self.outFilePath = QtWidgets.QLineEdit(buildHmmWin)
        self.outFilePath.setReadOnly(True)
        self.outFilePath.setObjectName("outFilePath")
        self.horizontalLayout_2.addWidget(self.outFilePath)
        self.selectOutfileBtn = QtWidgets.QPushButton(buildHmmWin)
        self.selectOutfileBtn.setObjectName("selectOutfileBtn")
        self.horizontalLayout_2.addWidget(self.selectOutfileBtn)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.buildHmmBtn = QtWidgets.QPushButton(buildHmmWin)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.buildHmmBtn.sizePolicy().hasHeightForWidth())
        self.buildHmmBtn.setSizePolicy(sizePolicy)
        self.buildHmmBtn.setObjectName("buildHmmBtn")
        self.horizontalLayout_3.addWidget(self.buildHmmBtn)
        self.verticalLayout.addLayout(self.horizontalLayout_3)

        self.retranslateUi(buildHmmWin)
        QtCore.QMetaObject.connectSlotsByName(buildHmmWin)

    def retranslateUi(self, buildHmmWin):
        _translate = QtCore.QCoreApplication.translate
        buildHmmWin.setWindowTitle(_translate("buildHmmWin", "Build HMM from Multiple Sequence Alignment"))
        self.label.setText(_translate("buildHmmWin", "Input File:"))
        self.selectInFileBtn.setText(_translate("buildHmmWin", "Select File"))
        self.label_2.setText(_translate("buildHmmWin", "Output File: "))
        self.selectOutfileBtn.setText(_translate("buildHmmWin", "Save As"))
        self.buildHmmBtn.setText(_translate("buildHmmWin", "Create HMM"))

