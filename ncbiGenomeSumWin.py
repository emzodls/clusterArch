# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ncbiGenomeSumWin.ui'
#
# Created by: PyQt5 UI code generator 5.7.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_ncbiGenomeSummaryWin(object):
    def setupUi(self, ncbiGenomeSummaryWin):
        ncbiGenomeSummaryWin.setObjectName("ncbiGenomeSummaryWin")
        ncbiGenomeSummaryWin.resize(381, 403)
        self.verticalLayout = QtWidgets.QVBoxLayout(ncbiGenomeSummaryWin)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtWidgets.QLabel(ncbiGenomeSummaryWin)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.dlList = QtWidgets.QListWidget(ncbiGenomeSummaryWin)
        self.dlList.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.dlList.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.dlList.setObjectName("dlList")
        self.verticalLayout.addWidget(self.dlList)
        self.label_2 = QtWidgets.QLabel(ncbiGenomeSummaryWin)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.dlFailedList = QtWidgets.QListWidget(ncbiGenomeSummaryWin)
        self.dlFailedList.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.dlFailedList.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.dlFailedList.setObjectName("dlFailedList")
        self.verticalLayout.addWidget(self.dlFailedList)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.closeBtn = QtWidgets.QPushButton(ncbiGenomeSummaryWin)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.closeBtn.sizePolicy().hasHeightForWidth())
        self.closeBtn.setSizePolicy(sizePolicy)
        self.closeBtn.setObjectName("closeBtn")
        self.horizontalLayout.addWidget(self.closeBtn)
        self.saveSummaryBtn = QtWidgets.QPushButton(ncbiGenomeSummaryWin)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.saveSummaryBtn.sizePolicy().hasHeightForWidth())
        self.saveSummaryBtn.setSizePolicy(sizePolicy)
        self.saveSummaryBtn.setObjectName("saveSummaryBtn")
        self.horizontalLayout.addWidget(self.saveSummaryBtn)
        self.verticalLayout.addLayout(self.horizontalLayout)

        self.retranslateUi(ncbiGenomeSummaryWin)
        self.closeBtn.clicked.connect(ncbiGenomeSummaryWin.close)
        QtCore.QMetaObject.connectSlotsByName(ncbiGenomeSummaryWin)

    def retranslateUi(self, ncbiGenomeSummaryWin):
        _translate = QtCore.QCoreApplication.translate
        ncbiGenomeSummaryWin.setWindowTitle(_translate("ncbiGenomeSummaryWin", "NCBI Genome Download Summary"))
        self.label.setText(_translate("ncbiGenomeSummaryWin", "Files to Download:"))
        self.label_2.setText(_translate("ncbiGenomeSummaryWin", "Failed to Download:"))
        self.closeBtn.setText(_translate("ncbiGenomeSummaryWin", "OK"))
        self.saveSummaryBtn.setText(_translate("ncbiGenomeSummaryWin", "Save Summary"))

