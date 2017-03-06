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
        self.dlList.setObjectName("dlList")
        self.verticalLayout.addWidget(self.dlList)
        self.label_2 = QtWidgets.QLabel(ncbiGenomeSummaryWin)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.dlFailedList = QtWidgets.QListWidget(ncbiGenomeSummaryWin)
        self.dlFailedList.setObjectName("dlFailedList")
        self.verticalLayout.addWidget(self.dlFailedList)
        self.closeBtn = QtWidgets.QPushButton(ncbiGenomeSummaryWin)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.closeBtn.sizePolicy().hasHeightForWidth())
        self.closeBtn.setSizePolicy(sizePolicy)
        self.closeBtn.setObjectName("closeBtn")
        self.verticalLayout.addWidget(self.closeBtn, 0, QtCore.Qt.AlignHCenter)

        self.retranslateUi(ncbiGenomeSummaryWin)
        self.closeBtn.clicked.connect(ncbiGenomeSummaryWin.close)
        QtCore.QMetaObject.connectSlotsByName(ncbiGenomeSummaryWin)

    def retranslateUi(self, ncbiGenomeSummaryWin):
        _translate = QtCore.QCoreApplication.translate
        ncbiGenomeSummaryWin.setWindowTitle(_translate("ncbiGenomeSummaryWin", "NCBI Genome Download Summary"))
        self.label.setText(_translate("ncbiGenomeSummaryWin", "Files to Download:"))
        self.label_2.setText(_translate("ncbiGenomeSummaryWin", "Failed to Download:"))
        self.closeBtn.setText(_translate("ncbiGenomeSummaryWin", "OK"))

