# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gbDivSumWin.ui'
#
# Created by: PyQt5 UI code generator 5.7.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_DatabaseSummaryWin(object):
    def setupUi(self, DatabaseSummaryWin):
        DatabaseSummaryWin.setObjectName("DatabaseSummaryWin")
        DatabaseSummaryWin.resize(381, 403)
        self.verticalLayout = QtWidgets.QVBoxLayout(DatabaseSummaryWin)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtWidgets.QLabel(DatabaseSummaryWin)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.dlList = QtWidgets.QListWidget(DatabaseSummaryWin)
        self.dlList.setObjectName("dlList")
        self.verticalLayout.addWidget(self.dlList)
        self.label_2 = QtWidgets.QLabel(DatabaseSummaryWin)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.dlFailedList = QtWidgets.QListWidget(DatabaseSummaryWin)
        self.dlFailedList.setObjectName("dlFailedList")
        self.verticalLayout.addWidget(self.dlFailedList)
        self.label_3 = QtWidgets.QLabel(DatabaseSummaryWin)
        self.label_3.setObjectName("label_3")
        self.verticalLayout.addWidget(self.label_3)
        self.processFailedList = QtWidgets.QListWidget(DatabaseSummaryWin)
        self.processFailedList.setObjectName("processFailedList")
        self.verticalLayout.addWidget(self.processFailedList)
        self.closeBtn = QtWidgets.QPushButton(DatabaseSummaryWin)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.closeBtn.sizePolicy().hasHeightForWidth())
        self.closeBtn.setSizePolicy(sizePolicy)
        self.closeBtn.setObjectName("closeBtn")
        self.verticalLayout.addWidget(self.closeBtn, 0, QtCore.Qt.AlignHCenter)

        self.retranslateUi(DatabaseSummaryWin)
        self.closeBtn.clicked.connect(DatabaseSummaryWin.close)
        QtCore.QMetaObject.connectSlotsByName(DatabaseSummaryWin)

    def retranslateUi(self, DatabaseSummaryWin):
        _translate = QtCore.QCoreApplication.translate
        DatabaseSummaryWin.setWindowTitle(_translate("DatabaseSummaryWin", "Genbank Division Download Summary"))
        self.label.setText(_translate("DatabaseSummaryWin", "Files to Download:"))
        self.label_2.setText(_translate("DatabaseSummaryWin", "Failed to Download:"))
        self.label_3.setText(_translate("DatabaseSummaryWin", "Failed to Read:"))
        self.closeBtn.setText(_translate("DatabaseSummaryWin", "OK"))

