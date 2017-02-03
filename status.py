# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'status.ui'
#
# Created by: PyQt5 UI code generator 5.7.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_runSearch(object):
    def setupUi(self, runSearch):
        runSearch.setObjectName("runSearch")
        runSearch.resize(367, 139)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(runSearch.sizePolicy().hasHeightForWidth())
        runSearch.setSizePolicy(sizePolicy)
        self.verticalLayout = QtWidgets.QVBoxLayout(runSearch)
        self.verticalLayout.setObjectName("verticalLayout")
        self.title = QtWidgets.QLabel(runSearch)
        self.title.setObjectName("title")
        self.verticalLayout.addWidget(self.title)
        self.currentTask = QtWidgets.QLabel(runSearch)
        self.currentTask.setAlignment(QtCore.Qt.AlignCenter)
        self.currentTask.setObjectName("currentTask")
        self.verticalLayout.addWidget(self.currentTask)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.percentCmpBar = QtWidgets.QProgressBar(runSearch)
        self.percentCmpBar.setProperty("value", 0)
        self.percentCmpBar.setObjectName("percentCmpBar")
        self.horizontalLayout.addWidget(self.percentCmpBar)
        self.percentCmpLabel = QtWidgets.QLabel(runSearch)
        self.percentCmpLabel.setObjectName("percentCmpLabel")
        self.horizontalLayout.addWidget(self.percentCmpLabel)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.viewResultsBtn = QtWidgets.QPushButton(runSearch)
        self.viewResultsBtn.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.viewResultsBtn.sizePolicy().hasHeightForWidth())
        self.viewResultsBtn.setSizePolicy(sizePolicy)
        self.viewResultsBtn.setObjectName("viewResultsBtn")
        self.verticalLayout.addWidget(self.viewResultsBtn, 0, QtCore.Qt.AlignHCenter)

        self.retranslateUi(runSearch)
        QtCore.QMetaObject.connectSlotsByName(runSearch)

    def retranslateUi(self, runSearch):
        _translate = QtCore.QCoreApplication.translate
        runSearch.setWindowTitle(_translate("runSearch", "Search Status"))
        self.title.setText(_translate("runSearch", "Running Architecture Search:"))
        self.currentTask.setText(_translate("runSearch", "Building BLAST Database"))
        self.percentCmpLabel.setText(_translate("runSearch", "0%"))
        self.viewResultsBtn.setText(_translate("runSearch", "View Results"))

