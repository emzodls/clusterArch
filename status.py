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
        runSearch.resize(342, 155)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(runSearch.sizePolicy().hasHeightForWidth())
        runSearch.setSizePolicy(sizePolicy)
        self.verticalLayout = QtWidgets.QVBoxLayout(runSearch)
        self.verticalLayout.setObjectName("verticalLayout")
        self.currentTask = QtWidgets.QLabel(runSearch)
        self.currentTask.setAlignment(QtCore.Qt.AlignCenter)
        self.currentTask.setWordWrap(True)
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
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.cancelBtn = QtWidgets.QPushButton(runSearch)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.cancelBtn.sizePolicy().hasHeightForWidth())
        self.cancelBtn.setSizePolicy(sizePolicy)
        self.cancelBtn.setObjectName("cancelBtn")
        self.horizontalLayout_2.addWidget(self.cancelBtn)
        self.viewResultsBtn = QtWidgets.QPushButton(runSearch)
        self.viewResultsBtn.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.viewResultsBtn.sizePolicy().hasHeightForWidth())
        self.viewResultsBtn.setSizePolicy(sizePolicy)
        self.viewResultsBtn.setObjectName("viewResultsBtn")
        self.horizontalLayout_2.addWidget(self.viewResultsBtn)
        self.verticalLayout.addLayout(self.horizontalLayout_2)

        self.retranslateUi(runSearch)
        QtCore.QMetaObject.connectSlotsByName(runSearch)

    def retranslateUi(self, runSearch):
        _translate = QtCore.QCoreApplication.translate
        runSearch.setWindowTitle(_translate("runSearch", "Search Status"))
        self.currentTask.setText(_translate("runSearch", "Running Checks"))
        self.percentCmpLabel.setText(_translate("runSearch", "1/6"))
        self.cancelBtn.setText(_translate("runSearch", "Cancel "))
        self.viewResultsBtn.setText(_translate("runSearch", "View Results"))

