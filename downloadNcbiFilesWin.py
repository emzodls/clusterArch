# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'downloadNcbiFilesWin.ui'
#
# Created by: PyQt5 UI code generator 5.7.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_downloadNcbiWin(object):
    def setupUi(self, downloadNcbiWin):
        downloadNcbiWin.setObjectName("downloadNcbiWin")
        downloadNcbiWin.resize(316, 134)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(downloadNcbiWin.sizePolicy().hasHeightForWidth())
        downloadNcbiWin.setSizePolicy(sizePolicy)
        self.verticalLayout = QtWidgets.QVBoxLayout(downloadNcbiWin)
        self.verticalLayout.setObjectName("verticalLayout")
        self.title = QtWidgets.QLabel(downloadNcbiWin)
        self.title.setObjectName("title")
        self.verticalLayout.addWidget(self.title)
        self.currentTask = QtWidgets.QLabel(downloadNcbiWin)
        self.currentTask.setText("")
        self.currentTask.setAlignment(QtCore.Qt.AlignCenter)
        self.currentTask.setObjectName("currentTask")
        self.verticalLayout.addWidget(self.currentTask)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.progressBar = QtWidgets.QProgressBar(downloadNcbiWin)
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName("progressBar")
        self.horizontalLayout.addWidget(self.progressBar)
        self.progressLabel = QtWidgets.QLabel(downloadNcbiWin)
        self.progressLabel.setTextInteractionFlags(QtCore.Qt.NoTextInteraction)
        self.progressLabel.setObjectName("progressLabel")
        self.horizontalLayout.addWidget(self.progressLabel, 0, QtCore.Qt.AlignRight)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.cancelBtn = QtWidgets.QPushButton(downloadNcbiWin)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.cancelBtn.sizePolicy().hasHeightForWidth())
        self.cancelBtn.setSizePolicy(sizePolicy)
        self.cancelBtn.setObjectName("cancelBtn")
        self.horizontalLayout_2.addWidget(self.cancelBtn)
        self.doneBtn = QtWidgets.QPushButton(downloadNcbiWin)
        self.doneBtn.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.doneBtn.sizePolicy().hasHeightForWidth())
        self.doneBtn.setSizePolicy(sizePolicy)
        self.doneBtn.setObjectName("doneBtn")
        self.horizontalLayout_2.addWidget(self.doneBtn)
        self.verticalLayout.addLayout(self.horizontalLayout_2)

        self.retranslateUi(downloadNcbiWin)
        QtCore.QMetaObject.connectSlotsByName(downloadNcbiWin)

    def retranslateUi(self, downloadNcbiWin):
        _translate = QtCore.QCoreApplication.translate
        downloadNcbiWin.setWindowTitle(_translate("downloadNcbiWin", "Downloading Files From NCBI"))
        self.title.setText(_translate("downloadNcbiWin", "Currently Downloading:"))
        self.progressLabel.setText(_translate("downloadNcbiWin", "0/0"))
        self.cancelBtn.setText(_translate("downloadNcbiWin", "Cancel "))
        self.doneBtn.setText(_translate("downloadNcbiWin", "OK"))

