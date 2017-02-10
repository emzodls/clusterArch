# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'createDbStatus.ui'
#
# Created by: PyQt5 UI code generator 5.7.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_createDbWin(object):
    def setupUi(self, createDbWin):
        createDbWin.setObjectName("createDbWin")
        createDbWin.resize(316, 134)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(createDbWin.sizePolicy().hasHeightForWidth())
        createDbWin.setSizePolicy(sizePolicy)
        self.verticalLayout = QtWidgets.QVBoxLayout(createDbWin)
        self.verticalLayout.setObjectName("verticalLayout")
        self.title = QtWidgets.QLabel(createDbWin)
        self.title.setObjectName("title")
        self.verticalLayout.addWidget(self.title)
        self.currentTask = QtWidgets.QLabel(createDbWin)
        self.currentTask.setText("")
        self.currentTask.setAlignment(QtCore.Qt.AlignCenter)
        self.currentTask.setObjectName("currentTask")
        self.verticalLayout.addWidget(self.currentTask)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.progressBar = QtWidgets.QProgressBar(createDbWin)
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName("progressBar")
        self.horizontalLayout.addWidget(self.progressBar)
        self.percentLabel = QtWidgets.QLabel(createDbWin)
        self.percentLabel.setTextInteractionFlags(QtCore.Qt.NoTextInteraction)
        self.percentLabel.setObjectName("percentLabel")
        self.horizontalLayout.addWidget(self.percentLabel, 0, QtCore.Qt.AlignRight)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.cancelBtn = QtWidgets.QPushButton(createDbWin)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.cancelBtn.sizePolicy().hasHeightForWidth())
        self.cancelBtn.setSizePolicy(sizePolicy)
        self.cancelBtn.setObjectName("cancelBtn")
        self.horizontalLayout_2.addWidget(self.cancelBtn)
        self.doneBtn = QtWidgets.QPushButton(createDbWin)
        self.doneBtn.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.doneBtn.sizePolicy().hasHeightForWidth())
        self.doneBtn.setSizePolicy(sizePolicy)
        self.doneBtn.setObjectName("doneBtn")
        self.horizontalLayout_2.addWidget(self.doneBtn)
        self.verticalLayout.addLayout(self.horizontalLayout_2)

        self.retranslateUi(createDbWin)
        QtCore.QMetaObject.connectSlotsByName(createDbWin)

    def retranslateUi(self, createDbWin):
        _translate = QtCore.QCoreApplication.translate
        createDbWin.setWindowTitle(_translate("createDbWin", "Creating Database"))
        self.title.setText(_translate("createDbWin", "Currently Processing:"))
        self.percentLabel.setText(_translate("createDbWin", "0%"))
        self.cancelBtn.setText(_translate("createDbWin", "Cancel "))
        self.doneBtn.setText(_translate("createDbWin", "OK"))

