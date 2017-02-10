# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'resultsWindow.ui'
#
# Created by: PyQt5 UI code generator 5.7.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Results(object):
    def setupUi(self, Results):
        Results.setObjectName("Results")
        Results.resize(929, 335)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(1)
        sizePolicy.setVerticalStretch(1)
        sizePolicy.setHeightForWidth(Results.sizePolicy().hasHeightForWidth())
        Results.setSizePolicy(sizePolicy)
        self.verticalLayout = QtWidgets.QVBoxLayout(Results)
        self.verticalLayout.setObjectName("verticalLayout")
        self.resultsList = QtWidgets.QTableWidget(Results)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(1)
        sizePolicy.setVerticalStretch(1)
        sizePolicy.setHeightForWidth(self.resultsList.sizePolicy().hasHeightForWidth())
        self.resultsList.setSizePolicy(sizePolicy)
        self.resultsList.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.resultsList.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.resultsList.setTextElideMode(QtCore.Qt.ElideMiddle)
        self.resultsList.setObjectName("resultsList")
        self.resultsList.setColumnCount(0)
        self.resultsList.setRowCount(0)
        self.resultsList.horizontalHeader().setCascadingSectionResizes(True)
        self.resultsList.horizontalHeader().setStretchLastSection(True)
        self.verticalLayout.addWidget(self.resultsList)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.exportResults = QtWidgets.QPushButton(Results)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.exportResults.sizePolicy().hasHeightForWidth())
        self.exportResults.setSizePolicy(sizePolicy)
        self.exportResults.setObjectName("exportResults")
        self.horizontalLayout_2.addWidget(self.exportResults)
        self.closeBtn = QtWidgets.QPushButton(Results)
        self.closeBtn.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.closeBtn.sizePolicy().hasHeightForWidth())
        self.closeBtn.setSizePolicy(sizePolicy)
        self.closeBtn.setObjectName("closeBtn")
        self.horizontalLayout_2.addWidget(self.closeBtn, 0, QtCore.Qt.AlignRight)
        self.verticalLayout.addLayout(self.horizontalLayout_2)

        self.retranslateUi(Results)
        self.closeBtn.pressed.connect(Results.close)
        QtCore.QMetaObject.connectSlotsByName(Results)

    def retranslateUi(self, Results):
        _translate = QtCore.QCoreApplication.translate
        Results.setWindowTitle(_translate("Results", "Results"))
        self.exportResults.setText(_translate("Results", "Save Results"))
        self.closeBtn.setText(_translate("Results", "Close"))

