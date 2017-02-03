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
        self.gridLayout = QtWidgets.QGridLayout(Results)
        self.gridLayout.setObjectName("gridLayout")
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
        self.gridLayout.addWidget(self.resultsList, 0, 0, 1, 1)
        self.closeBtn = QtWidgets.QPushButton(Results)
        self.closeBtn.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.closeBtn.sizePolicy().hasHeightForWidth())
        self.closeBtn.setSizePolicy(sizePolicy)
        self.closeBtn.setObjectName("closeBtn")
        self.gridLayout.addWidget(self.closeBtn, 1, 0, 1, 1, QtCore.Qt.AlignRight)

        self.retranslateUi(Results)
        self.closeBtn.pressed.connect(Results.close)
        QtCore.QMetaObject.connectSlotsByName(Results)

    def retranslateUi(self, Results):
        _translate = QtCore.QCoreApplication.translate
        Results.setWindowTitle(_translate("Results", "Results"))
        self.closeBtn.setText(_translate("Results", "Close"))

