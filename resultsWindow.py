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
        self.selectAllBtn = QtWidgets.QPushButton(Results)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.selectAllBtn.sizePolicy().hasHeightForWidth())
        self.selectAllBtn.setSizePolicy(sizePolicy)
        self.selectAllBtn.setObjectName("selectAllBtn")
        self.horizontalLayout_2.addWidget(self.selectAllBtn)
        self.exportSummaryResultsBtn = QtWidgets.QPushButton(Results)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.exportSummaryResultsBtn.sizePolicy().hasHeightForWidth())
        self.exportSummaryResultsBtn.setSizePolicy(sizePolicy)
        self.exportSummaryResultsBtn.setObjectName("exportSummaryResultsBtn")
        self.horizontalLayout_2.addWidget(self.exportSummaryResultsBtn)
        self.viewHTMLvisBtn = QtWidgets.QPushButton(Results)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.viewHTMLvisBtn.sizePolicy().hasHeightForWidth())
        self.viewHTMLvisBtn.setSizePolicy(sizePolicy)
        self.viewHTMLvisBtn.setObjectName("viewHTMLvisBtn")
        self.horizontalLayout_2.addWidget(self.viewHTMLvisBtn)
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
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.label = QtWidgets.QLabel(Results)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setObjectName("label")
        self.verticalLayout_2.addWidget(self.label)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.gbDlWindowSize = QtWidgets.QLineEdit(Results)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.gbDlWindowSize.sizePolicy().hasHeightForWidth())
        self.gbDlWindowSize.setSizePolicy(sizePolicy)
        self.gbDlWindowSize.setInputMethodHints(QtCore.Qt.ImhFormattedNumbersOnly)
        self.gbDlWindowSize.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.gbDlWindowSize.setObjectName("gbDlWindowSize")
        self.horizontalLayout_3.addWidget(self.gbDlWindowSize)
        self.label_3 = QtWidgets.QLabel(Results)
        self.label_3.setObjectName("label_3")
        self.horizontalLayout_3.addWidget(self.label_3)
        self.verticalLayout_2.addLayout(self.horizontalLayout_3)
        self.horizontalLayout.addLayout(self.verticalLayout_2)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.label_2 = QtWidgets.QLabel(Results)
        self.label_2.setObjectName("label_2")
        self.verticalLayout_3.addWidget(self.label_2)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.gbkExportDlDir = QtWidgets.QLineEdit(Results)
        self.gbkExportDlDir.setReadOnly(True)
        self.gbkExportDlDir.setObjectName("gbkExportDlDir")
        self.horizontalLayout_4.addWidget(self.gbkExportDlDir)
        self.selectGbkExportDirBtn = QtWidgets.QPushButton(Results)
        self.selectGbkExportDirBtn.setObjectName("selectGbkExportDirBtn")
        self.horizontalLayout_4.addWidget(self.selectGbkExportDirBtn)
        self.verticalLayout_3.addLayout(self.horizontalLayout_4)
        self.horizontalLayout.addLayout(self.verticalLayout_3)
        self.exportSelectedGbkBtn = QtWidgets.QPushButton(Results)
        self.exportSelectedGbkBtn.setEnabled(False)
        self.exportSelectedGbkBtn.setObjectName("exportSelectedGbkBtn")
        self.horizontalLayout.addWidget(self.exportSelectedGbkBtn)
        self.verticalLayout.addLayout(self.horizontalLayout)

        self.retranslateUi(Results)
        self.closeBtn.pressed.connect(Results.close)
        QtCore.QMetaObject.connectSlotsByName(Results)

    def retranslateUi(self, Results):
        _translate = QtCore.QCoreApplication.translate
        Results.setWindowTitle(_translate("Results", "Results"))
        self.resultsList.setSortingEnabled(True)
        self.selectAllBtn.setText(_translate("Results", "Select All"))
        self.exportSummaryResultsBtn.setText(_translate("Results", "Save Results Summary"))
        self.viewHTMLvisBtn.setText(_translate("Results", "View HTML Visualization"))
        self.closeBtn.setText(_translate("Results", "Close"))
        self.label.setText(_translate("Results", "Window Size:"))
        self.gbDlWindowSize.setText(_translate("Results", "100000"))
        self.label_3.setText(_translate("Results", " (bp from midpoint)"))
        self.label_2.setText(_translate("Results", "Genbank Download Directory:"))
        self.selectGbkExportDirBtn.setText(_translate("Results", "Select Directory"))
        self.exportSelectedGbkBtn.setText(_translate("Results", "Download Selected"))

