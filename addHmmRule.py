# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'addHmmRule.ui'
#
# Created by: PyQt5 UI code generator 5.7.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_addHmmRuleWindow(object):
    def setupUi(self, addHmmRuleWindow):
        addHmmRuleWindow.setObjectName("addHmmRuleWindow")
        addHmmRuleWindow.resize(413, 89)
        self.verticalLayout = QtWidgets.QVBoxLayout(addHmmRuleWindow)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(addHmmRuleWindow)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        self.hmmRule = QtWidgets.QLineEdit(addHmmRuleWindow)
        self.hmmRule.setObjectName("hmmRule")
        self.horizontalLayout.addWidget(self.hmmRule)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.addHmmRuleBtn = QtWidgets.QPushButton(addHmmRuleWindow)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.addHmmRuleBtn.sizePolicy().hasHeightForWidth())
        self.addHmmRuleBtn.setSizePolicy(sizePolicy)
        self.addHmmRuleBtn.setObjectName("addHmmRuleBtn")
        self.horizontalLayout_3.addWidget(self.addHmmRuleBtn)
        self.closeBtn = QtWidgets.QPushButton(addHmmRuleWindow)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.closeBtn.sizePolicy().hasHeightForWidth())
        self.closeBtn.setSizePolicy(sizePolicy)
        self.closeBtn.setObjectName("closeBtn")
        self.horizontalLayout_3.addWidget(self.closeBtn)
        self.verticalLayout.addLayout(self.horizontalLayout_3)

        self.retranslateUi(addHmmRuleWindow)
        self.closeBtn.clicked.connect(addHmmRuleWindow.close)
        QtCore.QMetaObject.connectSlotsByName(addHmmRuleWindow)

    def retranslateUi(self, addHmmRuleWindow):
        _translate = QtCore.QCoreApplication.translate
        addHmmRuleWindow.setWindowTitle(_translate("addHmmRuleWindow", "Add HMM Rule"))
        self.label.setText(_translate("addHmmRuleWindow", "Rule:"))
        self.addHmmRuleBtn.setText(_translate("addHmmRuleWindow", "Add HMM Rule"))
        self.closeBtn.setText(_translate("addHmmRuleWindow", "Close"))

