# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'firstgui.ui'
#
# Created by: PyQt5 UI code generator 5.7.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_myfirstgui(object):
    def setupUi(self, myfirstgui):
        myfirstgui.setObjectName("myfirstgui")
        myfirstgui.resize(784, 505)
        self.buttonBox = QtWidgets.QDialogButtonBox(myfirstgui)
        self.buttonBox.setGeometry(QtCore.QRect(390, 460, 381, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Close)
        self.buttonBox.setObjectName("buttonBox")
        self.myTextInput = QtWidgets.QLineEdit(myfirstgui)
        self.myTextInput.setGeometry(QtCore.QRect(10, 10, 281, 21))
        self.myTextInput.setObjectName("myTextInput")
        self.listWidget = QtWidgets.QListWidget(myfirstgui)
        self.listWidget.setGeometry(QtCore.QRect(10, 40, 281, 192))
        self.listWidget.setObjectName("listWidget")
        self.clearBtn = QtWidgets.QPushButton(myfirstgui)
        self.clearBtn.setGeometry(QtCore.QRect(190, 230, 101, 23))
        self.clearBtn.setObjectName("clearBtn")
        self.addBtn = QtWidgets.QPushButton(myfirstgui)
        self.addBtn.setGeometry(QtCore.QRect(0, 230, 101, 23))
        self.addBtn.setObjectName("addBtn")
        self.listWidget_2 = QtWidgets.QListWidget(myfirstgui)
        self.listWidget_2.setGeometry(QtCore.QRect(10, 290, 281, 192))
        self.listWidget_2.setObjectName("listWidget_2")
        self.myTextInput_2 = QtWidgets.QLineEdit(myfirstgui)
        self.myTextInput_2.setGeometry(QtCore.QRect(10, 260, 281, 21))
        self.myTextInput_2.setObjectName("myTextInput_2")
        self.addBtn_2 = QtWidgets.QPushButton(myfirstgui)
        self.addBtn_2.setGeometry(QtCore.QRect(10, 480, 101, 23))
        self.addBtn_2.setObjectName("addBtn_2")
        self.clearBtn_2 = QtWidgets.QPushButton(myfirstgui)
        self.clearBtn_2.setGeometry(QtCore.QRect(200, 480, 101, 23))
        self.clearBtn_2.setObjectName("clearBtn_2")
        self.tableWidget = QtWidgets.QTableWidget(myfirstgui)
        self.tableWidget.setGeometry(QtCore.QRect(470, 210, 256, 192))
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(0)
        self.tableWidget.setRowCount(0)

        self.retranslateUi(myfirstgui)
        self.buttonBox.accepted.connect(myfirstgui.accept)
        self.buttonBox.rejected.connect(myfirstgui.reject)
        self.clearBtn.clicked.connect(self.listWidget.clear)
        QtCore.QMetaObject.connectSlotsByName(myfirstgui)

    def retranslateUi(self, myfirstgui):
        _translate = QtCore.QCoreApplication.translate
        myfirstgui.setWindowTitle(_translate("myfirstgui", "My First Gui!"))
        self.clearBtn.setText(_translate("myfirstgui", "clear"))
        self.addBtn.setText(_translate("myfirstgui", "add"))
        self.addBtn_2.setText(_translate("myfirstgui", "add"))
        self.clearBtn_2.setText(_translate("myfirstgui", "clear"))

