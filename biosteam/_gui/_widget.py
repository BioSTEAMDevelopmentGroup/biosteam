# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 21:59:13 2020

@author: yoelr
"""
import os
from PyQt5.QtWidgets import (QAction, QApplication, QFileDialog,
                             QLabel, QWidget, QMainWindow, QMenu,
                             QMessageBox, QScrollArea, QSizePolicy, 
                             QGridLayout, QSizePolicy)
from PyQt5.QtCore import QSize, QTimer, Qt
from PyQt5.QtGui import QPixmap, QPalette, QImage
from ._digraph import new_digraph, get_connections, update_digraph

__all__ = ('FlowsheetWidget',)

# image_path = os.path.join(os.path.dirname(__file__), "flowsheet_digraph.png") 
image_format = 'png'

class FlowsheetWidget(QMainWindow):

    def __init__(self, flowsheet, autorefresh):
        super().__init__()
        self._autorefresh = False
        self.qtimer = QTimer(self)
        self.scaleFactor = 1
        self.flowsheet = flowsheet
        self.connections = ()
        self.refresh_rate = 3 * 1000
        
        # Set main window
        minsize = QSize(400, 200)
        self.setMinimumSize(minsize) 
        self.imageLabel = imageLabel = QLabel()
        imageLabel.setAlignment(Qt.AlignCenter)
        imageLabel.setScaledContents(True)
        
        # Scroll Area
        scrollArea = QScrollArea(self)
        scrollArea.setWidgetResizable(True)
        scrollArea.setFocusPolicy(Qt.NoFocus)
        
        # Flowsheet content
        self.content = content = QWidget(self)
        scrollArea.ensureWidgetVisible(imageLabel)
        content.setFocusPolicy(Qt.NoFocus)
        
        # Layout which will contain the flowsheet
        gridLayout = QGridLayout(self)
        gridLayout.addWidget(imageLabel, 0, 0)
        content.setLayout(gridLayout)
        
        scrollArea.setWidget(content)
        self.scrollArea = scrollArea
        self.setCentralWidget(scrollArea)
        self.createActions()
        self.createMenus()
        self.refresh()
        self.autorefresh = autorefresh

    def moveUp(self):
        bar = self.scrollArea.verticalScrollBar()
        bar.setValue(bar.value() - bar.singleStep())
    
    def moveDown(self):
        bar = self.scrollArea.verticalScrollBar()
        bar.setValue(bar.value() + bar.singleStep())
    
    def moveLeft(self):
        bar = self.scrollArea.horizontalScrollBar()
        bar.setValue(bar.value() - bar.singleStep())
    
    def moveRight(self):
        bar = self.scrollArea.horizontalScrollBar()
        bar.setValue(bar.value() + bar.singleStep())

    def zoomIn(self):
        self.scaleImage(1.10)

    def zoomOut(self):
        self.scaleImage(0.9)
        
    def normalSize(self):
        self.scaleFactor = 1.0
        imageLabel = self.imageLabel
        imageLabel.adjustSize()
        
        width = self.width()
        height = self.height()
        label_width = imageLabel.width()
        label_height = imageLabel.height()
        min_width = label_width + 50
        min_height = label_height + 50
        self.resize(min_width, min_height)
        
    def scaleImage(self, factor):
        self.scaleFactor *= factor
        self.imageLabel.resize(self.scaleFactor * self.imageLabel.pixmap().size())

        self.adjustScrollBar(self.scrollArea.horizontalScrollBar(), factor)
        self.adjustScrollBar(self.scrollArea.verticalScrollBar(), factor)
        
        self.zoomInAct.setEnabled(self.scaleFactor < 4.0)
        self.zoomOutAct.setEnabled(self.scaleFactor > 0.25)
        
    def adjustScrollBar(self, scrollBar, factor):
        value = int(factor * scrollBar.value()
                    + ((factor - 1) * scrollBar.pageStep()/2))
        scrollBar.setValue(value)

    def simulate(self):
        self.refresh()
        system = self.flowsheet.create_system()
        system.simulate()
        print("SIMULATION WAS SUCCESSFUL\n")

    def createActions(self):
        self.exitAct = QAction("E&xit", self, shortcut="Ctrl+Q",
                               triggered=self.close)
        self.refreshAct = QAction("&Refresh", self, shortcut="Ctrl+R",
                                  triggered=self.refresh)
        self.zoomInAct = QAction("Zoom &in (10%)", self, shortcut="Ctrl++",
                                 enabled=True, triggered=self.zoomIn)
        self.zoomOutAct = QAction("Zoom &out (10%)", self, shortcut="Ctrl+-",
                                  enabled=True, triggered=self.zoomOut)
        self.moveUpAct = QAction("&Move up", self, shortcut="Up",
                                 triggered=self.moveUp)
        self.moveDownAct = QAction("&Move down", self, shortcut="Down",
                                   triggered=self.moveDown)
        self.moveLeftAct = QAction("&Move down", self, shortcut="Left",
                                   triggered=self.moveLeft)
        self.moveRightAct = QAction("&Move down", self, shortcut="Right",
                                   triggered=self.moveRight)
        self.normalSizeAct = QAction("&Normal size", self, shortcut="Ctrl+N",
                                     triggered=self.normalSize)
        self.simulateAct = QAction("&Simulate", self, shortcut="Shift+Return",
                                   triggered=self.simulate)

    def createMenus(self):
        self.fileMenu = fileMenu = QMenu("&File", self)
        fileMenu.addSeparator()
        fileMenu.addAction(self.exitAct)

        self.viewMenu = viewMenu = QMenu("&View", self)
        viewMenu.addAction(self.refreshAct)
        viewMenu.addSeparator()
        viewMenu.addAction(self.normalSizeAct)
        viewMenu.addAction(self.zoomInAct)
        viewMenu.addAction(self.zoomOutAct)
        viewMenu.addSeparator()
        viewMenu.addAction(self.moveUpAct)
        viewMenu.addAction(self.moveDownAct)
        viewMenu.addAction(self.moveLeftAct)
        viewMenu.addAction(self.moveRightAct)
        
        self.simulationMenu = simulationMenu = QMenu("&Simulation", self)
        simulationMenu.addAction(self.simulateAct)
        
        menuBar = self.menuBar()
        menuBar.addMenu(fileMenu)
        menuBar.addMenu(viewMenu)
        menuBar.addMenu(simulationMenu)

    @property
    def title(self):
        ID = self.flowsheet.ID
        ID = ID.replace('_', ' ').capitalize()
        return f"Main flowsheet - {ID}"

    @property
    def autorefresh(self):
        return self._autorefresh
    @autorefresh.setter
    def autorefresh(self, autorefresh):
        autorefresh = bool(autorefresh)
        qtimer = self.qtimer
        if autorefresh and not self._autorefresh:
            qtimer.timeout.connect(self.refresh)
            qtimer.start(self.refresh_rate)
        elif not autorefresh and self._autorefresh:
            qtimer.stop()
        self._autorefresh = autorefresh    

    def refresh(self):
        self.setWindowTitle(self.title)
        flowsheet = self.flowsheet
        connections = get_connections(flowsheet.stream)
        if self.connections != connections:
            self.connections = connections
            digraph = new_digraph()
            update_digraph(digraph, flowsheet.unit, connections)
            img_data = digraph.pipe(image_format)
            img = QImage.fromData(img_data)
            pixmap = QPixmap.fromImage(img)
            if not pixmap.isNull():
                self.imageLabel.setPixmap(pixmap)
                self.normalSize()


# def moveUpLeft(self):
#     self.moveUp()
#     self.moveLeft()

# def moveUpRight(self):
#     self.moveUp()
#     self.moveRight()

# def moveDownLeft(self):
#     self.moveDown()
#     self.moveLeft()
    
# def moveDownRight(self):
#     self.moveDown()
#     self.moveRight()

# viewMenu.addAction(self.moveUpRightAct)
# viewMenu.addAction(self.moveUpLeftAct)
# viewMenu.addAction(self.moveDownLeftAct)
# viewMenu.addAction(self.moveDownRightAct)
# self.moveUpLeftAct = QAction("&Move up left", self, shortcut="Up+Left",
#                          triggered=self.moveUpLeft)
# self.moveUpRightAct = QAction("&Move up right", self, shortcut="Up+Right",
#                            triggered=self.moveUpRight)
# self.moveDownLeftAct = QAction("&Move down left", self, shortcut="Down+Left",
#                            triggered=self.moveDownLeft)
# self.moveDownRightAct = QAction("&Move down right", self, shortcut="Down+Right",
#                            triggered=self.moveDownRight)