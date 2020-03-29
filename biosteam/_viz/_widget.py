# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 21:59:13 2020

@author: yoelr
"""
from PyQt5.QtWidgets import (QAction, QApplication, QFileDialog,
                             QLabel, QWidget, QMainWindow, QMenu,
                             QMessageBox, QScrollArea, QSizePolicy, 
                             QGridLayout, QSizePolicy, QFrame)
from PyQt5.QtCore import QSize, QTimer, Qt
from PyQt5.QtGui import QPixmap, QPalette, QImage, QFont
from ._digraph import new_digraph, get_connections, update_digraph

__all__ = ('FlowsheetWidget',)

class FlowsheetWidget(QMainWindow):

    def __init__(self, flowsheet, autorefresh):
        super().__init__()
        self._autorefresh = False
        self.moveScale = 1
        self.scaleFactor = 1
        self.flowsheet = flowsheet
        self.connections = ()
        self.refresh_rate = 3 * 1000
        
        # Set main window
        minsize = QSize(400, 200)
        self.setMinimumSize(minsize) 
        self.flowsheetLabel = flowsheetLabel = QLabel()
        flowsheetLabel.setAlignment(Qt.AlignCenter)
        flowsheetLabel.setScaledContents(True)
        flowsheetLabel.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        
        self.notificationLabel = notificationLabel = QLabel()
        self.lastScrollPosition = None
        
        # Scroll Area
        scrollArea = QScrollArea(self)
        scrollArea.setWidgetResizable(True)
        scrollArea.setFocusPolicy(Qt.NoFocus)
        
        # Flowsheet content
        self.content = content = QWidget(self)
        scrollArea.ensureWidgetVisible(flowsheetLabel)
        content.setFocusPolicy(Qt.NoFocus)
        
        # Layout which will contain the flowsheet
        self.gridLayout = gridLayout = QGridLayout(self)
        gridLayout.addWidget(flowsheetLabel, 0, 0)
        content.setLayout(gridLayout)
        
        scrollArea.setWidget(content)
        self.scrollArea = scrollArea
        self.setCentralWidget(scrollArea)
        self.createActions()
        self.createMenus()
        self.refresh()
        self.autorefresh = autorefresh

        from biosteam import Unit
        Unit.IPYTHON_DISPLAY_UNIT_OPERATIONS = False

    def close(self):
        super().close()
        from biosteam import Unit
        Unit.IPYTHON_DISPLAY_UNIT_OPERATIONS = True

    def moveUp(self):
        bar = self.scrollArea.verticalScrollBar()
        bar.setValue(bar.value() - self.moveScale * bar.singleStep())
    
    def moveDown(self):
        bar = self.scrollArea.verticalScrollBar()
        bar.setValue(bar.value() + self.moveScale * bar.singleStep())
    
    def moveLeft(self):
        bar = self.scrollArea.horizontalScrollBar()
        bar.setValue(bar.value() - self.moveScale * bar.singleStep())
    
    def moveRight(self):
        bar = self.scrollArea.horizontalScrollBar()
        bar.setValue(bar.value() + self.moveScale * bar.singleStep())

    def zoomIn(self):
        self.scaleImage(1.10)

    def zoomOut(self):
        self.scaleImage(0.9)
        
    def normalSize(self):
        self.scaleFactor = 1.0
        flowsheetLabel = self.flowsheetLabel
        flowsheetLabel.adjustSize()
        
        width = self.width()
        height = self.height()
        label_width = flowsheetLabel.width()
        label_height = flowsheetLabel.height()
        min_width = label_width + 50
        min_height = label_height + 50
        self.resize(min_width, min_height)
        
    def scaleImage(self, factor):
        self.scaleFactor *= factor
        self.flowsheetLabel.resize(self.scaleFactor * self.flowsheetLabel.pixmap().size())

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
        self.successfulSimulation()

    def successfulSimulation(self):
        notificationLabel = self.notificationLabel
        notificationLabel.setParent(self)
        notificationLabel.setAlignment(Qt.AlignCenter)
        notificationLabel.setAutoFillBackground(True)
        notificationLabel.setText("SIMULATION WAS SUCCESSFUL")
        font = QFont("Times", 18, QFont.Bold)
        notificationLabel.setFont(font)
        self.gridLayout.addWidget(notificationLabel, 0, 0)
        self.centerScrollArea()
        qtimer = QTimer()
        qtimer.singleShot(2000, self.removeNotification)
        
    def centerScrollArea(self):
        scrollArea = self.scrollArea
        hbar = scrollArea.horizontalScrollBar()
        vbar = scrollArea.verticalScrollBar()
        self.lastScrollPosition = (vbar.value(), hbar.value())
        vbar = scrollArea.verticalScrollBar()
        hbar = scrollArea.horizontalScrollBar()
        vbar.setValue((vbar.minimum() + vbar.maximum()) / 2)
        hbar.setValue((hbar.minimum() + hbar.maximum()) / 2)
        
    def removeNotification(self):
        notificationLabel = self.notificationLabel
        notificationLabel.setParent(None)
        lastScrollPosition = self.lastScrollPosition
        if lastScrollPosition:
            vval, hval = lastScrollPosition
            scrollArea = self.scrollArea
            vbar = scrollArea.verticalScrollBar()
            vbar.setValue(vval)
            hbar = scrollArea.horizontalScrollBar()
            hbar.setValue(hval)
            self.lastScrollPosition = None
        
    def createActions(self):
        self.exitAct = QAction("E&xit", self, shortcut="Ctrl+Q",
                               triggered=self.close)
        self.refreshAct = QAction("&Refresh", self, shortcut="Ctrl+R",
                                  triggered=self.refresh)
        self.zoomInAct = QAction("Zoom &in (10%)", self, shortcut="Ctrl++",
                                 enabled=True, triggered=self.zoomIn)
        self.zoomOutAct = QAction("Zoom &out (10%)", self, shortcut="Ctrl+-",
                                  enabled=True, triggered=self.zoomOut)
        self.moveUpAct = QAction("Move &up", self, shortcut="Up",
                                 triggered=self.moveUp)
        self.moveDownAct = QAction("Move &down", self, shortcut="Down",
                                   triggered=self.moveDown)
        self.moveLeftAct = QAction("Move &down", self, shortcut="Left",
                                   triggered=self.moveLeft)
        self.moveRightAct = QAction("Move &down", self, shortcut="Right",
                                   triggered=self.moveRight)

        self.normalSizeAct = QAction("Normal &size", self, shortcut="Ctrl+N",
                                     triggered=self.normalSize)
        self.simulateAct = QAction("&Simulate", self, shortcut="Shift+Return",
                                   triggered=self.simulate)
        self.removeNotificationAct = QAction("&Remove notification", self, shortcut="Space",
                                   triggered=self.removeNotification)

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
        flowsheet = self.flowsheet
        ID = flowsheet.ID
        ID = ID.replace('_', ' ').capitalize()
        return f"{flowsheet.line} - {ID}"

    @property
    def autorefresh(self):
        return self._autorefresh
    @autorefresh.setter
    def autorefresh(self, autorefresh):
        autorefresh = bool(autorefresh)
        if autorefresh and not self._autorefresh:
            self.qtimer = qtimer = QTimer(self)
            qtimer.timeout.connect(self.refresh)
            qtimer.start(self.refresh_rate)
        elif not autorefresh and self._autorefresh:
            qtimer = self.qtimer
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
            img_data = digraph.pipe('png')
            img = QImage.fromData(img_data)
            pixmap = QPixmap.fromImage(img)
            if not pixmap.isNull():
                self.flowsheetLabel.setPixmap(pixmap)
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