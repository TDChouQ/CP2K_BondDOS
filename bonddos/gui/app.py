"""Main window."""
from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import Qt, QSize
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QListWidget,
                               QListWidgetItem, QStackedWidget, QPlainTextEdit, QTextBrowser, QLabel, QStatusBar)

from .. import __version__, APP_TITLE
from ..session import BondSession
from .common import app_icon, resource_path
from .dos_page import DosPage
from .cohp_page import CohpPage
from .pop_page import PopPage

QSS = """
QMainWindow, QWidget { background: #f6f7f9; color: #1f2328; font-size: 10pt; }
QListWidget#nav { background: #14324f; color: #dfe8f1; border: none; font-size: 11pt; outline: 0; }
QListWidget#nav::item { padding: 12px 14px; border-left: 3px solid transparent; }
QListWidget#nav::item:selected { background: #1d4670; color: white; border-left: 3px solid #35c2a1; }
QListWidget#nav::item:hover:!selected { background: #1a3d60; }
QGroupBox { background: white; border: 1px solid #dfe3e8; border-radius: 6px; margin-top: 14px; padding: 8px; }
QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 4px; color: #14324f; font-weight: bold; }
QPushButton { background: white; border: 1px solid #c8cfd8; border-radius: 4px; padding: 5px 10px; }
QPushButton:hover { background: #eef3f8; }
QPushButton:disabled { color: #9aa3ad; }
QPushButton#primary { background: #1f6fd1; color: white; border: 1px solid #1f6fd1; font-weight: bold; }
QPushButton#primary:hover { background: #2a7de0; }
QPushButton#primary:disabled { background: #9db9de; border-color: #9db9de; }
QLineEdit, QComboBox, QDoubleSpinBox, QSpinBox, QPlainTextEdit { background: white; border: 1px solid #c8cfd8;
    border-radius: 3px; padding: 3px; }
QTableWidget { background: white; gridline-color: #e6e9ee; alternate-background-color: #f7f9fb; }
QHeaderView::section { background: #eef1f5; border: none; border-right: 1px solid #dde2e8; padding: 4px; }
QTabWidget::pane { border: 1px solid #dfe3e8; background: white; }
QTabBar::tab { padding: 6px 12px; }
QTabBar::tab:selected { background: white; font-weight: bold; }
QTextBrowser { background: white; border: 1px solid #dfe3e8; }
QProgressBar { border: none; background: #e3e7ec; border-radius: 3px; }
QProgressBar::chunk { background: #35c2a1; border-radius: 3px; }
QScrollArea { border: none; }
"""


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"{APP_TITLE}  v{__version__}")
        self.setWindowIcon(app_icon())
        self.resize(1480, 920)
        self.session = BondSession()

        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setFont(QFont("Consolas", 9))
        self.setStatusBar(QStatusBar())

        self.dos = DosPage(self.log)
        self.cohp = CohpPage(self.session, self.log, self.status)
        self.pop = PopPage(self.session, self.log, self.status)
        self.cohp.loaded.connect(self.pop.on_loaded)
        self.methods = QTextBrowser()
        self.methods.setOpenExternalLinks(True)
        doc = Path(resource_path("docs/方法说明.md"))
        if doc.exists():
            self.methods.setMarkdown(doc.read_text(encoding="utf-8"))
        else:
            self.methods.setPlainText("缺少 docs/方法说明.md")

        self.cp2k_doc = QTextBrowser()
        self.cp2k_doc.setOpenExternalLinks(True)
        doc = Path(resource_path("docs/CP2K输入设置.md"))
        if doc.exists():
            self.cp2k_doc.setMarkdown(doc.read_text(encoding="utf-8"))
        else:
            self.cp2k_doc.setPlainText("缺少 docs/CP2K输入设置.md")

        self.stack = QStackedWidget()
        self.nav = QListWidget()
        self.nav.setObjectName("nav")
        self.nav.setFixedWidth(190)
        self.nav.setIconSize(QSize(20, 20))
        pages = [("态密度 DOS / PDOS", self.dos), ("k 点成键 COHP", self.cohp), ("电荷与原子 PDOS", self.pop),
                 ("CP2K 输入设置", self.cp2k_doc), ("方法说明", self.methods), ("运行日志", self.log_view)]
        for name, w in pages:
            QListWidgetItem(name, self.nav)
            self.stack.addWidget(w)
        self.nav.currentRowChanged.connect(self.stack.setCurrentIndex)
        self.nav.setCurrentRow(0)

        side = QWidget()
        side.setObjectName("side")
        side.setStyleSheet("#side { background: #14324f; }")
        sl = QVBoxLayout(side)
        sl.setContentsMargins(0, 0, 0, 0)
        brand = QLabel()
        brand.setAlignment(Qt.AlignmentFlag.AlignCenter)
        pm = app_icon().pixmap(64, 64)
        if not pm.isNull():
            brand.setPixmap(pm)
        brand.setStyleSheet("background:#14324f; padding-top:16px;")
        title = QLabel("CP2K BondDOS")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title.setStyleSheet("background:#14324f; color:white; font-size:13pt; font-weight:bold; padding:6px;")
        ver = QLabel(f"v{__version__} · CP2K 2026.2+")
        ver.setAlignment(Qt.AlignmentFlag.AlignCenter)
        ver.setStyleSheet("background:#14324f; color:#9fb4c8; padding-bottom:12px;")
        sl.addWidget(brand)
        sl.addWidget(title)
        sl.addWidget(ver)
        sl.addWidget(self.nav, 1)

        central = QWidget()
        lay = QHBoxLayout(central)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        lay.addWidget(side)
        lay.addWidget(self.stack, 1)
        self.setCentralWidget(central)
        self.log(f"CP2K BondDOS {__version__} 启动")

    def log(self, text):
        self.log_view.appendPlainText(str(text))

    def status(self, text):
        self.statusBar().showMessage(str(text))


def run_gui():
    if sys.platform == "win32":
        try:
            import ctypes
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("CP2K.BondDOS.1")
        except Exception:
            pass
    app = QApplication.instance() or QApplication(sys.argv)
    app.setWindowIcon(app_icon())
    app.setFont(QFont("Microsoft YaHei UI", 9))
    app.setStyleSheet(QSS)
    w = MainWindow()
    w.show()
    return app.exec()
