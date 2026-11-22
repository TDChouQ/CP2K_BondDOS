"""Shared Qt widgets: background worker, plot panel, small helpers."""
from __future__ import annotations

import os
import sys
import traceback
from pathlib import Path

os.environ.setdefault("QT_API", "pyside6")

from PySide6.QtCore import QThread, Signal, Qt, QSettings  # noqa: E402
from PySide6.QtGui import QIcon  # noqa: E402
from PySide6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QDoubleSpinBox,  # noqa: E402
                               QSpinBox, QCheckBox, QFileDialog, QMessageBox, QGroupBox, QFormLayout,
                               QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView, QSizePolicy,
                               QToolButton, QComboBox, QLineEdit)

import matplotlib  # noqa: E402
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from ..plotting import PlotStyle, apply_rc  # noqa: E402

apply_rc()


def resource_path(rel):
    base = getattr(sys, "_MEIPASS", None)
    if base:
        return str(Path(base) / rel)
    return str(Path(__file__).resolve().parents[2] / rel)


def app_icon():
    for name in ("assets/app_icon.ico", "assets/app_icon.png"):
        p = resource_path(name)
        if Path(p).exists():
            return QIcon(p)
    return QIcon()


def settings():
    return QSettings("CP2K-BondDOS", "CP2K-BondDOS")


def last_dir(key="last_dir"):
    v = settings().value(key, "")
    return v if v and Path(v).exists() else str(Path.home())


def remember_dir(path, key="last_dir"):
    p = Path(path)
    settings().setValue(key, str(p if p.is_dir() else p.parent))


def error_box(parent, title, text, detail=None):
    box = QMessageBox(parent)
    box.setIcon(QMessageBox.Icon.Warning)
    box.setWindowTitle(title)
    box.setText(text)
    if detail:
        box.setDetailedText(detail)
    box.exec()


def dspin(value, lo, hi, step=0.1, decimals=3, suffix=""):
    w = QDoubleSpinBox()
    w.setRange(lo, hi)
    w.setDecimals(decimals)
    w.setSingleStep(step)
    w.setValue(value)
    if suffix:
        w.setSuffix(suffix)
    w.setKeyboardTracking(False)
    return w


def ispin(value, lo, hi):
    w = QSpinBox()
    w.setRange(lo, hi)
    w.setValue(value)
    w.setKeyboardTracking(False)
    return w


def combo(items, current=0):
    c = QComboBox()
    for it in items:
        if isinstance(it, tuple):
            c.addItem(it[0], it[1])
        else:
            c.addItem(it, it)
    c.setCurrentIndex(current)
    return c


def make_table(headers, stretch_last=True):
    t = QTableWidget(0, len(headers))
    t.setHorizontalHeaderLabels(headers)
    t.verticalHeader().setVisible(False)
    t.setAlternatingRowColors(True)
    t.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    t.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    h = t.horizontalHeader()
    h.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
    if stretch_last:
        h.setStretchLastSection(True)
    return t


def set_row(table, row, values, align_right_from=1):
    for j, v in enumerate(values):
        it = QTableWidgetItem("" if v is None else str(v))
        if j >= align_right_from:
            it.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        table.setItem(row, j, it)


def fmt(x, nd=4):
    try:
        if x is None or x != x:
            return "—"
        return f"{x:.{nd}f}"
    except Exception:
        return str(x)


class Worker(QThread):
    """Run fn(progress=cb, cancel=cb) in a thread; emits progress / done / failed."""
    progressed = Signal(float, str)
    done = Signal(object)
    failed = Signal(str, str)

    def __init__(self, fn, *args, **kwargs):
        super().__init__()
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self._cancel = False

    def cancel(self):
        self._cancel = True

    def run(self):
        try:
            res = self.fn(*self.args, progress=lambda f, t: self.progressed.emit(float(f), str(t)),
                          cancel=lambda: self._cancel, **self.kwargs)
            self.done.emit(res)
        except Exception as e:
            self.failed.emit(str(e), traceback.format_exc())


class PlotPanel(QWidget):
    """Matplotlib canvas + toolbar + style controls + export."""

    def __init__(self, parent=None, vertical_default=False):
        super().__init__(parent)
        self.figure = Figure(figsize=(7, 4.8), dpi=100, facecolor="white")
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        self.redraw_cb = None

        self.font = dspin(11, 6, 30, 1, 1)
        self.lw = dspin(1.6, 0.3, 6, 0.2, 1)
        self.fill = QCheckBox("填充")
        self.vertical = QCheckBox("能量竖排")
        self.vertical.setChecked(vertical_default)
        self.legend = QCheckBox("图例")
        self.legend.setChecked(True)
        self.grid = QCheckBox("网格")
        self.grid.setChecked(True)
        self.ef_line = QCheckBox("E_F 线")
        self.ef_line.setChecked(True)
        self.ymin = QLineEdit()
        self.ymin.setPlaceholderText("自动")
        self.ymax = QLineEdit()
        self.ymax.setPlaceholderText("自动")
        for w in (self.ymin, self.ymax):
            w.setMaximumWidth(70)
        self.title = QLineEdit()
        self.title.setPlaceholderText("图标题（可选）")
        self.exp_w = dspin(3.5, 1, 20, 0.1, 2, " in")
        self.exp_h = dspin(3.0, 1, 20, 0.1, 2, " in")
        self.exp_dpi = ispin(600, 72, 1200)
        save = QPushButton("保存图片…")
        save.clicked.connect(self.save_figure)

        style_box = QHBoxLayout()
        style_box.setSpacing(6)
        for lab, w in (("字号", self.font), ("线宽", self.lw)):
            style_box.addWidget(QLabel(lab))
            style_box.addWidget(w)
        for w in (self.fill, self.vertical, self.legend, self.grid, self.ef_line):
            style_box.addWidget(w)
        style_box.addWidget(QLabel("纵轴"))
        style_box.addWidget(self.ymin)
        style_box.addWidget(QLabel("~"))
        style_box.addWidget(self.ymax)
        style_box.addStretch(1)

        exp_box = QHBoxLayout()
        exp_box.addWidget(self.title, 2)
        exp_box.addWidget(QLabel("导出尺寸"))
        exp_box.addWidget(self.exp_w)
        exp_box.addWidget(QLabel("×"))
        exp_box.addWidget(self.exp_h)
        exp_box.addWidget(QLabel("DPI"))
        exp_box.addWidget(self.exp_dpi)
        exp_box.addWidget(save)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.toolbar)
        lay.addWidget(self.canvas, 1)
        lay.addLayout(style_box)
        lay.addLayout(exp_box)
        for w in (self.font, self.lw):
            w.valueChanged.connect(self._changed)
        for w in (self.fill, self.vertical, self.legend, self.grid, self.ef_line):
            w.toggled.connect(self._changed)
        for w in (self.ymin, self.ymax, self.title):
            w.editingFinished.connect(self._changed)

    def _changed(self, *a):
        if self.redraw_cb:
            self.redraw_cb()

    @staticmethod
    def _num(edit):
        t = edit.text().strip()
        try:
            return float(t) if t else None
        except ValueError:
            return None

    def style(self, xmin=None, xmax=None, mirror_beta=True):
        return PlotStyle(font_size=self.font.value(), line_width=self.lw.value(), fill=self.fill.isChecked(),
                         legend=self.legend.isChecked(), grid=self.grid.isChecked(),
                         fermi_line=self.ef_line.isChecked(), vertical=self.vertical.isChecked(),
                         mirror_beta=mirror_beta, xmin=xmin, xmax=xmax, ymin=self._num(self.ymin),
                         ymax=self._num(self.ymax), title=self.title.text().strip())

    def draw(self):
        self.canvas.draw_idle()

    def save_figure(self):
        path, _ = QFileDialog.getSaveFileName(self, "保存图片", str(Path(last_dir()) / "figure.png"),
                                              "PNG (*.png);;PDF (*.pdf);;SVG (*.svg);;TIFF (*.tif)")
        if not path:
            return
        remember_dir(path)
        old = self.figure.get_size_inches()
        try:
            self.figure.set_size_inches(self.exp_w.value(), self.exp_h.value())
            self.figure.tight_layout()
            self.figure.savefig(path, dpi=self.exp_dpi.value(), bbox_inches="tight")
        except Exception as e:
            error_box(self, "保存失败", str(e))
        finally:
            self.figure.set_size_inches(*old)
            self.figure.tight_layout()
            self.canvas.draw_idle()


def group(title, layout):
    g = QGroupBox(title)
    g.setLayout(layout)
    return g


def form(*rows):
    f = QFormLayout()
    f.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
    for lab, w in rows:
        if isinstance(w, (list, tuple)):
            h = QHBoxLayout()
            for x in w:
                if isinstance(x, str):
                    h.addWidget(QLabel(x))
                else:
                    h.addWidget(x)
            f.addRow(lab, h)
        else:
            f.addRow(lab, w)
    return f


def tool_button(text, tip=""):
    b = QToolButton()
    b.setText(text)
    b.setToolTip(tip)
    return b
