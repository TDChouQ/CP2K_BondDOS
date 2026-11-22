"""Page 3: charges, gross populations and atom/orbital-resolved PDOS from MO_KP."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QSplitter, QScrollArea, QLabel,
                               QPlainTextEdit, QCheckBox, QFileDialog, QTabWidget)

from ..dosanalysis import Series, DosProject, write_series_csv, write_stats_csv
from ..kspace import parse_atom_selection
from ..plotting import draw_lines
from ..population import PDOSGroup, gross_populations
from .common import (PlotPanel, Worker, combo, make_table, set_row, fmt, group, form, last_dir, remember_dir,
                     error_box, dspin)

RES = [("总和", "total"), ("s/p/d/f", "l"), ("磁量子数 (px, dxy…)", "lm"), ("基组壳层", "shell")]
METHODS = [("Mulliken", "mulliken"), ("Löwdin", "lowdin")]


class PopPage(QWidget):
    def __init__(self, session, log, status, parent=None):
        super().__init__(parent)
        self.session = session
        self.log = log
        self.status = status
        self.worker = None
        self.series = []

        self.groups = QPlainTextEdit()
        self.groups.setPlaceholderText("每行一组：名称: 原子选择\n例如\nW: W\nO_桥: 3,5,7-9\nFe1: 1")
        self.groups.setMaximumHeight(130)
        b_each = QPushButton("每种元素一组")
        b_each.clicked.connect(self.fill_elements)
        self.res = combo(RES, 1)
        self.method = combo(METHODS)
        self.tdos = QCheckBox("同时绘制总 DOS")
        self.tdos.setChecked(True)
        self.occ = QCheckBox("按占据数加权（电子/eV）")
        self.mirror = QCheckBox("β 自旋向下镜像")
        self.mirror.setChecked(True)
        self.smin = dspin(-10, -500, 500, 1, 2, " eV")
        self.smax = dspin(0, -500, 500, 1, 2, " eV")
        self.b_pdos = QPushButton("计算原子 PDOS")
        self.b_pdos.setObjectName("primary")
        self.b_pdos.clicked.connect(self.run_pdos)
        self.b_charge = QPushButton("计算电荷与布居")
        self.b_charge.clicked.connect(self.run_charges)
        b_csv = QPushButton("导出 PDOS CSV…")
        b_csv.clicked.connect(self.export_pdos)
        b_stat = QPushButton("导出统计 CSV…")
        b_stat.clicked.connect(self.export_stats)
        b_chg = QPushButton("导出电荷 CSV…")
        b_chg.clicked.connect(self.export_charges)

        gl = QVBoxLayout()
        gl.addWidget(self.groups)
        gl.addWidget(b_each)
        gl.addLayout(form(("轨道分解", self.res), ("投影方法", self.method), ("统计能区", [self.smin, "~", self.smax])))
        for w in (self.tdos, self.occ, self.mirror):
            gl.addWidget(w)
        gl.addWidget(self.b_pdos)
        note = QLabel("能量网格、展宽和能量零点沿用“k 点成键”页的设置。原子分辨 PDOS 直接由 MO_KP 的系数与 S(k) "
                      "投影得到，无需在 CP2K 中把目标原子拆成单独 KIND。Mulliken 单轨道权重可略为负。")
        note.setWordWrap(True)
        note.setStyleSheet("color:#666;")
        gl.addWidget(note)
        cl = QVBoxLayout()
        cl.addWidget(self.b_charge)
        h = QHBoxLayout()
        h.addWidget(b_chg)
        cl.addLayout(h)
        left = QWidget()
        ll = QVBoxLayout(left)
        ll.addWidget(group("原子 / 轨道分辨 PDOS（k 点）", gl))
        ll.addWidget(group("电荷与布居", cl))
        h = QHBoxLayout()
        h.addWidget(b_csv)
        h.addWidget(b_stat)
        ll.addLayout(h)
        ll.addStretch(1)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(left)
        scroll.setMinimumWidth(360)

        self.plot = PlotPanel()
        self.plot.redraw_cb = self.redraw
        self.stats = make_table(["名称", "通道", "自旋", "积分", "中心 ε (eV)", "宽度", "偏度", "峰度", "填充率",
                                 "E≤0 中心", "上边缘 (eV)"])
        self.charges = make_table(["原子", "元素", "Z_eff", "Mulliken 布居", "Mulliken 电荷", "Löwdin 布居",
                                   "Löwdin 电荷", "Mulliken 磁矩", "Löwdin 磁矩"])
        self.gross = make_table(["原子", "元素", "壳层", "l", "轨道", "Mulliken (up/总)", "Mulliken (down)", "合计"])
        tabs = QTabWidget()
        tabs.addTab(self.stats, "PDOS 统计（能带中心等）")
        tabs.addTab(self.charges, "原子电荷")
        tabs.addTab(self.gross, "轨道布居 (GROSSPOP)")
        self.tabs = tabs
        right = QSplitter(Qt.Orientation.Vertical)
        right.addWidget(self.plot)
        right.addWidget(tabs)
        right.setStretchFactor(0, 3)
        right.setStretchFactor(1, 2)
        split = QSplitter(Qt.Orientation.Horizontal)
        split.addWidget(scroll)
        split.addWidget(right)
        split.setStretchFactor(1, 1)
        lay = QHBoxLayout(self)
        lay.addWidget(split)

    def on_loaded(self):
        if not self.groups.toPlainText().strip():
            self.fill_elements()
        self.charges.setRowCount(0)
        self.gross.setRowCount(0)
        self.series = []
        self.stats.setRowCount(0)

    def fill_elements(self):
        if self.session.system is None:
            return
        els = []
        for a in self.session.system.d.atoms:
            if a.element not in els:
                els.append(a.element)
        self.groups.setPlainText("\n".join(f"{e}: {e}" for e in els))

    def parse_groups(self):
        sysm = self.session.system
        out = []
        for line in self.groups.toPlainText().splitlines():
            line = line.strip()
            if not line:
                continue
            if ":" in line or "：" in line:
                name, sel = line.replace("：", ":").split(":", 1)
            else:
                name, sel = line, line
            atoms = parse_atom_selection(sel, len(sysm.d.atoms), sysm.structure)
            out.append(PDOSGroup(name.strip(), atoms, self.res.currentData()))
        if not out:
            raise ValueError("请至少定义一组原子")
        return out

    def _start(self, fn, done, text):
        if self.worker is not None and self.worker.isRunning():
            return
        self.worker = Worker(fn)
        self.worker.progressed.connect(self._on_progress)
        self.worker.done.connect(done)
        self.worker.failed.connect(self._failed)
        self.status(text)
        self.worker.start()

    def _on_progress(self, f, t):
        self.status(t)

    def _failed(self, m, tb):
        self.log(f"[PDOS] 错误 {m}\n{tb}")
        error_box(self, "失败", m, tb)

    def run_pdos(self):
        s = self.session
        if s.system is None:
            error_box(self, "尚未读取数据", "请先在“k 点成键”页读取 MO_KP 文件")
            return
        if not s.system.has_S:
            error_box(self, "缺少重叠矩阵", "MO_KP 需使用 AO_EXPORT_TYPE OVERLAP_MATRIX")
            return
        try:
            groups = self.parse_groups()
        except Exception as e:
            error_box(self, "原子组有误", str(e))
            return
        method = self.method.currentData()
        tdos = self.tdos.isChecked()
        occ = self.occ.isChecked()

        def job(progress, cancel):
            s.update_reference()
            return s.pdos(groups, method, tdos, occ, progress)
        self._start(job, self._pdos_done, "正在计算原子 PDOS…")

    def _pdos_done(self, res):
        g, cur = res
        s = self.session
        spin = {0: "none"} if s.system.ns == 1 else {0: "alpha", 1: "beta"}
        self.series = []
        for sp, d in cur.items():
            for name, y in d.items():
                ch = name.split("·", 1)[1] if "·" in name else "total"
                base = name.split("·", 1)[0]
                self.series.append(Series(name, spin[sp], g, y, "MO_KP", ch, base))
        self.status("原子 PDOS 完成")
        self.redraw()
        rows = DosProject.statistics(self.series, self.smin.value(), self.smax.value())
        self._rows = rows
        self.stats.setRowCount(len(rows))
        for i, r in enumerate(rows):
            set_row(self.stats, i, [r["name"], r["channel"], r["spin"], fmt(r["integral"], 3), fmt(r["center"]),
                                    fmt(r["width"]), fmt(r["skewness"], 3), fmt(r["kurtosis"], 3),
                                    fmt(r["filling"], 3), fmt(r["occ_center"]), fmt(r["upper_edge"])],
                    align_right_from=3)
        self.tabs.setCurrentIndex(0)

    def redraw(self):
        if not self.series:
            return
        s = self.session.settings
        st = self.plot.style(s.emin, s.emax, self.mirror.isChecked())
        data = [dict(label=x.name, x=x.x, y=x.y, spin=x.spin, color_key=x.name) for x in self.series]
        draw_lines(self.plot.figure, data, st,
                   ylabel="PDOS (e/eV)" if self.occ.isChecked() else "PDOS (states/eV)", fermi=0.0)
        self.plot.draw()

    def run_charges(self):
        s = self.session
        if s.system is None:
            error_box(self, "尚未读取数据", "请先在“k 点成键”页读取 MO_KP 文件")
            return

        def job(progress, cancel):
            return s.compute_populations(progress)
        self._start(job, self._charges_done, "正在计算布居…")

    def _charges_done(self, res):
        pm, pl = res
        sysm = self.session.system
        self.charges.setRowCount(len(sysm.d.atoms))
        for i, a in enumerate(sysm.d.atoms):
            vals = [a.index, a.element, fmt(a.z_eff, 1), fmt(pm.atom_pop[:, i].sum()), fmt(pm.charges[i]),
                    fmt(pl.atom_pop[:, i].sum()), fmt(pl.charges[i]),
                    fmt(pm.moments[i]) if pm.moments is not None else "—",
                    fmt(pl.moments[i]) if pl.moments is not None else "—"]
            set_row(self.charges, i, vals, align_right_from=2)
        rows = gross_populations(sysm, pm)
        self.gross.setRowCount(len(rows))
        for i, r in enumerate(rows):
            set_row(self.gross, i, [r[0], r[1], r[2], r[3], r[4], fmt(r[5]), fmt(r[6]) if r[6] is not None else "—",
                                    fmt(r[7])], align_right_from=5)
        self.tabs.setCurrentIndex(1)
        self.status("布居分析完成")
        self.log(f"[电荷] Mulliken 电荷和 = {pm.charges.sum():.5f}，Löwdin 电荷和 = {pl.charges.sum():.5f}")

    def export_pdos(self):
        if not self.series:
            return
        p, _ = QFileDialog.getSaveFileName(self, "导出 PDOS", str(Path(last_dir()) / "atom_pdos.csv"), "CSV (*.csv)")
        if p:
            remember_dir(p)
            write_series_csv(p, self.series)

    def export_stats(self):
        if not getattr(self, "_rows", None):
            return
        p, _ = QFileDialog.getSaveFileName(self, "导出统计", str(Path(last_dir()) / "atom_pdos_stats.csv"), "CSV (*.csv)")
        if p:
            remember_dir(p)
            write_stats_csv(p, self._rows)

    def export_charges(self):
        s = self.session
        if s.pop_m is None:
            return
        p, _ = QFileDialog.getSaveFileName(self, "导出电荷", str(Path(last_dir()) / "charges.csv"), "CSV (*.csv)")
        if not p:
            return
        remember_dir(p)
        import csv
        with open(p, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.writer(fh)
            w.writerow([self.charges.horizontalHeaderItem(j).text() for j in range(self.charges.columnCount())])
            for i in range(self.charges.rowCount()):
                w.writerow([self.charges.item(i, j).text() for j in range(self.charges.columnCount())])
