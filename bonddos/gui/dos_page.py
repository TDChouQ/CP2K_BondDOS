"""Page 1: native CP2K .dos/.pdos files (Γ or k-point)."""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QListWidget, QListWidgetItem,
                               QCheckBox, QLineEdit, QSplitter, QScrollArea, QFileDialog, QLabel, QGridLayout)

from ..dosanalysis import DosProject, DosOptions, write_series_csv, write_stats_csv
from ..dosfiles import scan_folder
from ..plotting import draw_lines
from .common import (PlotPanel, dspin, combo, make_table, set_row, fmt, group, form, last_dir, remember_dir,
                     error_box, resource_path)

REFS = [("费米能 E_F = 0", "fermi"), ("最高占据态 HOCO = 0", "hoco"), ("文件原零点", "file"), ("绝对能量", "absolute")]
SHAPES = [("Gaussian", "gaussian"), ("Lorentzian", "lorentzian"), ("Pseudo-Voigt", "voigt")]


class DosPage(QWidget):
    def __init__(self, log, parent=None):
        super().__init__(parent)
        self.log = log
        self.project = DosProject()
        self.series = []
        self.timer = QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.setInterval(250)
        self.timer.timeout.connect(self.update_plot)

        # ---- file list
        b_files = QPushButton("导入文件…")
        b_folder = QPushButton("导入文件夹…")
        b_clear = QPushButton("清空")
        b_demo = QPushButton("WO₃ 示例")
        b_files.clicked.connect(self.add_files)
        b_folder.clicked.connect(self.add_folder)
        b_clear.clicked.connect(self.clear)
        b_demo.clicked.connect(self.load_example)
        self.list = QListWidget()
        self.list.setMinimumHeight(150)
        self.list.itemChanged.connect(self.schedule)
        b_all = QPushButton("全选")
        b_none = QPushButton("全不选")
        b_all.clicked.connect(lambda: self._check_all(True))
        b_none.clicked.connect(lambda: self._check_all(False))
        fl = QVBoxLayout()
        r1 = QHBoxLayout()
        for b in (b_files, b_folder, b_demo, b_clear):
            r1.addWidget(b)
        fl.addLayout(r1)
        fl.addWidget(self.list)
        r2 = QHBoxLayout()
        r2.addWidget(b_all)
        r2.addWidget(b_none)
        fl.addLayout(r2)

        # ---- channels
        self.ch_boxes = {}
        grid = QGridLayout()
        for i, ch in enumerate(["total", "s", "p", "d", "f"]):
            cb = QCheckBox({"total": "总投影"}.get(ch, ch))
            cb.setChecked(ch == "total")
            cb.toggled.connect(self.schedule)
            self.ch_boxes[ch] = cb
            grid.addWidget(cb, i // 3, i % 3)
        self.components = QLineEdit()
        self.components.setPlaceholderText("磁量子数分量，如 dxy,dz2,px（需 COMPONENTS T）")
        self.components.editingFinished.connect(self.schedule)
        chl = QVBoxLayout()
        chl.addLayout(grid)
        chl.addWidget(self.components)

        # ---- options
        self.ref = combo(REFS)
        self.shift = dspin(0.0, -100, 100, 0.1, 3, " eV")
        self.fwhm = dspin(0.10, 0.0, 5, 0.02, 3, " eV")
        self.shape = combo(SHAPES)
        self.merge = QCheckBox("合并（求和）所选 KIND/文件")
        self.occupied = QCheckBox("按占据数加权（仅原始逐轨道数据）")
        self.idos = QCheckBox("显示积分 IDOS")
        self.mirror = QCheckBox("β 自旋向下镜像")
        self.mirror.setChecked(True)
        self.rebroaden = QCheckBox("对已展宽曲线再展宽")
        self.xmin = dspin(-10, -1000, 1000, 1, 2, " eV")
        self.xmax = dspin(8, -1000, 1000, 1, 2, " eV")
        self.smin = QLineEdit()
        self.smin.setPlaceholderText("全部")
        self.smax = QLineEdit()
        self.smax.setPlaceholderText("全部")
        for w in (self.ref, self.shape):
            w.currentIndexChanged.connect(self.schedule)
        for w in (self.shift, self.fwhm, self.xmin, self.xmax):
            w.valueChanged.connect(self.schedule)
        for w in (self.merge, self.occupied, self.idos, self.mirror, self.rebroaden):
            w.toggled.connect(self.schedule)
        for w in (self.smin, self.smax):
            w.editingFinished.connect(self.schedule)
        opts = form(("能量零点", self.ref), ("额外平移", self.shift), ("展宽 FWHM", self.fwhm), ("线型", self.shape),
                    ("显示能区", [self.xmin, "~", self.xmax]), ("统计能区", [self.smin, "~", self.smax]))
        ol = QVBoxLayout()
        ol.addLayout(opts)
        for w in (self.merge, self.occupied, self.idos, self.mirror, self.rebroaden):
            ol.addWidget(w)
        note = QLabel("FWHM 只用于原始逐轨道（MO）数据；CP2K 的 _curve 文件已按输入的 WIDTH 展宽。"
                      "能带中心等统计量对原始数据用离散权重精确计算，不受展宽影响。")
        note.setWordWrap(True)
        note.setStyleSheet("color:#666;")
        ol.addWidget(note)

        left = QWidget()
        ll = QVBoxLayout(left)
        ll.addWidget(group("数据文件 (.pdos / .dos)", fl))
        ll.addWidget(group("轨道通道", chl))
        ll.addWidget(group("能量与展宽", ol))
        ll.addStretch(1)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(left)
        scroll.setMinimumWidth(360)

        # ---- right side
        self.plot = PlotPanel()
        self.plot.redraw_cb = self.redraw
        self.stats = make_table(["名称", "通道", "自旋", "依据", "积分", "中心 ε (eV)", "宽度 (eV)", "偏度",
                                 "峰度", "填充率", "E≤0 积分", "E≤0 中心", "上边缘 (eV)"])
        b_csv = QPushButton("导出曲线 CSV…")
        b_stat = QPushButton("导出统计 CSV…")
        b_csv.clicked.connect(self.export_curves)
        b_stat.clicked.connect(self.export_stats)
        self.info = QLabel("")
        self.info.setStyleSheet("color:#555;")
        bl = QHBoxLayout()
        bl.addWidget(self.info, 1)
        bl.addWidget(b_csv)
        bl.addWidget(b_stat)
        right = QSplitter(Qt.Orientation.Vertical)
        right.addWidget(self.plot)
        w = QWidget()
        wl = QVBoxLayout(w)
        wl.setContentsMargins(0, 0, 0, 0)
        wl.addWidget(self.stats)
        wl.addLayout(bl)
        right.addWidget(w)
        right.setStretchFactor(0, 3)
        right.setStretchFactor(1, 1)

        split = QSplitter(Qt.Orientation.Horizontal)
        split.addWidget(scroll)
        split.addWidget(right)
        split.setStretchFactor(1, 1)
        lay = QHBoxLayout(self)
        lay.addWidget(split)

    # ------------------------------------------------------------------ files
    def add_files(self):
        files, _ = QFileDialog.getOpenFileNames(self, "选择 CP2K DOS/PDOS 文件", last_dir(),
                                                "CP2K DOS (*.pdos *.dos);;所有文件 (*)")
        if files:
            remember_dir(files[0])
            self._add(files)

    def add_folder(self):
        d = QFileDialog.getExistingDirectory(self, "选择包含 .pdos/.dos 的文件夹", last_dir())
        if d:
            remember_dir(d)
            f = scan_folder(d)
            files = f["pdos"] + f["dos"]
            if not files:
                error_box(self, "没有找到文件", "该文件夹中没有 .pdos 或 .dos 文件")
                return
            self._add(files)

    def load_example(self):
        folder = Path(resource_path("examples/WO3_gamma_pdos"))
        self._add(sorted(folder.glob("*.pdos")))
        self.ref.setCurrentIndex(0)
        for ch in ("total",):
            self.ch_boxes[ch].setChecked(True)

    def _add(self, files):
        self.project.errors.clear()
        new = self.project.add_files(files)
        self.list.blockSignals(True)
        base = len(self.project.data) - len(new)
        for j, d in enumerate(new):
            mode = {"sticks": "逐轨道", "curve": "曲线", "bins": "直方图"}[d.mode]
            kp = " · k点" if d.kpoint else ""
            it = QListWidgetItem(f"{d.name}  |  {d.label}  |  {mode}{kp}")
            it.setFlags(it.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            # raw files and their _curve partners should not be summed together
            it.setCheckState(Qt.CheckState.Checked)
            it.setData(Qt.ItemDataRole.UserRole, base + j)
            it.setToolTip("\n".join(d.header[:8]))
            self.list.addItem(it)
        self.list.blockSignals(False)
        for e in self.project.errors:
            self.log(f"[DOS] 读取失败 {e}")
        if self.project.errors:
            error_box(self, "部分文件无法读取", "\n".join(self.project.errors[:10]))
        if new:
            self.log(f"[DOS] 导入 {len(new)} 组数据")
            has_curve = any(d.mode == "curve" for d in self.project.data)
            has_raw = any(d.mode == "sticks" for d in self.project.data)
            if has_curve and has_raw:
                self.info.setText("提示：同时导入了原始文件和 _curve 文件，合并时不要同时勾选同一 KIND 的两种文件")
        self.update_plot()

    def clear(self):
        self.project.clear()
        self.list.clear()
        self.series = []
        self.stats.setRowCount(0)
        self.plot.figure.clear()
        self.plot.draw()

    def _check_all(self, on):
        self.list.blockSignals(True)
        for i in range(self.list.count()):
            self.list.item(i).setCheckState(Qt.CheckState.Checked if on else Qt.CheckState.Unchecked)
        self.list.blockSignals(False)
        self.update_plot()

    def selected(self):
        out = []
        for i in range(self.list.count()):
            it = self.list.item(i)
            if it.checkState() == Qt.CheckState.Checked:
                out.append(self.project.data[it.data(Qt.ItemDataRole.UserRole)])
        return out

    # ------------------------------------------------------------------ analysis
    def schedule(self, *a):
        self.timer.start()

    def options(self):
        chans = [c for c, cb in self.ch_boxes.items() if cb.isChecked()]
        extra = [x.strip() for x in self.components.text().replace("，", ",").split(",") if x.strip()]
        alias = {"d-2": "dxy", "d-1": "dyz", "d0": "dz2", "d+1": "dxz", "d+2": "dx2-y2", "dx2": "dx2-y2"}
        chans += [alias.get(x.lower(), x.lower()) for x in extra]
        return DosOptions(reference=self.ref.currentData(), shift=self.shift.value(), fwhm=self.fwhm.value(),
                          shape=self.shape.currentData(), channels=tuple(chans), merge=self.merge.isChecked(),
                          occupied=self.occupied.isChecked(), integrate=self.idos.isChecked(),
                          rebroaden_curves=self.rebroaden.isChecked())

    @staticmethod
    def _num(edit):
        try:
            t = edit.text().strip()
            return float(t) if t else None
        except ValueError:
            return None

    def update_plot(self):
        items = self.selected()
        if not items:
            self.series = []
            self.plot.figure.clear()
            self.plot.draw()
            self.stats.setRowCount(0)
            return
        try:
            self.series, warns = self.project.build_series(items, self.options())
        except Exception as e:
            self.info.setText(f"错误: {e}")
            return
        msg = []
        if warns:
            msg.append(warns[0])
        present = set(self.project.available_channels(items))
        missing = [c for c in self.options().channels if c not in present]
        if missing:
            msg.append("文件中没有通道: " + ", ".join(missing))
        self.info.setText("；".join(msg))
        self.redraw()
        self.fill_stats()

    def redraw(self):
        st = self.plot.style(self.xmin.value(), self.xmax.value(), self.mirror.isChecked())
        data = []
        for s in self.series:
            data.append(dict(label=f"{s.name}·{s.channel}" if s.channel != "total" else s.name,
                             color_key=f"{s.kind}·{s.channel}", x=s.x, y=s.y, spin=s.spin))
        ylabel = "IDOS (states)" if self.idos.isChecked() else (
            "Occupied DOS (e/eV)" if self.occupied.isChecked() else "DOS (states/eV)")
        ref = self.ref.currentData()
        st.energy_label = {"fermi": "E − E$_F$ (eV)", "hoco": "E − E$_{HOCO}$ (eV)", "file": "E (eV)",
                           "absolute": "E (eV)"}[ref]
        try:
            draw_lines(self.plot.figure, data, st, ylabel=ylabel,
                       fermi=(0.0 + self.shift.value()) if ref in ("fermi", "hoco") else None)
        except Exception as e:
            self.info.setText(f"绘图失败: {e}")
        self.plot.draw()

    def fill_stats(self):
        if self.idos.isChecked():
            self.stats.setRowCount(0)
            return
        rows = self.project.statistics(self.series, self._num(self.smin), self._num(self.smax))
        self._rows = rows
        self.stats.setRowCount(len(rows))
        for i, r in enumerate(rows):
            set_row(self.stats, i, [r["name"], r["channel"], r["spin"], r["basis"], fmt(r["integral"], 3),
                                    fmt(r["center"]), fmt(r["width"]), fmt(r["skewness"], 3), fmt(r["kurtosis"], 3),
                                    fmt(r["filling"], 3), fmt(r["occ_integral"], 3), fmt(r["occ_center"]),
                                    fmt(r["upper_edge"])], align_right_from=4)

    # ------------------------------------------------------------------ export
    def export_curves(self):
        if not self.series:
            return
        p, _ = QFileDialog.getSaveFileName(self, "导出曲线", str(Path(last_dir()) / "dos_curves.csv"), "CSV (*.csv)")
        if p:
            remember_dir(p)
            try:
                write_series_csv(p, self.series)
                self.log(f"[DOS] 曲线已导出 {p}")
            except Exception as e:
                error_box(self, "导出失败", str(e))

    def export_stats(self):
        if not getattr(self, "_rows", None):
            return
        p, _ = QFileDialog.getSaveFileName(self, "导出统计", str(Path(last_dir()) / "band_center.csv"), "CSV (*.csv)")
        if p:
            remember_dir(p)
            try:
                write_stats_csv(p, self._rows)
                self.log(f"[DOS] 统计已导出 {p}")
            except Exception as e:
                error_box(self, "导出失败", str(e))
