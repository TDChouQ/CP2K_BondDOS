"""Page 2: k-point COHP / COOP / COBI from MO_KP + KS_CSR_WRITE."""
from __future__ import annotations

import html
from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLineEdit, QSplitter, QScrollArea,
                               QFileDialog, QLabel, QCheckBox, QTextBrowser, QTableWidgetItem, QComboBox,
                               QProgressBar, QTabWidget, QAbstractItemView, QMessageBox)

from ..bonding import QLABEL
from ..mokp import find_ks_sets
from ..session import BondSession, BondRequest, write_config, parse_config
from ..constants import FWHM_PER_SIGMA
from ..plotting import draw_cohp
from .common import (PlotPanel, Worker, dspin, combo, make_table, set_row, fmt, group, form, last_dir,
                     remember_dir, error_box)

BASES = [("直接 Gaussian AO（默认）", "ao"), ("Löwdin 正交化（类 LOBSTER）", "lowdin")]
RESOLUTIONS = [("不分解", "none"), ("按角动量 s/p/d/f", "l"), ("按磁量子数 (px, dxy…)", "lm"),
               ("按基组壳层", "shell"), ("逐个 AO", "ao")]
SPINS = [("自动", "auto"), ("只算 α/up", "up"), ("只算 β/down", "down")]
EREFS = [("自动（绝缘体 HOCO；金属用占据数估计 E_F）", "auto"), ("最高占据态 HOCO/VBM", "hoco"),
         ("带隙中点", "midgap"), ("占据数估计的费米能", "fermi"), ("手动输入绝对能量 (eV)", "manual")]
PERIODIC = ["XYZ", "XY", "XZ", "YZ", "X", "Y", "Z", "NONE"]


class CohpPage(QWidget):
    loaded = Signal()

    def __init__(self, session: BondSession, log, status, parent=None):
        super().__init__(parent)
        self.session = session
        self.log = log
        self.status = status
        self.worker = None

        # ---------------- inputs
        self.mokp = QLineEdit()
        self.mokp.setPlaceholderText("cohp_mo-1_0.mokp")
        b_mokp = QPushButton("浏览…")
        b_mokp.clicked.connect(self.pick_mokp)
        self.ksdir = QLineEdit()
        self.ksdir.setPlaceholderText("KS_CSR_WRITE 输出的 *.csr 所在文件夹（可选）")
        b_ks = QPushButton("浏览…")
        b_ks.clicked.connect(self.pick_ksdir)
        self.ksset = QComboBox()
        self.ksdir.editingFinished.connect(self.refresh_sets)
        self.periodic = combo(PERIODIC)
        self.triangle = combo([("完整矩阵 (UPPER_TRIANGULAR F)", "FULL"), ("上三角 (UPPER_TRIANGULAR T)", "UPPER")])
        self.electrons = QLineEdit()
        self.electrons.setPlaceholderText("可选，用于核对")
        b_load = QPushButton("读取并检查数据")
        b_load.setObjectName("primary")
        b_load.clicked.connect(self.load)
        b_demo = QPushButton("生成演示数据…")
        b_demo.clicked.connect(self.make_demo)
        b_imp = QPushButton("导入 cohp.in…")
        b_imp.clicked.connect(self.import_config)
        inp = form(("MO_KP 文件", [self.mokp, b_mokp]), ("KS 文件夹", [self.ksdir, b_ks]), ("KS 文件组", self.ksset),
                   ("周期方向", self.periodic), ("CSR 存储", self.triangle), ("预期电子数", self.electrons))
        il = QVBoxLayout()
        il.addLayout(inp)
        h = QHBoxLayout()
        h.addWidget(b_load)
        h.addWidget(b_imp)
        h.addWidget(b_demo)
        il.addLayout(h)

        # ---------------- settings
        self.basis = combo(BASES)
        self.q_cohp = QCheckBox("COHP")
        self.q_coop = QCheckBox("COOP")
        self.q_cobi = QCheckBox("COBI")
        for c in (self.q_cohp, self.q_coop, self.q_cobi):
            c.setChecked(True)
        self.resolution = combo(RESOLUTIONS, 1)
        self.spins = combo(SPINS)
        self.fwhm = dspin(0.2355, 0.001, 5, 0.05, 4, " eV")
        self.sigma_lbl = QLabel()
        self.fwhm.valueChanged.connect(lambda v: self.sigma_lbl.setText(f"σ = {v / FWHM_PER_SIGMA:.4f} eV"))
        self.fwhm.valueChanged.emit(self.fwhm.value())
        self.emin = dspin(-15, -500, 500, 1, 2, " eV")
        self.emax = dspin(10, -500, 500, 1, 2, " eV")
        self.de = dspin(0.01, 0.0005, 1, 0.005, 4, " eV")
        self.eref = combo(EREFS)
        self.eref_val = QLineEdit()
        self.eref_val.setPlaceholderText("绝对能量 eV，如 CP2K 输出的 Fermi energy")
        st = form(("COHP 表象", self.basis), ("计算的量", [self.q_cohp, self.q_coop, self.q_cobi]),
                  ("轨道分解", self.resolution), ("自旋", self.spins), ("展宽 FWHM", [self.fwhm, self.sigma_lbl]),
                  ("能窗", [self.emin, "~", self.emax]), ("能量步长", self.de), ("能量零点", self.eref),
                  ("手动零点", self.eref_val))

        # ---------------- bonds
        self.bond_edit = QLineEdit()
        self.bond_edit.setPlaceholderText("例: 1 2 | 1 2 ALL | 1 2 -1 0 0 | 1 2 PAIR，末尾加 ORB 做轨道分解")
        b_add = QPushButton("添加")
        b_add.clicked.connect(self.add_bond)
        self.bond_edit.returnPressed.connect(self.add_bond)
        self.bonds = make_table(["键定义"])
        self.bonds.setMaximumHeight(130)
        b_del = QPushButton("删除所选")
        b_del.clicked.connect(lambda: [self.bonds.removeRow(r.row()) for r in
                                       sorted(self.bonds.selectionModel().selectedRows(), key=lambda x: -x.row())])
        self.auto = QCheckBox("自动搜索距离范围内的所有键（类 LOBSTER cohpGenerator）")
        self.dmin = dspin(1.0, 0, 20, 0.1, 2, " Å")
        self.dmax = dspin(3.0, 0, 20, 0.1, 2, " Å")
        self.el_a = QLineEdit()
        self.el_a.setPlaceholderText("元素A 如 W（可空）")
        self.el_b = QLineEdit()
        self.el_b.setPlaceholderText("元素B 如 O（可空）")
        self.auto_orb = QCheckBox("自动键也做轨道分解")
        bl = QVBoxLayout()
        h = QHBoxLayout()
        h.addWidget(self.bond_edit, 1)
        h.addWidget(b_add)
        h.addWidget(b_del)
        bl.addLayout(h)
        bl.addWidget(self.bonds)
        bl.addWidget(self.auto)
        bl.addLayout(form(("距离", [self.dmin, "~", self.dmax]), ("元素", [self.el_a, "–", self.el_b])))
        bl.addWidget(self.auto_orb)
        hint = QLabel("原子编号与 CP2K &COORD 顺序一致（从 1 开始）。平移 R 以 CP2K 打印的原始坐标为准；"
                      "“ALL” 列出所有等距最近镜像；“PAIR” 把 A–B 的全部周期镜像合计。")
        hint.setWordWrap(True)
        hint.setStyleSheet("color:#666;")
        bl.addWidget(hint)

        self.b_run = QPushButton("开始计算")
        self.b_run.setObjectName("primary")
        self.b_run.clicked.connect(self.run)
        self.b_cancel = QPushButton("取消")
        self.b_cancel.clicked.connect(self.cancel)
        self.b_cancel.setEnabled(False)
        b_export = QPushButton("导出全部结果…")
        b_export.clicked.connect(self.export)
        b_save = QPushButton("另存 cohp.in…")
        b_save.clicked.connect(self.save_config)
        self.progress = QProgressBar()
        self.progress.setRange(0, 1000)
        self.progress.setTextVisible(False)
        self.progress.setMaximumHeight(8)
        rl = QVBoxLayout()
        h = QHBoxLayout()
        h.addWidget(self.b_run)
        h.addWidget(self.b_cancel)
        rl.addLayout(h)
        h = QHBoxLayout()
        h.addWidget(b_export)
        h.addWidget(b_save)
        rl.addLayout(h)
        rl.addWidget(self.progress)

        left = QWidget()
        ll = QVBoxLayout(left)
        ll.addWidget(group("1 · 数据", il))
        ll.addWidget(group("2 · 分析设置", st))
        ll.addWidget(group("3 · 目标键", bl))
        ll.addLayout(rl)
        ll.addStretch(1)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(left)
        scroll.setMinimumWidth(430)

        # ---------------- right
        self.summary = QTextBrowser()
        self.summary.setOpenExternalLinks(True)
        self.results = make_table(["#", "键", "长度 (Å)", "R", "模式", "−ICOHP (eV)", "ICOOP", "ICOBI", "状态"])
        self.results.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.results.itemSelectionChanged.connect(self.on_select)
        self.results.setSortingEnabled(False)
        self.orbtable = make_table(["轨道 A", "轨道 B", "积分"])

        self.plot = PlotPanel(vertical_default=True)
        self.plot.redraw_cb = self.redraw
        self.q_plot = combo([("COHP", "cohp"), ("COOP", "coop"), ("COBI", "cobi")])
        self.negate = QCheckBox("−COHP（成键为正）")
        self.negate.setChecked(True)
        self.show_int = QCheckBox("积分曲线")
        self.show_int.setChecked(True)
        self.sum_sel = QCheckBox("所选键求和")
        self.spin_plot = combo([("全部自旋", "all"), ("up/α", 0), ("down/β", 1), ("up+down", "sum")])
        self.channel = QComboBox()
        self.channel.addItem("总贡献", None)
        for w in (self.q_plot, self.spin_plot, self.channel):
            w.currentIndexChanged.connect(self.redraw)
        for w in (self.negate, self.show_int, self.sum_sel):
            w.toggled.connect(self.redraw)
        pc = QHBoxLayout()
        for w in (QLabel("绘制"), self.q_plot, self.negate, self.show_int, self.sum_sel, self.spin_plot,
                  QLabel("轨道对"), self.channel):
            pc.addWidget(w)
        pc.addStretch(1)
        pw = QWidget()
        pl = QVBoxLayout(pw)
        pl.setContentsMargins(0, 0, 0, 0)
        pl.addLayout(pc)
        pl.addWidget(self.plot, 1)

        tabs = QTabWidget()
        tabs.addTab(self.summary, "数据概况与检查")
        tabs.addTab(self.results, "成键结果")
        tabs.addTab(self.orbtable, "轨道分解积分")
        self.tabs = tabs
        right = QSplitter(Qt.Orientation.Vertical)
        right.addWidget(pw)
        right.addWidget(tabs)
        right.setStretchFactor(0, 3)
        right.setStretchFactor(1, 2)

        split = QSplitter(Qt.Orientation.Horizontal)
        split.addWidget(scroll)
        split.addWidget(right)
        split.setStretchFactor(1, 1)
        lay = QHBoxLayout(self)
        lay.addWidget(split)
        self.summary.setHtml(WELCOME)

    # ------------------------------------------------------------------ inputs
    def pick_mokp(self):
        p, _ = QFileDialog.getOpenFileName(self, "选择 MO_KP 文件", last_dir(), "MO_KP (*.mokp);;所有文件 (*)")
        if p:
            remember_dir(p)
            self.mokp.setText(p)
            if not self.ksdir.text():
                self.ksdir.setText(str(Path(p).parent))
            self.refresh_sets()

    def pick_ksdir(self):
        d = QFileDialog.getExistingDirectory(self, "选择 KS CSR 文件夹", last_dir())
        if d:
            remember_dir(d)
            self.ksdir.setText(d)
            self.refresh_sets()

    def refresh_sets(self):
        self.ksset.clear()
        d = self.ksdir.text().strip()
        if not d:
            return
        sets = find_ks_sets(d)
        for lab, files in sets.items():
            self.ksset.addItem(f"{lab}  ({len(files)} 个文件)", lab)
        if not sets:
            self.ksset.addItem("（未找到 KS_SPIN_*_K_*.csr）", None)

    def _apply_settings(self):
        s = self.session.settings
        s.basis = self.basis.currentData()
        qs = [q for q, c in (("cohp", self.q_cohp), ("coop", self.q_coop), ("cobi", self.q_cobi)) if c.isChecked()]
        if not qs:
            raise ValueError("请至少选择一种量（COHP/COOP/COBI）")
        s.quantities = tuple(qs)
        s.resolution = self.resolution.currentData()
        s.spins = self.spins.currentData()
        s.fwhm = self.fwhm.value()
        s.emin, s.emax, s.de = self.emin.value(), self.emax.value(), self.de.value()
        if s.emax <= s.emin:
            raise ValueError("能窗上限必须大于下限")
        s.eref_mode = self.eref.currentData()
        if s.eref_mode == "manual":
            try:
                s.eref_value = float(self.eref_val.text())
            except ValueError:
                raise ValueError("手动零点需要填写数值（绝对能量，eV）")

    def _start(self, fn, done, busy_text):
        if self.worker is not None and self.worker.isRunning():
            return
        self.worker = Worker(fn)
        self.worker.progressed.connect(self._on_progress)
        self.worker.done.connect(done)
        self.worker.failed.connect(self._failed)
        self.worker.finished.connect(self._on_finished)
        self.b_run.setEnabled(False)
        self.b_cancel.setEnabled(True)
        self.status(busy_text)
        self.worker.start()

    def _on_progress(self, f, t):
        self.progress.setValue(int(f * 1000))
        self.status(t)

    def _on_finished(self):
        self.b_run.setEnabled(True)
        self.b_cancel.setEnabled(False)

    def _failed(self, msg, tb):
        self.progress.setValue(0)
        self.status("出错")
        self.log(f"[COHP] 错误: {msg}\n{tb}")
        error_box(self, "处理失败", msg, tb)

    def cancel(self):
        if self.worker is not None:
            self.worker.cancel()

    def load(self):
        mokp = self.mokp.text().strip()
        if not mokp or not Path(mokp).is_file():
            error_box(self, "缺少文件", "请选择 MO_KP (.mokp) 文件")
            return
        try:
            self._apply_settings()
        except ValueError as e:
            error_box(self, "设置有误", str(e))
            return
        ks = self.ksdir.text().strip() or None
        kset = self.ksset.currentData()
        el = self.electrons.text().strip()
        try:
            electrons = float(el) if el else None
        except ValueError:
            electrons = None
        per = self.periodic.currentText()
        tri = self.triangle.currentData()

        def job(progress, cancel):
            self.session.load(mokp, ks if kset else None, kset, tri, per, electrons, progress)
            return True
        self._start(job, self._loaded, "正在读取 MO_KP / KS 文件…")

    def _loaded(self, _):
        self.progress.setValue(1000)
        self.status("数据已读取")
        self.show_summary()
        self.results.setRowCount(0)
        self.log(f"[COHP] 已读取 {self.mokp.text()}")
        self.loaded.emit()

    def show_summary(self):
        s = self.session
        info = s.summary()
        rows = []
        rows.append(f"<h3>{html.escape(info['formula'])} · {info['atoms']} 原子</h3>")
        gap = "—" if info["gap"] is None else f"{info['gap']:.4f} eV"
        rows.append("<table cellpadding=3>")
        for k, v in [("AO / MO", f"{info['nao']} / {info['nmo']}"), ("k 点", f"{info['nkp']}（网格 {info['mesh']}，"
                     f"{'完整等权网格' if info['full_mesh'] else '<b style=color:#c0392b>非完整网格</b>'}）"),
                     ("自旋", "限制性 (RKS)" if info["nspins"] == 1 else "非限制性 (UKS)"),
                     ("哈密顿量", info["H"]), ("重叠矩阵 S(k)", "有" if info["S"] else "无"),
                     ("电子数 Σwf", f"{info['electrons']:.6f}"),
                     ("VBM / CBM / 带隙", f"{fmt(info['vbm'])} / {fmt(info['cbm'])} eV / {gap}"),
                     ("部分占据（金属/展宽）", "是" if info["partial"] else "否"),
                     ("能量零点", f"{info['eref']:.6f} eV — {html.escape(info['eref_note'])}")]:
            rows.append(f"<tr><td><b>{k}</b></td><td>{v}</td></tr>")
        rows.append("</table><h4>数据一致性检查</h4><table cellpadding=3>")
        for c in s.system.checks:
            mark = "<span style='color:#1a7f37'>✔</span>" if c.ok else "<span style='color:#c0392b'>✘</span>"
            note = "" if c.ok else f"<br><span style='color:#c0392b'>{html.escape(c.note)}</span>"
            rows.append(f"<tr><td>{mark}</td><td>{html.escape(c.name)}</td><td>{c.value:.3e}{note}</td></tr>")
        rows.append("</table>")
        g = s.gauge
        if g.get("available"):
            ok = g["consistent"]
            rows.append(f"<p>{'✔' if ok else '✘'} 周期镜像相位约定检查（重叠矩阵随距离衰减）："
                        f"score(+)={g['score_plus']:.3g}，score(−)={g['score_minus']:.3g}"
                        f"{'' if ok else '<br><b style=color:#c0392b>与预期相反，请勿使用 IMAGE 结果并联系开发者</b>'}</p>")
        else:
            rows.append(f"<p>周期镜像相位检查：{html.escape(g.get('reason', '不可用'))}</p>")
        for w in info["warnings"] + info["mesh_messages"]:
            rows.append(f"<p style='color:#b35900'>⚠ {html.escape(w)}</p>")
        self.summary.setHtml("\n".join(rows))
        self.tabs.setCurrentIndex(0)

    # ------------------------------------------------------------------ bonds
    def add_bond(self):
        t = self.bond_edit.text().strip()
        if not t:
            return
        try:
            BondRequest.parse(t)
        except Exception as e:
            error_box(self, "键定义有误", str(e))
            return
        r = self.bonds.rowCount()
        self.bonds.insertRow(r)
        self.bonds.setItem(r, 0, QTableWidgetItem(t))
        self.bond_edit.clear()

    def bond_texts(self):
        return [self.bonds.item(i, 0).text() for i in range(self.bonds.rowCount()) if self.bonds.item(i, 0)]

    def auto_spec(self):
        if not self.auto.isChecked():
            return None
        if self.dmax.value() <= self.dmin.value():
            raise ValueError("自动搜索的最大距离必须大于最小距离")
        ea = [x.strip() for x in self.el_a.text().replace("，", ",").split(",") if x.strip()]
        eb = [x.strip() for x in self.el_b.text().replace("，", ",").split(",") if x.strip()]
        return {"dmin": self.dmin.value(), "dmax": self.dmax.value(), "elements_a": ea or None,
                "elements_b": eb or None, "resolve": self.auto_orb.isChecked()}

    def run(self):
        if self.session.system is None:
            error_box(self, "尚未读取数据", "请先点击“读取并检查数据”")
            return
        try:
            self._apply_settings()
            reqs = [BondRequest.parse(t) for t in self.bond_texts()]
            auto = self.auto_spec()
        except Exception as e:
            error_box(self, "设置有误", str(e))
            return
        if not reqs and not auto:
            error_box(self, "没有目标键", "请添加至少一根键，或勾选自动搜索")
            return

        def job(progress, cancel):
            return self.session.compute(reqs, auto, progress, cancel)
        self._start(job, self._computed, "正在计算成键贡献…")

    def _computed(self, results):
        self.progress.setValue(1000)
        self.status(f"完成：{len(results)} 根键")
        self.log(f"[COHP] 计算完成 {len(results)} 根键；能量零点 {self.session.eref:.4f} eV")
        self.fill_results()
        self.tabs.setCurrentIndex(1)
        if results:
            self.results.selectRow(0)

    def fill_results(self):
        res = self.session.results
        self.results.setRowCount(len(res))
        for i, r in enumerate(res):
            b = r.bond
            I = r.integrals
            vals = [i + 1, b.label, f"{b.length:.4f}", " ".join(str(x) for x in b.r_raw) if b.mode == "IMAGE" else "—",
                    b.mode, fmt(-I["cohp"]["total"]) if "cohp" in I else "—",
                    fmt(I["coop"]["total"]) if "coop" in I else "—", fmt(I["cobi"]["total"]) if "cobi" in I else "—",
                    {"ok": "✔", "nyquist": "⚠ Nyquist", "gamma": "⚠ Γ叠加", "aliased": "✘ 混叠",
                     "nonperiodic": "✘ 非周期"}.get(r.status, r.status)]
            set_row(self.results, i, vals, align_right_from=2)
            tip = "\n".join([f"{QLABEL[q]}: " + ", ".join(f"{k}={v:.5f}" for k, v in d.items()) for q, d in I.items()]
                            + r.notes)
            for j in range(self.results.columnCount()):
                it = self.results.item(i, j)
                if it:
                    it.setToolTip(tip)

    def selected_results(self):
        rows = sorted({i.row() for i in self.results.selectedIndexes()})
        return [self.session.results[r] for r in rows if r < len(self.session.results)]

    def on_select(self):
        sel = self.selected_results()
        self.channel.blockSignals(True)
        self.channel.clear()
        self.channel.addItem("总贡献", None)
        if len(sel) == 1 and sel[0].channels_a:
            r = sel[0]
            q = self.q_plot.currentData()
            for ia, (na, _) in enumerate(r.channels_a):
                for ib, (nb, _) in enumerate(r.channels_b):
                    self.channel.addItem(f"{na} – {nb}", (ia, ib))
            self.fill_orbitals(r)
        else:
            self.orbtable.setRowCount(0)
        self.channel.blockSignals(False)
        self.redraw()

    def fill_orbitals(self, r):
        q = self.q_plot.currentData()
        if q not in r.resolved_integrals:
            q = next(iter(r.resolved_integrals))
        rows = []
        for sp, M in r.resolved_integrals[q].items():
            for ia, (na, _) in enumerate(r.channels_a):
                for ib, (nb, _) in enumerate(r.channels_b):
                    rows.append((na, nb, sp, M[ia, ib]))
        agg = {}
        for na, nb, sp, v in rows:
            agg[(na, nb)] = agg.get((na, nb), 0.0) + v
        items = sorted(agg.items(), key=lambda kv: kv[1] if q == "cohp" else -kv[1])
        self.orbtable.setRowCount(len(items))
        self.orbtable.setHorizontalHeaderLabels(["轨道 A", "轨道 B", f"I{QLABEL[q]}（自旋求和）"])
        for i, ((na, nb), v) in enumerate(items):
            set_row(self.orbtable, i, [na, nb, fmt(v, 5)], align_right_from=2)

    def redraw(self):
        sel = self.selected_results()
        s = self.session
        if not sel or s.engine is None:
            self.plot.figure.clear()
            self.plot.draw()
            return
        q = self.q_plot.currentData()
        if q not in s.engine.quantities:
            self.plot.figure.clear()
            ax = self.plot.figure.add_subplot(111)
            ax.text(0.5, 0.5, f"{QLABEL[q]} 未计算（COHP 需要 KS 文件；COOP/COBI 需要 S(k)）", ha="center",
                    va="center", transform=ax.transAxes)
            self.plot.draw()
            return
        ch = self.channel.currentData()
        spin_sel = self.spin_plot.currentData()
        spin_name = {0: "up", 1: "down"} if s.system.ns == 2 else {0: None}
        curves, cumul = [], []
        try:
            if self.sum_sel.isChecked() and len(sel) > 1:
                g, c = s.sum_curves([r for r in sel if r.status not in ("aliased", "nonperiodic")], q)
                _, ci = s.sum_curves([r for r in sel if r.status not in ("aliased", "nonperiodic")], q, cumulative=True)
                groups = [("Σ selected", c, ci)]
            else:
                groups = []
                for r in sel[:12]:
                    if r.status in ("aliased", "nonperiodic"):
                        continue
                    use_ch = ch if (len(sel) == 1 and r.channels_a) else None
                    g, c = s.curves(r, q, channel=use_ch)
                    _, ci = s.curves(r, q, channel=use_ch, cumulative=True)
                    lab = r.bond.label + ("" if use_ch is None else f" {self.channel.currentText()}")
                    groups.append((lab, c, ci))
            for lab, c, ci in groups:
                if spin_sel == "sum" or s.system.ns == 1:
                    curves.append(dict(label=lab, x=g, y=sum(c.values()), spin="total"))
                    cumul.append(dict(label=lab, x=g, y=sum(ci.values()), spin="total"))
                else:
                    for sp in c:
                        if spin_sel != "all" and sp != spin_sel:
                            continue
                        curves.append(dict(label=lab, x=g, y=c[sp], spin=spin_name[sp]))
                        cumul.append(dict(label=lab, x=g, y=ci[sp], spin=spin_name[sp]))
            st = self.plot.style(s.settings.emin, s.settings.emax)
            draw_cohp(self.plot.figure, curves, cumul, st, QLABEL[q], self.negate.isChecked(),
                      self.show_int.isChecked())
        except Exception as e:
            self.log(f"[COHP] 绘图失败: {e}")
        self.plot.draw()

    # ------------------------------------------------------------------ export / config
    def export(self):
        if self.session.system is None:
            return
        d = QFileDialog.getExistingDirectory(self, "选择导出文件夹", last_dir())
        if not d:
            return
        remember_dir(d)

        self._export_dir = d

        def job(progress, cancel):
            return self.session.export(d, "CP2K", True, True, progress)
        self._start(job, self._exported, "正在导出…")

    def _exported(self, files):
        d = self._export_dir
        self.progress.setValue(1000)
        self.status(f"已导出 {len(files)} 个文件")
        self.log("[COHP] 导出:\n  " + "\n  ".join(str(f) for f in files))
        QMessageBox.information(self, "导出完成", f"已写出 {len(files)} 个文件到\n{d}\n\n"
                                "包含 CSV 表格、JSON 报告，以及 COHPCAR/ICOHPLIST/COBICAR/ICOBILIST/"
                                "COOPCAR/ICOOPLIST/CHARGE/GROSSPOP/DOSCAR.lobster 与 POSCAR。")

    def save_config(self):
        if self.session.system is None:
            return
        p, _ = QFileDialog.getSaveFileName(self, "保存 cohp.in", str(Path(last_dir()) / "cohp.in"), "cohp.in (*.in);;所有文件 (*)")
        if p:
            try:
                self._apply_settings()
                auto = self.auto_spec()
                write_config(p, self.session, self.bond_texts(), auto)
                self.log(f"[COHP] 配置已保存 {p}")
            except Exception as e:
                error_box(self, "保存失败", str(e))

    def import_config(self):
        p, _ = QFileDialog.getOpenFileName(self, "导入 cohp.in", last_dir(), "cohp.in (*.in);;所有文件 (*)")
        if not p:
            return
        try:
            cfg = parse_config(p)
        except Exception as e:
            error_box(self, "读取失败", str(e))
            return
        base = Path(p).parent

        def rel(x):
            q = Path(x)
            return str(q if q.is_absolute() else (base / q).resolve())
        if "MOKP" in cfg:
            self.mokp.setText(rel(cfg["MOKP"]))
        if "KS_DIR" in cfg:
            self.ksdir.setText(rel(cfg["KS_DIR"]))
            self.refresh_sets()
            if cfg.get("KS_PATTERN"):
                i = self.ksset.findData(cfg["KS_PATTERN"])
                if i >= 0:
                    self.ksset.setCurrentIndex(i)
        if "PERIODIC" in cfg:
            self.periodic.setCurrentText(cfg["PERIODIC"].upper())
        self.basis.setCurrentIndex(1 if cfg.get("COHP_BASIS", "AO").upper().startswith("L") else 0)
        res = cfg.get("ORBITAL_RESOLUTION", "L").lower()
        i = self.resolution.findData(res)
        if i >= 0:
            self.resolution.setCurrentIndex(i)
        if "FWHM_EV" in cfg:
            self.fwhm.setValue(float(cfg["FWHM_EV"]))
        elif "SIGMA_EV" in cfg:
            self.fwhm.setValue(float(cfg["SIGMA_EV"]) * FWHM_PER_SIGMA)
        for key, w in (("EMIN_EV", self.emin), ("EMAX_EV", self.emax), ("DE_EV", self.de)):
            if key in cfg:
                w.setValue(float(cfg[key]))
        ef = cfg.get("EF_EV", "auto").lower()
        i = self.eref.findData(ef)
        if i >= 0:
            self.eref.setCurrentIndex(i)
        else:
            self.eref.setCurrentIndex(self.eref.findData("manual"))
            self.eref_val.setText(ef)
        self.bonds.setRowCount(0)
        for b in cfg["BOND"]:
            r = self.bonds.rowCount()
            self.bonds.insertRow(r)
            self.bonds.setItem(r, 0, QTableWidgetItem(b))
        if "AUTO_BONDS" in cfg:
            t = cfg["AUTO_BONDS"].split()
            self.auto.setChecked(True)
            self.dmin.setValue(float(t[0]))
            self.dmax.setValue(float(t[1]))
        self.log(f"[COHP] 已导入配置 {p}（请点击“读取并检查数据”）")

    def make_demo(self):
        d = QFileDialog.getExistingDirectory(self, "选择一个空文件夹存放演示数据", last_dir())
        if not d:
            return
        from ..demo import generate_demo
        folder = Path(d) / "BondDOS_demo"
        try:
            generate_demo(folder, uks=False)
        except Exception as e:
            error_box(self, "生成失败", str(e))
            return
        self.mokp.setText(str(folder / "cohp_mo-1_0.mokp"))
        self.ksdir.setText(str(folder))
        self.refresh_sets()
        self.bonds.setRowCount(0)
        for t in ("1 2 ORB", "2 3 ALL", "1 2 PAIR"):
            r = self.bonds.rowCount()
            self.bonds.insertRow(r)
            self.bonds.setItem(r, 0, QTableWidgetItem(t))
        self.auto.setChecked(True)
        self.dmin.setValue(1.0)
        self.dmax.setValue(2.6)
        self.emin.setValue(-20)
        self.log(f"[COHP] 已生成合成演示数据 {folder}（不是 CP2K 计算结果）")
        self.load()


WELCOME = """
<h3>k 点 COHP / COOP / COBI 分析</h3>
<p>需要 CP2K 2026.2+ 单点计算的两类输出：</p>
<ol>
<li><b>&amp;MO_KP</b>（<code>AO_EXPORT_TYPE OVERLAP_MATRIX</code>，<code>NDIGITS 14</code>）→ <code>*.mokp</code></li>
<li><b>&amp;KS_CSR_WRITE</b>（<code>REAL_SPACE F</code>，<code>BINARY F</code>）→ <code>*KS_SPIN_s_K_k*.csr</code></li>
</ol>
<p>k 点需完整网格：<code>FULL_GRID T</code>、<code>SYMMETRY F</code>。没有 KS 文件时仍可算 COOP、COBI、电荷与原子 PDOS。</p>
<p>步骤：选择文件 → 读取并检查 → 添加键或勾选自动搜索 → 开始计算 → 在结果表中选键查看曲线 → 导出。</p>
<p>可点“生成演示数据”体验完整流程（合成模型，不是 CP2K 结果）。</p>
"""
