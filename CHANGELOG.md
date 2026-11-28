# 更新记录

## v1.1.0（2026-11-28）

用 WSL 下的 CP2K 2026.2 做了真实算例端到端测试后的修正与改进。

跨平台 / 无界面版

- **Linux、macOS 无界面命令行版**：`pip install git+https://github.com/TDChouQ/CP2K_BondDOS` 得到命令 `cp2k-bonddos`（只依赖 numpy、matplotlib）；`build_linux.sh` 生成单文件可执行程序（约 35 MB）；没有 PySide6 时不带参数运行会打印帮助而不是报错。
- 新命令 `pdos`：由 `.mokp` 算任意原子/原子组的 PDOS、能带中心（d 带中心等）和电荷，原来只有图形界面能做。
- `cohp --plot` 直接输出 −COHP/COOP/COBI 曲线图；`dos --scale`；`--version`。
- 附带 GitHub Actions 配置 `ci/github-workflow-build.yml`（三个平台跑测试，自动构建 Linux 命令行程序和 Windows exe；在 GitHub 网页上把它放到 `.github/workflows/build.yml` 即可启用）。
- 新增 `docs/使用教程.md`（从 CP2K 输入到结果解读的完整流程，数字均来自附带的 MnO 例子）。

修复

- 读取 `.mokp` 时，Fortran 输出的三位指数（如 Fermi 展宽占据数 `2.91755220991368-102`，没有字母 E）导致读取失败：用了 `&SMEAR` 且导出了高能空轨道的计算无法载入。
- `numpy < 2.0` 下 `np.trapezoid` 不存在，`dos` / `band_center` 命令和曲线统计直接报错（requirements 声明支持 numpy ≥ 1.24）。
- 费米能估计：费米能落在带隙里时线性插值可偏差 1 eV 以上，改为 Fermi–Dirac 拟合（与 CP2K 输出一致到 10⁻⁵ eV）；UKS 两个自旋分别估计，固定磁矩计算（两个费米能）给出提示。
- `DOSCAR.lobster`、`GROSSPOP.lobster` 原来不能被 pymatgen 读取，已修正格式；全部 9 种 LOBSTER 文件现在用 pymatgen 2026.10 和 LobsterPy 0.6.1 实测可读（自旋极化与非极化），LobsterPy 的自动成键分析可直接使用。
- `ENERGY_ZERO ABSOLUTE` 时 CP2K 在总 DOS 文件头里把费米能写成 0，原来会被当真；现在改从同一次计算的 PDOS 文件取。
- NMO = NAO 时不再因“最高导出轨道仍有占据”误报。
- `.gitignore` 的 `!tests/fixtures/*` 不覆盖子目录，补充 `!tests/fixtures/**`。

新增

- **AO-COHP 默认改为共价键能规范**：`COHP_GAUGE = ECOV | EREF | RAW`（界面：“AO-COHP 规范”）。原始 AO-COHP 依赖周期体系的势能零点，异核键连符号都可能出错（FM-MnO 的 Mn–O：RAW +1.02 eV，ECOV −1.34 eV）。**同一数据的 ICOHP 数值因此与 v1.0.0 不同**；要复现 v1.0.0 设 `COHP_GAUGE = RAW`。详见方法说明 2.6。
- **CP2K 输入说明**：新增 `docs/CP2K输入设置.md`（界面新增“CP2K 输入设置”页）、命令 `template`、`examples/cp2k_input/`（`bonddos_print.inc` 和两个可直接运行的完整输入）。每个关键字是否必须、不写的后果都用 CP2K 2026.2 实测。
- 输入不全时的处理：找不到 `.mokp`、KS 文件没有 k 点编号（缺 `&KPOINTS`）、实空间 KS 矩阵时给出对应的修改提示；`AO_EXPORT_TYPE GTO_BASIS` 且 NMO = NAO 时由 S = (CC†)⁻¹ 重建 S(k)；上三角 KS 矩阵自动识别。
- 固定磁矩计算的 `.dos/.pdos`：CP2K 的展宽曲线统一以较高的费米能为零点，两个自旋现在统一使用它（原来原始文件与曲线文件的 β 通道零点相差两费米能之差）。
- 界面里 Löwdin 表象不再标注“类 LOBSTER”（MOLOPT 基组下 ICOHP 数值偏大，只宜定性）。
- 时间反演约化的 k 网格（CP2K 默认 `FULL_GRID F`、`SYMMETRY F`）自动展开为完整网格，IMAGE 模式可直接使用。
- 偶数 MP 网格（不含 Γ）给出对称性破缺提示；命令行输出 k 网格相关警告。
- 命令行输出重定向到文件时不再刷进度条。
- 测试：新增 `tests/test_v101.py` 和两个真实 CP2K 2026.2 回归样例（石墨烯、铁磁 MnO），参考值来自等价 Γ 点超胞。
- PDOS 实测：由 `.mokp` 算的原子 PDOS 与 CP2K 自带的 k 点 PDOS 曲线逐点一致（10⁻⁹，闭壳层差一个自旋简并因子 2）；`&DOS` 的各种设置（`COMPONENTS`、`&LDOS`、三种展宽、四种能量零点、两种单位、`MP_GRID`、`APPEND`、OT、UKS、Γ 点原始文件）都用真实输出测过，结论见 `docs/CP2K输入设置.md`。
- 一个文件夹里有多组 KS 文件（结构优化每步一组）时，按 HC = SCε 自动选出与 `.mokp` 匹配的一组。
- 另用真实算例确认：DFT+U、实波函数 Γ 点（`WAVEFUNCTIONS REAL`）、`UPPER_TRIANGULAR T`、MO_KP 默认 `NDIGITS 9`、`GAMMA_CENTERED T` 偶数网格都能正确读入；

## v1.0.0（2026-11-22）

第一版。

- 态密度：读取 CP2K 2026.2 `&DOS` 输出的 `.dos/.pdos`（Γ 点逐轨道、`_curve` 曲线、k 点 PDOS、自旋极化），合并、重新展宽、IDOS，能带中心/宽度/偏度/峰度/填充率/上边缘。
- k 点成键：基于 `&MO_KP` + `&KS_CSR_WRITE` 的 COHP/COOP/COBI，AO 与 Löwdin 两种表象，周期镜像（IMAGE）与原子对合计（PAIR），轨道分解，自动键搜索。
- 电荷与布居：Mulliken/Löwdin 电荷、磁矩、轨道布居，原子/轨道分辨的 k 点 PDOS。
- 导出：CSV、JSON 报告，LOBSTER 格式文件（COHPCAR、ICOHPLIST、COOPCAR、ICOOPLIST、COBICAR、ICOBILIST、CHARGE、GROSSPOP、DOSCAR）与 POSCAR。
- 数据检查：S 正定性、归一化、广义本征方程、电子数、k 网格混叠指标、周期镜像相位约定检查。
- Windows 图形界面（PySide6）与跨平台命令行；一键打包脚本。
