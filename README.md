# CP2K BondDOS

CP2K BondDOS 是面向 CP2K 2026.2+ 的电子结构后处理软件（Windows 图形界面 + 跨平台命令行），把 DOS/PDOS 与 k 点 COHP/COOP/COBI 成键分析放在一个程序里。主要功能：

| 功能 | 说明 |
|---|---|
| 自动搜索全部键 | 按距离范围（可限定元素）列出所有原子对与周期镜像，一次算出 ICOHP/ICOOP/ICOBI 总表（类 LOBSTER cohpGenerator） |
| COBI / ICOBI 键级 | Löwdin 表象为 Wiberg 型（与 LOBSTER 定义一致），AO 表象为 Mayer 型 |
| 电荷与布居 | Mulliken / Löwdin 原子电荷、UKS 磁矩、逐轨道布居（GROSSPOP） |
| 原子分辨 k 点 PDOS | 由 MO_KP 直接投影到任意原子/原子组及 s/p/d/f、px/dxy…，**无需拆分 KIND** |
| LOBSTER 格式导出（pymatgen、LobsterPy 实测可读） | COHPCAR / ICOHPLIST / COOPCAR / ICOOPLIST / COBICAR / ICOBILIST / CHARGE / GROSSPOP / DOSCAR.lobster + POSCAR，可直接给 pymatgen、LobsterPy |
| 能带描述符 | 中心、宽度、偏度、峰度、填充率、费米能以下中心、上边缘（Hilbert 变换） |
| 自动识别 k 网格 | 无需手填网格尺寸；支持偏移 MP 网格 |
| 相位/混叠自检 | 重叠衰减判据检查周期镜像相位约定、±0.5 边界原子的回卷方向；给出 k 网格混叠指标 |
| 出版级作图 | 能量横排/竖排（LOBSTER 风格）、−COHP、积分曲线单独面板、β 镜像、填充、字号线宽、按英寸和 DPI 导出 PNG/PDF/SVG/TIFF |

**第一次用请看 `docs/使用教程.md`**：从 CP2K 输入到 COHP 曲线、d 带中心、电荷的完整流程，附可直接运行的例子。

| 平台 | 形态 |
|---|---|
| Windows | 图形界面 `CP2K_BondDOS.exe`（双击 `build_windows.bat` 生成） |
| Linux / macOS / 服务器 | 无界面命令行：`pip install git+https://github.com/TDChouQ/CP2K_BondDOS` 后用 `cp2k-bonddos`；或 `bash build_linux.sh` 生成单文件程序 |

## 〇、CP2K 输入里要加什么

在结构优化后的单点能输入（`RUN_TYPE ENERGY`）的 `&DFT` 里改三处，详见 **`docs/CP2K输入设置.md`**（软件内“CP2K 输入设置”页；逐项用 CP2K 2026.2 实测）：

| 位置 | 加什么 | 不加的后果 |
|---|---|---|
| `&DFT/&KPOINTS` | `SCHEME MONKHORST-PACK n n n`（只算 Γ 点写 `SCHEME GAMMA`）；偶数网格加 `GAMMA_CENTERED T` | 没有 `&KPOINTS` 就不输出 `.mokp` |
| `&DFT/&SCF` | 对角化（不能用 OT）+ `ADDED_MOS N` | 空轨道不够，费米能以上的反键区画不出来 |
| `&DFT/&PRINT` | `&MO_KP ON` 且 `AO_EXPORT_TYPE OVERLAP_MATRIX`；`&KS_CSR_WRITE ON` | 没有 S(k) / H(k)，COOP、COBI、电荷或 COHP 算不了 |

```bash
python main.py template --uks --smear --mesh 3 3 3   # 打印可直接粘贴的输入段
```

`examples/cp2k_input/` 里有两个可直接运行的完整例子（Si，铁磁 MnO）和 `bonddos_print.inc`（在 `&PRINT` 里 `@INCLUDE` 即可）。

## 一、生成 exe（只需做一次）

1. 安装 Python 3.12（https://www.python.org ，安装时勾选 **Add python.exe to PATH**）。
2. 双击 **`build_windows.bat`**。脚本会建立 `.venv`、安装 PySide6/Matplotlib/NumPy/PyInstaller、运行自检测试，然后打包。
3. 完成后得到 **`dist\CP2K_BondDOS.exe`**（单文件，带图标，可复制到任何 Windows 10/11 电脑使用，不需要 Python）。

不想打包时，也可以在第 2 步完成后双击 `run_from_source.bat` 直接运行。

> 国内网络若 pip 很慢，可先执行 `pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple` 再双击打包脚本。

## 二、三个工作页

### 1. 态密度 DOS / PDOS（CP2K 原生 `.dos/.pdos`）

导入文件或整个文件夹 → 勾选数据与轨道（总投影、s/p/d/f，或填 `dxy,dz2` 等分量）→ 选择能量零点（E_F / HOCO / 原零点 / 绝对能量）→ 图和统计表自动更新。

* “合并（求和）所选”：多个 KIND 逐点求和（命令行 `merge_pdos`）。
* 统计表给出能带中心、宽度、偏度、峰度、填充率和上边缘（命令行 `band_center`）。
* 原始逐轨道数据可设 FWHM 重新展宽；k 点 `_curve` 文件已按 CP2K 的 WIDTH 展宽。
* “WO₃ 示例”按钮载入随包的真实 WO₃ Γ 点 PDOS。

### 2. k 点成键 COHP

选择 `.mokp` 与 KS 文件夹 → **读取并检查数据** → 设置 → 添加键或勾选自动搜索 → **开始计算** → 在结果表选键看曲线 → **导出全部结果**。

键的写法：

| 写法 | 含义 |
|---|---|
| `1 2` | 1 号与 2 号原子的最近一根键 |
| `1 2 ALL` | 所有等距的最近周期镜像 |
| `1 2 -1 0 0` | 指定晶格平移（以 `.mokp` 打印的原始坐标为准） |
| `1 2 PAIR` | 1–2 的全部周期镜像合计（Γ-only 也可用） |
| 末尾加 `ORB` | 对这根键做轨道分解（分解方式在“轨道分解”中选） |

结果表的“状态”：✔ 正常；⚠ Nyquist（R 正好是网格一半，R 与 −R 混叠）；⚠ Γ叠加（该方向只有一个 k 点）；✘ 混叠（|R| 超过网格一半，不计算）。

“生成演示数据”会创建一个合成紧束缚模型（**不是 CP2K 结果**）用于体验完整流程。

### 3. 电荷与原子 PDOS

需先在第 2 页读取 `.mokp`（KS 文件可无）。原子组每行写 `名称: 选择`，选择可以是元素符号或编号区间，如 `Fe_surf: 1-4,9`。能量网格、展宽、零点沿用第 2 页设置。

## 三、命令行（Linux 服务器或批处理）

```bash
python main.py cohp cohp.in --check      # 只检查
python main.py cohp cohp.in --plot       # 计算并导出（配置示例见 cohp.in.example），同时出图
python main.py pdos cohp_mo-1_0.mokp --group "Mn: Mn" --group "O: 2,5-7" --stats bc.csv --charges q.csv --plot pdos.png
python main.py dos *.pdos --channels total,d --emin -10 --emax 5 --stats bc.csv --plot pdos.png
python main.py merge_pdos *.pdos         # -> merged_pdos.dat
python main.py band_center *.pdos        # -> band_center.dat
python main.py demo demo_folder          # 合成测试数据
python main.py template --uks --smear    # CP2K 输入中需要加入的内容
```

Linux 上 `pip install -r requirements-cli.txt`（只有 numpy、matplotlib，不需要 PySide6）；pip 安装后命令名为 `cp2k-bonddos`。

## 四、正确性

* 自带 34 项自动测试（`python -m pytest tests -q`），其中 k 点 COHP/COOP/COBI 与独立的 Γ 点超胞实空间计算一致到 10⁻⁸，覆盖 AO/Löwdin、RKS/UKS、偏移网格、坐标平移不变性、轨道分解求和、LOBSTER 文件格式。
* CP2K 文件格式按 v2026.2 源码（`kpoint_mo_dump.F`、`qs_scf_csr_write.F`、DBCSR `csr_write`、`qs_dos.F`、`qs_pdos.F`）核对；Fourier 相位与回卷约定按 `kpoint_methods.F`、`qs_neighbor_lists.F` 核对。
* **已用真实 CP2K 2026.2 输出做端到端验证（v1.1.0）**：Si（RKS）、石墨烯（Fermi 展宽、二维 3×3×1）、铁磁 MnO（UKS、固定磁矩展宽、原子在晶胞外、O 恰在 0.5 边界）三个体系中，原胞 k 点 IMAGE 结果与等价 Γ 点超胞的 PAIR 结果一致到 10⁻⁷–10⁻⁸；Mulliken 电荷/磁矩、费米能与 CP2K 主输出一致。其中两个小算例（SZV 基组）作为回归样例放在 `tests/fixtures/real_cp2k_2026.2/`。

方法细节见软件内“方法说明”页或 `docs/方法说明.md`；CP2K 输入设置与完整计算流程见 `docs/DOS和COHP计算流程.md`。

## 五、文件结构

```text
CP2K_BondDOS/
├─ build_windows.bat        一键打包
├─ run_from_source.bat      不打包直接运行
├─ build_linux.sh           Linux 无界面单文件程序
├─ pyproject.toml           pip 安装（命令 cp2k-bonddos）
├─ main.py                  入口（无参数=图形界面）
├─ CHANGELOG.md             更新记录
├─ CP2K_BondDOS.spec        PyInstaller 配置
├─ bonddos/                 源码（mokp/kspace/bonding/population/dosfiles/dosanalysis/lobster/plotting/gui）
├─ assets/                  图标 app_icon.svg/png/ico
├─ docs/                    使用教程、CP2K输入设置、方法说明、计算流程
├─ examples/                WO₃ Γ 点 PDOS 示例；cp2k_input/ 为 CP2K 输入模板与完整例子
└─ tests/                   自动测试与 2026.2 格式样例
```
