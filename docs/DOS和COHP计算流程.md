# CP2K 中 DOS、PDOS 与 COHP 的计算流程

## 通过 CP2K 输出 DOS 与 PDOS

**适用于 CP2K 2026.2 及以后版本（新版 DOS/PDOS 接口）**

CP2K 2026.2 通过 `&DFT / &PRINT / &DOS` 输出 DOS 和 PDOS。计算前需要保证有足够的空轨道，否则 DOS/PDOS 无法覆盖未占据态。

使用 Multiwfn 辅助创建单点能输入文件时，可设置计算额外空轨道：

```text
-9                // 其它设置
14                // 输出轨道能级
12                // 计算额外空轨道
```

在正常单点能输入文件的 `&DFT` 部分加入：

```text
    &PRINT
      &DOS
        DELTA_E 0.0003674932
        NLUMO 100
        NDIGITS 8

        &CURVE
          ENERGY_UNIT EV
          ENERGY_ZERO FERMI

          &BROADEN
            TYPE GAUSSIAN
            WIDTH [eV] 0.10
          &END BROADEN
        &END CURVE

        &PDOS
          COMPONENTS F
        &END PDOS

      &END DOS
    &END PRINT
```

其中：

- `DELTA_E 0.0003674932`：能量网格间隔，单位 Hartree（≈ 0.01 eV）。它总以 Hartree 解释，与 `ENERGY_UNIT` 无关。
- `NLUMO 100`：DOS/PDOS 中考虑的未占据轨道数。标准对角化时 CP2K 会按需增加 `ADDED_MOS`，OT 时会在 SCF 后另行计算空轨道；仍应检查输出覆盖到所需的能量范围。
- `ENERGY_UNIT EV`：输出能量单位，可选 `EV` 或 `HARTREE`。
- `ENERGY_ZERO`：能量零点，**2026.2 只有四个选项**：
  - `FERMI`：费米能为 0，适合金属、掺杂体系；
  - `HOCO`：最高占据轨道为 0，适合半导体、绝缘体；
  - `AUTO`：有展宽/分数占据时取 FERMI，否则取 HOCO（文件头写 `AUTO -> FERMI` 等）；
  - `ABSOLUTE`：不移动能量。
  - 注意：没有 `LUMO`、`VACUUM` 选项。
- `WIDTH [eV] 0.10`：展宽宽度，是 **FWHM**（半高全宽）。建议显式写单位。常用 0.05（尖锐）/ 0.10（常用）/ 0.20（平滑）。`TYPE` 可为 `GAUSSIAN`、`LORENTZIAN`、`PSEUDO_VOIGT`。
- `COMPONENTS`：`F` 按 s/p/d/f 输出；`T` 进一步分解到 py/pz/px、d-2…d+2（分别对应 dxy、dyz、dz²、dxz、dx²−y²）。

### 输出文件

| 文件 | 内容 |
|---|---|
| `*-1.dos` | 原始 total DOS（未展宽直方图；Density 列是**归一化的态分数**，Occupation 列是电子数） |
| `*-curve-1.dos` | 展宽后的 total DOS（自旋极化为 5 列：E、α-Density、α-Occ、β-Density、β-Occ） |
| `*-k<n>-1.pdos` | **仅 Γ 点**：KIND n 的逐轨道原始 PDOS（MO、能量、占据、Total、各轨道） |
| `*-k<n>_curve-1.pdos` | KIND n 的展宽 PDOS 曲线；**k 点计算只有这一种 PDOS 文件**（标题为 “K-point projected DOS”） |
| `*-ALPHA_k<n>_curve-1.pdos` / `BETA` | 自旋极化时的两个自旋 |

文件名中的 `n` 是 **KIND 编号**，不是 k 点编号。

### 单个原子的 PDOS

**方法一（推荐，k 点计算）**：在 COHP 所需的 `&MO_KP` 输出基础上，用 CP2K BondDOS 的“电荷与原子 PDOS”页直接计算任意原子/原子组的 PDOS（可按 s/p/d/f、px/dxy… 或基组壳层分解）。不需要修改 KIND。

**方法二（CP2K 原生）**：把目标原子单独定义为一个 KIND，例如 N 掺杂碳中与 N 相邻的某个 C：

```text
&KIND C
  ELEMENT C
  BASIS_SET DZVP-MOLOPT-SR-GTH
  POTENTIAL GTH-PBE
&END KIND

&KIND C_nearN
  ELEMENT C
  BASIS_SET DZVP-MOLOPT-SR-GTH
  POTENTIAL GTH-PBE
&END KIND

&KIND N
  ELEMENT N
  BASIS_SET DZVP-MOLOPT-SR-GTH
  POTENTIAL GTH-PBE
&END KIND
```

坐标中将目标 C 的标签改为 `C_nearN`。`C_nearN` 只是 KIND 标签，元素由 `ELEMENT C` 指定；必须保证基组、赝势以及磁矩初猜、DFT+U 等 KIND 相关设置与普通 C 完全一致。

### 合并 PDOS 与计算 Band Center

用 **CP2K BondDOS** 完成：

- 图形界面：“态密度 DOS / PDOS”页 → 导入文件/文件夹 → 勾选“合并（求和）所选 KIND/文件” → 下方表格直接给出积分、能带中心、宽度、偏度、峰度、填充率、费米能以下中心、上边缘（Hilbert 变换）。
- 命令行：

```bash
CP2K_BondDOS merge_pdos *.pdos            # -> merged_pdos.dat
CP2K_BondDOS band_center 1.pdos 2.pdos     # -> band_center.dat
```

（源码运行时把 `CP2K_BondDOS` 换成 `python main.py`。）

能带中心

$$
\varepsilon_\alpha=\frac{\int E\rho_\alpha(E)\,dE}{\int\rho_\alpha(E)\,dE}
$$

- **All**：整个窗口的能量重心。包含大量空态的轨道会偏向正能量，物理意义需谨慎。
- **Occupied(E≤0)**：只积分零点以下。仅当零点为费米能时它才是“占据部分”；有分数占据或展宽尾巴时与真实占据并不完全等价。
- 逐轨道原始数据用离散权重精确计算，不受展宽影响；曲线数据用梯形积分。
- 不同文件只有在同一计算、同一能量零点、同一网格时才能逐点相加；BondDOS 会按能量插值到公共网格。原始文件与对应的 `_curve` 文件不能同时相加。

## 通过 Multiwfn 计算 DOS 和 PDOS

### DOS 的计算

在 Multiwfn 的 CP2K 输入文件创建界面中完成基本参数设置后，依次选择：

```text
-9                  // 其它设置
14                  // 在 SCF 后输出轨道信息（包括所有 k 点）
10                  // 设置 PRINT_LEVEL
2                   // MEDIUM
12                  // 设置计算的空轨道数
20
```

为了获得较平滑、可靠的 DOS，建议采用较密的 k 点网格。CP2K 计算完成后，启动 Multiwfn，载入 `.out`，依次输入：

```text
cp
5                   // Obtain exact DOS based on energy levels at all k-points
1                   // Generate DOS curve and export to DOS.txt
```

生成 `DOS.txt`（第 1 列能量，第 2 列 DOS）。

### PDOS 的计算（Γ 点超胞法）

```text
-2                  // 输出 Molden 文件
-9                  // 其它设置
3                   // 设置超胞倍数
5,5,5
12                  // 设置空轨道数
500                 // 超胞越大，通常需要设置更多空轨道
```

计算完成后在 `.molden` 中补充 `[Cell]` 与 `[Nval]`，载入 Multiwfn：

```text
10                  // DOS 绘制功能
-1                  // 定义片段
1  cond  O  [Enter]  P  q      // O-2p
2  cond  W  [Enter]  D  q      // W-5d
0
-6  H                // HOMO = 0 eV
-0                   // 绘制 TDOS、O-2p、W-5d
s                    // 保存设置；以后输入 1 载入
```

## COHP 的计算

**适用于 CP2K 2026.2 及以后版本，使用 `MO_KP` 与 `KS_CSR_WRITE` 输出，用 CP2K BondDOS 后处理。**

本方法直接利用 CP2K 的局域 Gaussian AO 表象构造带 k 点的 COHP/ICOHP，不依赖 LOBSTER。同一套数据还给出 COOP、COBI（键指数）、Mulliken/Löwdin 电荷和原子分辨 PDOS。

### CP2K 输入文件设置

> 逐项核对过的最新版说明（每个关键字是否必须、不写的后果、可直接运行的例子）见 `docs/CP2K输入设置.md`；`python main.py template --uks --smear` 可直接打印需要加入的输入段。下面是简述。

**增加未占据轨道**（需要分析费米能以上反键态时）：

```text
&SCF
  ...
  ADDED_MOS 50
  ...
&END SCF
```

`MO_KP` 没有 `NLUMO` 关键字，只输出 SCF 中已经存在的轨道；所以 `ADDED_MOS` 要保证最高导出能级覆盖所需能窗。BondDOS 会检查最后一个导出轨道是否为空，若仍有占据会报警（说明可能漏掉占据态）。必须使用标准对角化（不能用 OT）。

**设置 k 点**

```text
&KPOINTS
  SCHEME MONKHORST-PACK 3 2 4
  FULL_GRID T
  SYMMETRY F
  WAVEFUNCTIONS COMPLEX
&END KPOINTS
```

- `FULL_GRID T`、`SYMMETRY F`：周期镜像分解（指定 A–B(R) 的单根键）需要完整、等权的网格；BondDOS 会从 k 点列表自动识别网格尺寸并检查。
- 只用 Γ 点时也要显式写 `SCHEME GAMMA` 以走 k 点输出路径：

```text
&KPOINTS
  SCHEME GAMMA
  FULL_GRID T
  SYMMETRY F
  WAVEFUNCTIONS COMPLEX
&END KPOINTS
```

> **注意**：只有一个 Γ 点时无法由 H(k) 分离不同的周期镜像，此时请用 PAIR 模式（A–B 全部镜像合计）或用大超胞。小晶胞中具体单根周期键的 COHP 需要多 k 点，并做 k 网格收敛。

**输出 COHP 后处理所需数据**（在 `&DFT` 下、`&SCF` 外）：

```text
&PRINT

  &MO_KP ON
    FILENAME ./cohp_mo
    AO_EXPORT_TYPE OVERLAP_MATRIX
    NDIGITS 14
    UNIT ANGSTROM
  &END MO_KP

  &KS_CSR_WRITE ON
    FILENAME ./cohp_ks
    REAL_SPACE F
    BINARY F
    THRESHOLD 0.0
    UPPER_TRIANGULAR F
    ADD_LAST SYMBOLIC

    &EACH
      QS_SCF 0
    &END EACH
  &END KS_CSR_WRITE

&END PRINT
```

- `AO_EXPORT_TYPE OVERLAP_MATRIX`：同时输出各 k 点的重叠矩阵 S(k)（默认 `GTO_BASIS` 不含 S，不能用于本分析）。
- `NDIGITS 14`：同时决定稀疏截断阈值 10⁻¹⁴，太小会使归一化检查失败。
- `REAL_SPACE F`：输出 k 空间的 H(k)；`BINARY F`：文本格式（不支持二进制）。
- 能量在 `.mokp` 中始终为 Hartree；`UNIT` 只影响长度。

日常计算中，除 k 点网格外通常只需调整 `ADDED_MOS`。

### CP2K 输出文件

```text
cohp_mo-1_0.mokp

cohp_ks-KS_SPIN_1_K_1-1_0.csr
cohp_ks-KS_SPIN_1_K_2-1_0.csr
...
cohp_ks-KS_SPIN_2_K_n-1_0.csr      # 仅 UKS
```

把一个 `.mokp` 和同一次最终电子步的全部 `.csr` 放在同一文件夹。若文件夹中有多组 KS 文件（不同前缀/后缀），BondDOS 会让你选择一组；它会用 HC = SCε 检查 KS 文件与 MO_KP 是否匹配。

### 运行 CP2K BondDOS

**图形界面（Windows）**：双击 `CP2K_BondDOS.exe` → “k 点成键 COHP” 页：

1. 选择 `.mokp` 与 KS 文件夹 → **读取并检查数据**（显示网格、电子数、带隙、各项一致性检查）；
2. 选择表象（直接 AO / Löwdin）、计算的量（COHP/COOP/COBI）、轨道分解、展宽 FWHM、能窗、能量零点；
3. 添加键：`1 2`（最近镜像）、`1 2 ALL`（全部等距最近镜像）、`1 2 -1 0 0`（指定平移）、`1 2 PAIR`（全部镜像合计），末尾加 `ORB` 做轨道分解；或勾选**自动搜索**（例如 1.5–3.0 Å，可限定元素 W–O）；
4. **开始计算** → 结果表中选择一根或多根键查看 −COHP / COOP / COBI 曲线与积分；可对所选键求和、选择轨道对；
5. **导出全部结果** → CSV、JSON 报告，以及 `COHPCAR/ICOHPLIST/COBICAR/ICOBILIST/COOPCAR/ICOOPLIST/CHARGE/GROSSPOP/DOSCAR.lobster` 与 `POSCAR`（可直接给 pymatgen / LobsterPy 使用）。

**命令行（如在 Linux 服务器上用源码运行）**：

```bash
python main.py cohp cohp.in --check   # 只检查数据
python main.py cohp cohp.in           # 计算并导出
```

`cohp.in` 示例：

```text
MOKP = cohp_mo-1_0.mokp
KS_DIR = .
PERIODIC = XYZ
COHP_BASIS = AO          # AO 或 LOWDIN
COHP_GAUGE = ECOV        # 默认；ECOV | EREF | RAW，见方法说明 2.6
COOP = T
COBI = T
ORBITAL_RESOLUTION = L   # NONE | L | LM | SHELL | AO
SPIN = AUTO
FWHM_EV = 0.2355         # 或 SIGMA_EV = 0.10（高斯标准差），FWHM = 2.3548 σ
EMIN_EV = -10.0
EMAX_EV = 5.0
DE_EV = 0.01
EF_EV = auto             # auto | hoco | midgap | fermi | 绝对能量(eV)
BOND = 1 2
BOND = 1 3 ALL
BOND = 2 4 -1 0 0
AUTO_BONDS = 1.5 3.0     # 可选：自动搜索全部键，可加元素 AUTO_BONDS = 1.5 3.0 W O
OUTPUT_PREFIX = results/COHP
```

k 点网格尺寸由程序从 `.mokp` 的 k 点列表自动识别，无需在 `cohp.in` 中填写。完整参数见 `cohp.in.example`。

### 符号约定与单位

- raw COHP：**负值成键、正值反键**；图中默认画 −COHP（与 LOBSTER 图一致，成键为正），表中给出 −ICOHP。
- COOP、COBI：正值成键。
- ICOHP（eV）、ICOOP（电子）、ICOBI（无量纲）由 CP2K 实际占据数对离散态精确求和，**与展宽、绘图能窗无关**，支持金属的分数占据。
- 自旋极化时分别给出 up、down 与 total。

### 方法原理

LOBSTER 面向平面波/PAW 计算，需要把 Kohn–Sham 波函数投影到局域原子轨道基组后再做 Hamilton population 分析，称为 projected COHP（pCOHP）。

CP2K 的 Quickstep（GPW/GAPW）本身用局域 Gaussian AO 展开 Kohn–Sham 轨道，因此可以直接得到 AO 表象中的 MO 系数 C(k)、Kohn–Sham 矩阵 H(k) 和重叠矩阵 S(k)，不需要“平面波 → 局域轨道”的投影。

对每个 k 点与能带 n，原子 A、B 之间的 Hamilton population 为

$$
h_{AB,nk}=2\,\mathrm{Re}\left[C_{A,nk}^{\dagger}H_{AB}(k)C_{B,nk}\right]
$$

按 k 点权重求和并用 Gaussian 展宽得到能量分辨 COHP；ICOHP 是对占据态（实际占据数）的求和。

需要区分具体周期镜像时，先对 H(k) 做离散 Fourier 变换（CP2K 约定 H(k)=Σ_R H(R)e^{+2πik·R}）

$$
\tilde H_{AB}(R)=\sum_k w_k e^{-2\pi i k\cdot R}H_{AB}(k),\qquad
h_{ABR,nk}=2\,\mathrm{Re}\left[C_A^{\dagger}\tilde H_{AB}(R)e^{2\pi i k\cdot R}C_B\right]
$$

得到 A(0)–B(R) 的 bond-resolved COHP。CP2K 内部按回卷坐标定义 R，程序自动换算为你输入的原始坐标平移。把 H 换成 S 得到 COOP；把 H 换成密度矩阵（Löwdin 表象，Wiberg 型）或 PS（AO 表象，Mayer 型）得到 COBI。完整公式与验证见软件内“方法说明”页。

因此本方法的准确名称是：

**direct AO-based COHP in the CP2K Gaussian atomic-orbital representation**

而不是 LOBSTER pCOHP。由于具体周期键的 H(R) 来自有限 k 网格的离散 Fourier 变换，正式计算应做 **k 点网格收敛**；等价化学键的 ICOHP 是否随 k 点加密趋于一致也是有用的检查。BondDOS 报告中的“k 网格混叠指标” max|S̃_μμ(0)−1| 也可帮助判断网格是否足够。

**参考文献**

- R. Dronskowski, P. E. Blöchl, *J. Phys. Chem.* **97**, 8617–8624 (1993). DOI: 10.1021/j100135a014 — COHP
- S. Maintz et al., *J. Comput. Chem.* **37**, 1030–1035 (2016). DOI: 10.1002/jcc.24300 — LOBSTER
- P. C. Müller et al., *J. Phys. Chem. C* **125**, 7959–7970 (2021). DOI: 10.1021/acs.jpcc.1c00718 — COBI
