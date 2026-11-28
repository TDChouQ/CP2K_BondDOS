# CP2K 输入里要加什么

适用于 CP2K 2026.2 及以后版本。完整的操作流程见 `docs/使用教程.md`。下面每一条都用 CP2K 2026.2 实际跑过（包括“不写会怎样”）。

做法：在**结构优化完成后的单点能输入**（`RUN_TYPE ENERGY`）里改三处，全部在 `&FORCE_EVAL / &DFT` 里面。命令行 `python main.py template --uks --smear --mesh 3 3 3` 可以直接打印这一段；`examples/cp2k_input/` 里有两个能直接运行的完整例子。

## 一、必须加的三处

### ① `&KPOINTS`（只算 Γ 点也要写）

```text
&KPOINTS
  SCHEME MONKHORST-PACK 3 3 3      ! 只算 Γ 点的大超胞写 SCHEME GAMMA
  SYMMETRY F                       ! 默认值；不要开点群对称性约化
&END KPOINTS
```

* **不写 `&KPOINTS` 就没有 `.mokp` 文件**，KS 矩阵也变成没有 k 点编号的 `…KS_SPIN_1-1_0.csr`，BondDOS 无法使用。
* 偶数网格（如 4 4 4）请加一行 `GAMMA_CENTERED T`，否则网格不含 Γ 点，fcc、六方晶格里对称等价的键会算出不同的值（Si 4×4×4 实测相差约 3%，加上后完全相等）。
* `FULL_GRID T` 可写可不写：CP2K 默认只保留一半 k 点（时间反演约化），BondDOS 会自动展开。
* `WAVEFUNCTIONS` 用默认的 `COMPLEX`。Γ 点用 `REAL` 也能读。
* **conda-forge 版 CP2K 2026.2 的已知问题**：实数对角化（`WAVEFUNCTIONS REAL`，或完全不写 `&KPOINTS` 的普通 Γ 点计算，包括 OT 的 `FULL_ALL` 预条件）在 8 个 Si 原子这样的体系上就会在 ELPA 里段错误（`elpa_autotune`）。在 `&GLOBAL` 里加 `PREFERRED_DIAG_LIBRARY SCALAPACK` 即可；复数 k 点路径不受影响。

### ② `&SCF`：对角化 + 足够的空轨道

```text
&SCF
  ADDED_MOS 50                     ! 按体系调整，见下
  &DIAGONALIZATION
  &END DIAGONALIZATION
  &MIXING
    METHOD BROYDEN_MIXING
    ALPHA 0.2
  &END MIXING
  ! 金属、窄带隙磁性氧化物再加：
  ! &SMEAR ON
  !   METHOD FERMI_DIRAC
  !   ELECTRONIC_TEMPERATURE [K] 300
  ! &END SMEAR
&END SCF
```

* k 点计算不能用 OT，必须用对角化。
* `MO_KP` **只导出 SCF 里实际存在的轨道**，没有 `NLUMO` 关键字。`ADDED_MOS` 决定 COHP/PDOS 曲线在费米能以上能画到多高：想看到 E_F 以上约 5 eV 的反键态，一般每个原子给 3–5 个空轨道；`ADDED_MOS -1` 导出全部轨道（小体系可用，文件最大）。
* BondDOS 读入后会检查“最高导出轨道占据数”，不为 0 说明 `ADDED_MOS` 太少、可能漏掉占据态。
* UKS、Fermi 展宽、固定磁矩（`FIXED_MAGNETIC_MOMENT`）都支持。

### ③ `&DFT / &PRINT` 里加两个输出段

```text
&PRINT                             ! 已有 &PRINT 时把下面两段放进去，不要再写一个 &PRINT
  &MO_KP ON
    FILENAME ./cohp_mo
    AO_EXPORT_TYPE OVERLAP_MATRIX  ! 必须
    NDIGITS 12
    UNIT ANGSTROM
  &END MO_KP
  &KS_CSR_WRITE ON
    FILENAME ./cohp_ks
  &END KS_CSR_WRITE
&END PRINT
```

**注释不要写在 `FILENAME` 的同一行**：CP2K 会把 `!` 后面的内容也当成文件名（实测得到 `cohp_mo!->…mokp`）。

也可以把 `examples/cp2k_input/bonddos_print.inc` 复制到计算目录，在 `&PRINT` 里写一行 `@INCLUDE 'bonddos_print.inc'`。

| 关键字 | 是否必须 | 说明（不写的后果） |
|---|---|---|
| `&MO_KP ON` | **必须** | 输出 `.mokp`：k 点、本征值、占据数、MO 系数（`ON` 可省，建议写上） |
| `AO_EXPORT_TYPE OVERLAP_MATRIX` | **必须** | 默认是 `GTO_BASIS`，文件里没有 S(k)，COOP、COBI、电荷、原子 PDOS 和 E_cov 规范都算不了（仅当 `ADDED_MOS -1` 导出了全部轨道时 BondDOS 能自行重建 S） |
| `&KS_CSR_WRITE ON` | **算 COHP 必须** | 输出每个自旋、每个 k 点的 H(k)（实测写不写 `ON` 都会输出，建议写上）。没有这一段仍可算 COOP、COBI、电荷、PDOS，但算不了 COHP |
| `NDIGITS 12` | 建议 | 默认 9 也能用（已测）；12–14 时一致性检查更严格，文件稍大 |
| `UNIT ANGSTROM` | 可选 | 默认 Bohr，BondDOS 两种都认 |
| `FILENAME` | 可选 | 不写时文件名是 `<PROJECT>-1_0.mokp`、`<PROJECT>-KS_SPIN_1_K_1-1_0.csr` |
| `REAL_SPACE F`、`BINARY F` | 保持默认 | `REAL_SPACE T`（实空间矩阵）和 `BINARY T`（二进制）不支持 |
| `UPPER_TRIANGULAR` | 随意 | `F`（默认）或 `T` 都能读，自动识别 |
| `&EACH` / `ADD_LAST` | 不需要 | 单点计算默认只在 SCF 结束后写一次 |

### 可选：CP2K 自带的 DOS/PDOS（BondDOS 第一页用）

```text
&DOS
  DELTA_E 0.0003674932             ! Hartree，约 0.01 eV
  NLUMO 50
  &CURVE
    ENERGY_UNIT EV
    ENERGY_ZERO FERMI
    &BROADEN
      TYPE GAUSSIAN
      WIDTH [eV] 0.10              ! FWHM
    &END BROADEN
  &END CURVE
  &PDOS
    COMPONENTS F                   ! T：分解到 px/py/pz、dxy…
  &END PDOS
&END DOS
```

按 KIND 输出；要看单个原子的 PDOS 不必拆 KIND，用 BondDOS 的“电荷与原子 PDOS”页或 `pdos` 命令即可（只需要 `.mokp`）。

`&DOS` 各选项的实测结果（CP2K 2026.2）：

| 设置 | 实测 |
|---|---|
| k 点计算 | 只输出展宽曲线 `*-k<n>_curve-1.pdos`，没有逐轨道的原始文件 |
| `&LDOS`、`&R_LDOS`（按原子列表投影） | **k 点计算时被 CP2K 忽略**（输出里有警告）；只在不写 `&KPOINTS` 的 Γ 点计算里生效（`*-list1-1.pdos`）。k 点下要单个原子的 PDOS 用 BondDOS |
| 闭壳层（RKS）曲线的强度 | 每个轨道计 1 个态（不含自旋简并）。BondDOS 由 `.mokp` 算的原子 PDOS 含两个自旋，正好是它的 2 倍（逐点核对一致到 10⁻⁹）；`dos` 命令加 `--scale 2` 可换算 |
| `ENERGY_ZERO` | `FERMI` / `HOCO` / `ABSOLUTE` / `AUTO` 都能读，BondDOS 可互相换算。UKS 两个自旋费米能不同时，CP2K 所有曲线都以较高者为零点 |
| `ENERGY_UNIT` | 默认 `HARTREE`，建议写 `EV`；两种都能读 |
| `TYPE` | `GAUSSIAN`、`LORENTZIAN`、`PSEUDO_VOIGT` 都能读（Lorentzian 尾巴长，窗口内积分略小于态数） |
| `MP_GRID n n n` | 按指定的 MP 网格计算 DOS/PDOS（仅 k 点计算；实测可用）；不影响 `.mokp` |
| `NLUMO` | 对角化计算会自动把 `ADDED_MOS` 增大到 `NLUMO`；OT（仅 Γ 点）在 SCF 后另算空轨道 |
| `APPEND T` | 结构优化每一步追加一段，BondDOS 读最后一段 |

## 二、算完后交给 BondDOS 的文件

```text
cohp_mo-1_0.mokp                      1 个
cohp_ks-KS_SPIN_1_K_1-1_0.csr …       每个 k 点 1 个；UKS 还有 KS_SPIN_2_…
*.pdos / *.dos                        可选，第一页用
```

把 `.mokp` 和全部 `.csr` 放在同一个文件夹。界面：“k 点成键 COHP”页选 `.mokp` 和该文件夹 →“读取并检查数据”。命令行：`python main.py cohp cohp.in`（见 `cohp.in.example`）。

文件大小参考（文本格式，NDIGITS 14）：

| 体系 | `.mokp` | 全部 `.csr` |
|---|---|---|
| Si 原胞 2 原子，3×3×3，26 AO | 2 MB | 1.6 MB |
| MnO 原胞 2 原子 UKS，3×3×3，39 AO | 6 MB | 7 MB |
| Si 54 原子超胞，Γ 点，702 AO / 228 MO | 19 MB | 28 MB |
| MnO 54 原子超胞 UKS，Γ 点，1053 AO / 501 MO | 62 MB | 134 MB |

大小大致正比于 k 点数 × 自旋数 × AO 数 × 导出的 MO 数。上百原子的体系用 Γ 点（PAIR 模式）并控制 `ADDED_MOS`。

## 三、读入后先看这几项

“读取并检查数据”（或 `python main.py cohp cohp.in --check`）应当全部通过：

| 检查 | 不通过时 |
|---|---|
| MO 归一化 `max\|C†SC−I\|` | `NDIGITS` 太小，或 `.mokp` 不完整 |
| 广义本征方程 `HC = SCε` | `.csr` 和 `.mokp` 不是同一次计算（旧文件没删干净） |
| 最高导出轨道占据数 | `ADDED_MOS` 太少 |
| k 网格 | 提示“点群对称性约化”→ 改 `SYMMETRY F`；提示“不含 Γ 点”→ 加 `GAMMA_CENTERED T` 或用奇数网格 |
| 相位约定检查 | 应显示“一致” |

## 四、常见情况

* **结构优化（GEO_OPT / CELL_OPT）时不要开这两个输出**：先优化，再用最终结构单独做一次 `RUN_TYPE ENERGY`。实测优化时每一步都会各写一套文件（`cohp_mo-1_0.mokp`、`-1_1.mokp`…，`.csr` 同理），又大又乱。真遇到这种文件夹，选最后一个 `.mokp` 即可，BondDOS 会按 HC = SCε 自动找出与它匹配的那组 `.csr`。
* **小晶胞里具体某一根键**（IMAGE 模式）需要多 k 点，至少 3×3×3，并做 k 网格收敛；大超胞 Γ 点只能得到“A–B 全部周期镜像合计”（PAIR 模式），超胞足够大时两者相同。
* **表面、二维体系**：非周期方向取 1 个 k 点（如 3 3 1），BondDOS 的“周期方向”保持 XYZ 即可。
* **DFT+U、杂化泛函**：输出的 H(k) 已包含这些项，流程不变（未单独测试杂化泛函）。
