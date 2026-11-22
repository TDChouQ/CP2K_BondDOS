# 更新记录

## v1.0.0（2026-11-22）

第一版。

- 态密度：读取 CP2K 2026.2 `&DOS` 输出的 `.dos/.pdos`（Γ 点逐轨道、`_curve` 曲线、k 点 PDOS、自旋极化），合并、重新展宽、IDOS，能带中心/宽度/偏度/峰度/填充率/上边缘。
- k 点成键：基于 `&MO_KP` + `&KS_CSR_WRITE` 的 COHP/COOP/COBI，AO 与 Löwdin 两种表象，周期镜像（IMAGE）与原子对合计（PAIR），轨道分解，自动键搜索。
- 电荷与布居：Mulliken/Löwdin 电荷、磁矩、轨道布居，原子/轨道分辨的 k 点 PDOS。
- 导出：CSV、JSON 报告，LOBSTER 格式文件（COHPCAR、ICOHPLIST、COOPCAR、ICOOPLIST、COBICAR、ICOBILIST、CHARGE、GROSSPOP、DOSCAR）与 POSCAR。
- 数据检查：S 正定性、归一化、广义本征方程、电子数、k 网格混叠指标、周期镜像相位约定检查。
- Windows 图形界面（PySide6）与跨平台命令行；一键打包脚本。
