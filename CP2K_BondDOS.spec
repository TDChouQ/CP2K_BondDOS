# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec: python -m PyInstaller --noconfirm CP2K_BondDOS.spec
block_cipher = None

datas = [
    ("assets/app_icon.ico", "assets"),
    ("assets/app_icon.png", "assets"),
    ("docs/方法说明.md", "docs"),
    ("docs/CP2K输入设置.md", "docs"),
    ("docs/使用教程.md", "docs"),
    ("examples/cp2k_input", "examples/cp2k_input"),
    ("examples/WO3_gamma_pdos", "examples/WO3_gamma_pdos"),
]

excludes = ["tkinter", "scipy", "pandas", "IPython", "jupyter", "pytest", "PyQt5", "PyQt6", "PySide2",
            "PySide6.QtWebEngineCore", "PySide6.QtWebEngineWidgets", "PySide6.QtQml", "PySide6.QtQuick",
            "PySide6.Qt3DCore", "PySide6.QtMultimedia", "PySide6.QtCharts", "PySide6.QtDataVisualization"]

a = Analysis(["main.py"], pathex=["."], binaries=[], datas=datas,
             hiddenimports=["matplotlib.backends.backend_qtagg", "matplotlib.backends.backend_pdf",
                            "matplotlib.backends.backend_svg"],
             hookspath=[], runtime_hooks=[], excludes=excludes, noarchive=False)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)
exe = EXE(pyz, a.scripts, a.binaries, a.zipfiles, a.datas, [],
          name="CP2K_BondDOS", debug=False, strip=False, upx=False, runtime_tmpdir=None,
          console=False, icon="assets/app_icon.ico", version=None)
