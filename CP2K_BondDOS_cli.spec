# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec for the command-line-only build (Linux / macOS, also works on Windows):
#   python -m PyInstaller --noconfirm CP2K_BondDOS_cli.spec   ->  dist/cp2k-bonddos
block_cipher = None

excludes = ["tkinter", "scipy", "pandas", "IPython", "jupyter", "pytest", "PyQt5", "PyQt6", "PySide2", "PySide6",
            "bonddos.gui"]

a = Analysis(["main.py"], pathex=["."], binaries=[], datas=[], hiddenimports=["matplotlib.backends.backend_agg"],
             hookspath=[], runtime_hooks=[], excludes=excludes, cipher=block_cipher, noarchive=False)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)
exe = EXE(pyz, a.scripts, a.binaries, a.zipfiles, a.datas, [], name="cp2k-bonddos", debug=False,
          bootloader_ignore_signals=False, strip=False, upx=False, runtime_tmpdir=None, console=True)
