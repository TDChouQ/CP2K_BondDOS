#!/usr/bin/env bash
# Build a single-file, GUI-less Linux executable:  dist/cp2k-bonddos
# Needs Python 3.10+ with venv.  Usage:  bash build_linux.sh
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
if [ ! -d .venv-linux ]; then "$PY" -m venv .venv-linux; fi
. .venv-linux/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-cli.txt pyinstaller pytest scipy
python -m pytest tests -q || echo "[warning] some self tests failed - please report the output above"
python -m PyInstaller --noconfirm --clean CP2K_BondDOS_cli.spec
echo
echo "Done: $(pwd)/dist/cp2k-bonddos"
./dist/cp2k-bonddos --version
