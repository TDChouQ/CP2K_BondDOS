@echo off
chcp 65001 >nul
cd /d "%~dp0"
if exist .venv\Scripts\activate.bat (call .venv\Scripts\activate.bat) else (
  echo 首次运行请先双击 build_windows.bat（会创建 .venv 并安装依赖）
  pause & exit /b 1
)
python main.py %*
