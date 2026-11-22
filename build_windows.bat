@echo off
chcp 65001 >nul
setlocal
cd /d "%~dp0"
echo ============================================================
echo   CP2K BondDOS - 一键打包 Windows 程序
echo   需要已安装 Python 3.10 - 3.13 (python.org, 勾选 Add to PATH)
echo ============================================================

set PY=
where py >nul 2>nul && set PY=py -3
if "%PY%"=="" (where python >nul 2>nul && set PY=python)
if "%PY%"=="" (
  echo [错误] 没有找到 Python。请先从 https://www.python.org 安装 Python 3.12，安装时勾选 "Add python.exe to PATH"。
  pause & exit /b 1
)

if not exist .venv (
  echo [1/4] 创建虚拟环境 .venv ...
  %PY% -m venv .venv || (echo 创建虚拟环境失败 & pause & exit /b 1)
)
call .venv\Scripts\activate.bat

echo [2/4] 安装依赖（首次约需几分钟）...
python -m pip install --upgrade pip
python -m pip install -r requirements-dev.txt || (echo 依赖安装失败，请检查网络或 pip 镜像 & pause & exit /b 1)

echo [3/4] 运行自检测试 ...
python -m pytest tests -q
if errorlevel 1 (
  echo [警告] 自检有失败项，请把上面的输出发给开发者。仍继续打包。
)

echo [4/4] 打包 exe ...
python -m PyInstaller --noconfirm --clean CP2K_BondDOS.spec || (echo 打包失败 & pause & exit /b 1)

echo.
echo 完成！程序位于: %cd%\dist\CP2K_BondDOS.exe
echo 首次启动需解包运行库，可能要等几秒。
explorer dist
pause
