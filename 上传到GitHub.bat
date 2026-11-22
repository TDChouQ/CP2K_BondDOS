@echo off
chcp 65001 >nul
setlocal EnableExtensions
cd /d "%~dp0"
title CP2K BondDOS - 第一次上传到 GitHub
echo ============================================================
echo   CP2K BondDOS - 第一次上传到 GitHub
echo ============================================================
echo.

where git >nul 2>nul
if errorlevel 1 (
  echo [缺少 Git] 请先安装 Git for Windows（一路默认选项即可），装好后重新双击本脚本。
  echo 正在打开下载页面: https://git-scm.com/download/win
  start "" https://git-scm.com/download/win
  pause & exit /b 1
)

rem ---- 清理两个过时文件（内容已改名为 README.md 和 docs\DOS和COHP计算流程.md）
if exist "使用说明.md" del /q "使用说明.md"
if exist "docs\DOS和COHP（修订版）.md" del /q "docs\DOS和COHP（修订版）.md"

rem ---- 提交者身份（只需设置一次）
for /f "delims=" %%i in ('git config --global user.name') do set GN=%%i
if "%GN%"=="" (
  set /p GN=请输入提交者名字（如你的 GitHub 用户名）: 
  call git config --global user.name "%%GN%%"
)
for /f "delims=" %%i in ('git config --global user.email') do set GE=%%i
if "%GE%"=="" (
  set /p GE=请输入提交者邮箱（GitHub 注册邮箱或 用户名@users.noreply.github.com）: 
  call git config --global user.email "%%GE%%"
)

rem ---- 建立本地仓库并提交
if not exist ".git" (
  git init -b main || (echo git init 失败 & pause & exit /b 1)
)
git add -A
git diff --cached --quiet
if errorlevel 1 (
  git commit -q -m "CP2K BondDOS v1.0.0" -m "DOS/PDOS and k-point COHP/COOP/COBI analysis for CP2K 2026.2+ outputs, with charges, atom-resolved PDOS, LOBSTER-format export, a PySide6 GUI and a command line interface." -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" -m "Claude-Session: https://claude.ai/code/session_01DueWTEkSU3ei8o8LfpPzSv" || (echo 提交失败 & pause & exit /b 1)
  echo 已在本地提交 v1.0.0
) else (
  echo 没有新的改动需要提交
)

rem ---- 远程仓库地址
echo.
echo 请先在 https://github.com/new 新建一个空仓库：
echo   Repository name 填 CP2K_BondDOS，选 Private，不要勾选 "Add a README file"
echo 建好后复制页面上的地址，形如 https://github.com/你的用户名/CP2K_BondDOS.git
start "" https://github.com/new
echo.
set /p URL=把仓库地址粘贴到这里后回车: 
if "%URL%"=="" (echo 没有输入地址 & pause & exit /b 1)
git remote get-url origin >nul 2>nul
if errorlevel 1 (git remote add origin "%URL%") else (git remote set-url origin "%URL%")

rem ---- 推送（第一次会弹出浏览器窗口让你登录 GitHub 授权）
echo.
echo 正在上传……如果弹出 GitHub 登录窗口，请选择 "Sign in with your browser" 并授权。
git push -u origin main
if errorlevel 1 (
  echo.
  echo [上传失败] 常见原因：
  echo   1. 仓库不是空的（建仓库时勾选了 README）——删掉重建一个空仓库再运行本脚本；
  echo   2. 登录没有完成或网络不通——重新运行本脚本即可；
  echo   3. 地址填错——重新运行，粘贴正确地址。
  pause & exit /b 1
)
echo.
echo ============================================================
echo   上传成功！以后有改动时双击 "同步更新到GitHub.bat" 即可。
echo ============================================================
start "" "%URL:.git=%"
pause
