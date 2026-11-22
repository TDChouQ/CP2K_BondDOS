@echo off
chcp 65001 >nul
setlocal EnableExtensions
cd /d "%~dp0"
title CP2K BondDOS - 同步更新到 GitHub
where git >nul 2>nul || (echo 未安装 Git，请先运行 "上传到GitHub.bat" & pause & exit /b 1)
if not exist ".git" (echo 还没有初始化仓库，请先运行 "上传到GitHub.bat" & pause & exit /b 1)

git add -A
git status --short
git diff --cached --quiet
if not errorlevel 1 (
  echo 没有新的改动，尝试同步远程……
  git push
  pause & exit /b 0
)
set MSG=
set /p MSG=请输入这次更新的说明（直接回车则用 "更新"）: 
if "%MSG%"=="" set MSG=更新
git commit -q -m "%MSG%" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" || (echo 提交失败 & pause & exit /b 1)
git push || (echo 上传失败：检查网络后重新运行本脚本 & pause & exit /b 1)
echo.
echo 已同步到 GitHub。
pause
