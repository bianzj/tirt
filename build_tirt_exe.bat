@echo off
setlocal
cd /d "%~dp0"

set "CONDA_ROOT=C:\work\miniconda"
set "PYTHON_EXE=%CONDA_ROOT%\python.exe"
if not exist "%PYTHON_EXE%" set "PYTHON_EXE=python"
if exist "%CONDA_ROOT%" set "PATH=%CONDA_ROOT%;%CONDA_ROOT%\Library\bin;%CONDA_ROOT%\Scripts;%PATH%"
set "MPLBACKEND=Agg"

"%PYTHON_EXE%" -m PyInstaller ^
  --noconfirm ^
  --clean ^
  --onefile ^
  --name TiRT ^
  --exclude-module PyQt6 ^
  --exclude-module PySide6 ^
  --add-data "gui;gui" ^
  --add-data "data;data" ^
  --add-data "input.csv;." ^
  gui\launcher.py

if errorlevel 1 (
    echo Build failed.
    pause
    exit /b 1
)

echo Built: %CD%\dist\TiRT.exe
pause
