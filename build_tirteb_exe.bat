@echo off
setlocal
cd /d "%~dp0"

set "PYTHON_EXE=C:\work\miniconda\envs\python311\python.exe"
set "PYTHON_DIR=C:\work\miniconda\envs\python311"
set "PATH=%PYTHON_DIR%;%PYTHON_DIR%\Library\bin;%PYTHON_DIR%\Scripts;%PATH%"

if not exist "%PYTHON_EXE%" (
    echo Python not found: %PYTHON_EXE%
    pause
    exit /b 1
)

"%PYTHON_EXE%" -m PyInstaller ^
  --noconfirm ^
  --clean ^
  --onefile ^
  --windowed ^
  --paths "." ^
  --hidden-import tirteb_run ^
  --hidden-import tirteb_server ^
  --collect-submodules base ^
  --name TIRTEB ^
  --exclude-module matplotlib ^
  --exclude-module PyQt6 ^
  --exclude-module PySide6 ^
  --exclude-module tkinter ^
  --add-data "gui;gui" ^
  --add-data "data;data" ^
  --add-data "input_tirteb.csv;." ^
  gui\tirteb_launcher.py

if errorlevel 1 (
    echo Build failed.
    pause
    exit /b 1
)

echo Built: %CD%\dist\TIRTEB.exe
pause
