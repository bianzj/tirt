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

"%PYTHON_EXE%" gui\tirteb_launcher.py --host 127.0.0.1 --port 8901
if errorlevel 1 pause
