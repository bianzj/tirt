@echo off
setlocal
cd /d "%~dp0"

rem Override with: set "PYTHON_EXE=<path to python.exe with PyInstaller>"
if not defined PYTHON_EXE set "PYTHON_EXE=C:\work\miniconda\envs\python311\python.exe"
if not exist "%PYTHON_EXE%" set "PYTHON_EXE=python"
for %%P in ("%PYTHON_EXE%") do if exist "%%~dpPLibrary\bin" set "PATH=%%~dpP;%%~dpPLibrary\bin;%%~dpPScripts;%PATH%"
set "MPLBACKEND=Agg"

if not exist "node_modules\electron\dist\electron.exe" (
    echo Installing desktop dependencies...
    call npm install --no-audit --no-fund
    if errorlevel 1 (
        echo npm install failed.
        pause
        exit /b 1
    )
)

"%PYTHON_EXE%" -m PyInstaller ^
  --noconfirm ^
  --clean ^
  --onefile ^
  --windowed ^
  --icon "%CD%\tirt.ico" ^
  --name TiRT-service ^
  --distpath "dist" ^
  --workpath "build\pyinstaller" ^
  --specpath "build" ^
  --exclude-module matplotlib ^
  --exclude-module PIL ^
  --exclude-module PyQt6 ^
  --exclude-module PySide6 ^
  --add-data "%CD%\gui;gui" ^
  --add-data "%CD%\data;data" ^
  --add-data "%CD%\input.csv;." ^
  --add-data "%CD%\input_tirteb.csv;." ^
  gui\launcher.py

if errorlevel 1 (
    echo Build failed.
    pause
    exit /b 1
)

call npm run pack:win -- --config.electronDist=node_modules/electron/dist
if errorlevel 1 (
    echo Desktop packaging failed.
    pause
    exit /b 1
)

for /f "delims=" %%v in ('node -p "require('./package.json').version"') do set "APP_VERSION=%%v"
echo Built installer: %CD%\release\TiRT-Setup-%APP_VERSION%-x64.exe
pause
