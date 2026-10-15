@echo off
setlocal
cd /d "%~dp0"

set "CONDA_BAT="
for /f "delims=" %%i in ('where conda.bat 2^>nul') do if not defined CONDA_BAT set "CONDA_BAT=%%i"
if not defined CONDA_BAT if exist "%USERPROFILE%\miniconda3\condabin\conda.bat" set "CONDA_BAT=%USERPROFILE%\miniconda3\condabin\conda.bat"
if not defined CONDA_BAT if exist "%USERPROFILE%\anaconda3\condabin\conda.bat" set "CONDA_BAT=%USERPROFILE%\anaconda3\condabin\conda.bat"

if not defined CONDA_BAT (
    echo 未找到 Miniconda 或 Anaconda。
    echo 请先安装 Miniconda，然后重新运行此文件。
    pause
    exit /b 1
)

call "%CONDA_BAT%" run -n tirt python --version >nul 2>&1
if errorlevel 1 (
    echo 首次运行，正在创建 TiRT 环境...
    call "%CONDA_BAT%" env create -f environment.yml
    if errorlevel 1 (
        echo 环境创建失败。
        pause
        exit /b 1
    )
)

call "%CONDA_BAT%" activate tirt
if errorlevel 1 (
    echo TiRT 环境激活失败。
    pause
    exit /b 1
)

echo TiRT 正在启动: http://127.0.0.1:8765/
start "TiRT" http://127.0.0.1:8765/
python gui\server.py --host 127.0.0.1 --port 8765
pause
