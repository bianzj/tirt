#!/bin/bash

set -e
cd "$(dirname "$0")"

if command -v conda >/dev/null 2>&1; then
    CONDA_BIN="$(command -v conda)"
elif [ -x "$HOME/miniconda3/bin/conda" ]; then
    CONDA_BIN="$HOME/miniconda3/bin/conda"
elif [ -x "$HOME/anaconda3/bin/conda" ]; then
    CONDA_BIN="$HOME/anaconda3/bin/conda"
else
    echo "未找到 Miniconda 或 Anaconda。"
    echo "请先安装 Miniconda，然后重新运行此文件。"
    read -r -p "按回车键退出..."
    exit 1
fi

CONDA_ROOT="$(dirname "$(dirname "$CONDA_BIN")")"
source "$CONDA_ROOT/etc/profile.d/conda.sh"

if ! conda env list | awk '{print $1}' | grep -qx tirt; then
    echo "首次运行，正在创建 TiRT 环境..."
    conda env create -f environment.yml
fi

conda activate tirt
export TIRT_PYTHON="$(command -v python)"

if [ ! -x "node_modules/electron/dist/Electron.app/Contents/MacOS/Electron" ]; then
    npm install --no-audit --no-fund
fi
npm run desktop
