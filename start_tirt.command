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
echo "TiRT 正在启动: http://127.0.0.1:8765/"

python gui/server.py --host 127.0.0.1 --port 8765 &
SERVER_PID=$!
trap 'kill "$SERVER_PID" 2>/dev/null || true' EXIT

sleep 1
open "http://127.0.0.1:8765/"
wait "$SERVER_PID"
