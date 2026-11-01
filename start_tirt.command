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
# launcher.py searches for a free port from 8765 and opens the browser once the server is listening.
echo "TiRT 正在启动，浏览器会在服务就绪后自动打开；按 Ctrl+C 或关闭本窗口退出。"
set +e
python gui/launcher.py --host 127.0.0.1 --port 8765
STATUS=$?
if [ "$STATUS" -ne 0 ]; then
    echo "TiRT 启动失败，请查看上方错误信息。"
    read -r -p "按回车键退出..."
fi
exit "$STATUS"
