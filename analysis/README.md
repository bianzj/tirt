# `analysis/` 结果分析

这里放面向实验结果的离线分析脚本，不参与模型核心计算，也不作为日常运行入口。

当前脚本：

```powershell
python analysis/plot_observations.py --result-dir cases/<case_name>
```

分析结果建议写入对应案例的 `plots/` 或 `analysis/` 子目录，避免把中间文件放到项目根目录。通用方向性绘图请使用 `plot/`。
