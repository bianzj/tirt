# `cases/` 案例与结果

每个案例目录保存一组可复现输入及其结果，推荐结构如下：

```text
cases/my_case/
├── input.csv          # TiRT 输入，或 TIRTEB 运行时保存的输入快照
├── output.csv         # 主结果
├── observations.csv   # TIRTEB 观测结果（如有）
└── plots/             # 可重新生成的图像
```

`input.csv` 和必要的外部数据可以纳入版本管理；`output.csv`、`observations.csv`、`intermediate.csv`、`summary.csv` 和 `plots/` 都是可重新生成的产物，已由 `.gitignore` 排除。

目录名以下划线开头的目录（如 `_verify_run`）由验证脚本生成，不作为正式案例，整体由 `.gitignore` 排除。
