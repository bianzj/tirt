# 固定参考数据

本目录是 TIRT 快速验证的默认数据来源，冻结自 `docs/validation/current_model_20261021/`，复制时核对所有文件字节及 `manifest.json` 中的 SHA-256。

- `inputs/`：10 个复现专用场景 JSON，包含光学与热红外配置。
- `histream_reference.csv`：已计算的 HiSTREAM 三波段参考值。
- `all_results.csv`：TIRT 的已保存输出，作为数值复现基准。
- `accuracy.csv`：已保存的分场景误差统计。
- `manifest.json`：场景、波段、源代码记录及固定文件校验值。
- `reference_manifest.json`、`reference_protocols.json`：历史来源和运行参数记录，内部路径按原始档案保留。

在 TIRT 根目录执行 `python validation/validate_current_model.py`，或在 `validation/` 直接运行 `python validate_current_model.py`。结果写入 `validation/results/`。
