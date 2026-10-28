# 分场景案例验证记录（2026-10-24）

- 入口：[案例说明](../../examples/README.md)、[全部运行](../../examples/run_all.py)。
- 10 个光学、16 个热红外案例成功；2,138 条有效记录；52 张 PNG。
- 52 张图已通过六张总览逐一查看，中文标题、图例、波长/角度和纵轴单位可读。
- 9 项针对性测试通过，覆盖裸地极限、等温黑体、角度符号、光谱参数检查、
  垄行分派、部分/全部坡面背面缺测、失败验证清除旧成功回执。未修改核心模型。
- 源文件哈希、参数快照、CSV 条目数、温度换算、缺测字段、PNG 可读性检查通过。
- 独立复核发现“热红外整个光谱方向不可见”会中断案例，已修复并加入回归测试。
- 所有 26 个案例在带空格的干净解压路径、不同工作目录下重新运行成功，
  52 个 CSV 的所有字段与原结果完全相同；另直接调用光学、热红外独立入口各一次。

当前公开光学接口不支持坡面/复合地形的 bare、hom、row 共 6 种组合；
这些组合列于覆盖表与批处理报告，未替换成其他植被类型。

## 原有模型复算

按用户指示，将 `scripts/validate_current_model.py` 中数值和辐亮度复算的
绝对误差阈值统一为 `1e-6`，相对误差阈值仍为 0。
原有 Planck 正反转换的内部检查未改。

- 10 个场景、3 个波段；1,830 条请求、1,806 条有效预测/参考配对。
- 最大绝对预测复算差异：`7.555229331046576e-9`。
- 当前阈值：`1e-6`；检查通过；输入/参考文件校验通过。
- 这是已有参考结果复算，并非重新运行 HiSTREAM 或增加独立精度验证。

## 本地产物

- `examples/results/batch_all.json`：实际批处理状态与未支持列表。
- `examples/results/outputs_verified.json`：输出完整性记录。
- `cases/_scenario_examples_regression_20261024/verification.json`：模型复算记录。
- `cases/_scenario_examples_regression_20261024/contact_1.png` 至 `contact_6.png`：图面检查总览。
- `release/scenario_examples_20261024/TiRT_scenario_examples_20261024.zip`：
  案例脚本、必要模型源码、依赖清单、参考结果和测试（约 5.8 MB）。
- 同目录 `extraction_verification.json`、`SHA256SUMS.txt`：干净解压复算记录和包校验值。

压缩包 SHA-256：`1ec637481406bd9e1b5f7eeda312ea29f9daeccc44a367a292d024ab0ed17b7f`。
