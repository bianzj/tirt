# TIRT 复现与模型对比

## 直接复现 TIRT

在项目根目录运行：

```powershell
python validation/validate_current_model.py
# 可选：同时生成逐场景对比图
python validation/validate_current_model.py --plots
# 修正物理算法后，保留旧冻结基准，显式记录变化并与同一 HiSTREAM 结果比较
python validation/validate_current_model.py --allow-model-changes --output validation/results/physics_corrections_20261026
```

也可以在本目录运行 `python validate_current_model.py`。默认输入与输出路径根据代码位置解析，不依赖当前工作目录。

- 主程序：`validate_current_model.py`。
- 输入：`reference/inputs/` 中的 10 个场景 JSON。
- 已保存的 HiSTREAM 结果：`reference/histream_reference.csv`。
- 已保存的 TIRT 结果：同一目录的 `all_results.csv`。
- 默认输出：`results/`；可用 `--output <目录>` 指定其他位置。

此入口只重算 TIRT，读取已有 HiSTREAM 结果作对比，不启动 HiSTREAM 或 StreamSim。覆盖 670 nm、860 nm 和 10.5 μm，共 1,830 条请求；排除 24 条坡面背向观测后有 1,806 条有效配对。复现检查使用绝对容差 `1e-6`、相对容差 `0`，它不是 TIRT 相对 HiSTREAM 的物理误差阈值。

2026-10-26 的物理修正有意改变标准短波结果。默认严格复现模式仍会报告相对旧预测的变化；`--allow-model-changes` 保留输入/参考哈希检查并输出当前误差，`verification.json` 明确区分 `prediction_replay_passed` 与当前计算完成，不会覆盖冻结参考。修正内容、共同收敛及精度限制见 [修正报告](../docs/validation/physics_corrections_20261026/report.md)。

输出包括 `all_results.csv`（TIRT 结果）、`reference_comparison.csv`（逐方向对比）、`accuracy.csv`（误差统计）与 `verification.json`（复现回执）。加 `--plots` 时另生成 `figures/`。运行需要 NumPy、pandas；绘图另需 matplotlib。

`streamsim/` 保留原来的参考数据和引擎案例；`scope/` 保留 SCOPE 对比代码。原命令 `python scripts/validate_current_model.py` 作为兼容入口转发到此处。

## TIRTEB 与 TiRT 方向观测一致性

`compare_tirt_tirteb.py` 在相同场景、组分温度、观测几何和波段光学下，逐方向比较 TIRTEB 的短波反射率、热红外辐亮度和亮温与 TiRT 的 `run._surface_result` / `run._optical_band`。覆盖 plane、single slope、compound terrain、urban 与 bare/hom/row/crown 的组合；默认容差 `1e-9`，`--full` 增加天顶角和坡面背面方向。

```powershell
python validation/compare_tirt_tirteb.py --full
```

## 多山脊 DEM 对 HiSTREAM

`ridge_20261029/` 用 StreamStudio 的横向沙丘 DEM（201×201，100.5 m × 100.5 m，起伏 6.1 m）做多山脊对照：HiSTREAM eRaytracing 生成参考，TiRT 用 `surface_model=dem` 的逐坡面/坡向面元集合近似。10.5 µm 逐方向 RMSE 为 0.29–0.57 K，形态一致、整体偏低约 0.3–0.6 K。说明和复现命令见 [ridge_20261029/README.md](ridge_20261029/README.md)。

## 目录结构

```text
validation/
  validate_current_model.py  # TIRT 模拟与验证入口
  reference/                # 固定输入和已计算参考结果
    inputs/                 # 10 个场景 JSON
    histream_reference.csv  # 光学、热红外 HiSTREAM 参考值
    all_results.csv         # TIRT 已保存的复现基准
    accuracy.csv            # 已保存的误差统计
    manifest.json           # 场景清单和文件校验值
  results/                  # 每次运行产生的 CSV 和验证回执
  streamsim/                # 既有 HiSTREAM 案例与结果归档
  ridge_20261029/           # 多山脊 DEM 对 HiSTREAM 的对照
  scope/                    # 既有 SCOPE 对比
  compare_tirt_tirteb.py    # TIRTEB 方向观测与 TiRT 的逐方向一致性检查
```

`reference/` 从 `docs/validation/current_model_20261021/` 数据快照逐文件复制并核对哈希；原快照保留用于已发布报告追溯。日常验证读取本目录的 `reference/`，新输出写入 `results/`，不覆盖基准。在 IDE 中直接运行本入口即可验证全部场景，无需先启动界面或外部引擎。
