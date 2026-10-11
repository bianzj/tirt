# `tests/` 测试

测试覆盖核心物理和 GUI 服务接口（方向性、能量平衡和多角度）。

```powershell
python tests/test_model_invariants.py
python tests/test_urban_wall_radiation.py
python tests/test_tirteb_tirt_consistency.py
python tests/test_canopy_radiation.py
python tests/test_tirt_gui.py --mode source --port 8765
python tests/test_tirteb_gui.py --mode source --port 8901
```

`test_model_invariants.py` 对 16 种地表 × 植被组合检查：全角度结果有限；等温黑体在各波段返回自身温度；地形植被结果随发射率变化。

`test_tirteb_tirt_consistency.py` 用同一个 `input_tirteb.csv` 和相同组分温度分别运行 TiRT（`run.py`）和 TIRTEB 的热红外部分，要求各地表 × 植被组合的辐亮度一致。

`test_canopy_radiation.py` 检查 TIRTEB 冠层辐射的能量守恒：等温冠层的叶片和土壤长波吸收等于自身发射，黑色冠层吸收未透过的短波，通量权重覆盖整块地面和全部叶面积，场景 Rn = LE + H + G。

GUI 测试会自动启动 `gui/server.py`，只依赖 Python 标准库；`--mode exe` 改为测试打包后的 `dist/TiRT-service.exe`。完整回归优先使用 `scripts/verify_run.py` 和 `scripts/verify_tirteb.py`。
