# TiRT 测试汇总

## 已执行测试

- 源码 GUI：`tests/test_tirt_gui.py --mode source`
- 单文件 exe GUI：`tests/test_tirt_gui.py --mode exe`
- `run.py` 与 `input.csv` 组合验证：`scripts/verify_run.py`

## GUI 覆盖项

- `/api/health` 健康检查。
- GUI 首页访问。
- `app.js`、`styles.css`、Three.js、OrbitControls 静态资源访问。
- 默认模拟运行。
- 输出字段完整性检查。
- 多个地表/植被组合。
- 手动观测角度。
- 上传观测几何文本。
- 上传光谱文本。
- 非法输入错误处理。

## 结果文件

- `docs/gui_test_report_source.md`
- `docs/gui_test_report_exe.md`


## Firefox 兼容测试

- 本机 Firefox 页面加载：通过。
- 本机 Firefox 点击 `Run simulation`：通过。
- 结果区显示：通过。
- 截图：`docs/firefox_results.png`。
