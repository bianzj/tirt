# TIRTEB 测试与打包报告

## 测试环境

- Windows 10/11，C:\work\miniconda\envs\python311\python.exe
- Python 3.11
- NumPy、SciPy、Matplotlib、Selenium
- Firefox：C:\Program Files\Mozilla Firefox\firefox.exe

## 核心运行测试

scripts/verify_tirteb.py 已通过：

- 默认案例和 12 个 surface × vegetation 组合。
- 单坡与复合地形。
- 方向文件与波长文件模式。
- 指定 DOY、日内小时筛选和 intermediate.csv。
- soil_heat_method × stomatal_method × temperature_method 的 12 组比较，汇总 48 行。
- 所有 output.csv 和 observations.csv 均有数据，主要数值均为有限值。

运行命令：

~~~bash
python scripts/verify_tirteb.py
~~~

中断后可以复用已经通过的案例：

~~~bash
python scripts/verify_tirteb.py --resume
~~~

## GUI 测试

tests/test_tirteb_gui.py 覆盖页面、健康接口、默认运行、平面/地形/城市场景、方向/波长/叶片/土壤/城市文件上传和指定观测时间：

~~~bash
python tests/test_tirteb_gui.py --mode source --port 8901
python tests/test_tirteb_gui.py --mode exe --port 8902
~~~

tests/test_tirteb_firefox.py 使用已安装 Firefox 进行页面加载、运行和截图验证。

## 打包结果

~~~text
dist/TIRTEB.exe    约 48 MB，单文件
~~~

默认 `gui/`、`data/` 和 `input_tirteb.csv` 已嵌入 EXE。首次启动会解压到临时目录，退出后由系统清理。测试确认 EXE 能提供 `/time.html`，并能读取内置默认数据。

## 视频

scripts/make_tirteb_videos.py 生成 4 段：理论概览、GUI 快速上手、场景和文件输入、结果导出与分发。视频由实际 GUI 截图和结果页面生成，文件位于 videos。
