# TIRTEB 视频介绍脚本

## 01_theory_overview.mp4

展示模型输入、太阳几何、辐射传输、能量平衡迭代和方向性亮温输出。适合理论和方法介绍。

## 02_gui_quick_start.mp4

展示启动 EXE、打开本地 GUI、设置时间范围、运行时间序列、读取图表并导出 CSV。

## 03_scenes_and_files.mp4

展示平面、地形、城市、裸土、均质、垄行和树冠场景，以及方向文件、波长文件和 measured 光谱上传。

## 04_results_and_distribution.mp4

展示时间序列指标、方向性观测结果、绘图脚本、CSV 导出和 dist/TIRTEB.exe 分发方式。

视频生成：

~~~bash
python scripts/capture_tirteb_screenshots.py
python scripts/make_tirteb_videos.py
~~~
