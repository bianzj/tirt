# TiRT 能量平衡操作说明

## 直接使用 EXE

双击 `dist/TiRT.exe`，程序会启动统一的 TiRT 本地服务并打开：

~~~text
http://127.0.0.1:8901/time.html
~~~

给别人使用时，直接分发 `dist/TiRT.exe` 即可。单文件程序已嵌入 `gui/`、`data/` 和默认能量平衡输入，接收方不需要安装 Python、NumPy 或 SciPy。若需要修改默认数据，改用源代码方式运行，并保留项目根目录的输入文件、`data/`、`base/` 和 `gui/`。

## 源代码运行

~~~bat
start_tirt.bat
~~~

或：

~~~bash
python gui/server.py --port 8765
~~~

浏览器打开 http://127.0.0.1:8765/time.html。命令行运行单个案例：

~~~bash
python tirteb_run.py --input input_tirteb.csv --output cases/my_case/output.csv
~~~

## GUI 操作顺序

1. 在 Site and period 设置经纬度、DOY 范围和气象文件。
2. 在 Observation 设置方向、波长和观测时间筛选。
3. 在 Scene structure 选择平面、地形或城市，以及裸土、均质、垄行或树冠。
4. 在 Spectral information 设置 PROSPECT、BSM、热红外参数，或上传 measured 光谱。
5. 在 Plant physiology 和 Thermal/energy balance 中调整方法和参数。
6. 点击 Run time series，等待状态变为 Complete。
7. 使用图表选择 Tv_C、LE、H 或 G，再保存 `output.csv`、`observations.csv` 和当前 `input_tirteb.csv`。

## 输入文件

方向文件每行至少包含：

~~~text
vza vaa
~~~

波长文件每行一个微米值，例如：

~~~text
0.55
0.86
10.5
~~~

叶片 measured 文件可使用三列 wavelength、leaf_reflectance、leaf_transmittance；土壤 measured 文件使用两列 wavelength、soil_reflectance。城市 measured 文件使用七列：波长、屋顶/墙面/街道反射率、屋顶/墙面/街道发射率。文件波长可以用 um，也可以使用 nm；程序会自动识别小于等于 20 的波长为 um。

## 绘图

安装 matplotlib 后：

~~~bash
python analysis/plot_observations.py --result-dir cases/my_case --time 214.5 --bands 0p55 0p86 10p5
~~~

图件写入案例目录的 `plots/`。0p55、0p86 和 10p5 是结果列名中的波长标签。

## 常见问题

- 端口被占用：使用 --port 8890，再打开对应地址。
- 没有观测记录：检查 observation_time_mode 和气象文件中的实际时间是否匹配。
- measured 文件报错：至少提供两行有效数值，并保证波长列和光谱列顺序正确。
- 需要调试能量平衡：设置 output_intermediate=yes，检查 intermediate.csv。
