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
3. 在 Scene structure 选择平面、单坡、复合地形或城市，以及裸土、均质、垄行或树冠（复合地形和城市上的植被按树冠计算）。
4. 在 Spectral information 设置 PROSPECT、BSM、热红外发射率，或上传叶片/土壤可见光近红外光谱与城市光谱。
5. 在 Plant physiology 和 Thermal/energy balance 中调整方法和参数。
6. 点击 Run time series，等待状态变为 Complete。
7. 使用图表选择 Tv_C、LE、H 或 G，再保存 `output.csv`、`observations.csv` 和当前 `input_tirteb.csv`。

## 输入文件

`input_tirteb.csv` 与方向性模型的 `input.csv` 使用同一套参数名称和含义（地表、坡面、地形、建筑、垄行、树冠、热红外发射率、观测几何和波段），参数说明见 [TIRTEB 理论与方法说明](TIRTEB_theory_and_methods.md#参数体系)。

方向文件（`geometry_mode=2` 时的 `geometry_file`）每行至少包含：

~~~text
vza vaa
~~~

波长文件每行一个微米值，例如：

~~~text
0.55
0.86
10.5
~~~

叶片可见光近红外文件（`leaf_optics_source=file`）可使用三列 wavelength、reflectance、transmittance；土壤文件（`soil_optics_source=file`）使用两列 wavelength、reflectance。热红外发射率光谱文件与 TiRT 相同：`vegetation_spectrum_file`、`soil_spectrum_file` 为 wavelength、emissivity 两列；`urban_spectrum_file` 为波长加屋顶/墙面/街道发射率四列，或七列（波长、屋顶/墙面/街道可见光近红外反射率、屋顶/墙面/街道发射率）。文件波长可以用 µm，也可以使用 nm；程序会自动识别小于等于 20 的波长为 µm。

## 绘图

安装 matplotlib 后：

~~~bash
python analysis/plot_observations.py --result-dir cases/my_case --time 214.5 --bands 0p55 0p86 10p5
~~~

图件写入案例目录的 `plots/`。0p55、0p86 和 10p5 是结果列名中的波长标签。

## 常见问题

- 端口被占用：使用 --port 8890，再打开对应地址。
- 没有观测记录：检查 observation_time_mode 和气象文件中的实际时间是否匹配。
- 叶片或土壤光谱文件报错：至少提供两行有效数值，并保证波长列和光谱列顺序正确。
- 复合地形报错“总占据率必须小于 1”：`terrain_density` 是每平方米的圆锥数，需满足 `terrain_density·π·terrain_radius² < 1`。
- 需要调试能量平衡：设置 output_intermediate=yes，检查 intermediate.csv。
