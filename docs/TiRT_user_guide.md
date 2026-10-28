# TiRT 操作文档

## 1. 运行方式

### 1.1 安装版

运行安装包 `TiRT-Setup-1.2.0-x64.exe`，可以选择安装目录，安装后桌面和开始菜单出现 `TiRT` 快捷方式。程序不需要另外安装 Python。

启动 `TiRT` 后打开独立桌面窗口：默认是方向性模型（TIRT-Static），点击 `TIRT-Dynamic` 进入能量平衡时间序列。关闭窗口会同时结束后台服务。

安装目录中另有：

- `docs/`：操作说明、理论说明、测试文档（与 HiSTREAM 对比）和项目代码说明，均为 PDF；开始菜单中的 `TiRT 文档` 快捷方式打开这个文件夹。
- `TiRT-1.2.0-source.zip`：代码压缩包（模型源代码、输入模板、数据和文档），用于命令行运行或二次开发。

### 1.2 源码运行

解压源代码后，在项目根目录安装依赖（`pip install -r requirements.txt`），然后：

```bat
python gui\server.py --port 8765
```

浏览器打开：

```text
http://127.0.0.1:8765/
```

### 1.3 命令行运行模型

```bat
python run.py --output cases\demo\output.csv
```

## 2. GUI 页面区域

GUI 由三部分组成：

- 顶部工具栏：运行模拟、保存输出、保存输入。
- 左侧参数面板：地表结构、观测角度、光谱、温度。
- 右侧结果区：三维场景预览、方向性结果图和指标。

## 3. 场景结构设置

### 3.1 Surface mode

选择地表类型：

- `Plane`：平面地表。
- `Slope`：单坡地表。
- `Terrain`：复合地形。
- `Urban`：城市建筑与街道。

### 3.2 Vegetation type

选择植被类型：

- `Bare soil`：裸土或无植被表面。
- `Turbid veg`：均质植被。
- `Crop rows`：垄行作物。
- `Forest crowns`：离散树冠。

不同组合会自动显示对应参数。

## 4. 结构参数

### 4.1 坡面参数

- `Slope angle`：坡度。
- `Slope aspect`：坡向。

### 4.2 地形参数

- `Terrain height`：地形单元高度。
- `Terrain radius`：地形单元半径。
- `Terrain density`：地形单元密度。

### 4.3 城市参数

- `Building density`：建筑密度。
- `Building height`：建筑高度。
- `Building length`：建筑长度。
- `Building width`：建筑宽度。
- `Building azimuth`：建筑方位角。

### 4.4 植被参数

- `LAI`：叶面积指数。可输入一个值，也可输入分号分隔的多个值，例如 `1;3;5`；多值时结果区会出现 LAI 下拉，主平面曲线叠加显示各 LAI。
- `Hotspot`：热点参数。
- `Row width/gap/height/azimuth`：垄行结构。
- `Crown rad-a/rad-b/height/density`：树冠结构。

## 5. 角度设置

### 5.1 太阳角度

- `Solar zenith`：太阳天顶角。
- `Solar azimuth`：太阳方位角。

### 5.2 自动方向

勾选：

- `Principal plane`：生成太阳主平面和垂直主平面方向。
- `Hemisphere`：追加半球观测方向。

`Principal-plane step` 控制天顶角间隔，`Maximum VZA` 控制最大观测天顶角。

### 5.3 手动方向

选择 `Manual` 后输入：

```text
0/0; 10/90; 20/180
```

每组为：

```text
VZA/VAA
```

## 6. 光谱设置

### 6.1 直接输入

`Spectra` 分组是一张波段参数表，光学和热红外共用，每行用分号给出每个波段一个值：

```text
Wavelength µm: 0.67;0.86;10.5
Leaf ρ:        0.05;0.45;0.015
Leaf τ:        0.03;0.40;0.0
Soil ρ:        0.12;0.25;0.045
Street ρ:      0.12;0.18;0.045
Wall ρ:        0.20;0.25;0.07
Roof ρ:        0.15;0.20;0.05
```

小于 2.5 µm 的波段输出方向反射率（默认 80% 直射 + 20% 漫射；`optical,shortwave_diffuse_fraction,0` 可改为纯直射 BRF）；大于等于 2.5 µm 的波段输出亮温，发射率由同一列参数换算：叶片 `ε = 1 − ρ − τ`，土壤、街道、墙壁、屋顶 `ε = 1 − ρ`。表格下方实时显示各热红外波段换算出的发射率，便于与光学参数相互核对。每行的值数量必须与波段数量一致；叶片反射率 + 透过率不能超过 1。

### 6.2 文件输入

波段行右侧的 `...` 可读取 txt/csv 文件，文件按同样顺序共 7 行：波段、叶片反射率、叶片透过率、土壤反射率、街道反射率、墙壁反射率、屋顶反射率。每行可以带文字标签，`#` 之后为注释，例如：

```text
波段       0.55;0.67;0.86;8.6;10.5
叶片反射率 0.10;0.05;0.45;0.02;0.015
叶片透过率 0.08;0.03;0.40;0;0
土壤反射率 0.10;0.12;0.25;0.08;0.045
街道反射率 0.10;0.12;0.18;0.06;0.045
墙壁反射率 0.20;0.20;0.25;0.08;0.07
屋顶反射率 0.15;0.15;0.20;0.07;0.05
```

读取后各行填入表格，可以继续手动修改。示例文件见 `data/spectrals/optical_properties.txt`。

## 7. 温度设置

每个组分分别设置：

- `Sunlit °C`：日照温度。
- `Shaded °C`：阴影温度。

城市场景会显示屋顶、墙体、街道温度；地形场景会显示地形温度；植被场景会显示叶片温度。

## 8. 运行模拟

点击：

```text
Run simulation
```

运行完成后顶部状态变为 `Complete`，结果区显示：

- 平均反射率（短波）或平均亮温（热红外）。
- 方向性变化范围。
- 当前波段及类型（Reflectance BRF / Brightness temperature °C）。
- 当前波段下观测方向数量。

## 9. 查看结果

结果图包括：

- `Polar directional map`：半球方向图。
- `Parallel response`：太阳主平面方向曲线。
- `Perpendicular response`：垂直太阳主平面方向曲线。

只保留波段选择（多 LAI 时另有 LAI 切换）。波段决定显示量：

- 小于 2.5 µm：显示反射率 BRF。
- 大于等于 2.5 µm：显示亮温 °C。

不再提供 Quantity 下拉，也不显示开尔文亮温或辐亮度。

## 10. 保存文件

### 10.1 保存输出

点击：

```text
Save observations.csv
```

保存当前模拟得到的方向性观测表。

### 10.2 保存输入

点击：

```text
Save input.csv
```

保存当前 GUI 参数为输入文件，便于复现实验。

## 11. 常见问题

### 11.1 窗口提示服务启动失败

检查是否被防火墙或安全软件拦截。后台服务只监听本机地址 `127.0.0.1`，不对外提供服务。

### 11.2 输入后报错

优先检查：

- 波段参数表每行的值数量是否与波段数量一致。
- 角度是否为数字。
- 地形占据率 `terrain_density·π·terrain_radius²`、建筑占地比例 `building_density·length·width` 是否小于 1。
- `surface_model` 和 `vegetation_model` 是否为支持项。
- 单坡场景中位于坡面背面的观测方向不会输出，计算结果会提示被省略的方向数。
