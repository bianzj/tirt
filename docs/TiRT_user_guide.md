# TiRT 操作文档

## 1. 运行方式

### 1.1 直接运行发布版

给别人使用时，只需要提供：

```text
C:\work\tirt\dist\TiRT.exe
```

`TiRT.exe` 已内置默认 `gui/`、`data/` 和 `input.csv`。双击后会启动本地服务并打开浏览器界面。

默认地址：

```text
http://127.0.0.1:8765/
```

支持 Firefox、Chrome、Edge 等现代浏览器；已用本机 Firefox 完成运行测试。

### 1.2 源码运行

```bat
cd C:\work\tirt
C:\work\miniconda\envs\python311\python.exe gui\server.py --port 8765
```

浏览器打开：

```text
http://127.0.0.1:8765/
```

### 1.3 命令行运行模型

```bat
cd C:\work\tirt
C:\work\miniconda\envs\python311\python.exe run.py --output cases\demo\output.csv
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

- `LAI`：叶面积指数。
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

```text
Wavelength bands: 8;10.5;12
Leaf: 0.975;0.985;0.980
Soil: 0.920;0.955;0.965
```

每个组分的发射率数量必须与波段数量一致。

### 6.2 文件输入

点击对应行右侧 `...` 选择文件。支持 txt/csv 文本。文件中的波段会被插值到当前目标波段。

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

- 平均亮温或平均输出量。
- 方向性变化范围。
- 平均辐亮度。
- 当前波段下观测方向数量。

## 9. 查看结果

结果图包括：

- `Polar directional map`：半球方向图。
- `Parallel response`：太阳主平面方向曲线。
- `Perpendicular response`：垂直太阳主平面方向曲线。

可切换：

- `Temperature (°C)`
- `Temperature (K)`
- `Radiance`

可切换场景编号和波段。

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

### 11.1 双击后没有页面

检查是否被防火墙或安全软件拦截。也可以手动打开：

```text
http://127.0.0.1:8765/
```

### 11.2 端口被占用

从命令行指定端口：

```bat
TiRT.exe --port 8766
```

### 11.3 输入后报错

优先检查：

- 波段数量与发射率数量是否一致。
- 角度是否为数字。
- 地形占据率是否过大。
- `surface_model` 和 `vegetation_model` 是否为支持项。

## 12. 分发方式

把下面目录压缩后发给别人：

```text
C:\work\tirt\dist\
```

当前 `dist` 目录内只有：

```text
TiRT.exe
```

默认数据已内置，无需额外携带 `data/` 或 `input.csv`。
