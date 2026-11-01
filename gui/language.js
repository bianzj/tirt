const STORAGE_KEY = "tirt-language";
const language = localStorage.getItem(STORAGE_KEY) === "zh-CN" ? "zh-CN" : "en-US";
const sourceText = new WeakMap();
const sourceAttributes = new WeakMap();

const ZH = Object.freeze({
  // Page frame
  "TIRT-Static 1.2.0 — Static scene": "TIRT-Static 1.2.0 — 静态场景",
  "TIRT-Dynamic 1.2.0 — Dynamic scene": "TIRT-Dynamic 1.2.0 — 动态场景",
  "Static scene": "静态场景",
  "Dynamic scene": "动态场景",
  "TERRESTRIAL INHOMOGENEOUS RADIATIVE TRANSFER": "陆表非均质辐射传输",
  "LAND SURFACE ENERGY BALANCE": "陆表能量平衡",
  "Version 1.2.0": "版本 1.2.0",
  "Produced by TIRTLab · Email:": "出品：TIRTLab　邮箱：",
  "Switch to Chinese": "切换为中文",
  "Switch to English": "切换为英文",
  "Ready": "就绪",
  "Running": "运行中",
  "Running...": "运行中...",
  "Running {seconds} s": "运行中 {seconds} 秒",
  "Complete": "完成",
  "Error": "错误",
  "Input error": "输入错误",
  "Exiting": "正在退出",
  "Scene unavailable": "场景不可用",
  "Run simulation": "运行模拟",
  "Run time series": "运行时间序列",
  "Save observations.csv": "保存 observations.csv",
  "Save processes.csv": "保存 processes.csv",
  "Save input.csv": "保存 input.csv",
  "Save inputs (.zip)": "保存输入 (.zip)",
  "Exit": "退出",
  "Exit TiRT": "退出 TiRT",
  "CONTROL DECK": "控制面板",
  "TIME DRIVER": "时间驱动",
  "Input settings": "输入设置",
  "Reset settings": "重置设置",

  // Scene structure
  "Scene structure": "场景结构",
  "Surface mode": "地表类型",
  "Vegetation type": "植被类型",
  "Plane": "平地",
  "Slope": "坡地",
  "Terrain": "地形",
  "DEM / ridges": "DEM / 多山脊",
  "Urban": "城市",
  "Bare soil": "裸土",
  "Turbid veg": "均质植被",
  "Crop rows": "垄行作物",
  "Crown forest": "离散森林",
  "Surface structure": "地表结构",
  "Slope angle": "坡度",
  "Slope aspect": "坡向",
  "Slope angle (°)": "坡度 (°)",
  "Slope aspect (°)": "坡向 (°)",
  "Slope (°)": "坡度 (°)",
  "Aspect (°)": "坡向 (°)",
  "Terrain height (m)": "地形高度 (m)",
  "Terrain radius (m)": "地形半径 (m)",
  "Terrain density": "地形密度",
  "Terrain density (m⁻²)": "地形密度 (m⁻²)",
  "DEM file (.asc)": "DEM 文件（.asc）",
  "Slope bin (°)": "坡度分箱 (°)",
  "Aspect bin (°)": "坡向分箱 (°)",
  "Cast shadow": "投影阴影",
  "Building density": "建筑密度",
  "Building height (m)": "建筑高度 (m)",
  "Building length (m)": "建筑长度 (m)",
  "Building width (m)": "建筑宽度 (m)",
  "Building azimuth (°)": "建筑方位角 (°)",
  "Vegetation structure": "植被结构",
  "Hotspot": "热点参数",
  "Canopy height (m)": "冠层高度 (m)",
  "Leaf width (m)": "叶片宽度 (m)",
  "Row width": "垄宽",
  "Row gap": "垄间距",
  "Row height": "垄高",
  "Row azimuth": "垄方位角",
  "Row width (m)": "垄宽 (m)",
  "Row gap (m)": "垄间距 (m)",
  "Row height (m)": "垄高 (m)",
  "Row azimuth (°)": "垄方位角 (°)",
  "Crown rad-a": "树冠半径 a",
  "Crown rad-b": "树冠半径 b",
  "Crown rad-a (m)": "树冠半径 a (m)",
  "Crown rad-b (m)": "树冠半径 b (m)",
  "Trunk height (m)": "树干高度 (m)",
  "Crown height": "树干高度",
  "Crown density": "树冠密度",
  "Crown density (m⁻²)": "树冠密度 (m⁻²)",
  "Crown base height above the ground; sets how crown shadows overlap the visible soil.": "树冠底部离地高度，决定树冠阴影与可见土壤的重叠方式。",
  "Bare soil selected: vegetation structure is not used.": "已选择裸土：不使用植被结构参数。",
  "Thapke": "Hapke 方向性",
  "Off": "关闭",
  "On": "开启",
  "Anisotropy K": "各向异性 K",

  // Angles and observation
  "Angles": "角度",
  "Observation": "观测",
  "Solar zenith": "太阳天顶角",
  "Solar azimuth": "太阳方位角",
  "SZA": "太阳天顶角",
  "SAA": "太阳方位角",
  "Principal-plane step (°)": "主平面步长 (°)",
  "Maximum VZA (°)": "最大观测天顶角 (°)",
  "Principal plane": "主平面",
  "Hemisphere": "半球",
  "Additional observation angles": "附加观测角度",
  "File": "文件",
  "Manual": "手动",
  "Manual VZA/VAA pairs": "手动 VZA/VAA 角度对",
  "Geometry file": "几何文件",
  "Select geometry file": "选择几何文件",
  "Observation angle file": "观测角度文件",
  "Select observation angle file": "选择观测角度文件",
  "No file selected": "未选择文件",
  "Using entered values": "使用输入值",
  "Selected {name}": "已选择 {name}",
  "Time": "时间",
  "All times": "全部时间",
  "Specific times": "指定时间",
  "Specific hours": "指定小时",
  "Values": "数值",
  "214.5 or 12": "214.5 或 12",
  "12 or 12.5": "12 或 12.5",

  // Site and meteorology
  "Site and period": "地点和时段",
  "Meteorological driver": "气象驱动",
  "START-END (DOY)": "起止日序 (DOY)",
  "UTC offset (h)": "UTC 时差 (h)",
  "Latitude (°)": "纬度 (°)",
  "Longitude (°)": "经度 (°)",
  "Default: data/defined/meteo.txt": "默认：data/defined/meteo.txt",
  "Meteorological file": "气象文件",
  "Select meteorological file": "选择气象文件",
  "Using default meteorological file": "使用默认气象文件",
  "SW diffuse fraction": "短波散射比例",
  "-1: automatic Erbs partition; 0-1: measured diffuse fraction": "-1：按 Erbs 方法自动分离；0-1：实测散射比例",
  "0.2 = 80% direct + 20% diffuse (default); -1: automatic Erbs partition; 0-1: measured diffuse fraction": "0.2 = 80% 直射 + 20% 散射（默认）；-1：按 Erbs 方法自动分离；0-1：实测散射比例",
  "0.2 = 80% direct + 20% diffuse (default); 0-1 allowed": "0.2 = 80% 直射 + 20% 散射（默认）；可填 0-1",
  "Sky temperature (K)": "天空亮温 (K)",
  "Equivalent 8-14 um sky brightness temperature; none = black sky (emission only)": "8-14 µm 天空等效亮温；填 none 为黑天空（仅发射）",
  "Carotenoids ug cm-2 (SCOPE default 10)": "类胡萝卜素 µg cm⁻²（SCOPE 默认 10）",
  "TiRT is already running another simulation; wait for it to finish and try again.": "TiRT 正在运行另一个模拟，请等待完成后再试。",

  // Spectra
  "Spectra": "光谱",
  "Spectral information": "光谱信息",
  "Component": "组分",
  "Emissivity": "发射率",
  "Wavelength bands": "波段",
  "Wavelengths (µm)": "波长 (µm)",
  "Thermal infrared": "热红外",
  "Leaf": "叶片",
  "Soil": "土壤",
  "Roof": "屋顶",
  "Wall": "墙面",
  "Street": "街道",
  "Street vegetation": "街道植被",
  "Leaf emissivity": "叶片发射率",
  "Soil emissivity": "土壤发射率",
  "Terrain emissivity": "地形发射率",
  "Roof emissivity": "屋顶发射率",
  "Wall emissivity": "墙面发射率",
  "Street emissivity": "街道发射率",
  "Property": "参数",
  "One value per band (;)": "每个波段一个值（;）",
  "Table file": "表格文件",
  "Wavelength µm": "波段 µm",
  "Band table file": "波段参数表文件",
  "Observation band table": "观测波段参数表",
  "Table": "表格",
  "One value per band (;). \"-\" uses the PROSPECT/BSM or measured spectrum (shortwave) or the default (thermal). Thermal emissivity = 1 − ρ − τ, also used by the energy balance.": "每个波段一个值（;）。\"-\" 表示短波取 PROSPECT/BSM 或实测光谱，热红外取默认值。热红外发射率 = 1 − ρ − τ，能量平衡也使用该值。",
  "Load a 7-row band table": "读取 7 行波段参数表",
  "Band table contains a non-numeric value": "波段参数表包含非数字内容",
  "Band table needs {count} rows, found {found}": "波段参数表需要 {count} 行，实际 {found} 行",
  "Thermal emissivity ε = 1 − ρ − τ": "热红外发射率 ε = 1 − ρ − τ",
  "No thermal band (≥ 2.5 µm): all bands give reflectance BRF": "没有热红外波段（≥ 2.5 µm）：所有波段输出反射率 BRF",
  "leaf": "叶片",
  "soil": "土壤",
  "street": "街道",
  "wall": "墙面",
  "roof": "屋顶",
  "Reflectance": "反射率",
  "Leaf ρ": "叶片 ρ",
  "Leaf τ": "叶片 τ",
  "Soil ρ": "土壤 ρ",
  "Roof ρ": "屋顶 ρ",
  "Wall ρ": "墙面 ρ",
  "Street ρ": "街道 ρ",
  "Leaf reflectance": "叶片反射率",
  "Leaf transmittance": "叶片透射率",
  "Soil reflectance": "土壤反射率",
  "Roof reflectance": "屋顶反射率",
  "Wall reflectance": "墙面反射率",
  "Street reflectance": "街道反射率",
  "Wavelength file": "波段文件",
  "Select wavelength file": "选择波段文件",
  "Leaf spectrum file": "叶片光谱文件",
  "Soil spectrum file": "土壤光谱文件",
  "Terrain spectrum file": "地形光谱文件",
  "Roof spectrum file": "屋顶光谱文件",
  "Wall spectrum file": "墙面光谱文件",
  "Street spectrum file": "街道光谱文件",
  "Urban spectrum file": "城市光谱文件",
  "Select leaf spectrum file": "选择叶片光谱文件",
  "Select soil spectrum file": "选择土壤光谱文件",
  "Select terrain spectrum file": "选择地形光谱文件",
  "Select roof spectrum file": "选择屋顶光谱文件",
  "Select wall spectrum file": "选择墙面光谱文件",
  "Select street spectrum file": "选择街道光谱文件",
  "Select urban spectrum file": "选择城市光谱文件",
  "Leaf VNIR": "叶片可见光-近红外",
  "Soil VNIR": "土壤可见光-近红外",
  "Spectrum files": "光谱文件",
  "Leaf structure N": "叶片结构参数 N",
  "BSM brightness": "BSM 亮度",
  "BSM latitude": "BSM 纬度",
  "BSM longitude": "BSM 经度",
  "Soil moisture (m³ m⁻³)": "土壤含水量 (m³ m⁻³)",
  "Volumetric water content; values above 1 are read as percent": "体积含水量；大于 1 的值按百分数读取",

  // Temperatures
  "Temperatures": "温度",
  "Sunlit": "光照",
  "Shaded": "阴影",
  "Sunlit °C": "光照 °C",
  "Shaded °C": "阴影 °C",
  "Leaf sunlit temperature": "叶片光照温度",
  "Leaf shaded temperature": "叶片阴影温度",
  "Soil sunlit temperature": "土壤光照温度",
  "Soil shaded temperature": "土壤阴影温度",
  "Terrain sunlit temperature": "地形光照温度",
  "Terrain shaded temperature": "地形阴影温度",
  "Roof sunlit temperature": "屋顶光照温度",
  "Roof shaded temperature": "屋顶阴影温度",
  "Wall sunlit temperature": "墙面光照温度",
  "Wall shaded temperature": "墙面阴影温度",
  "Street sunlit temperature": "街道光照温度",
  "Street shaded temperature": "街道阴影温度",
  "Temperature parameters": "温度参数",
  "Annual temperature": "年均温度",

  // Plant physiology
  "Plant physiology": "植物生理",
  "Leaf type": "叶片类型",
  "Stomatal method": "气孔导度模型",
  "Stomatal g0": "气孔导度 g0",
  "Stress factor": "胁迫因子",
  "Photoprotection qls": "光保护 qls",
  "Vertical profile kv": "垂直廓线 kv",

  // Physical properties and options
  "Physical properties": "物理属性",
  "Soil conductivity": "土壤导热率",
  "Soil heat capacity": "土壤热容",
  "Soil density": "土壤密度",
  "Soil boundary resistance": "土壤边界层阻抗",
  "Soil surface resistance (s m⁻¹)": "土壤表面阻抗 (s m⁻¹)",
  "Urban conductivity": "城市材料导热率",
  "Urban heat capacity": "城市材料热容",
  "Urban density": "城市材料密度",
  "Urban boundary resistance": "城市边界层阻抗",
  "Urban surface resistance": "城市表面阻抗",
  "Wall thickness (m)": "墙体厚度 (m)",
  "Wall layers": "墙体层数",
  "Indoor wall boundary (C)": "室内墙面边界温度 (C)",
  "Wall diffuse factor": "墙面散射因子",
  "Urban aerodynamic transfer": "城市空气动力传输",
  "Urban roughness (m, 0=auto)": "城市粗糙度 (m，0=自动)",
  "Canyon roughness ratio": "街谷粗糙度比",
  "Minimum canyon / roof ratio": "最小街谷/屋顶比",
  "Leaf wind decay": "叶片风速衰减",
  "Leaf boundary factor": "叶片边界层因子",
  "Options": "选项",
  "Soil/street heat method": "土壤/街道热传导方法",
  "Building heat method": "建筑热传导方法",
  "Direct": "直接",
  "Depth": "深度",
  "Surface": "表面",
  "Temperature method": "温度求解方法",
  "Robust": "稳健",

  // Scene preview
  "SCENE PREVIEW": "场景预览",
  "Scene preview": "场景预览",
  "Sun": "太阳",
  "View rays": "观测方向",
  "View": "观测",
  "Vegetation": "植被",
  "Canopy": "冠层",
  "canopy": "冠层",
  "Loading Three.js scene": "正在加载 Three.js 场景",
  "Loading scene": "正在加载场景",
  "Three.js could not be loaded: {message}": "无法加载 Three.js：{message}",
  "Three.js load failed": "Three.js 加载失败",
  "North is azimuth 0 degrees": "北向为方位角 0°",
  "Solar": "太阳",
  "Scene": "场景",
  "Domain": "范围",
  "100 m × 100 m · North 0°": "100 m × 100 m · 北向 0°",
  "Domain 100 m × 100 m · North 0°": "范围 100 m × 100 m · 北向 0°",
  "Time driver · input.csv": "时间驱动 · input.csv",
  "{count} principal": "主平面 {count}",
  "hemisphere": "半球",
  "file": "文件",
  "{count} manual": "手动 {count}",
  "No directions": "无观测方向",
  "{count} directions": "{count} 个方向",
  "{zenith} zenith × {azimuth} azimuth": "{zenith} 个天顶角 × {azimuth} 个方位角",

  // Directional results
  "MODEL OUTPUT": "模型输出",
  "Directional results": "方向性结果",
  "Polar": "极坐标",
  "Parallel": "平行",
  "Vertical": "垂直",
  "Band": "波段",
  "default": "默认",
  "Run a simulation to inspect the directional field": "运行模拟以查看方向性结果",
  "Choose a geometry mode and thermal band, then run the current input through TiRT.": "选择几何模式和热红外波段，然后运行当前 TiRT 输入。",
  "Mean temperature": "平均温度",
  "Mean {quantity}": "平均{quantity}",
  "selected band": "当前波段",
  "Directional span": "方向性跨度",
  "max - min": "最大值 - 最小值",
  "Observations": "观测数",
  "VIEW HEMISPHERE": "观测半球",
  "Polar directional map": "极坐标方向图",
  "Color scale": "色标",
  "SOLAR PRINCIPAL PLANE": "太阳主平面",
  "Parallel response": "平行方向响应",
  "ORTHOGONAL PLANE": "正交平面",
  "Perpendicular response": "垂直方向响应",
  "Reflectance (BRF)": "反射率 (BRF)",
  "Brightness temperature (°C)": "亮温 (°C)",
  "VZA radius · {band} µm · {quantity}": "径向为 VZA · {band} µm · {quantity}",
  "RAA 0° / 180° · both sides": "RAA 0° / 180° · 双侧",
  "RAA 90° / 270° · both sides": "RAA 90° / 270° · 双侧",
  "No matching directions": "没有匹配的方向",

  // Static-page messages
  "Invalid value for {field}": "{field} 的数值无效",
  "Wavelengths must be numbers of at least 0.4 µm": "波段必须是不小于 0.4 µm 的数值",
  "LAI must contain numbers in 0..20, for example 1;3;5": "LAI 必须是 0 到 20 之间的数值，例如 1;3;5",
  "Hotspot must contain nonnegative numbers, for example 0.1;0.2": "热点参数必须为非负数，例如 0.1;0.2",
  "{label} must be between 0 and 1": "{label}必须在 0 到 1 之间",
  "{label}: enter {count} values, one per wavelength": "{label}：请输入 {count} 个值，每个波段一个",
  "Leaf reflectance + transmittance must not exceed 1": "叶片反射率 + 透射率不能超过 1",
  "Manual angles must be VZA/VAA pairs with VZA 0–90° and VAA 0–360°": "手动角度必须为 VZA/VAA 对，VZA 为 0–90°，VAA 为 0–360°",
  "Select an observation plane or provide observation angles": "请选择观测平面或提供观测角度",
  "One value or a semicolon-separated list, e.g. 1;3;5": "单个值或以分号分隔的列表，例如 1;3;5",
  "One value or a semicolon-separated list, e.g. 0.1;0.2": "单个值或以分号分隔的列表，例如 0.1;0.2",
  "Simulation failed": "模拟失败",
  "Simulation complete - {count} output records": "模拟完成，共 {count} 条输出记录",
  "Run a simulation first": "请先运行模拟",
  "Saved observations.csv · {count} directions · BRF + brightness temperature (°C)": "已保存 observations.csv · {count} 个方向 · BRF + 亮温 (°C)",
  "Input export failed": "输入导出失败",
  "Saved input ZIP with spectra and observation files": "已保存包含光谱和观测文件的输入 ZIP",
  "Loaded {name}": "已加载 {name}",
  "Saved {name}": "已保存 {name}",
  "TiRT has exited.": "TiRT 已退出。",
  "TiRT could not exit": "TiRT 无法退出",
  "Exit failed (HTTP {status})": "退出失败 (HTTP {status})",

  // Temporal results
  "TEMPORAL OUTPUT": "时间序列输出",
  "Time evolution": "时间变化",
  "Time output": "输出量",
  "TIME": "时间",
  "Energy balance": "能量平衡",
  "Brightness temperature": "亮温",
  "Air temperature": "气温",
  "Sunlit soil temperature": "光照土壤温度",
  "Shaded soil temperature": "阴影土壤温度",
  "Sunlit leaf temperature": "光照叶片温度",
  "Shaded leaf temperature": "阴影叶片温度",
  "Sunlit roof temperature": "光照屋顶温度",
  "Shaded roof temperature": "阴影屋顶温度",
  "Sunlit wall temperature": "光照墙面温度",
  "Shaded wall temperature": "阴影墙面温度",
  "Sunlit street temperature": "光照街道温度",
  "Shaded street temperature": "阴影街道温度",
  "Latent heat flux": "潜热通量",
  "Sensible heat flux": "感热通量",
  "Ground heat flux": "土壤热通量",
  "Brightness temperature over time": "亮温随时间变化",
  "{quantity} over time": "{quantity}随时间变化",
  "No time series": "暂无时间序列",
  "No time series data": "暂无时间序列数据",
  "Run a time series to inspect the temporal response": "运行时间序列以查看时间响应",
  "Use a meteorological file or the default driver.": "使用气象文件或默认驱动。",
  "Radiance": "辐亮度",
  "Mean": "平均",
  "Minimum": "最小值",
  "Maximum": "最大值",
  "Time steps": "时间步数",
  "Selected quantity": "当前量",
  "time series": "时间序列",
  "Selected interval": "选定区间",
  "selected interval": "选定区间",
  "No finite values": "没有有效数值",
  "URBAN COMPONENTS": "城市组分",
  "Resistance and energy diagnostics": "阻抗与能量诊断",
  "LATEST STEP": "最新时间步",
  "State": "状态",
  "SW direct": "短波直射",
  "SW diffuse": "短波散射",
  "LW input": "长波输入",
  "LW out": "长波输出",
  "{steps} time steps · {records} observation records": "{steps} 个时间步 · {records} 条观测记录",
  "Energy balance: {closed}/{total} closed · max residual {residual} W/m²": "能量平衡：{closed}/{total} 闭合 · 最大残差 {residual} W/m²",
  "Invalid server response (HTTP {status})": "服务器响应无效 (HTTP {status})",
  "Time series failed (HTTP {status})": "时间序列运行失败 (HTTP {status})",
  "Simulation completed; display warning: {details}": "模拟已完成；显示警告：{details}",
  "Time series complete · {count} unclosed time steps; see energy balance residuals": "时间序列完成 · {count} 个时间步未闭合，请查看能量平衡残差",
  "Time series complete · {count} time steps": "时间序列完成 · {count} 个时间步",
});

function translated(value) {
  return language === "zh-CN" ? (ZH[value] || value) : value;
}

function t(key, params = {}) {
  return translated(key).replace(/\{(\w+)\}/g, (match, name) => (name in params ? String(params[name]) : match));
}

window.tirtI18n = Object.freeze({ language, t });

function splitText(value) {
  const match = String(value).match(/^(\s*)(.*?)(\s*)$/s);
  return match ? `${match[1]}${translated(match[2])}${match[3]}` : value;
}

const sourceTitle = document.title;

function applyLanguage() {
  const root = document.body;
  if (!root) return;
  document.title = translated(sourceTitle);
  const walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);
  const nodes = [];
  while (walker.nextNode()) nodes.push(walker.currentNode);
  for (const node of nodes) {
    if (node.parentElement?.closest("script,style,[data-i18n-ignore]")) continue;
    const current = node.nodeValue || "";
    const previous = sourceText.get(node);
    const previousLocalized = previous == null ? "" : splitText(previous);
    const source = previous == null || (current !== previous && current !== previousLocalized) ? current : previous;
    sourceText.set(node, source);
    node.nodeValue = splitText(source);
  }
  const attributes = ["title", "aria-label", "placeholder"];
  for (const element of root.querySelectorAll("*[title], *[aria-label], *[placeholder]")) {
    if (element.closest("[data-i18n-ignore]")) continue;
    let originals = sourceAttributes.get(element);
    if (!originals) {
      originals = {};
      sourceAttributes.set(element, originals);
    }
    for (const name of attributes) {
      if (!element.hasAttribute(name)) continue;
      const current = element.getAttribute(name) || "";
      const previous = originals[name];
      const previousLocalized = previous == null ? "" : translated(previous);
      const source = previous == null || (current !== previous && current !== previousLocalized) ? current : previous;
      originals[name] = source;
      element.setAttribute(name, translated(source));
    }
  }
  document.documentElement.lang = language;
  const button = document.querySelector("#language-toggle");
  if (button) {
    button.querySelector("[data-language-label]").textContent = language === "zh-CN" ? "EN" : "中文";
    button.title = t(language === "zh-CN" ? "Switch to English" : "Switch to Chinese");
    button.setAttribute("aria-label", button.title);
  }
}

function initializeLanguage() {
  document.querySelector("#language-toggle")?.addEventListener("click", () => {
    localStorage.setItem(STORAGE_KEY, language === "zh-CN" ? "en-US" : "zh-CN");
    window.location.reload();
  });
  applyLanguage();
  const observer = new MutationObserver(() => {
    observer.disconnect();
    applyLanguage();
    observer.observe(document.body, { childList: true, characterData: true, subtree: true, attributes: true, attributeFilter: ["title", "aria-label", "placeholder"] });
  });
  observer.observe(document.body, { childList: true, characterData: true, subtree: true, attributes: true, attributeFilter: ["title", "aria-label", "placeholder"] });
}

if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", initializeLanguage, { once: true });
else initializeLanguage();
