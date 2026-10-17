const STORAGE_KEY = "tirt-language";
const language = localStorage.getItem(STORAGE_KEY) === "zh-CN" ? "zh-CN" : "en-US";
const sourceText = new WeakMap();
const sourceAttributes = new WeakMap();

const ZH = Object.freeze({
  "Static scene": "静态场景",
  "Dynamic scene": "动态场景",
  "THERMAL INFRARED RADIATIVE TRANSFER": "热红外辐射传输",
  "THERMAL INFRARED ENERGY BALANCE": "热红外能量平衡",
  "Directional Simulator": "方向性模拟器",
  "Temporal Simulator": "时间序列模拟器",
  "Energy balance": "能量平衡",
  "Directional model": "方向性模型",
  "Ready": "就绪",
  "Running...": "运行中...",
  "Run simulation": "运行模拟",
  "Run time series": "运行时间序列",
  "Save observations.csv": "保存 observations.csv",
  "Save processes.csv": "保存 processes.csv",
  "Save input.csv": "保存 input.csv",
  "Exit": "退出",
  "Input settings": "输入设置",
  "Scene structure": "场景结构",
  "Site and period": "地点和时段",
  "Meteorological driver": "气象驱动",
  "Observation": "观测",
  "Angles": "角度",
  "Spectra": "光谱",
  "Temperatures": "温度",
  "Scene preview": "场景预览",
  "Directional results": "方向性结果",
  "Time evolution": "时间变化",
  "Brightness temperature": "亮温",
  "No time series": "暂无时间序列",
  "No time series data": "暂无时间序列数据",
  "Run a simulation to inspect the directional field": "运行模拟以查看方向性结果",
  "Run a time series to inspect the temporal response": "运行时间序列以查看时间响应",
  "Choose a geometry mode and thermal band, then run the current input through TiRT.": "选择几何模式和热红外波段，然后运行当前 TiRT 输入。",
  "Use a meteorological file or the default driver.": "使用气象文件或默认驱动。",
  "Control deck": "控制面板",
  "Reset settings": "重置设置",
  "Solar": "太阳",
  "Observation": "观测",
  "Scene": "场景",
  "Domain": "范围",
  "Mean temperature": "平均温度",
  "Directional span": "方向性跨度",
  "Radiance": "辐亮度",
  "Observations": "观测记录",
  "Mean": "平均",
  "Minimum": "最小值",
  "Maximum": "最大值",
  "Time steps": "时间步数",
  "Selected quantity": "当前量",
  "time series": "时间序列",
  "Selected interval": "选定区间",
  "View hemisphere": "观测半球",
  "Polar directional map": "极坐标方向图",
  "Parallel response": "平行方向响应",
  "Perpendicular response": "垂直方向响应",
  "Run a time series to inspect the temporal response": "运行时间序列以查看时间响应",
  "No finite values": "没有有效数值",
  "Switch to Chinese": "切换为中文",
  "Switch to English": "切换为英文",
});

function translated(value) {
  return language === "zh-CN" ? (ZH[value] || value) : value;
}

function splitText(value) {
  const match = String(value).match(/^(\s*)(.*?)(\s*)$/s);
  return match ? `${match[1]}${translated(match[2])}${match[3]}` : value;
}

function applyLanguage() {
  const root = document.body;
  if (!root) return;
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
    button.title = language === "zh-CN" ? "Switch to English" : "Switch to Chinese";
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
    observer.observe(document.body, { childList: true, characterData: true, subtree: true });
  });
  observer.observe(document.body, { childList: true, characterData: true, subtree: true });
}

if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", initializeLanguage, { once: true });
else initializeLanguage();
