# Firefox 兼容测试

## 结论

TiRT GUI 可以在 Firefox 中使用。

## 测试环境

- 浏览器：`C:\Program Files\Mozilla Firefox\firefox.exe`
- 测试方式：Selenium 控制本机 Firefox headless 模式
- 测试页面：`http://127.0.0.1:8785/`

## 测试结果

- 首页加载：通过
- `Run simulation` 点击运行：通过
- 结果区显示：通过
- 状态显示 `Complete`：通过
- 观测数量显示：通过
- 截图输出：`docs/firefox_results.png`

## 说明

TiRT GUI 使用标准 HTML、CSS、JavaScript module、Fetch API 和 WebGL/Three.js。Firefox、Chrome、Edge 均可使用。若双击 exe 后浏览器未自动打开，可手动访问：

```text
http://127.0.0.1:8765/
```
