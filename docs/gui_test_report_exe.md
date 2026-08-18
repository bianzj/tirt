# TiRT GUI 测试报告

- 测试模式：`exe`
- 测试地址：`http://127.0.0.1:8782`
- 通过：`21/21`

| 测试项 | 结果 | 说明 |
|---|---:|---|
| health endpoint | PASS | `{'ok': True, 'project': 'C:\\Users\\jiank\\AppData\\Local\\Temp\\_MEI0000df882', 'engine': 'tirt'}` |
| index page | PASS | `HTTP 200` |
| static asset app.js | PASS | `HTTP 200, bytes~46732` |
| static asset styles.css | PASS | `HTTP 200, bytes~17026` |
| static asset vendor/three.module.js | PASS | `HTTP 200, bytes~1325994` |
| static asset vendor/controls/OrbitControls.js | PASS | `HTTP 200, bytes~31285` |
| default run | PASS | `HTTP 200, count=84` |
| output fields | PASS | `fields=['brightness_temperature_C', 'brightness_temperature_K', 'emissivity_leaf', 'emissivity_roof', 'emissivity_soil', 'emissivity_street', 'emissivity_wall', 'index', 'model', 'raa', 'radiance', 'saa', 'scenario_id', 'scenario_label', 'surface_model', 'sza', 'vaa', 'vegetation_model', 'vza', 'wavelength_um']` |
| model plane+bare | PASS | `HTTP 200, count=12` |
| model plane+hom | PASS | `HTTP 200, count=12` |
| model plane+row | PASS | `HTTP 200, count=12` |
| model plane+crown | PASS | `HTTP 200, count=12` |
| model slope+hom | PASS | `HTTP 200, count=12` |
| model terrain+bare | PASS | `HTTP 200, count=12` |
| model terrain+crown | PASS | `HTTP 200, count=12` |
| model urban+bare | PASS | `HTTP 200, count=12` |
| model urban+crown | PASS | `HTTP 200, count=12` |
| manual geometry | PASS | `HTTP 200, count=3` |
| uploaded geometry text | PASS | `HTTP 200, count=3` |
| uploaded spectral files | PASS | `HTTP 200, count=3` |
| invalid input error | PASS | `HTTP 400, error=surface_model 必须是 plane、slope、terrain 或 urban` |
