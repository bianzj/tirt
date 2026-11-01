const $ = (id) => document.getElementById(id);
const t = (key, params = {}) => window.tirtI18n
  ? window.tirtI18n.t(key, params)
  : key.replace(/\{(\w+)\}/g, (match, name) => (name in params ? String(params[name]) : match));
if (window.location.protocol === "file:") window.location.replace("http://127.0.0.1:8765/time.html");
let meteoText = "";
const spectrumFiles = { directions: "", vegetation_spectrum: "", soil_spectrum: "", urban_spectrum: "" };
const BAND_TABLE_IDS = ["wavelengths", "leaf-reflectance", "leaf-transmittance", "soil-reflectance",
  "street-reflectance", "wall-reflectance", "roof-reflectance"];
// Thermal values used for "-" cells; must match TIR_DEFAULTS in base/io.py.
const TIR_DEFAULTS = { "leaf-reflectance": 0.015, "leaf-transmittance": 0, "soil-reflectance": 0.045,
  "street-reflectance": 0.045, "wall-reflectance": 0.07, "roof-reflectance": 0.05 };
let seriesRows = [];
let observationRows = [];
let diagnosticRows = [];
let toastTimer;
window.temporalState = { get seriesRows() { return seriesRows; }, get observationRows() { return observationRows; } };

const COLORS = { cyan: "#d53b36", amber: "#d38a24", coral: "#d53b36", soft: "#35566b", muted: "#5d7383", grid: "#c9dae4" };
const quantityLabels = {
  Tv_C: "Brightness temperature", Tss_C: "Sunlit soil temperature", Tsh_C: "Shaded soil temperature",
  Tls_C: "Sunlit leaf temperature", Tlh_C: "Shaded leaf temperature", Ta_C: "Air temperature",
  Troof_sunlit_C: "Sunlit roof temperature", Troof_shaded_C: "Shaded roof temperature",
  Twall_sunlit_C: "Sunlit wall temperature", Twall_shaded_C: "Shaded wall temperature",
  Tstreet_sunlit_C: "Sunlit street temperature", Tstreet_shaded_C: "Shaded street temperature",
  Tleaf_sunlit_C: "Sunlit leaf temperature", Tleaf_shaded_C: "Shaded leaf temperature",
  LE: "Latent heat flux", H: "Sensible heat flux", G: "Ground heat flux",
};

const commonQuantities = ["Tv_C", "Ta_C"];
const fluxQuantities = ["LE", "H", "G"];
const soilQuantities = ["Tss_C", "Tsh_C"];
const leafQuantities = ["Tls_C", "Tlh_C"];
const urbanQuantities = [
  "Troof_sunlit_C", "Troof_shaded_C", "Twall_sunlit_C", "Twall_shaded_C",
  "Tstreet_sunlit_C", "Tstreet_shaded_C",
];
const urbanLeafQuantities = ["Tleaf_sunlit_C", "Tleaf_shaded_C"];

function diagnosticValue(row, key) {
  if (row?.[key] == null || row[key] === "") return null;
  const value = Number(row?.[key]);
  return Number.isFinite(value) ? value : null;
}

function diagnosticText(value, digits = 2) {
  return value === null ? "--" : value.toFixed(digits);
}

function updateDiagnostics() {
  const card = $("diagnostics-card");
  const urban = value("surface-model") === "urban";
  const hasVegetation = value("vegetation-model") !== "bare";
  const row = diagnosticRows[diagnosticRows.length - 1];
  const visible = urban && Boolean(row);
  card.classList.toggle("is-hidden", !visible);
  if (!visible) return;

  const components = [
    ["roof", "Roof"], ["wall", "Wall"], ["street", "Street"],
    ...(hasVegetation ? [["leaf", "Street vegetation"]] : []),
  ];
  const states = [["sunlit", "Sunlit"], ["shaded", "Shaded"]];
  const body = $("diagnostics-body");
  body.replaceChildren(...components.flatMap(([component, label]) => states.map(([state, stateLabel]) => {
    const suffix = `${component}_${state}`;
    const fluxSuffix = `${component}${state}`;
    const tr = document.createElement("tr");
    const values = [
      t(label), t(stateLabel),
      diagnosticText(diagnosticValue(row, `SW_direct_${suffix}`)),
      diagnosticText(diagnosticValue(row, `SW_diffuse_${suffix}`)),
      diagnosticText(diagnosticValue(row, `LW_input_${suffix}`)),
      diagnosticText(diagnosticValue(row, `LW_out_${fluxSuffix}`)),
      diagnosticText(diagnosticValue(row, `H_${fluxSuffix}`)),
      diagnosticText(diagnosticValue(row, `LE_${fluxSuffix}`)),
      diagnosticText(diagnosticValue(row, `G_${fluxSuffix}`)),
      diagnosticText(diagnosticValue(row, `raa_${suffix}`)),
      diagnosticText(diagnosticValue(row, `rss_${suffix}`)),
    ];
    tr.replaceChildren(...values.map((text, index) => {
      const cell = document.createElement(index < 2 ? "th" : "td");
      cell.textContent = text;
      if (index === 0) cell.scope = "rowgroup";
      return cell;
    }));
    return tr;
  })));
  const meta = [
    `DOY ${diagnosticText(diagnosticValue(row, "time"), 4)}`,
    `z ${diagnosticText(diagnosticValue(row, "z"), 2)} m`,
    `ustar ${diagnosticText(diagnosticValue(row, "ustar"), 4)} m/s`,
    `u canyon ${diagnosticText(diagnosticValue(row, "u_canopy"), 3)} m/s`,
    `u leaf ${diagnosticText(diagnosticValue(row, "u_leaf"), 3)} m/s`,
  ];
  $("diagnostics-meta").textContent = meta.join(" · ");
}

const DISPLAY_QUANTITY = "Tv_C";

function setStatus(text, kind = "idle") { const pill = $("status-pill"); pill.className = `status-pill ${kind}`; pill.querySelector("span").textContent = t(text); if (kind !== "error") pill.title = ""; }
function notify(message, error = false) { const toast = $("toast"); toast.textContent = message; toast.classList.toggle("error", error); toast.classList.add("show"); clearTimeout(toastTimer); toastTimer = setTimeout(() => toast.classList.remove("show"), 4200); }
function value(id) { return $(id).value; }
function anglePairs(raw) {
  return String(raw || "").split(/[;\n]+/).map((item) => item.trim()).filter(Boolean).map((item) => item.split(/\s*\/\s*|\s*,\s*|\s+/).map(Number)).filter((pair) => pair.length >= 2 && Number.isFinite(pair[0]) && Number.isFinite(pair[1]));
}
function timeRange(raw) {
  const values = String(raw || "").split(/\s*[-/]\s*/).map(Number).filter(Number.isFinite);
  if (values.length >= 2) return values.slice(0, 2);
  if (values.length === 1) return [values[0], values[0]];
  return [214, 215];
}
function fmt(value, digits = 3) { return Number.isFinite(Number(value)) ? Number(value).toFixed(digits) : "--"; }

// Same layout as base/inputio.py read_band_cells: optional text label, then ";"-separated cells; "-" or empty = model value.
function bandCells(text) {
  const raw = String(text ?? "").trim();
  if (!raw) return [];
  const parts = /[;,，；]/.test(raw) ? raw.split(/[;,，；]/) : raw.split(/\s+/);
  return parts.map((part) => {
    const item = part.trim();
    if (item === "" || item === "-") return null;
    const number = Number(item);
    if (!Number.isFinite(number)) throw new Error(t("Band table contains a non-numeric value"));
    return number;
  });
}

function parseBandTable(text) {
  const isLabel = (token) => token !== "-" && !/[;,，；]/.test(token) && !Number.isFinite(Number(token));
  const rows = [];
  for (const line of String(text).replace(/^\uFEFF/, "").split(/\r?\n/)) {
    let body = line.split("#")[0].trim();
    const colon = body.search(/[:：]/);
    if (colon >= 0) body = body.slice(colon + 1);
    const tokens = body.split(/\s+/).filter(Boolean);
    while (tokens.length && isLabel(tokens[0])) tokens.shift();
    if (tokens.length) rows.push(bandCells(tokens.join(" ")));
  }
  if (rows.length !== BAND_TABLE_IDS.length) {
    throw new Error(t("Band table needs {count} rows, found {found}", { count: BAND_TABLE_IDS.length, found: rows.length }));
  }
  return rows;
}

function updateEmissivityReadout() {
  const readout = $("thermal-emissivity-readout");
  let rows;
  try { rows = Object.fromEntries(BAND_TABLE_IDS.map((id) => [id, bandCells(value(id))])); }
  catch (error) { readout.textContent = error.message; return; }
  const cell = (id, index) => rows[id][index] ?? TIR_DEFAULTS[id];
  const parts = rows.wavelengths.map((wave, index) => {
    if (!(wave >= 2.5)) return null;
    const leaf = 1 - cell("leaf-reflectance", index) - cell("leaf-transmittance", index);
    const others = ["soil", "street", "wall", "roof"].map((name) => `${t(name)} ${fmt(1 - cell(`${name}-reflectance`, index), 3)}`);
    return `${wave} µm: ${t("leaf")} ${fmt(leaf, 3)} · ${others.join(" · ")}`;
  }).filter(Boolean);
  readout.textContent = parts.length
    ? `${t("Thermal emissivity ε = 1 − ρ − τ")} — ${parts.join("; ")}`
    : t("No thermal band (≥ 2.5 µm): all bands give reflectance BRF");
}

function currentOverrides() {
  const angles = anglePairs(value("observation-angles"));
  const observationVza = angles.map(([vza]) => vza).join(";");
  const observationVaa = angles.map(([, vaa]) => vaa).join(";");
  const observationTimeMode = value("observation-time-mode");
  const observationTimeValues = value("observation-time-values");
  const directionsFromFile = Boolean(spectrumFiles.directions);
  const [startTime, endTime] = timeRange(value("time-range"));
  return {
    surface_model: value("surface-model"), vegetation_model: value("vegetation-model"),
    lai: value("lai"), canopy_height: value("canopy-height"), hspot: value("hspot"), leaf_width: value("leaf-width"),
    row_width: value("row-width"), row_blank: value("row-blank"), row_height: value("row-height"), row_azimuth: value("row-azimuth"),
    crown_rad_a: value("crown-rad-a"), crown_rad_b: value("crown-rad-b"), crown_height: value("crown-height"), crown_density: value("crown-density"),
    slope_angle: value("slope-angle"), slope_aspect: value("slope-aspect"),
    terrain_height: value("terrain-height"), terrain_radius: value("terrain-radius"), terrain_density: value("terrain-density"),
    building_density: value("building-density"), building_length: value("building-length"), building_width: value("building-width"), building_height: value("building-height"), building_azimuth: value("building-azimuth"),
    urban_roughness_length: value("urban-roughness-length"), urban_canyon_roughness_ratio: value("urban-canyon-roughness-ratio"), urban_canyon_min_ratio: value("urban-canyon-min-ratio"), urban_leaf_wind_decay: value("urban-leaf-wind-decay"), urban_leaf_boundary_factor: value("urban-leaf-boundary-factor"),
    latitude: value("latitude"), longitude: value("longitude"), utc_offset_hours: value("utc-offset"), start_time: startTime, end_time: endTime,
    geometry_mode: directionsFromFile ? "2" : "0", geometry_pair_mode: "zip", direction_source: directionsFromFile ? "file" : "direct", vza: observationVza, vaa: observationVaa,
    spectral_source: "direct", wavelengths: value("wavelengths"),
    leaf_reflectance: value("leaf-reflectance"), leaf_transmittance: value("leaf-transmittance"),
    soil_reflectance: value("soil-reflectance"), street_reflectance: value("street-reflectance"),
    wall_reflectance: value("wall-reflectance"), roof_reflectance: value("roof-reflectance"),
    cab: value("cab"), cw: value("cw"), cdm: value("cdm"), cs: value("cs"), n: value("leaf-n"), cca: value("cca"),
    bsm_brightness: value("soil-brightness"), bsm_latitude: value("soil-latitude"), bsm_longitude: value("soil-longitude"),
    vcmax: value("vcmax"), ballberry_m: value("ballberry-m"), leaf_type: value("leaf-type"), stomatal_method: value("stomatal-method"), stomatal_g0: value("stomatal-g0"),
    rdparam: value("rdparam"), tyear: value("tyear"), beta: value("beta"), knpqs: value("knpqs"), qls: value("qls"), kv: value("kv"), stressfactor: value("stressfactor"),
    tparam1: value("tparam1"), tparam2: value("tparam2"), tparam3: value("tparam3"), tparam4: value("tparam4"), tparam5: value("tparam5"),
    shortwave_diffuse_fraction: value("shortwave-diffuse-fraction"), soil_moisture: value("soil-moisture"), soil_rss: value("soil-rss"), soil_cs: value("soil-cs"), soil_rhos: value("soil-rhos"), soil_lambda: value("soil-lambda"), soil_rbs: value("soil-rbs"),
    urban_rss: value("urban-rss"), urban_cs: value("urban-cs"), urban_rhos: value("urban-rhos"), urban_lambda: value("urban-lambda"), urban_rbs: value("urban-rbs"), urban_wall_thickness: value("urban-wall-thickness"), urban_wall_layers: value("urban-wall-layers"), urban_internal_temperature_c: value("urban-internal-temperature"), urban_wall_diffuse_factor: value("urban-wall-diffuse-factor"),
    soil_heat_method: value("soil-heat-method"), urban_heat_method: value("urban-heat-method"), temperature_method: value("temperature-method"),
    observation_time_mode: observationTimeMode, observation_times: observationTimeMode === "time" ? observationTimeValues : "", observation_hours: observationTimeMode === "hour" ? observationTimeValues : "", observation_quantity: "temperature", output_intermediate: "yes",
  };
}

function updateConditionalFields() {
  const surface = value("surface-model");
  const vegetation = value("vegetation-model");
  const crowns = vegetation === "crown";
  $("slope-fields").classList.toggle("is-hidden", surface !== "slope");
  $("terrain-fields").classList.toggle("is-hidden", surface !== "terrain");
  $("urban-fields").classList.toggle("is-hidden", surface !== "urban");
  $("row-fields").classList.toggle("is-hidden", vegetation !== "row" || crowns);
  $("crown-fields").classList.toggle("is-hidden", !crowns);
}

function updateObservationTimeFields() {
  const mode = value("observation-time-mode");
  const input = $("observation-time-values");
  input.disabled = mode === "all";
  input.placeholder = t(mode === "hour" ? "12 or 12.5" : "214.5");
}

function timeValue(row) { return Number(row.time); }
function finiteRows(rows, key) { return rows.filter(row => row[key] != null && row[key] !== "").map((row) => ({ x: timeValue(row), y: Number(row[key]) })).filter((point) => Number.isFinite(point.x) && Number.isFinite(point.y)); }
function xLabel(x, first) { return `${((x - first) * 24).toFixed(1)} h`; }

function canvasContext(canvas) {
  const rect = canvas.getBoundingClientRect(); const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const width = Math.max(1, Math.floor(rect.width * dpr)); const height = Math.max(1, Math.floor(rect.height * dpr));
  canvas.width = width; canvas.height = height; const ctx = canvas.getContext("2d"); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, width: rect.width, height: rect.height };
}

function rangeFor(series) {
  const values = series.flatMap((item) => item.points.map((point) => point.y)).filter(Number.isFinite);
  if (!values.length) return [0, 1];
  let lo = Math.min(...values); let hi = Math.max(...values); if (Math.abs(hi - lo) < 1e-9) { lo -= 1; hi += 1; }
  const pad = (hi - lo) * .12; return [lo - pad, hi + pad];
}

function drawChart(canvasId, series, yUnit = "") {
  const { ctx, width, height } = canvasContext($(canvasId)); const pad = { left: 48, right: 16, top: 18, bottom: 29 };
  ctx.clearRect(0, 0, width, height); ctx.font = "10px system-ui, sans-serif"; ctx.lineWidth = 1; ctx.strokeStyle = COLORS.grid; ctx.fillStyle = COLORS.muted;
  const points = series.flatMap((item) => item.points); if (!points.length) { ctx.textAlign = "center"; ctx.fillText(t("No finite values"), width / 2, height / 2); return; }
  const xs = points.map((point) => point.x); const first = Math.min(...xs); let xlo = first; let xhi = Math.max(...xs); if (xhi === xlo) xhi = xlo + 1 / 24;
  const [ylo, yhi] = rangeFor(series); const plotW = width - pad.left - pad.right; const plotH = height - pad.top - pad.bottom;
  const px = (x) => pad.left + (x - xlo) / (xhi - xlo) * plotW; const py = (y) => pad.top + (yhi - y) / (yhi - ylo) * plotH;
  for (let i = 0; i <= 4; i++) { const y = pad.top + plotH * i / 4; ctx.beginPath(); ctx.moveTo(pad.left, y); ctx.lineTo(pad.left + plotW, y); ctx.stroke(); const val = yhi - (yhi - ylo) * i / 4; ctx.textAlign = "right"; ctx.fillText(fmt(val, 3), pad.left - 7, y + 3); }
  [0, .5, 1].forEach((fraction) => { const x = pad.left + plotW * fraction; ctx.beginPath(); ctx.moveTo(x, pad.top); ctx.lineTo(x, pad.top + plotH); ctx.stroke(); ctx.textAlign = fraction === 0 ? "left" : fraction === 1 ? "right" : "center"; ctx.fillText(xLabel(xlo + (xhi - xlo) * fraction, first), x, height - 9); });
  series.forEach((item) => { const valid = item.points.filter((point) => Number.isFinite(point.y)); if (!valid.length) return; ctx.strokeStyle = item.color; ctx.lineWidth = 2; ctx.beginPath(); valid.forEach((point, index) => { if (index) ctx.lineTo(px(point.x), py(point.y)); else ctx.moveTo(px(point.x), py(point.y)); }); ctx.stroke(); ctx.fillStyle = item.color; valid.forEach((point) => { ctx.beginPath(); ctx.arc(px(point.x), py(point.y), 2.4, 0, Math.PI * 2); ctx.fill(); }); });
  ctx.fillStyle = COLORS.soft; ctx.textAlign = "left"; ctx.fillText(yUnit, pad.left, 11);
  let legendX = Math.max(pad.left + 80, width - series.length * 118); series.forEach((item) => { ctx.fillStyle = item.color; ctx.fillRect(legendX, 7, 14, 3); ctx.fillStyle = COLORS.soft; ctx.fillText(item.label, legendX + 19, 11); legendX += 118; });
}

function drawTimeCharts() {
  const quantity = DISPLAY_QUANTITY; const current = [{ label: t(quantityLabels[quantity] || quantity), color: COLORS.cyan, points: finiteRows(seriesRows, quantity) }];
  drawChart("time-chart", current, "°C");
}

function updateMetrics() {
  const quantity = DISPLAY_QUANTITY; const points = finiteRows(seriesRows, quantity); const label = t(quantityLabels[quantity] || quantity); $("metric-label").textContent = label; $("chart-title").textContent = t("{quantity} over time", { quantity: label }); $("metric-count").textContent = seriesRows.length || "--"; $("series-count").textContent = seriesRows.length ? t("{steps} time steps · {records} observation records", { steps: seriesRows.length, records: observationRows.length }) : t("No time series");
  if (!points.length) { $("metric-min").textContent = "--"; $("metric-max").textContent = "--"; $("metric-min-time").textContent = "--"; $("metric-max-time").textContent = "--"; return; }
  const min = points.reduce((a, b) => b.y < a.y ? b : a); const max = points.reduce((a, b) => b.y > a.y ? b : a);
  $("metric-min").textContent = `${fmt(min.y, 3)} °C`; $("metric-max").textContent = `${fmt(max.y, 3)} °C`; $("metric-min-time").textContent = `DOY ${fmt(min.x, 4)}`; $("metric-max-time").textContent = `DOY ${fmt(max.x, 4)}`;
}

function csvCell(value) { const text = String(value ?? ""); return /[",\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text; }
function outputValue(value) { if (value == null || value === "") return ""; const number = Number(value); return Number.isFinite(number) ? number.toFixed(3) : value; }
function saveRowsCsv(rows, filename) { if (!rows.length) return; const fields = [...new Set(rows.flatMap((row) => Object.keys(row)))]; const body = [fields.map(csvCell).join(","), ...rows.map((row) => fields.map((field) => csvCell(outputValue(row[field]))).join(","))].join("\n"); const link = document.createElement("a"); link.href = URL.createObjectURL(new Blob([body], { type: "text/csv;charset=utf-8" })); link.download = filename; link.click(); URL.revokeObjectURL(link.href); notify(t("Saved {name}", { name: filename })); }
function saveCsv() { saveRowsCsv(seriesRows, "processes.csv"); }
function saveText(filename, text) { const link = document.createElement("a"); link.href = URL.createObjectURL(new Blob([text], { type: "text/plain;charset=utf-8" })); link.download = filename; link.click(); URL.revokeObjectURL(link.href); notify(t("Saved {name}", { name: filename })); }
function saveObservations() { saveRowsCsv(observationRows, "observations.csv"); }
function saveInput() {
  const config = currentOverrides();
  const lines = ["# TiRT energy-balance GUI input", ",,,", "[gui],,,"];
  Object.entries(config).forEach(([key, value]) => { if (value !== "" && value !== undefined && value !== null) lines.push(["gui", key, value, ""].map(csvCell).join(",")); });
  saveText("input.csv", `${lines.join("\n")}\n`);
}

async function runSimulation() {
  setStatus("Running", "busy"); $("run-button").disabled = true;
  try {
    const response = await fetch("/api/run_time", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ overrides: currentOverrides(), meteo_text: meteoText, files: spectrumFiles }) });
    const responseText = await response.text();
    let payload;
    try { payload = JSON.parse(responseText); }
    catch (_) { throw new Error(t("Invalid server response (HTTP {status})", { status: response.status })); }
    if (!response.ok || !payload.ok) throw new Error(payload.error || t("Time series failed (HTTP {status})", { status: response.status }));

    seriesRows = payload.series || [];
    observationRows = payload.observations || [];
    diagnosticRows = payload.intermediate || [];
    const unclosed = seriesRows.filter(row => row.closed === false || String(row.closed).toLowerCase() === "false");
    const maxResidual = Math.max(0, ...seriesRows.map(row => Number(row.energy_balance_residual_W_m2)).filter(Number.isFinite));
    $("convergence-summary").textContent = t("Energy balance: {closed}/{total} closed · max residual {residual} W/m²", { closed: seriesRows.length - unclosed.length, total: seriesRows.length, residual: maxResidual.toFixed(3) });
    $("save-button").disabled = !seriesRows.length;
    $("save-observations-button").disabled = !observationRows.length;
    $("empty-state").classList.toggle("is-hidden", Boolean(seriesRows.length));

    const displayErrors = [];
    [["metrics", updateMetrics], ["chart", drawTimeCharts], ["diagnostics", updateDiagnostics], ["scene", () => window.refreshTemporalScene?.()]].forEach(([name, render]) => {
      try { render(); }
      catch (error) { console.error(`TiRT ${name} refresh failed`, error); displayErrors.push(`${name}: ${error.message}`); }
    });
    setStatus("Complete", "idle");
    if (displayErrors.length) notify(t("Simulation completed; display warning: {details}", { details: displayErrors.join("; ") }), true);
    else if (unclosed.length) notify(t("Time series complete · {count} unclosed time steps; see energy balance residuals", { count: unclosed.length }), true);
    else notify(t("Time series complete · {count} time steps", { count: payload.count }));
  } catch (error) {
    console.error("TiRT energy-balance run failed", error);
    setStatus("Error", "error");
    $("status-pill").title = error.message;
    notify(error.message, true);
  }
  finally { $("run-button").disabled = false; }
}

async function shutdownTiRT() {
  const button = document.getElementById("exit-button");
  button.disabled = true;
  setStatus("Exiting", "busy");
  try {
    const response = await fetch("/api/shutdown", { method: "POST" });
    if (!response.ok) throw new Error(t("Exit failed (HTTP {status})", { status: response.status }));
    document.body.textContent = t("TiRT has exited.");
  } catch (error) {
    button.disabled = false;
    setStatus("Error", "error");
    notify(error.message || t("TiRT could not exit"), true);
  }
}

function bindEvents() {
  document.getElementById("exit-button").addEventListener("click", shutdownTiRT);
  $("run-button").addEventListener("click", runSimulation); $("save-button").addEventListener("click", saveCsv); $("save-observations-button").addEventListener("click", saveObservations); $("save-input-button").addEventListener("click", saveInput); $("reset-button").addEventListener("click", () => window.location.reload()); $("surface-model").addEventListener("change", () => { updateConditionalFields(); updateDiagnostics(); window.refreshTemporalScene?.(); }); $("vegetation-model").addEventListener("change", () => { updateConditionalFields(); updateDiagnostics(); window.refreshTemporalScene?.(); });
  $("meteo-upload").addEventListener("change", async (event) => { const file = event.target.files[0]; if (!file) return; meteoText = await file.text(); $("meteo-status").textContent = "✓"; $("meteo-status").className = "file-status has-file"; $("meteo-status").setAttribute("aria-label", t("Selected {name}", { name: file.name })); notify(t("Loaded {name}", { name: file.name })); });
  [["angle-upload", "directions", "angle-status"], ["vegetation-spectrum-upload", "vegetation_spectrum", "vegetation-spectrum-status"], ["soil-spectrum-upload", "soil_spectrum", "soil-spectrum-status"], ["urban-spectrum-upload", "urban_spectrum", "urban-spectrum-status"]].forEach(([id, key, statusId]) => {
    $(id).addEventListener("change", async (event) => {
      const file = event.target.files[0]; if (!file) return;
      spectrumFiles[key] = await file.text(); const status = $(statusId); status.textContent = "✓"; status.className = "file-status has-file"; status.setAttribute("aria-label", t("Selected {name}", { name: file.name }));
      if (key === "directions") { const pair = spectrumFiles[key].split(/\r?\n/).map((line) => line.trim().split(/[;,\s]+/).map(Number)).find((items) => items.length >= 2 && items.slice(0, 2).every(Number.isFinite)); if (pair) $("observation-angles").value = `${pair[0]}/${pair[1]}`; }
      notify(t("Loaded {name}", { name: file.name }));
    });
  });
  const tableStatus = (text, className, label) => {
    const status = $("band-table-status");
    status.textContent = text; status.className = `file-status ${className}`; status.setAttribute("aria-label", label);
  };
  $("band-table-upload").addEventListener("change", async (event) => {
    const file = event.target.files[0]; if (!file) return;
    const contents = await file.text();
    try {
      parseBandTable(contents).forEach((row, index) => { $(BAND_TABLE_IDS[index]).value = row.map((cell) => cell ?? "-").join(";"); });
    } catch (error) {
      event.target.value = ""; tableStatus("•", "no-file", t("No file selected"));
      notify(`${file.name}: ${error.message}`, true); return;
    }
    tableStatus("✓", "has-file", t("Selected {name}", { name: file.name }));
    updateEmissivityReadout(); notify(t("Loaded {name}", { name: file.name }));
  });
  BAND_TABLE_IDS.forEach((id) => $(id).addEventListener("input", () => {
    $("band-table-upload").value = ""; tableStatus("•", "no-file", t("Using entered values")); updateEmissivityReadout();
  }));
  $("observation-time-mode").addEventListener("change", updateObservationTimeFields);
  const sceneInputs = ["lai", "canopy-height", "hspot", "slope-angle", "slope-aspect", "terrain-height", "terrain-radius", "terrain-density", "row-width", "row-blank", "row-height", "row-azimuth", "crown-rad-a", "crown-rad-b", "crown-height", "crown-density", "building-density", "building-height", "building-length", "building-width", "building-azimuth", "latitude", "longitude"];
  sceneInputs.push("urban-roughness-length", "urban-canyon-roughness-ratio", "urban-canyon-min-ratio", "urban-leaf-wind-decay", "urban-leaf-boundary-factor");
  sceneInputs.forEach((id) => $(id)?.addEventListener("input", () => window.refreshTemporalScene?.()));
  window.addEventListener("resize", () => { if (seriesRows.length) drawTimeCharts(); window.refreshTemporalScene?.(); });
}

bindEvents(); updateConditionalFields(); updateObservationTimeFields(); updateDiagnostics(); updateEmissivityReadout();
