const $ = (id) => document.getElementById(id);
if (window.location.protocol === "file:") window.location.replace("http://127.0.0.1:8765/time.html");
let meteoText = "";
const spectrumFiles = { directions: "", wavelengths: "", vegetation_spectrum: "", soil_spectrum: "", urban_spectrum: "" };
let seriesRows = [];
let observationRows = [];
let toastTimer;
window.temporalState = { get seriesRows() { return seriesRows; }, get observationRows() { return observationRows; } };

const COLORS = { cyan: "#57d3c6", amber: "#edb866", coral: "#ed8e7d", soft: "#b7c3c0", muted: "#687774", grid: "#293735" };
const quantityLabels = {
  Tv_C: "Brightness temperature", Tss_C: "Sunlit soil temperature", Tsh_C: "Shaded soil temperature",
  Tls_C: "Sunlit leaf temperature", Tlh_C: "Shaded leaf temperature", Ta_C: "Air temperature",
  LE: "Latent heat flux", H: "Sensible heat flux", G: "Ground heat flux",
};

function setStatus(text, kind = "idle") { const pill = $("status-pill"); pill.className = `status-pill ${kind}`; pill.querySelector("span").textContent = text; }
function notify(message, error = false) { const toast = $("toast"); toast.textContent = message; toast.classList.toggle("error", error); toast.classList.add("show"); clearTimeout(toastTimer); toastTimer = setTimeout(() => toast.classList.remove("show"), 4200); }
function value(id) { return $(id).value; }
function numberList(raw) { return String(raw || "").split(/[;,:\s]+/).map(Number).filter(Number.isFinite); }
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

function currentOverrides() {
  const wavelengths = numberList(value("wavelengths"));
  const angles = anglePairs(value("observation-angles"));
  const observationVza = angles.map(([vza]) => vza).join(";");
  const observationVaa = angles.map(([, vaa]) => vaa).join(";");
  const observationTimeMode = value("observation-time-mode");
  const observationTimeValues = value("observation-time-values");
  const directionsFromFile = Boolean(spectrumFiles.directions);
  const [startTime, endTime] = timeRange(value("time-range"));
  return {
    surface_model: value("surface-model"), vegetation_model: value("vegetation-model"),
    lai: value("lai"), height: value("height"), hspot: value("hspot"), leaf_width: value("leaf-width"),
    row_width: value("row-width"), row_blank: value("row-blank"), row_height: value("row-height"), row_azimuth: value("row-azimuth"),
    crown_radius: value("crown-radius"), crown_width: value("crown-width"), crown_height: value("crown-height"), crown_density: value("crown-density"),
    terrain_type: value("surface-model") === "slope" ? "1" : value("terrain-type"), terrain_slope: value("terrain-slope"), terrain_aspect: value("terrain-aspect"),
    terrain_height: value("terrain-height"), terrain_radius: value("terrain-radius"), terrain_density: value("terrain-density"),
    building_density: value("building-density"), building_length: value("building-length"), building_width: value("building-width"), building_height: value("building-height"),
    latitude: value("latitude"), longitude: value("longitude"), utc_offset_hours: value("utc-offset"), start_time: startTime, end_time: endTime,
    direction_source: directionsFromFile ? "file" : "direct", observation_direction_mode: directionsFromFile ? "file" : "direct", vza: observationVza, vaa: observationVaa, observation_vza: observationVza, observation_vaa: observationVaa,
    spectral_source: "direct", observation_spectral_mode: "direct", wavelengths: wavelengths.join(";"), observation_wavelengths: wavelengths.join(";"),
    vegetation_spectrum_source: "model", soil_spectrum_source: "model", urban_spectrum_source: "model",
    cab: value("cab"), cw: value("cw"), cdm: value("cdm"), cs: value("cs"), n: value("leaf-n"),
    leaf_reflectance: value("leaf-reflectance-tir"), leaf_transmittance: value("leaf-transmittance-tir"),
    soil_brightness: value("soil-brightness"), soil_latitude: value("soil-latitude"), soil_longitude: value("soil-longitude"), soil_reflectance: value("soil-reflectance-tir"),
    vcmax: value("vcmax"), ballberry_m: value("ballberry-m"), leaf_type: value("leaf-type"), stomatal_method: value("stomatal-method"), stomatal_g0: value("stomatal-g0"),
    rdparam: value("rdparam"), tyear: value("tyear"), beta: value("beta"), knpqs: value("knpqs"), qls: value("qls"), kv: value("kv"), stressfactor: value("stressfactor"),
    tparam1: value("tparam1"), tparam2: value("tparam2"), tparam3: value("tparam3"), tparam4: value("tparam4"), tparam5: value("tparam5"),
    soil_moisture: value("soil-moisture"), soil_rss: value("soil-rss"), soil_cs: value("soil-cs"), soil_rhos: value("soil-rhos"), soil_lambda: value("soil-lambda"), soil_rbs: value("soil-rbs"),
    urban_rss: value("urban-rss"), urban_cs: value("urban-cs"), urban_rhos: value("urban-rhos"), urban_lambda: value("urban-lambda"), urban_rbs: value("urban-rbs"),
    soil_heat_method: value("soil-heat-method"), urban_heat_method: value("urban-heat-method"), temperature_method: value("temperature-method"),
    observation_time_mode: observationTimeMode, observation_times: observationTimeMode === "time" ? observationTimeValues : "", observation_hours: observationTimeMode === "hour" ? observationTimeValues : "", observation_quantity: "temperature",
  };
}

function updateConditionalFields() {
  const surface = value("surface-model");
  const vegetation = value("vegetation-model");
  $("slope-fields").classList.toggle("is-hidden", surface !== "slope" && surface !== "terrain");
  $("terrain-fields").classList.toggle("is-hidden", surface !== "terrain");
  $("urban-fields").classList.toggle("is-hidden", surface !== "urban");
  $("row-fields").classList.toggle("is-hidden", vegetation !== "row");
  $("crown-fields").classList.toggle("is-hidden", vegetation !== "crown");
}

function updateObservationTimeFields() {
  const mode = value("observation-time-mode");
  const input = $("observation-time-values");
  input.disabled = mode === "all";
  input.placeholder = mode === "hour" ? "12 or 12.5" : "214.5";
}

function timeValue(row) { return Number(row.time); }
function finiteRows(rows, key) { return rows.map((row) => ({ x: timeValue(row), y: Number(row[key]) })).filter((point) => Number.isFinite(point.x) && Number.isFinite(point.y)); }
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
  const points = series.flatMap((item) => item.points); if (!points.length) { ctx.textAlign = "center"; ctx.fillText("No finite values", width / 2, height / 2); return; }
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
  const quantity = $("quantity").value; const current = [{ label: quantityLabels[quantity] || quantity, color: COLORS.cyan, points: finiteRows(seriesRows, quantity) }];
  drawChart("time-chart", current, quantity === "LE" || quantity === "H" || quantity === "G" ? "W m⁻²" : "°C");
}

function updateMetrics() {
  const quantity = $("quantity").value; const points = finiteRows(seriesRows, quantity); $("metric-label").textContent = quantityLabels[quantity] || quantity; $("chart-title").textContent = `${quantityLabels[quantity] || quantity} over time`; $("metric-count").textContent = seriesRows.length || "--"; $("series-count").textContent = seriesRows.length ? `${seriesRows.length} time steps · ${observationRows.length} observation records` : "No time series";
  if (!points.length) { $("metric-min").textContent = "--"; $("metric-max").textContent = "--"; $("metric-min-time").textContent = "--"; $("metric-max-time").textContent = "--"; return; }
  const min = points.reduce((a, b) => b.y < a.y ? b : a); const max = points.reduce((a, b) => b.y > a.y ? b : a); const unit = ["LE", "H", "G"].includes(quantity) ? " W m⁻²" : " °C";
  $("metric-min").textContent = `${fmt(min.y, 3)}${unit}`; $("metric-max").textContent = `${fmt(max.y, 3)}${unit}`; $("metric-min-time").textContent = `DOY ${fmt(min.x, 4)}`; $("metric-max-time").textContent = `DOY ${fmt(max.x, 4)}`;
}

function csvCell(value) { const text = String(value ?? ""); return /[",\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text; }
function outputValue(value) { const number = Number(value); return Number.isFinite(number) ? number.toFixed(3) : value; }
function saveRowsCsv(rows, filename) { if (!rows.length) return; const fields = [...new Set(rows.flatMap((row) => Object.keys(row)))]; const body = [fields.map(csvCell).join(","), ...rows.map((row) => fields.map((field) => csvCell(outputValue(row[field]))).join(","))].join("\n"); const link = document.createElement("a"); link.href = URL.createObjectURL(new Blob([body], { type: "text/csv;charset=utf-8" })); link.download = filename; link.click(); URL.revokeObjectURL(link.href); notify(`Saved ${filename}`); }
function saveCsv() { saveRowsCsv(seriesRows, "temporal_output.csv"); }
function saveText(filename, text) { const link = document.createElement("a"); link.href = URL.createObjectURL(new Blob([text], { type: "text/plain;charset=utf-8" })); link.download = filename; link.click(); URL.revokeObjectURL(link.href); notify(`Saved ${filename}`); }
function saveObservations() { saveRowsCsv(observationRows, "observations.csv"); }
function saveInput() {
  const config = currentOverrides();
  const lines = ["# TiRT temporal GUI input", ",,,", "[gui],,,"];
  Object.entries(config).forEach(([key, value]) => { if (value !== "" && value !== undefined && value !== null) lines.push(["gui", key, value, ""].map(csvCell).join(",")); });
  saveText("input.csv", `${lines.join("\n")}\n`);
}

async function runSimulation() {
  setStatus("Running", "busy"); $("run-button").disabled = true;
  try {
    const response = await fetch("/api/run_time", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ overrides: currentOverrides(), meteo_text: meteoText, files: spectrumFiles }) });
    const payload = await response.json(); if (!response.ok || !payload.ok) throw new Error(payload.error || "Time series failed");
    seriesRows = payload.series || []; observationRows = payload.observations || []; $("save-button").disabled = !seriesRows.length; $("save-observations-button").disabled = !observationRows.length; $("empty-state").classList.toggle("is-hidden", Boolean(seriesRows.length)); updateMetrics(); drawTimeCharts(); window.refreshTemporalScene?.(); setStatus("Complete", "idle"); notify(`Time series complete · ${payload.count} time steps`);
  } catch (error) { setStatus("Error", "error"); notify(error.message, true); }
  finally { $("run-button").disabled = false; }
}

function bindEvents() {
  $("run-button").addEventListener("click", runSimulation); $("save-button").addEventListener("click", saveCsv); $("save-observations-button").addEventListener("click", saveObservations); $("save-input-button").addEventListener("click", saveInput); $("reset-button").addEventListener("click", () => window.location.reload()); $("quantity").addEventListener("change", () => { updateMetrics(); drawTimeCharts(); }); $("surface-model").addEventListener("change", () => { updateConditionalFields(); window.refreshTemporalScene?.(); }); $("vegetation-model").addEventListener("change", () => { updateConditionalFields(); window.refreshTemporalScene?.(); });
  $("meteo-upload").addEventListener("change", async (event) => { const file = event.target.files[0]; if (!file) return; meteoText = await file.text(); $("meteo-status").textContent = "✓"; $("meteo-status").className = "file-status has-file"; $("meteo-status").setAttribute("aria-label", `Selected ${file.name}`); notify(`Loaded ${file.name}`); });
  [["angle-upload", "directions", "angle-status"], ["wavelength-upload", "wavelengths", "wavelength-status"], ["vegetation-spectrum-upload", "vegetation_spectrum", "vegetation-spectrum-status"], ["soil-spectrum-upload", "soil_spectrum", "soil-spectrum-status"], ["urban-spectrum-upload", "urban_spectrum", "urban-spectrum-status"]].forEach(([id, key, statusId]) => {
    $(id).addEventListener("change", async (event) => {
      const file = event.target.files[0]; if (!file) return;
      spectrumFiles[key] = await file.text(); const status = $(statusId); status.textContent = "✓"; status.className = "file-status has-file"; status.setAttribute("aria-label", `Selected ${file.name}`);
      if (key === "wavelengths") { const bands = spectrumFiles[key].split(/\r?\n/).map((line) => Number(line.trim().split(/[;,\s]+/)[0])).filter(Number.isFinite); if (bands.length) $("wavelengths").value = bands.join(";"); }
      notify(`Loaded ${file.name}`);
    });
  });
  $("observation-time-mode").addEventListener("change", updateObservationTimeFields);
  const sceneInputs = ["lai", "height", "hspot", "terrain-type", "terrain-slope", "terrain-aspect", "terrain-height", "terrain-radius", "terrain-density", "row-width", "row-blank", "row-height", "row-azimuth", "crown-radius", "crown-width", "crown-height", "crown-density", "building-density", "building-height", "building-length", "building-width", "latitude", "longitude"];
  sceneInputs.forEach((id) => $(id)?.addEventListener("input", () => window.refreshTemporalScene?.()));
  window.addEventListener("resize", () => { if (seriesRows.length) drawTimeCharts(); window.refreshTemporalScene?.(); });
}

bindEvents(); updateConditionalFields(); updateObservationTimeFields();
