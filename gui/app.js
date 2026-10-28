let THREE;
let OrbitControls;
let scene;
let camera;
let renderer;
let controls;
let sceneObjects;
let uploadedGeometryText = "";
let simulationRows = [];
let selectedWavelength = null;
let selectedSceneKey = null;
let toastTimer;

const $ = (id) => document.getElementById(id);
const t = (key, params) => window.tirtI18n.t(key, params);
if (window.location.protocol === "file:") {
  window.location.replace("http://127.0.0.1:8765/");
}
const root = document.documentElement;
const colors = { cyan: "#1683a3", amber: "#d38a24", coral: "#d53b36", muted: "#5d7383", grid: "#c9dae4", ink: "#173044" };
const SCENE_SIZE = 100;
const HALF_SCENE = SCENE_SIZE / 2;
const NORTH_AZIMUTH = 0;

function numberList(value) {
  return String(value || "").trim().split(/[;,\s]+/).filter(Boolean).map(Number);
}

function firstNumber(raw, fallback = 0) {
  const values = numberList(raw).filter(Number.isFinite);
  return values.length ? values[0] : fallback;
}

const BAND_TABLE_IDS = ["wavelengths", "leaf-reflectance", "leaf-transmittance", "soil-reflectance",
  "street-reflectance", "wall-reflectance", "roof-reflectance"];

// Same layout as run.py read_band_table: optional text label, then values; "#" starts a comment.
function parseBandTable(text) {
  const rows = [];
  for (const line of String(text).replace(/^\uFEFF/, "").split(/\r?\n/)) {
    const fields = line.split("#")[0].split(/[\s,;:，；：]+/).filter(Boolean);
    while (fields.length && !Number.isFinite(Number(fields[0]))) fields.shift();
    if (!fields.length) continue;
    const values = fields.map(Number);
    if (!values.every(Number.isFinite)) throw new Error(t("Band table contains a non-numeric value"));
    rows.push(values);
  }
  if (rows.length !== BAND_TABLE_IDS.length) {
    throw new Error(t("Band table needs {count} rows, found {found}", { count: BAND_TABLE_IDS.length, found: rows.length }));
  }
  return rows;
}

function updateEmissivityReadout() {
  const wavelengths = numberList(value("wavelengths"));
  const row = (id) => numberList(value(id));
  const [leafRho, leafTau, soil, street, wall, roof] = BAND_TABLE_IDS.slice(1).map(row);
  const parts = wavelengths.map((wave, index) => {
    if (!(wave >= 2.5)) return null;
    const show = (rho, tau = 0) => fmt(1 - rho - tau, 3);
    return `${wave} µm: ${t("leaf")} ${show(leafRho[index], leafTau[index])} · ${t("soil")} ${show(soil[index])} · ${t("street")} ${show(street[index])} · ${t("wall")} ${show(wall[index])} · ${t("roof")} ${show(roof[index])}`;
  }).filter(Boolean);
  $("thermal-emissivity-readout").textContent = parts.length
    ? `${t("Thermal emissivity ε = 1 − ρ − τ")} — ${parts.join("; ")}`
    : t("No thermal band (≥ 2.5 µm): all bands give reflectance BRF");
}

// Y is vertical; zenith angles are measured from the vertical axis.
function zenithDirection(zenithDegrees, azimuthDegrees, distance = 1, originY = 0) {
  const zenith = zenithDegrees * Math.PI / 180;
  const azimuth = azimuthDegrees * Math.PI / 180;
  return [Math.sin(azimuth) * Math.sin(zenith) * distance,
    originY + Math.cos(zenith) * distance, Math.cos(azimuth) * Math.sin(zenith) * distance];
}

function value(id) { return $(id).value; }
function selectedAdditionalSource() { return document.querySelector('input[name="additional-source"]:checked')?.value || "file"; }
function manualGeometryPairs() {
  return String($("manual-geometry")?.value || "").split(/[;\n]+/)
    .map((entry) => entry.trim()).filter(Boolean)
    .map((entry) => entry.split(/\s*\/\s*/).map((part) => part.trim() ? Number(part) : NaN));
}
function setStatus(text, kind = "idle") {
  const pill = $("status-pill");
  pill.className = `status-pill ${kind}`;
  pill.querySelector("span").textContent = t(text);
}
function notify(message, error = false) {
  const toast = $("toast");
  toast.textContent = message;
  toast.classList.toggle("error", error);
  toast.classList.add("show");
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => toast.classList.remove("show"), 4200);
}

function currentOverrides() {
  const overrides = {
    surface_model: value("surface-model"),
    vegetation_model: value("vegetation-model"),
    lai: value("lai"),
    hspot: value("hspot"),
    geometry_view: "both",
    include_principal: $("include-principal").checked ? "1" : "0",
    include_hemisphere: $("include-hemisphere").checked ? "1" : "0",
    sza: value("sza"),
    saa: value("saa"),
    spectral_source: "direct",
    wavelengths: value("wavelengths"),
    leaf_reflectance: value("leaf-reflectance"),
    leaf_transmittance: value("leaf-transmittance"),
    soil_reflectance: value("soil-reflectance"),
    street_reflectance: value("street-reflectance"),
    wall_reflectance: value("wall-reflectance"),
    roof_reflectance: value("roof-reflectance"),
    hapke_enable: value("hapke-enable"),
    hapke_K: value("hapke-k"),
    building_density: value("building-density"),
    building_height: value("building-height"),
    building_length: value("building-length"),
    building_width: value("building-width"),
    building_azimuth: value("building-azimuth"),
    slope_angle: value("slope-angle"),
    slope_aspect: value("slope-aspect"),
    terrain_height: value("terrain-height"),
    terrain_radius: value("terrain-radius"),
    terrain_density: value("terrain-density"),
    row_width: value("row-width"),
    row_blank: value("row-blank"),
    row_height: value("row-height"),
    row_azimuth: value("row-azimuth"),
    crown_rad_a: value("crown-rad-a"),
    crown_rad_b: value("crown-rad-b"),
    crown_height: value("crown-height"),
    crown_density: value("crown-density"),
    leaf_sunlit_c: value("leaf-sunlit-c"),
    leaf_shaded_c: value("leaf-shaded-c"),
    soil_sunlit_c: value("soil-sunlit-c"),
    soil_shaded_c: value("soil-shaded-c"),
    roof_sunlit_c: value("roof-sunlit-c"),
    roof_shaded_c: value("roof-shaded-c"),
    wall_sunlit_c: value("wall-sunlit-c"),
    wall_shaded_c: value("wall-shaded-c"),
    street_sunlit_c: value("street-sunlit-c"),
    street_shaded_c: value("street-shaded-c"),
    terrain_sunlit_c: value("terrain-sunlit-c"),
    terrain_shaded_c: value("terrain-shaded-c"),
  };
  const additionalSource = selectedAdditionalSource();
  overrides.geometry_mode = "3";
  overrides.principal_vza_step = value("principal-vza-step");
  overrides.principal_vza_max = value("principal-vza-max");
  const includeHemisphere = $("include-hemisphere").checked;
  if (includeHemisphere) overrides.hemisphere_geometry_file = "data/directions/hemisphere_directions.txt";
  overrides.additional_geometry_source = additionalSource;
  if (additionalSource === "file") {
    if (uploadedGeometryText) {
      overrides.additional_geometry_file = "";
    } else {
      overrides.additional_geometry_source = "none";
    }
  } else if (additionalSource === "manual") {
    const pairs = manualGeometryPairs();
    overrides.observation_vza = pairs.map(([vza]) => vza).join(";");
    overrides.observation_vaa = pairs.map(([, vaa]) => vaa).join(";");
  }
  return overrides;
}

function currentPayload() {
  return { overrides: currentOverrides(),
    geometry_text: selectedAdditionalSource() === "file" ? uploadedGeometryText : "" };
}

function updateGeometryFields() {
  const source = selectedAdditionalSource();
  document.querySelectorAll(".source-row").forEach((row) => row.classList.toggle("is-active", row.dataset.source === source));
  updateSceneReadout();
}

function updateConditionalFields() {
  const surface = value("surface-model");
  const vegetation = value("vegetation-model");
  $("slope-structure-fields").classList.toggle("is-hidden", surface !== "slope");
  $("terrain-structure-fields").classList.toggle("is-hidden", surface !== "terrain");
  $("urban-structure-fields").classList.toggle("is-hidden", surface !== "urban");
  $("vegetation-common-fields").classList.toggle("is-hidden", vegetation === "bare");
  $("row-structure-fields").classList.toggle("is-hidden", vegetation !== "row");
  $("crown-structure-fields").classList.toggle("is-hidden", vegetation !== "crown");
  $("bare-structure-note").classList.toggle("is-hidden", vegetation !== "bare");
  $("bare-soil-anisotropy-fields").classList.toggle("is-hidden", vegetation !== "bare");
  $("leaf-temperature-row").classList.toggle("is-hidden", vegetation === "bare");
  $("terrain-temperature-row").classList.toggle("is-hidden", surface !== "terrain");
  $("roof-temperature-row").classList.toggle("is-hidden", surface !== "urban");
  $("wall-temperature-row").classList.toggle("is-hidden", surface !== "urban");
  $("street-temperature-row").classList.toggle("is-hidden", surface !== "urban");
}

function updateSceneReadout() {
  const sza = firstNumber(value("sza"), 30);
  const saa = firstNumber(value("saa"), 0);
  $("solar-readout").textContent = `SZA ${fmt(sza, 0)}° · SAA ${fmt(saa, 0)}°`;
  const surfaceLabel = $("surface-model").selectedOptions[0].textContent;
  const vegetationLabel = $("vegetation-model").selectedOptions[0].textContent;
  $("scene-title").textContent = `${surfaceLabel} · ${vegetationLabel}`;
  $("scene-readout").textContent = `${surfaceLabel} / ${vegetationLabel}`;
  const principalCount = Math.ceil(Number(value("principal-vza-max")) / Math.max(Number(value("principal-vza-step")), 0.1)) + 1;
  const additionalCount = selectedAdditionalSource() === "manual" ? manualGeometryPairs().length : 0;
  const parts = [];
  if ($("include-principal").checked) parts.push(t("{count} principal", { count: principalCount * 4 }));
  if ($("include-hemisphere").checked) parts.push(t("hemisphere"));
  if (selectedAdditionalSource() === "file" && !$("include-hemisphere").checked) parts.push(t("file"));
  if (additionalCount) parts.push(t("{count} manual", { count: additionalCount }));
  let count = parts.join(" + ") || t("No directions");
  if (simulationRows.length) count = t("{count} directions", { count: new Set(simulationRows.map((row) => `${row.scenario_id || 1}|${row.vza}|${row.vaa}|${row.sza}|${row.saa}`)).size });
  $("observation-readout").textContent = count;
}

function fieldName(input) {
  const labelText = input.closest("label")?.firstChild?.textContent?.trim();
  return input.getAttribute("aria-label") || labelText || input.id;
}

function validateInputs() {
  const activeNumberInputs = [...document.querySelectorAll('input[type="number"]')]
    .filter((input) => !input.closest(".is-hidden"));
  const invalidNumber = activeNumberInputs.find((input) => {
    return input.value.trim() === "" || !Number.isFinite(input.valueAsNumber) || !input.checkValidity();
  });
  if (invalidNumber) {
    return { input: invalidNumber, message: t("Invalid value for {field}", { field: fieldName(invalidNumber) }) };
  }

  const wavelengths = numberList(value("wavelengths"));
  if (!wavelengths.length || wavelengths.some((item) => !Number.isFinite(item) || item < 0.4)) {
    return { input: $("wavelengths"), message: t("Wavelengths must be numbers of at least 0.4 µm") };
  }
  const laiValues = numberList(value("lai"));
  if (!laiValues.length || laiValues.some((item) => !Number.isFinite(item) || item < 0 || item > 20)) {
    return { input: $("lai"), message: t("LAI must contain numbers in 0..20, for example 1;3;5") };
  }
  const hotspotValues = numberList(value("hspot"));
  if (!$("hspot").closest(".is-hidden") && (!hotspotValues.length || hotspotValues.some((item) => !Number.isFinite(item) || item < 0))) {
    return { input: $("hspot"), message: t("Hotspot must contain nonnegative numbers, for example 0.1;0.2") };
  }

  const propertyIds = [
    ["leaf-reflectance", "Leaf reflectance"],
    ["leaf-transmittance", "Leaf transmittance"],
    ["soil-reflectance", "Soil reflectance"],
    ["street-reflectance", "Street reflectance"],
    ["wall-reflectance", "Wall reflectance"],
    ["roof-reflectance", "Roof reflectance"],
  ];
  for (const [id, label] of propertyIds) {
    const input = $(id);
    const values = numberList(input.value);
    if (!values.length || values.some((item) => !Number.isFinite(item) || item < 0 || item > 1)) {
      return { input, message: t("{label} must be between 0 and 1", { label: t(label) }) };
    }
    if (values.length !== wavelengths.length) {
      return { input, message: t("{label}: enter {count} values, one per wavelength", { label: t(label), count: wavelengths.length }) };
    }
  }
  const leafRho = numberList(value("leaf-reflectance"));
  const leafTau = numberList(value("leaf-transmittance"));
  if (leafRho.some((item, index) => item + leafTau[index] > 1)) {
    return { input: $("leaf-reflectance"), message: t("Leaf reflectance + transmittance must not exceed 1") };
  }

  const source = selectedAdditionalSource();
  if (source === "manual" && manualGeometryPairs().some((pair) => pair.length !== 2 ||
      !pair.every(Number.isFinite) || pair[0] < 0 || pair[0] > 90 || pair[1] < 0 || pair[1] > 360)) {
    return { input: $("manual-geometry"), message: t("Manual angles must be VZA/VAA pairs with VZA 0–90° and VAA 0–360°") };
  }
  const hasAdditional = source === "manual" ? manualGeometryPairs().length > 0 : Boolean(uploadedGeometryText);
  if (!$("include-principal").checked && !$("include-hemisphere").checked && !hasAdditional) {
    const input = source === "manual" ? $("manual-geometry") : $("geometry-upload");
    return { input, message: t("Select an observation plane or provide observation angles") };
  }
  return null;
}

function setRunBusy(busy) {
  const button = $("run-button");
  button.disabled = busy;
  button.setAttribute("aria-busy", String(busy));
  button.classList.toggle("is-running", busy);
  button.querySelector(".button-label").textContent = t(busy ? "Running..." : "Run simulation");
}

async function runSimulation() {
  const validation = validateInputs();
  if (validation) {
    setStatus("Input error", "error");
    notify(validation.message, true);
    validation.input?.focus();
    return;
  }

  const started = Date.now();
  const timer = setInterval(() => {
    setStatus(t("Running {seconds} s", { seconds: Math.floor((Date.now() - started) / 1000) }), "busy");
  }, 1000);
  setStatus("Running", "busy");
  setRunBusy(true);
  try {
    const response = await fetch("/api/run", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(currentPayload()),
    });
    const payload = await response.json();
    if (!response.ok || !payload.ok) throw new Error(payload.error || t("Simulation failed"));
    simulationRows = payload.rows;
    $("save-output-button").disabled = false;
    updateBandOptions();
    updateSceneOptions();
    renderResults();
    updateSceneReadout();
    updateScene();
    setStatus("Complete", "idle");
    notify(t("Simulation complete - {count} output records", { count: payload.count }));
  } catch (error) {
    setStatus("Error", "error");
    notify(error.message, true);
  } finally {
    clearInterval(timer);
    setRunBusy(false);
  }
}

function updateBandOptions() {
  const bands = [...new Set(simulationRows.map((row) => Number(row.wavelength_um)))].sort((a, b) => a - b);
  const select = $("band-select");
  const previous = Number(selectedWavelength);
  select.innerHTML = "";
  bands.forEach((band) => {
    const option = document.createElement("option");
    option.value = String(band);
    const kind = Number(band) < 2.5 ? "BRF" : "BT";
    option.textContent = `${fmt(band, Number(band) < 2 ? 2 : 1)} µm · ${kind}`;
    select.appendChild(option);
  });
  if (bands.length) {
    selectedWavelength = bands.includes(previous)
      ? previous
      : (bands.includes(10.5) ? 10.5 : bands[0]);
    select.value = String(selectedWavelength);
  }
}

function sceneKey(row) { return `${Number(row.lai)}|${Number(row.scenario_id) || 1}`; }

function sceneList(rows) {
  const scenes = new Map();
  rows.forEach((row) => {
    const key = sceneKey(row);
    if (!scenes.has(key)) scenes.set(key, { key, lai: Number(row.lai), scenarioId: Number(row.scenario_id) || 1, scenarioLabel: row.scenario_label || "default" });
  });
  const list = [...scenes.values()].sort((a, b) => a.lai - b.lai || a.scenarioId - b.scenarioId);
  const multipleLai = new Set(list.map((scene) => scene.lai)).size > 1;
  list.forEach((scene) => {
    const parts = [];
    if (multipleLai && Number.isFinite(scene.lai)) parts.push(`LAI=${scene.lai}`);
    if (scene.scenarioLabel !== "default") parts.push(...scene.scenarioLabel.split("; ").map((part) => {
      const [name, setting] = part.split("=");
      return setting === undefined ? t(part) : `${t(name)}=${setting}`;
    }));
    scene.label = parts.join("; ") || t("default");
  });
  return list;
}

function updateSceneOptions() {
  const scenes = sceneList(simulationRows);
  const select = $("scene-select");
  select.closest("label").classList.toggle("is-hidden", scenes.length <= 1);
  select.innerHTML = "";
  scenes.forEach((item) => {
    const option = document.createElement("option");
    option.value = item.key;
    option.textContent = item.label;
    select.appendChild(option);
  });
  if (scenes.length) {
    selectedSceneKey = scenes.some((item) => item.key === selectedSceneKey) ? selectedSceneKey : scenes[0].key;
    select.value = selectedSceneKey;
  }
}

function bandRows() {
  const wavelength = Number($("band-select").value);
  selectedWavelength = wavelength;
  return simulationRows.filter((row) => Math.abs(Number(row.wavelength_um) - wavelength) < 1e-7);
}

function selectedRows() {
  const rows = bandRows();
  if (selectedSceneKey === null) return rows;
  return rows.filter((row) => sceneKey(row) === selectedSceneKey);
}

function isShortwaveBand(band) { return Number(band) < 2.5; }
function bandQuantity(band = $("band-select")?.value) {
  return isShortwaveBand(band) ? "reflectance_brf" : "brightness_temperature_C";
}
function outputValue(row, quantity = bandQuantity()) {
  const raw = row[quantity];
  return raw === null || raw === undefined || raw === "" ? NaN : Number(raw);
}
function quantityValue(row) { return outputValue(row, bandQuantity()); }
function quantityLabel(quantity = bandQuantity()) {
  return t(quantity === "reflectance_brf" ? "Reflectance (BRF)" : "Brightness temperature (°C)");
}
function fmt(number, digits = 2) { return Number.isFinite(Number(number)) ? Number(number).toFixed(digits) : "--"; }
function fmtSignificant(number, significant = 3) {
  const numeric = Number(number);
  if (!Number.isFinite(numeric)) return "--";
  if (numeric === 0) return "0";
  const exponent = Math.floor(Math.log10(Math.abs(numeric)));
  return numeric.toFixed(Math.max(significant - 1 - exponent, 0));
}
function formatOutputNumber(number) { return fmtSignificant(number, 3); }
function formatOutputValue(number, quantity) {
  return `${formatOutputNumber(number)}${quantity === "brightness_temperature_C" ? " °C" : ""}`;
}

function renderResults() {
  if (!simulationRows.length) return;
  $("empty-state").classList.add("is-hidden");
  $("results-content").classList.remove("is-hidden");
  const rows = selectedRows();
  const values = rows.map(quantityValue).filter(Number.isFinite);
  const quantity = bandQuantity();
  const mean = values.length ? values.reduce((sum, item) => sum + item, 0) / values.length : NaN;
  const span = values.length ? Math.max(...values) - Math.min(...values) : NaN;
  $("mean-bt").closest(".metric").querySelector("span").textContent = t("Mean {quantity}", { quantity: quantityLabel(quantity) });
  $("mean-bt").textContent = formatOutputValue(mean, quantity);
  $("range-bt").textContent = fmtSignificant(span);
  $("metric-band").textContent = `${fmt(selectedWavelength, 2)} µm`;
  $("metric-band-kind").textContent = quantityLabel(quantity);
  $("observation-count").textContent = String(rows.length);
  renderCharts(rows);
  updateChartVisibility();
}

function observationTableRows(rows) {
  const bands = [...new Set(rows.map((row) => Number(row.wavelength_um)).filter(Number.isFinite))].sort((a, b) => a - b);
  const observations = new Map();
  rows.forEach((row) => {
    const scenarioId = Number(row.scenario_id) || 1; const scenarioLabel = row.scenario_label || "default";
    const lai = Number(row.lai);
    const vza = Number(row.vza); const vaa = Number(row.vaa); const raa = Number(row.raa);
    const sza = Number(row.sza); const saa = Number(row.saa); const wavelength = Number(row.wavelength_um);
    const output = outputValue(row, bandQuantity(wavelength));
    if (![lai, vza, vaa, raa, wavelength, output].every(Number.isFinite)) return;
    const key = [scenarioId, lai, vza, vaa, sza, saa].map((value) => Number.isFinite(value) ? value : "").join("|");
    let observation = observations.get(key);
    if (!observation) {
      observation = { scenarioId, scenarioLabel, lai, vza, vaa, raa, values: new Map() };
      observations.set(key, observation);
    }
    observation.values.set(wavelength, output);
  });
  return { bands, observationRows: [...observations.values()] };
}

function csvEscape(value) {
  const text = value == null ? "" : String(value);
  return /[",\n\r]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
}

function exportOutput() {
  if (!simulationRows.length) {
    notify(t("Run a simulation first"), true);
    return;
  }
  const { bands, observationRows } = observationTableRows(simulationRows);
  const fields = ["index", "scenario_id", "scenario_label", "lai", "vza", "vaa", "raa",
    ...bands.map((band) => `${fmt(band, 2)}um_${isShortwaveBand(band) ? "brf" : "temperature_C"}`)];
  const lines = [
    fields.join(","),
    ...observationRows.map((observation, index) => [
      index, observation.scenarioId, observation.scenarioLabel, observation.lai, observation.vza, observation.vaa, observation.raa,
      ...bands.map((band) => formatOutputNumber(observation.values.get(band))),
    ].map(csvEscape).join(",")),
  ];
  const blob = new Blob([`${lines.join("\n")}\n`], { type: "text/csv;charset=utf-8" });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = "observations.csv";
  link.click();
  setTimeout(() => URL.revokeObjectURL(link.href), 0);
  notify(t("Saved observations.csv · {count} directions · BRF + brightness temperature (°C)", { count: observationRows.length }));
}

function canvasContext(canvas) {
  const rect = canvas.getBoundingClientRect();
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  canvas.width = Math.max(1, Math.floor(rect.width * dpr));
  canvas.height = Math.max(1, Math.floor(rect.height * dpr));
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, rect.width, rect.height);
  return { ctx, width: rect.width, height: rect.height };
}

function colorFor(value, min, max) {
  const t = max === min ? .5 : (value - min) / (max - min);
  const mix = Math.max(0, Math.min(1, t));
  // A high-contrast blue-white-red scale keeps adjacent filled cells distinct.
  const stops = [[30, 74, 161], [245, 245, 245], [190, 28, 28]];
  const segment = mix < .5 ? 0 : 1;
  const local = segment === 0 ? mix * 2 : (mix - .5) * 2;
  const left = stops[segment];
  const right = stops[segment + 1];
  return `rgb(${left.map((item, index) => Math.round(item + (right[index] - item) * local)).join(",")})`;
}

function fillPolarCell(ctx, cx, cy, innerRadius, outerRadius, startAngle, endAngle, color) {
  ctx.fillStyle = color;
  ctx.beginPath();
  if (innerRadius < 1e-6) {
    ctx.moveTo(cx, cy);
    ctx.arc(cx, cy, outerRadius, startAngle, endAngle);
  } else {
    ctx.arc(cx, cy, outerRadius, startAngle, endAngle);
    ctx.arc(cx, cy, innerRadius, endAngle, startAngle, true);
  }
  ctx.closePath();
  ctx.fill();
}

function drawPolar(rows) {
  const { ctx, width, height } = canvasContext($("polar-chart"));
  const cx = width / 2, cy = height / 2 + 4, radius = Math.min(width, height) * .36;
  const finiteRows = rows.filter((row) => Number.isFinite(quantityValue(row)) && Number.isFinite(Number(row.vza)) && Number.isFinite(Number(row.vaa)));
  if (!finiteRows.length) return;
  const values = finiteRows.map(quantityValue); const min = Math.min(...values), max = Math.max(...values);
  const configuredMax = Math.max(Number(value("principal-vza-max")), 1);
  const dataMax = Math.max(0, ...finiteRows.map((row) => Number(row.vza)));
  const vzaMax = Math.max(configuredMax, dataMax);
  const normalizeAzimuth = (angle) => ((angle % 360) + 360) % 360;
  const keyFor = (vza, vaa) => `${vza.toFixed(8)}|${normalizeAzimuth(vaa).toFixed(8)}`;
  const buckets = new Map();
  finiteRows.forEach((row) => {
    const vza = Number(row.vza); const vaa = normalizeAzimuth(Number(row.vaa)); const key = keyFor(vza, vaa);
    const bucket = buckets.get(key) || { vza, vaa, sum: 0, count: 0 };
    bucket.sum += quantityValue(row); bucket.count += 1; buckets.set(key, bucket);
  });
  const levels = [...new Set([...buckets.values()].map((item) => item.vza))].sort((a, b) => a - b);
  const azimuths = [...new Set([...buckets.values()].map((item) => item.vaa))].sort((a, b) => a - b);
  if (levels.length < 2 || !azimuths.length) return;
  const radialMin = levels[0]; const radialMax = levels[levels.length - 1];
  const levelSamples = levels.map((level) => azimuths.map((azimuth) => {
    const bucket = buckets.get(keyFor(level, azimuth));
    return bucket ? bucket.sum / bucket.count : NaN;
  }));
  const circularInterpolate = (points, target) => {
    const valid = points.map((item, index) => ({ angle: azimuths[index], value: item })).filter((item) => Number.isFinite(item.value));
    if (!valid.length) return NaN;
    if (valid.length === 1) return valid[0].value;
    let angle = normalizeAzimuth(target);
    if (angle < valid[0].angle) angle += 360;
    let lower = valid[valid.length - 1]; let upper = { angle: valid[0].angle + 360, value: valid[0].value };
    for (let index = 0; index < valid.length - 1; index += 1) {
      if (angle <= valid[index + 1].angle) { lower = valid[index]; upper = valid[index + 1]; break; }
    }
    const fraction = (angle - lower.angle) / Math.max(upper.angle - lower.angle, 1e-9);
    return lower.value + (upper.value - lower.value) * fraction;
  };
  const sample = (vza, vaa) => {
    if (vza < radialMin - 1e-9 || vza > radialMax + 1e-9) return NaN;
    const radial = levelSamples.map((valuesAtLevel) => circularInterpolate(valuesAtLevel, vaa));
    const valid = radial.map((item, index) => ({ level: levels[index], value: item })).filter((item) => Number.isFinite(item.value));
    if (!valid.length) return NaN;
    if (valid.length === 1 || vza <= valid[0].level) return valid[0].value;
    for (let index = 0; index < valid.length - 1; index += 1) {
      if (vza <= valid[index + 1].level) {
        const fraction = (vza - valid[index].level) / Math.max(valid[index + 1].level - valid[index].level, 1e-9);
        return valid[index].value + (valid[index + 1].value - valid[index].value) * fraction;
      }
    }
    return valid[valid.length - 1].value;
  };
  // Render a dense interpolated polar field.
  const radialSteps = Math.max(1, Math.ceil((vzaMax - radialMin) / 1));
  const azimuthSteps = 180;
  for (let radialIndex = 0; radialIndex < radialSteps; radialIndex += 1) {
    const innerVza = radialMin + (vzaMax - radialMin) * radialIndex / radialSteps;
    const outerVza = radialMin + (vzaMax - radialMin) * (radialIndex + 1) / radialSteps;
    for (let azimuthIndex = 0; azimuthIndex < azimuthSteps; azimuthIndex += 1) {
      const startAzimuth = azimuthIndex * 360 / azimuthSteps;
      const cellValue = sample((innerVza + outerVza) / 2, startAzimuth + 180 / azimuthSteps);
      if (!Number.isFinite(cellValue)) continue;
      const valueForCell = cellValue;
      const startAngle = -Math.PI / 2 + startAzimuth * Math.PI / 180;
      const endAngle = -Math.PI / 2 + (startAzimuth + 360 / azimuthSteps) * Math.PI / 180;
      fillPolarCell(ctx, cx, cy, radius * innerVza / vzaMax, radius * outerVza / vzaMax, startAngle, endAngle, colorFor(valueForCell, min, max));
    }
  }
  ctx.font = "10px system-ui"; ctx.textAlign = "center"; ctx.textBaseline = "middle";
  [vzaMax / 3, vzaMax * 2 / 3, vzaMax].forEach((ring) => { const r = radius * ring / vzaMax; ctx.strokeStyle = colors.grid; ctx.lineWidth = 1; ctx.beginPath(); ctx.arc(cx, cy, r, 0, Math.PI * 2); ctx.stroke(); ctx.fillStyle = colors.muted; ctx.fillText(`${fmt(ring, 0)}°`, cx + 9, cy - r + 1); });
  [0, 90, 180, 270].forEach((angle) => { const rad = angle * Math.PI / 180; ctx.strokeStyle = "#b8cbd7"; ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(cx + Math.sin(rad) * radius, cy - Math.cos(rad) * radius); ctx.stroke(); ctx.fillStyle = colors.muted; ctx.fillText(`${angle}°`, cx + Math.sin(rad) * (radius + 16), cy - Math.cos(rad) * (radius + 16)); });
  $("polar-colorbar-min").textContent = fmt(min);
  $("polar-colorbar-max").textContent = fmt(max);
  ctx.fillStyle = colors.muted; ctx.textAlign = "left"; ctx.fillText(t("VZA radius · {band} µm · {quantity}", { band: fmt(selectedWavelength, 2), quantity: quantityLabel() }), 8, height - 10);
}

function planeRows(rows, plane) {
  return rows.filter((row) => {
    const delta = ((Number(row.vaa) - Number(row.saa)) % 360 + 360) % 360;
    if (plane === "parallel") return Math.min(delta, 360 - delta) < 1e-5 || Math.abs(delta - 180) < 1e-5;
    return Math.abs(delta - 90) < 1e-5 || Math.abs(delta - 270) < 1e-5;
  }).map((row) => {
    const delta = ((Number(row.vaa) - Number(row.saa)) % 360 + 360) % 360;
    const negative = plane === "parallel"
      ? Math.abs(delta - 180) < 1e-5
      : Math.abs(delta - 270) < 1e-5;
    return { row, x: Number(row.vza) * (negative ? -1 : 1), negative };
  }).sort((a, b) => a.x - b.x);
}

function lineChartRange(rows) {
  const values = ["parallel", "perpendicular"].flatMap((plane) => planeRows(rows, plane).map((point) => quantityValue(point.row))).filter(Number.isFinite);
  const yMin = values.length ? Math.min(...values) : 0;
  const yMax = values.length ? Math.max(...values) : 1;
  const margin = Math.max((yMax - yMin) * .15, Math.abs(yMax) > 10 ? .5 : .01);
  return { lo: yMin - margin, hi: yMax + margin };
}

function drawLineChart(rows, canvasId, plane, scale, color) {
  const { ctx, width, height } = canvasContext($(canvasId));
  const points = planeRows(rows, plane);
  const pad = { left: 38, right: 12, top: 26, bottom: 29 };
  const xMax = Math.max(30, ...points.map((point) => Math.abs(point.x)), 80);
  const { lo, hi } = scale; const plotW = width - pad.left - pad.right, plotH = height - pad.top - pad.bottom;
  ctx.font = "10px system-ui"; ctx.strokeStyle = colors.grid; ctx.fillStyle = colors.muted; ctx.lineWidth = 1;
  [-xMax, 0, xMax].forEach((x) => { const px = pad.left + (x + xMax) / (xMax * 2) * plotW; ctx.beginPath(); ctx.moveTo(px, pad.top); ctx.lineTo(px, pad.top + plotH); ctx.stroke(); ctx.textAlign = "center"; ctx.fillText(`${fmt(x, 0)}°`, px, height - 11); });
  [lo, (lo + hi) / 2, hi].forEach((y) => { const py = pad.top + (hi - y) / (hi - lo) * plotH; ctx.beginPath(); ctx.moveTo(pad.left, py); ctx.lineTo(pad.left + plotW, py); ctx.stroke(); ctx.textAlign = "right"; ctx.fillText(fmt(y), pad.left - 7, py); });
  if (!points.length) { ctx.fillStyle = colors.muted; ctx.textAlign = "center"; ctx.fillText(t("No matching directions"), width / 2, height / 2); return; }
  const px = (x) => pad.left + (x + xMax) / (xMax * 2) * plotW; const py = (y) => pad.top + (hi - y) / (hi - lo) * plotH;
  const branches = [points.filter((point) => point.negative), points.filter((point) => !point.negative)];
  branches.forEach((branch) => {
    ctx.strokeStyle = color; ctx.lineWidth = 2; ctx.beginPath();
    branch.forEach((point, index) => { if (index) ctx.lineTo(px(point.x), py(quantityValue(point.row))); else ctx.moveTo(px(point.x), py(quantityValue(point.row))); });
    ctx.stroke();
  });
  points.forEach((point) => { ctx.fillStyle = color; ctx.beginPath(); ctx.arc(px(point.x), py(quantityValue(point.row)), 3.5, 0, Math.PI * 2); ctx.fill(); });
  ctx.textAlign = "center"; ctx.fillStyle = colors.muted;
  const sideText = t(plane === "perpendicular" ? "RAA 90° / 270° · both sides" : "RAA 0° / 180° · both sides");
  ctx.fillText(`${sideText} · ${t("{count} directions", { count: points.length })}`, width / 2, 12);
}

function renderCharts(rows) {
  drawPolar(rows);
  const scale = lineChartRange(rows);
  drawLineChart(rows, "parallel-chart", "parallel", scale, colors.amber);
  drawLineChart(rows, "perpendicular-chart", "perpendicular", scale, colors.coral);
}

function updateChartVisibility() {
  $("polar-card").classList.toggle("is-hidden", !$("show-polar").checked);
  $("parallel-card").classList.toggle("is-hidden", !$("show-parallel").checked);
  $("perpendicular-card").classList.toggle("is-hidden", !$("show-perpendicular").checked);
  if (simulationRows.length) renderCharts(selectedRows());
}

function makeMaterial(color, roughness = .8) { return new THREE.MeshStandardMaterial({ color, roughness, metalness: .05 }); }
function addMesh(geometry, material, position = [0, 0, 0], rotation = [0, 0, 0]) { const mesh = new THREE.Mesh(geometry, material); mesh.position.set(...position); mesh.rotation.set(...rotation); sceneObjects.add(mesh); return mesh; }
function addLine(start, end, color, opacity = 1) { const geometry = new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(...start), new THREE.Vector3(...end)]); const material = new THREE.LineBasicMaterial({ color, transparent: opacity < 1, opacity }); sceneObjects.add(new THREE.Line(geometry, material)); }
function addNorthMarker(baseHeight = 0) {
  const material = makeMaterial(0x57d3c6, .7);
  const north = NORTH_AZIMUTH * Math.PI / 180;
  const end = [Math.sin(north) * 10, baseHeight + .08, Math.cos(north) * 10];
  addLine([0, baseHeight + .08, 0], end, 0x57d3c6, .9);
  addMesh(new THREE.ConeGeometry(.55, 1.8, 8), material, [0, baseHeight + .08, 10.7], [Math.PI / 2, 0, 0]);
}
function fitSceneView() {
  if (!camera || !controls || !sceneObjects) return;
  const bounds = new THREE.Box3().setFromObject(sceneObjects);
  if (bounds.isEmpty()) return;
  const center = bounds.getCenter(new THREE.Vector3());
  const size = bounds.getSize(new THREE.Vector3());
  const maxDimension = Math.max(size.x, size.y, size.z, 1);
  const distance = maxDimension / (2 * Math.tan(camera.fov * Math.PI / 360)) * 1.25;
  const direction = camera.position.clone().sub(controls.target);
  if (direction.lengthSq() < 1.0e-6) direction.set(1, .75, 1);
  direction.normalize();
  controls.target.copy(center);
  camera.position.copy(center).add(direction.multiplyScalar(Math.max(distance, 80)));
  camera.near = Math.max(distance / 1000, .1);
  camera.far = Math.max(distance * 4, 250);
  camera.updateProjectionMatrix();
  controls.update();
}

function updateScene() {
  if (!scene || !THREE) return;
  sceneObjects.clear();
  const surface = value("surface-model"); const vegetation = value("vegetation-model");
  const ground = makeMaterial(surface === "urban" ? 0x3a4642 : surface === "slope" || surface === "terrain" ? 0x5c6552 : 0x384d42);
  const slopeAngle = firstNumber(value("slope-angle"), 30) * Math.PI / 180;
  const slopeAspect = firstNumber(value("slope-aspect"), 0) * Math.PI / 180;
  const seeded = (index) => { const x = Math.sin(index * 12.9898 + 78.233) * 43758.5453; return x - Math.floor(x); };
  const terrainHeight = Math.max(Number(value("terrain-height")) || 20, 0);
  const terrainRadius = Math.max(Number(value("terrain-radius")) || 10, 5);
  const terrainDensity = Math.max(Number(value("terrain-density")) || .001, .0001);
  const terrainMounds = [];
  if (surface === "terrain") {
    const moundCount = Math.min(24, Math.max(4, Math.round(terrainDensity * SCENE_SIZE * SCENE_SIZE)));
    for (let index = 0; index < moundCount; index++) {
      terrainMounds.push({
        x: seeded(index * 3 + 20) * (SCENE_SIZE - terrainRadius * 2) - (HALF_SCENE - terrainRadius),
        z: seeded(index * 3 + 21) * (SCENE_SIZE - terrainRadius * 2) - (HALF_SCENE - terrainRadius),
        height: terrainHeight * (.45 + seeded(index * 3 + 22) * .55),
        radius: terrainRadius * (.7 + seeded(index * 3 + 23) * .6),
      });
    }
  }
  const surfaceHeightAt = (x, z) => {
    if (surface === "slope") return -Math.tan(slopeAngle) * (x * Math.sin(slopeAspect) + z * Math.cos(slopeAspect));
    if (surface !== "terrain") return 0;
    return terrainMounds.reduce((height, mound) => {
      const distance = Math.hypot(x - mound.x, z - mound.z) / mound.radius;
      return distance < 1 ? Math.max(height, mound.height * (1 - distance * distance)) : height;
    }, 0);
  };
  const surfaceNormalAt = (x, z) => {
    const delta = .5;
    const dx = (surfaceHeightAt(x + delta, z) - surfaceHeightAt(x - delta, z)) / (2 * delta);
    const dz = (surfaceHeightAt(x, z + delta) - surfaceHeightAt(x, z - delta)) / (2 * delta);
    return new THREE.Vector3(-dx, 1, -dz).normalize();
  };
  const surfaceGeometry = new THREE.PlaneGeometry(SCENE_SIZE, SCENE_SIZE, surface === "plane" || surface === "urban" ? 1 : 30, surface === "plane" || surface === "urban" ? 1 : 30);
  surfaceGeometry.rotateX(-Math.PI / 2);
  const positions = surfaceGeometry.attributes.position;
  for (let index = 0; index < positions.count; index++) positions.setY(index, surfaceHeightAt(positions.getX(index), positions.getZ(index)));
  positions.needsUpdate = true; surfaceGeometry.computeVertexNormals();
  addMesh(surfaceGeometry, ground);
  if (surface === "plane" || surface === "urban") {
    const grid = new THREE.GridHelper(SCENE_SIZE, 20, 0x45615b, 0x263732); grid.position.y = .02; sceneObjects.add(grid);
  } else {
    const wire = new THREE.LineSegments(new THREE.WireframeGeometry(surfaceGeometry), new THREE.LineBasicMaterial({ color: 0x45615b, transparent: true, opacity: .35 })); sceneObjects.add(wire);
  }
  if (surface === "urban") {
    const building = makeMaterial(0x79827b); const buildingHeight = Math.max(Number(value("building-height")), 1); const buildingLength = Math.max(Number(value("building-length")), 1); const buildingWidth = Math.max(Number(value("building-width")), 1);
    [[-20, .08, -12], [4, .14, -8], [24, .1, 12], [-12, .06, 20]].forEach(([x, scale, z]) => { const height = buildingHeight * scale; addMesh(new THREE.BoxGeometry(buildingLength * (z > 0 ? .8 : 1), height, buildingWidth * (z > 0 ? 1.2 : .8)), building, [x, height / 2, z], [0, Number(value("building-azimuth")) * Math.PI / 180, 0]); });
  }
  if (vegetation === "row") {
    const rowWidth = Math.max(firstNumber(value("row-width"), .5), .1); const rowGap = Math.max(firstNumber(value("row-blank"), .5), .05); const pitch = rowWidth + rowGap;
    const rowAngle = firstNumber(value("row-azimuth"), 0) * Math.PI / 180; const along = [Math.cos(rowAngle), 0, Math.sin(rowAngle)]; const across = [-Math.sin(rowAngle), 0, Math.cos(rowAngle)];
    const rowFacet = new THREE.ConeGeometry(.18, .9, 6); const rowMaterial = makeMaterial(0x6dbb86); const rowOffsetStep = Math.max(rowWidth / 2, .65);
    let rowIndex = 0;
    for (let cross = -HALF_SCENE; cross <= HALF_SCENE; cross += pitch) {
      for (let offset = -rowWidth / 2; offset <= rowWidth / 2; offset += rowOffsetStep) {
        for (let distance = -HALF_SCENE; distance <= HALF_SCENE; distance += 1.35) {
          const x = along[0] * distance + across[0] * (cross + offset); const z = along[2] * distance + across[2] * (cross + offset);
          const rowMesh = addMesh(rowFacet, rowMaterial, [x, surfaceHeightAt(x, z) + .42 + (rowIndex % 3) * .06, z]); rowMesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), surfaceNormalAt(x, z)); rowIndex += 1;
        }
      }
    }
  } else if (vegetation === "crown") {
    const crownRadA = Math.max(firstNumber(value("crown-rad-a"), 3), .35); const crownRadB = Math.max(firstNumber(value("crown-rad-b"), 3), .35); const crownHeight = Math.max(firstNumber(value("crown-height"), 3), 1);
    for (let index = 0; index < 24; index++) {
      const x = seeded(index * 2) * (SCENE_SIZE - 8) - (HALF_SCENE - 4); const z = seeded(index * 2 + 1) * (SCENE_SIZE - 8) - (HALF_SCENE - 4); const height = crownHeight * (.75 + seeded(index + 100) * .5);
      const baseHeight = surfaceHeightAt(x, z); const trunk = addMesh(new THREE.CylinderGeometry(.13, .2, height, 6), makeMaterial(0x725a45), [x, baseHeight + height / 2, z]);
      trunk.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), surfaceNormalAt(x, z));
      const crownScale = .75 + seeded(index + 200) * .35;
      const crown = addMesh(new THREE.SphereGeometry(crownRadA * crownScale, 12, 8), makeMaterial(0x4d9872), [x, baseHeight + height + crownRadB * crownScale * .65, z]);
      crown.scale.set(1, crownRadB / crownRadA, 1);
      crown.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), surfaceNormalAt(x, z));
    }
  } else if (vegetation === "hom") {
    const facet = new THREE.CylinderGeometry(.24, .24, .18, 10); const facetMesh = new THREE.InstancedMesh(facet, makeMaterial(0x65ab80), 1150); const dummy = new THREE.Object3D();
    const facetCount = 900 + Math.round(Math.min(Math.max(firstNumber(value("lai"), 1.5), .1), 5) * 50);
    for (let index = 0; index < facetCount; index++) {
      const x = seeded(index * 2 + 1) * (SCENE_SIZE - 4) - (HALF_SCENE - 2);
      const z = seeded(index * 2 + 2) * (SCENE_SIZE - 4) - (HALF_SCENE - 2);
      const y = surfaceHeightAt(x, z) + .35 + seeded(index * 2 + 3) * 1.55;
      const tiltX = (seeded(index * 2 + 4) - .5) * .35;
      const tiltZ = (seeded(index * 2 + 5) - .5) * .35;
      dummy.position.set(x, y, z); dummy.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), surfaceNormalAt(x, z)); dummy.rotateX(tiltX); dummy.rotateY(seeded(index * 2 + 6) * Math.PI); dummy.rotateZ(tiltZ); dummy.updateMatrix(); facetMesh.setMatrixAt(index, dummy.matrix);
    }
    facetMesh.instanceMatrix.needsUpdate = true; facetMesh.frustumCulled = false; sceneObjects.add(facetMesh);
  }
  const originY = surfaceHeightAt(0, 0) + 1;
  addNorthMarker(surfaceHeightAt(0, 0));
  const sun = zenithDirection(firstNumber(value("sza"), 30), firstNumber(value("saa"), 0), 45, originY);
  addMesh(new THREE.SphereGeometry(.35, 16, 12), makeMaterial(0xedb866, .45), sun); addLine([0, originY, 0], sun, 0xedb866, .55);
  const rows = simulationRows.length ? selectedRows().slice(0, 28) : [];
  rows.forEach((row) => { const end = zenithDirection(Number(row.vza), Number(row.vaa), 30, originY); addLine([0, originY, 0], end, 0xed8e7d, .32); });
  fitSceneView();
  $("scene-loading").classList.add("is-hidden");
}

function initScene() {
  const container = $("scene-container");
  scene = new THREE.Scene(); scene.background = new THREE.Color(0xe8f1f5);
  camera = new THREE.PerspectiveCamera(42, container.clientWidth / container.clientHeight, .1, 250); camera.position.set(70, 55, 75);
  renderer = new THREE.WebGLRenderer({ antialias: true }); renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2)); renderer.setSize(container.clientWidth, container.clientHeight); container.appendChild(renderer.domElement);
  controls = new OrbitControls(camera, renderer.domElement); controls.enableDamping = true; controls.target.set(0, 2, 0);
  scene.add(new THREE.HemisphereLight(0xbed5ce, 0x1c2623, 1.8)); const light = new THREE.DirectionalLight(0xffe1b0, 2.2); light.position.set(8, 13, 6); scene.add(light);
  sceneObjects = new THREE.Group(); scene.add(sceneObjects); updateScene();
  const resize = () => { const width = container.clientWidth; const height = container.clientHeight; camera.aspect = width / height; camera.updateProjectionMatrix(); renderer.setSize(width, height); if (simulationRows.length) renderResults(); };
  window.addEventListener("resize", resize);
  const animate = () => { requestAnimationFrame(animate); controls.update(); renderer.render(scene, camera); }; animate();
}

function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url; link.download = filename; link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

async function exportInput() {
  const validation = validateInputs();
  if (validation) {
    notify(validation.message, true); validation.input?.focus(); return;
  }
  const button = $("export-button");
  button.disabled = true;
  try {
    const response = await fetch("/api/export", { method: "POST",
      headers: { "Content-Type": "application/json" }, body: JSON.stringify(currentPayload()) });
    if (!response.ok) {
      const payload = await response.json();
      throw new Error(payload.error || t("Input export failed"));
    }
    downloadBlob(await response.blob(), "tirt_input.zip");
    notify(t("Saved input ZIP with spectra and observation files"));
  } catch (error) {
    notify(error.message, true);
  } finally {
    button.disabled = false;
  }
}

async function shutdownTiRT() {
  const button = document.getElementById('exit-button');
  button.disabled = true;
  try {
    await fetch("/api/shutdown", { method: "POST" });
    document.body.textContent = t("TiRT has exited.");
  } catch (error) {
    button.disabled = false;
    notify(error.message || t("TiRT could not exit"), true);
  }
}

function bindEvents() {
  $("run-button").addEventListener("click", runSimulation); $("save-output-button").addEventListener("click", exportOutput); $("export-button").addEventListener("click", exportInput); document.getElementById('exit-button').addEventListener('click', shutdownTiRT); $("reset-button").addEventListener("click", () => window.location.reload());
  document.querySelectorAll('input[name="additional-source"]').forEach((input) => input.addEventListener("change", updateGeometryFields)); $("band-select").addEventListener("change", renderResults); $("scene-select").addEventListener("change", () => { selectedSceneKey = $("scene-select").value; renderResults(); });
  ["show-polar", "show-parallel", "show-perpendicular"].forEach((id) => $(id).addEventListener("change", updateChartVisibility));
  ["surface-model", "vegetation-model"].forEach((id) => $(id).addEventListener("change", () => { updateConditionalFields(); updateSceneReadout(); updateScene(); }));
  const sceneInputIds = ["sza", "saa", "principal-vza-step", "principal-vza-max", "manual-geometry", "lai", "hspot", "slope-angle", "slope-aspect", "terrain-height", "terrain-radius", "terrain-density", "row-width", "row-blank", "row-height", "row-azimuth", "crown-rad-a", "crown-rad-b", "crown-height", "crown-density", "building-density", "building-height", "building-length", "building-width", "building-azimuth"];
  const refreshScene = () => { updateSceneReadout(); updateScene(); };
  sceneInputIds.forEach((id) => { const input = $(id); input?.addEventListener("input", refreshScene); input?.addEventListener("change", refreshScene); });
  ["include-principal", "include-hemisphere"].forEach((id) => $(id).addEventListener("change", () => { updateSceneReadout(); updateScene(); }));
  $("geometry-upload").addEventListener("change", async (event) => {
    const file = event.target.files[0];
    const status = $("geometry-file-status");
    if (!file) {
      status.textContent = "•";
      status.className = "file-status no-file";
      status.setAttribute("aria-label", t("No file selected"));
      uploadedGeometryText = "";
      return;
    }
    const contents = await file.text();
    if (event.target.files[0] !== file) return;
    uploadedGeometryText = contents;
    status.textContent = "✓";
    status.className = "file-status has-file";
    status.setAttribute("aria-label", t("Selected {name}", { name: file.name }));
    notify(t("Loaded {name}", { name: file.name }));
  });
  const tableStatus = (text, className, label) => {
    const status = $("band-table-file-status");
    status.textContent = text; status.className = `file-status ${className}`;
    status.setAttribute("aria-label", label);
  };
  $("band-table-file").addEventListener("change", async (event) => {
    const file = event.target.files[0];
    if (!file) return;
    const contents = await file.text();
    if (event.target.files[0] !== file) return;
    try {
      parseBandTable(contents).forEach((row, index) => { $(BAND_TABLE_IDS[index]).value = row.join(";"); });
    } catch (error) {
      event.target.value = "";
      tableStatus("•", "no-file", t("No file selected"));
      notify(`${file.name}: ${error.message}`, true);
      return;
    }
    tableStatus("✓", "has-file", t("Selected {name}", { name: file.name }));
    updateEmissivityReadout();
    notify(t("Loaded {name}", { name: file.name }));
  });
  BAND_TABLE_IDS.forEach((id) => $(id).addEventListener("input", () => {
    $("band-table-file").value = "";
    tableStatus("•", "no-file", t("Using entered values"));
    updateEmissivityReadout();
  }));
}

async function boot() {
  bindEvents(); updateGeometryFields(); updateConditionalFields(); updateSceneReadout(); updateChartVisibility(); updateEmissivityReadout();
  try {
    THREE = await import("three"); ({ OrbitControls } = await import("three/addons/controls/OrbitControls.js")); initScene();
  } catch (error) {
    console.error(error);
    $("scene-loading").textContent = t("Three.js could not be loaded: {message}", { message: error.message }); notify(t("Three.js load failed"), true); setStatus("Scene unavailable", "error");
  }
}

boot();
