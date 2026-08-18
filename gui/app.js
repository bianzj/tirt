let THREE;
let OrbitControls;
let scene;
let camera;
let renderer;
let controls;
let sceneObjects;
let uploadedGeometryText = "";
const uploadedFiles = { wavelengths: "", leaf_spectrum: "", soil_spectrum: "", terrain_spectrum: "", roof_spectrum: "", wall_spectrum: "", street_spectrum: "" };
let simulationRows = [];
let selectedWavelength = null;
let selectedScenarioId = "1";
let toastTimer;

const $ = (id) => document.getElementById(id);
if (window.location.protocol === "file:") {
  window.location.replace("http://127.0.0.1:8765/");
}
const root = document.documentElement;
const colors = { cyan: "#57d3c6", amber: "#edb866", coral: "#ed8e7d", muted: "#687774", grid: "#293735", ink: "#e9f0ee" };
const SCENE_SIZE = 100;
const HALF_SCENE = SCENE_SIZE / 2;
const NORTH_AZIMUTH = 0;

function numberList(value) {
  return String(value || "").split(/[;,\s]+/).map(Number).filter(Number.isFinite);
}

function firstNumber(raw, fallback = 0) {
  const values = numberList(raw);
  return values.length ? values[0] : fallback;
}

function value(id) { return $(id).value; }
function selectedAdditionalSource() { return document.querySelector('input[name="additional-source"]:checked')?.value || "file"; }
function manualGeometryPairs() {
  return String($("manual-geometry")?.value || "").split(/[;\n]+/).map((entry) => {
    const pair = entry.trim().split(/\s*\/\s*/).map(Number);
    return [pair[0], pair[1]];
  }).filter((pair) => Number.isFinite(pair[0]) && Number.isFinite(pair[1]));
}
function setStatus(text, kind = "idle") {
  const pill = $("status-pill");
  pill.className = `status-pill ${kind}`;
  pill.querySelector("span").textContent = text;
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
  const wavelengths = numberList(value("wavelengths"));
  const bandValues = (id, fallback) => {
    const values = numberList(value(id));
    if (!wavelengths.length) return value(id);
    if (values.length === wavelengths.length) return value(id);
    if (values.length === 1) return wavelengths.map(() => values[0]).join(";");
    if (wavelengths.length === 1 && values.length) return String(values[0]);
    if (values.length) return wavelengths.map((_, index) => values[index % values.length]).join(";");
    return fallback;
  };
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
    wavelengths: value("wavelengths"),
    emissivity_leaf: bandValues("emissivity-leaf", "0.985"),
    emissivity_soil: bandValues("emissivity-soil", "0.955"),
    hapke_enable: value("hapke-enable"),
    hapke_K: value("hapke-k"),
    emissivity_terrain: bandValues("emissivity-terrain", "0.975"),
    emissivity_roof: bandValues("emissivity-roof", "0.95"),
    emissivity_wall: bandValues("emissivity-wall", "0.92"),
    emissivity_street: bandValues("emissivity-street", "0.955"),
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
  overrides.additional_geometry_source = additionalSource;
  if (additionalSource === "file") {
    if (includeHemisphere) {
      overrides.additional_geometry_file = "data/directions/hemisphere_directions.txt";
    } else if (uploadedGeometryText) {
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
  $("leaf-component-row").classList.toggle("is-hidden", vegetation === "bare");
  $("bare-soil-anisotropy-fields").classList.toggle("is-hidden", vegetation !== "bare");
  $("terrain-component-row").classList.toggle("is-hidden", surface !== "terrain");
  $("roof-component-row").classList.toggle("is-hidden", surface !== "urban");
  $("wall-component-row").classList.toggle("is-hidden", surface !== "urban");
  $("street-component-row").classList.toggle("is-hidden", surface !== "urban");
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
  if ($("include-principal").checked) parts.push(`${principalCount * 4} principal`);
  if ($("include-hemisphere").checked) parts.push("hemisphere");
  if (selectedAdditionalSource() === "file" && !$("include-hemisphere").checked) parts.push("file");
  if (additionalCount) parts.push(`${additionalCount} manual`);
  let count = parts.join(" + ") || "No directions";
  if (simulationRows.length) count = `${new Set(simulationRows.map((row) => `${row.scenario_id || 1}|${row.vza}|${row.vaa}|${row.sza}|${row.saa}`)).size} directions`;
  $("observation-readout").textContent = count;
}

async function runSimulation() {
  setStatus("Running", "busy");
  $("run-button").disabled = true;
  try {
    const response = await fetch("/api/run", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ overrides: currentOverrides(), geometry_text: uploadedGeometryText, files: uploadedFiles }),
    });
    const payload = await response.json();
    if (!response.ok || !payload.ok) throw new Error(payload.error || "Simulation failed");
    simulationRows = payload.rows;
    $("save-output-button").disabled = false;
    updateBandOptions();
    updateScenarioOptions();
    renderResults();
    updateSceneReadout();
    updateScene();
    setStatus("Complete", "idle");
    notify(`Simulation complete · ${payload.count} output records`);
  } catch (error) {
    setStatus("Error", "error");
    notify(error.message, true);
  } finally {
    $("run-button").disabled = false;
  }
}

function updateBandOptions() {
  const bands = [...new Set(simulationRows.map((row) => Number(row.wavelength_um)))].sort((a, b) => a - b);
  const select = $("band-select");
  select.innerHTML = "";
  bands.forEach((band) => {
    const option = document.createElement("option");
    option.value = String(band);
    option.textContent = `${fmt(band, 1)} µm`;
    select.appendChild(option);
  });
  if (bands.length) {
    selectedWavelength = bands.includes(10.5) ? 10.5 : bands[0];
    select.value = String(selectedWavelength);
  }
}

function updateScenarioOptions() {
  const scenarios = [...new Map(simulationRows.map((row) => [String(row.scenario_id || 1), row.scenario_label || "default"])).entries()];
  const select = $("scenario-select");
  select.innerHTML = "";
  scenarios.forEach(([id, label]) => {
    const option = document.createElement("option");
    option.value = id;
    option.textContent = `S${String(id).padStart(2, "0")} · ${label}`;
    select.appendChild(option);
  });
  if (!scenarios.some(([id]) => id === selectedScenarioId)) selectedScenarioId = scenarios[0]?.[0] || "1";
  select.value = selectedScenarioId;
}

function selectedScenarioRows() {
  selectedScenarioId = $("scenario-select")?.value || selectedScenarioId;
  return simulationRows.filter((row) => String(row.scenario_id || 1) === selectedScenarioId);
}

function selectedRows() {
  const wavelength = Number($("band-select").value);
  selectedWavelength = wavelength;
  return selectedScenarioRows().filter((row) => Math.abs(Number(row.wavelength_um) - wavelength) < 1e-7);
}

function outputValue(row, quantity = $("quantity").value) { return Number(row[quantity]); }
function quantityValue(row) { return outputValue(row); }
function quantityLabel(quantity = $("quantity").value) {
  return { brightness_temperature_C: "Temperature (°C)", brightness_temperature_K: "Temperature (K)", radiance: "Radiance" }[quantity] || quantity;
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
  const suffix = quantity === "brightness_temperature_C" ? " °C" : quantity === "brightness_temperature_K" ? " K" : "";
  return `${formatOutputNumber(number)}${suffix}`;
}

function renderResults() {
  if (!simulationRows.length) return;
  $("empty-state").classList.add("is-hidden");
  $("results-content").classList.remove("is-hidden");
  const rows = selectedRows();
  const values = rows.map(quantityValue).filter(Number.isFinite);
  const radiances = rows.map((row) => Number(row.radiance)).filter(Number.isFinite);
  const mean = values.reduce((sum, item) => sum + item, 0) / Math.max(values.length, 1);
  const span = values.length ? Math.max(...values) - Math.min(...values) : NaN;
  const meanRadiance = radiances.reduce((sum, item) => sum + item, 0) / Math.max(radiances.length, 1);
  $("mean-bt").closest(".metric").querySelector("span").textContent = `Mean ${quantityLabel()}`;
  $("mean-bt").textContent = formatOutputValue(mean, $("quantity").value);
  $("range-bt").textContent = fmtSignificant(span);
  $("mean-radiance").textContent = fmtSignificant(meanRadiance);
  $("observation-count").textContent = String(rows.length);
  renderCharts(rows);
  updateChartVisibility();
}

function observationTableRows(rows, quantity) {
  const bands = [...new Set(rows.map((row) => Number(row.wavelength_um)).filter(Number.isFinite))].sort((a, b) => a - b);
  const observations = new Map();
  rows.forEach((row) => {
    const scenarioId = Number(row.scenario_id) || 1; const scenarioLabel = row.scenario_label || "default";
    const vza = Number(row.vza); const vaa = Number(row.vaa); const raa = Number(row.raa);
    const sza = Number(row.sza); const saa = Number(row.saa); const wavelength = Number(row.wavelength_um); const output = outputValue(row, quantity);
    if (![vza, vaa, raa, wavelength, output].every(Number.isFinite)) return;
    const key = [scenarioId, vza, vaa, sza, saa].map((value) => Number.isFinite(value) ? value : "").join("|");
    let observation = observations.get(key);
    if (!observation) {
      observation = { scenarioId, scenarioLabel, vza, vaa, raa, values: new Map() };
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
    notify("Run a simulation first", true);
    return;
  }
  const quantity = $("quantity").value;
  const { bands, observationRows } = observationTableRows(simulationRows, quantity);
  const valueName = quantity === "radiance" ? "radiance" : quantity === "brightness_temperature_K" ? "temperature_K" : "temperature_C";
  const fields = ["index", "scenario_id", "scenario_label", "vza", "vaa", "raa", ...bands.map((band) => `${fmt(band, 1)}um_${valueName}`)];
  const lines = [
    fields.join(","),
    ...observationRows.map((observation, index) => [
      index, observation.scenarioId, observation.scenarioLabel, observation.vza, observation.vaa, observation.raa,
      ...bands.map((band) => formatOutputNumber(observation.values.get(band))),
    ].map(csvEscape).join(",")),
  ];
  const blob = new Blob([`${lines.join("\n")}\n`], { type: "text/csv;charset=utf-8" });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = "observations.csv";
  link.click();
  setTimeout(() => URL.revokeObjectURL(link.href), 0);
  notify(`Saved observations.csv · ${observationRows.length} directions · all cases · ${quantityLabel(quantity)}`);
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
  [0, 90, 180, 270].forEach((angle) => { const rad = angle * Math.PI / 180; ctx.strokeStyle = "#22302e"; ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(cx + Math.sin(rad) * radius, cy - Math.cos(rad) * radius); ctx.stroke(); ctx.fillStyle = colors.muted; ctx.fillText(`${angle}°`, cx + Math.sin(rad) * (radius + 16), cy - Math.cos(rad) * (radius + 16)); });
  $("polar-colorbar-min").textContent = fmt(min);
  $("polar-colorbar-max").textContent = fmt(max);
  ctx.fillStyle = colors.muted; ctx.textAlign = "left"; ctx.fillText(`VZA radius · ${fmt(selectedWavelength, 1)} µm`, 8, height - 10);
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
  const margin = Math.max((yMax - yMin) * .15, .5);
  return { lo: yMin - margin, hi: yMax + margin };
}

function drawLineChart(rows, canvasId, plane, color, scale) {
  const { ctx, width, height } = canvasContext($(canvasId));
  const points = planeRows(rows, plane); const pad = { left: 38, right: 12, top: 26, bottom: 29 };
  const xMax = Math.max(30, ...points.map((point) => Math.abs(point.x)), 80); const { lo, hi } = scale; const plotW = width - pad.left - pad.right, plotH = height - pad.top - pad.bottom;
  ctx.font = "10px system-ui"; ctx.strokeStyle = colors.grid; ctx.fillStyle = colors.muted; ctx.lineWidth = 1;
  [-xMax, 0, xMax].forEach((x) => { const px = pad.left + (x + xMax) / (xMax * 2) * plotW; ctx.beginPath(); ctx.moveTo(px, pad.top); ctx.lineTo(px, pad.top + plotH); ctx.stroke(); ctx.textAlign = "center"; ctx.fillText(`${fmt(x, 0)}°`, px, height - 11); });
  [lo, (lo + hi) / 2, hi].forEach((y) => { const py = pad.top + (hi - y) / (hi - lo) * plotH; ctx.beginPath(); ctx.moveTo(pad.left, py); ctx.lineTo(pad.left + plotW, py); ctx.stroke(); ctx.textAlign = "right"; ctx.fillText(fmt(y), pad.left - 7, py); });
  if (!points.length) { ctx.fillStyle = colors.muted; ctx.textAlign = "center"; ctx.fillText("No matching directions", width / 2, height / 2); return; }
  const px = (x) => pad.left + (x + xMax) / (xMax * 2) * plotW; const py = (y) => pad.top + (hi - y) / (hi - lo) * plotH;
  // Keep the two sides separate so the perpendicular plane cannot look truncated
  // when one side is rendered after the other at the nadir point.
  const branches = [points.filter((point) => point.negative), points.filter((point) => !point.negative)];
  branches.forEach((branch) => {
    ctx.strokeStyle = color; ctx.lineWidth = 2; ctx.beginPath();
    branch.forEach((point, index) => { if (index) ctx.lineTo(px(point.x), py(quantityValue(point.row))); else ctx.moveTo(px(point.x), py(quantityValue(point.row))); });
    ctx.stroke();
  });
  points.forEach((point) => { ctx.fillStyle = color; ctx.beginPath(); ctx.arc(px(point.x), py(quantityValue(point.row)), 3.5, 0, Math.PI * 2); ctx.fill(); });
  ctx.fillStyle = colors.muted; ctx.textAlign = "center";
  const sideText = plane === "perpendicular" ? "RAA 90° / 270° · both sides" : "RAA 0° / 180° · both sides";
  ctx.fillText(`${sideText} · ${points.length} directions`, width / 2, 12);
}

function renderCharts(rows) {
  drawPolar(rows);
  const scale = lineChartRange(rows);
  drawLineChart(rows, "parallel-chart", "parallel", colors.amber, scale);
  drawLineChart(rows, "perpendicular-chart", "perpendicular", colors.coral, scale);
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
  const sza = firstNumber(value("sza"), 30) * Math.PI / 180; const saa = firstNumber(value("saa"), 0) * Math.PI / 180; const sunDistance = 45; const sun = [Math.sin(saa) * Math.cos(sza) * sunDistance, originY + Math.sin(sza) * sunDistance, Math.cos(saa) * Math.cos(sza) * sunDistance];
  addMesh(new THREE.SphereGeometry(.35, 16, 12), makeMaterial(0xedb866, .45), sun); addLine([0, originY, 0], sun, 0xedb866, .55);
  const rows = simulationRows.length ? selectedRows().slice(0, 28) : [];
  rows.forEach((row) => { const vza = Number(row.vza) * Math.PI / 180; const vaa = Number(row.vaa) * Math.PI / 180; const viewDistance = 30; const end = [Math.sin(vaa) * Math.cos(vza) * viewDistance, originY + Math.sin(vza) * viewDistance, Math.cos(vaa) * Math.cos(vza) * viewDistance]; addLine([0, originY, 0], end, 0xed8e7d, .32); });
  fitSceneView();
  $("scene-loading").classList.add("is-hidden");
}

function initScene() {
  const container = $("scene-container");
  scene = new THREE.Scene(); scene.background = new THREE.Color(0x0c1413);
  camera = new THREE.PerspectiveCamera(42, container.clientWidth / container.clientHeight, .1, 250); camera.position.set(70, 55, 75);
  renderer = new THREE.WebGLRenderer({ antialias: true }); renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2)); renderer.setSize(container.clientWidth, container.clientHeight); container.appendChild(renderer.domElement);
  controls = new OrbitControls(camera, renderer.domElement); controls.enableDamping = true; controls.target.set(0, 2, 0);
  scene.add(new THREE.HemisphereLight(0xbed5ce, 0x1c2623, 1.8)); const light = new THREE.DirectionalLight(0xffe1b0, 2.2); light.position.set(8, 13, 6); scene.add(light);
  sceneObjects = new THREE.Group(); scene.add(sceneObjects); updateScene();
  const resize = () => { const width = container.clientWidth; const height = container.clientHeight; camera.aspect = width / height; camera.updateProjectionMatrix(); renderer.setSize(width, height); if (simulationRows.length) renderResults(); };
  window.addEventListener("resize", resize);
  const animate = () => { requestAnimationFrame(animate); controls.update(); renderer.render(scene, camera); }; animate();
}

function csvCell(value) { const text = String(value ?? ""); return /[",\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text; }
function csvRow(section, key, val, description = "") { return [section, key, val, description].map(csvCell).join(","); }
function buildInputCsv() {
  const o = currentOverrides();
  const rows = [
    "# TiRT GUI exported input", ",,,", "[case],,,",
    csvRow("case", "case_name", "gui_case", "case name"),
    csvRow("case", "result_dir", "cases/gui_case", "result directory"),
    ",,,", "[surface],,,",
    csvRow("surface", "surface_model", o.surface_model), csvRow("surface", "vegetation_model", o.vegetation_model),
    ",,,", "[vegetation_structure],,,",
    csvRow("vegetation_structure", "lai", o.lai), csvRow("vegetation_structure", "hspot", o.hspot),
    csvRow("row_structure", "row_width", o.row_width), csvRow("row_structure", "row_blank", o.row_blank), csvRow("row_structure", "row_height", o.row_height), csvRow("row_structure", "row_azimuth", o.row_azimuth),
    csvRow("crown_structure", "crown_rad_a", o.crown_rad_a), csvRow("crown_structure", "crown_rad_b", o.crown_rad_b), csvRow("crown_structure", "crown_height", o.crown_height), csvRow("crown_structure", "crown_density", o.crown_density),
    ",,,", "[slope],,,", csvRow("slope", "slope_angle", o.slope_angle), csvRow("slope", "slope_aspect", o.slope_aspect),
    ",,,", "[terrain],,,", csvRow("terrain", "terrain_height", o.terrain_height), csvRow("terrain", "terrain_radius", o.terrain_radius), csvRow("terrain", "terrain_density", o.terrain_density),
    ",,,", "[urban],,,", csvRow("urban", "building_density", o.building_density), csvRow("urban", "building_height", o.building_height), csvRow("urban", "building_length", o.building_length), csvRow("urban", "building_width", o.building_width), csvRow("urban", "building_azimuth", o.building_azimuth),
    ",,,", "[spectral],,,", csvRow("spectral", "wavelengths", o.wavelengths), csvRow("spectral", "emissivity_leaf", o.emissivity_leaf), csvRow("spectral", "emissivity_soil", o.emissivity_soil), csvRow("spectral", "emissivity_terrain", o.emissivity_terrain),
    ",,,", "[soil_anisotropy],,,", csvRow("soil_anisotropy", "hapke_enable", o.hapke_enable), csvRow("soil_anisotropy", "hapke_K", o.hapke_K),
    csvRow("urban_spectrum", "emissivity_roof", o.emissivity_roof), csvRow("urban_spectrum", "emissivity_wall", o.emissivity_wall), csvRow("urban_spectrum", "emissivity_street", o.emissivity_street),
    ",,,", "[thermal],,,", csvRow("thermal", "leaf_sunlit_c", o.leaf_sunlit_c), csvRow("thermal", "leaf_shaded_c", o.leaf_shaded_c), csvRow("thermal", "soil_sunlit_c", o.soil_sunlit_c), csvRow("thermal", "soil_shaded_c", o.soil_shaded_c), csvRow("thermal", "terrain_sunlit_c", o.terrain_sunlit_c), csvRow("thermal", "terrain_shaded_c", o.terrain_shaded_c), csvRow("thermal", "roof_sunlit_c", o.roof_sunlit_c), csvRow("thermal", "roof_shaded_c", o.roof_shaded_c), csvRow("thermal", "wall_sunlit_c", o.wall_sunlit_c), csvRow("thermal", "wall_shaded_c", o.wall_shaded_c), csvRow("thermal", "street_sunlit_c", o.street_sunlit_c), csvRow("thermal", "street_shaded_c", o.street_shaded_c),
    ",,,", "[geometry],,,", csvRow("geometry", "geometry_mode", o.geometry_mode), csvRow("geometry", "geometry_view", o.geometry_view), csvRow("geometry", "include_principal", o.include_principal), csvRow("geometry", "include_hemisphere", o.include_hemisphere), csvRow("geometry", "sza", o.sza), csvRow("geometry", "saa", o.saa),
  ];
  rows.push(csvRow("geometry", "principal_vza_step", o.principal_vza_step), csvRow("geometry", "principal_vza_max", o.principal_vza_max), csvRow("geometry", "additional_geometry_source", o.additional_geometry_source));
  if (o.additional_geometry_source === "file") rows.push(csvRow("geometry", "additional_geometry_file", o.additional_geometry_file));
  if (o.additional_geometry_source === "manual") rows.push(csvRow("geometry", "observation_vza", o.observation_vza), csvRow("geometry", "observation_vaa", o.observation_vaa));
  return `${rows.join("\n")}\n`;
}
function exportInput() { const blob = new Blob([buildInputCsv()], { type: "text/csv;charset=utf-8" }); const link = document.createElement("a"); link.href = URL.createObjectURL(blob); link.download = "input.csv"; link.click(); URL.revokeObjectURL(link.href); notify("Saved input.csv"); }

function bindEvents() {
  $("run-button").addEventListener("click", runSimulation); $("save-output-button").addEventListener("click", exportOutput); $("export-button").addEventListener("click", exportInput); $("reset-button").addEventListener("click", () => window.location.reload());
  document.querySelectorAll('input[name="additional-source"]').forEach((input) => input.addEventListener("change", updateGeometryFields)); $("band-select").addEventListener("change", renderResults); $("quantity").addEventListener("change", renderResults); $("scenario-select").addEventListener("change", renderResults);
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
      status.setAttribute("aria-label", "No file selected");
      uploadedGeometryText = "";
      return;
    }
    uploadedGeometryText = await file.text();
    status.textContent = "✓";
    status.className = "file-status has-file";
    status.setAttribute("aria-label", `Selected ${file.name}`);
    notify(`Loaded ${file.name}`);
  });
  [["wavelength-file", "wavelengths"], ["leaf-spectrum-file", "leaf_spectrum"], ["soil-spectrum-file", "soil_spectrum"], ["terrain-spectrum-file", "terrain_spectrum"], ["roof-spectrum-file", "roof_spectrum"], ["wall-spectrum-file", "wall_spectrum"], ["street-spectrum-file", "street_spectrum"]].forEach(([id, key]) => {
    $(id)?.addEventListener("change", async (event) => {
      const file = event.target.files[0];
      const status = $(`${id}-status`);
      if (!file) return;
      uploadedFiles[key] = await file.text();
      if (status) {
        status.textContent = "✓";
        status.className = "file-status has-file";
        status.setAttribute("aria-label", `Selected ${file.name}`);
      }
      if (key === "wavelengths") {
        const fileWavelengths = uploadedFiles[key].split(/\r?\n/).map((line) => Number(line.trim().split(/[;,\s]+/)[0])).filter(Number.isFinite);
        if (fileWavelengths.length) $("wavelengths").value = fileWavelengths.join(";");
      }
      notify(`Loaded ${file.name}`);
    });
  });
}

async function boot() {
  bindEvents(); updateGeometryFields(); updateConditionalFields(); updateSceneReadout(); updateChartVisibility();
  try {
    THREE = await import("three"); ({ OrbitControls } = await import("three/addons/controls/OrbitControls.js")); initScene();
  } catch (error) {
    console.error(error);
    $("scene-loading").textContent = `Three.js could not be loaded: ${error.message}`; notify("Three.js load failed", true); setStatus("Scene unavailable", "error");
  }
}

boot();
