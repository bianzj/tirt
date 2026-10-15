// Headless regression checks for the actual GUI input parser and payload builder.
const assert = require("node:assert/strict"); // Explicit CommonJS in ESM repositories.
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");

const root = path.resolve(__dirname, "..");
const html = fs.readFileSync(path.join(root, "gui/index.html"), "utf8");
const source = fs.readFileSync(path.join(root, "gui/app.js"), "utf8").replace(/boot\(\);\s*$/, "");
const elements = new Map();
const hidden = new Set();
const parents = Object.fromEntries(["leaf", "terrain", "roof", "wall", "street"].map((name) =>
  [`emissivity-${name}`, `${name}-component-row`]));
function element(id) {
  if (!elements.has(id)) elements.set(id, { id, value: "", checked: false, disabled: false,
    textContent: "", files: [], handlers: {},
    classList: { toggle(name, enabled) { if (name === "is-hidden") enabled ? hidden.add(id) : hidden.delete(id); },
      add() {}, remove() {} },
    closest(selector) { return selector === ".is-hidden" && hidden.has(parents[id] || id) ? {} : null; },
    addEventListener(name, handler) { this.handlers[name] = handler; },
    setAttribute() {}, getAttribute() { return id; }, focus() {},
    querySelector() { return element(`${id}-child`); },
  });
  return elements.get(id);
}
for (const tag of html.matchAll(/<(input|select)[^>]*id="([^"]+)"[^>]*>/g)) {
  const item = element(tag[2]);
  item.value = (tag[0].match(/value="([^"]*)"/) || [])[1] || "";
  item.checked = /\schecked(?:\s|>)/.test(tag[0]);
  if (tag[1] === "select") {
    const options = html.slice(tag.index + tag[0].length).split("</select>")[0];
    item.value = (options.match(/<option[^>]*value="([^"]+)"[^>]*selected/) ||
      options.match(/<option[^>]*value="([^"]+)"/) || [])[1] || "";
    item.selectedOptions = [{ textContent: item.value }];
  }
}
let additionalSource = "file";
const context = { console, setTimeout: () => 0, clearTimeout() {},
  window: { location: { protocol: "http:" } },
  document: { documentElement: {}, getElementById: element,
    querySelector: () => ({ value: additionalSource }), querySelectorAll: () => [] } };
vm.createContext(context);
vm.runInContext(source, context);
const execute = (code) => vm.runInContext(code, context);
const plain = (value) => JSON.parse(JSON.stringify(value));
execute("updateConditionalFields(); bindEvents();");

assert.deepEqual(plain(execute('parseWavelengthFile("# um\\r\\n8\\r\\n\\r\\n10.5\\r\\n12\\r\\n")')), [8, 10.5, 12]);
assert.deepEqual(plain(execute('parseWavelengthFile("8\\n0\\n")')), [8, 0]);
assert.equal(execute("validateInputs()"), null);

element("emissivity-soil").value = "0.8;0.9";
assert.equal(execute("validateInputs().input.id"), "emissivity-soil");
assert.match(execute("validateInputs().message"), /3 values/);
assert.equal(execute("currentOverrides().emissivity_soil"), "0.8;0.9");
element("emissivity-soil").value = "0.8";
assert.equal(execute("validateInputs().input.id"), "emissivity-soil");
element("emissivity-soil").value = "0.8;bad;0.9";
assert.match(execute("validateInputs().message"), /between 0 and 1/);
execute('uploadedFiles.soil_spectrum="8 0.8\\n10.5 0.9\\n12 0.95\\n";');
assert.equal(execute("validateInputs()"), null); // Uploaded spectra supersede direct values.
execute('uploadedGeometryText="40 180\\n";');
additionalSource = "manual";
assert.equal(execute("currentPayload().geometry_text"), "");
assert.equal(execute("currentPayload().overrides.additional_geometry_source"), "manual");
assert.equal(execute("currentPayload().overrides.hemisphere_geometry_file"), "data/directions/hemisphere_directions.txt");
element("manual-geometry").value = "10/90; invalid";
assert.equal(execute("validateInputs().input.id"), "manual-geometry");
element("manual-geometry").value = "10/90";
additionalSource = "file";
assert.equal(execute("currentPayload().geometry_text"), "40 180\n");

const vertical = plain(execute("zenithDirection(0, 90, 30, 2)"));
assert.deepEqual(vertical, [0, 32, 0]);
const horizon = plain(execute("zenithDirection(90, 90, 30, 2)"));
assert.ok(Math.abs(horizon[0] - 30) < 1e-12 && Math.abs(horizon[1] - 2) < 1e-12 && Math.abs(horizon[2]) < 1e-12);
const oblique = plain(execute("zenithDirection(30, 0, 1)"));
assert.ok(Math.abs(oblique[1] - Math.cos(Math.PI / 6)) < 1e-12);

(async () => {
  const file = { name: "bands.txt", text: async () => "8\n10.5\n12\n" };
  const upload = element("wavelength-file"); upload.files = [file];
  await upload.handlers.change({ target: upload });
  assert.equal(element("wavelengths").value, "8;10.5;12");
  assert.equal(execute("validateInputs()"), null);
  element("wavelengths").handlers.input();
  assert.equal(execute("currentPayload().files.wavelengths"), "");
  element("emissivity-soil").value = "0.8;0.9;0.95";
  element("emissivity-soil").handlers.input();
  assert.equal(execute("currentPayload().files.soil_spectrum"), "");
  assert.equal(execute("validateInputs()"), null);
  console.log("PASS GUI input parsing, spectrum counts, active uploads, and zenith preview");
})().catch((error) => { console.error(error); process.exitCode = 1; });
