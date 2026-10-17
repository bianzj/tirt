const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const root = path.resolve(__dirname, '..');
const html = fs.readFileSync(path.join(root, 'gui/time.html'), 'utf8');
const source = fs.readFileSync(path.join(root, 'gui/time.js'), 'utf8')
  .replace(/\nbindEvents\(\); updateConditionalFields\(\);[\s\S]*$/, '\n');
const elements = new Map();
function element(id) {
  if (!elements.has(id)) elements.set(id, { value: '', handlers: {},
    addEventListener(name, handler) { this.handlers[name] = handler; },
    classList: { toggle() {}, add() {}, remove() {} }, setAttribute() {}, textContent: '' });
  return elements.get(id);
}
for (const tag of html.matchAll(/<(input|select)[^>]*id="([^"]+)"[^>]*>/g)) {
  const item = element(tag[2]);
  item.value = (tag[0].match(/value="([^"]*)"/) || [])[1] || '';
  if (tag[1] === 'select') {
    const options = html.slice(tag.index + tag[0].length).split('</select>')[0];
    item.value = (options.match(/<option[^>]*value="([^"]+)"[^>]*selected/) ||
      options.match(/<option[^>]*value="([^"]+)"/) || [])[1] || '';
  }
}
const context = { console, setTimeout() {}, clearTimeout() {},
  window: { location: { protocol: 'http:' }, addEventListener() {} },
  document: { getElementById: element } };
vm.createContext(context);
vm.runInContext(source, context);
const run = code => vm.runInContext(code, context);
assert.equal(run('currentOverrides().direction_source'), 'direct');
assert.equal(run('currentOverrides().observation_time_mode'), 'all');
assert.equal(run('currentOverrides().wavelengths'), '10.5');
element('observation-angles').value = '15/77;30/160';
assert.equal(run('currentOverrides().vza'), '15;30');
assert.equal(run('currentOverrides().vaa'), '77;160');
assert.equal(run('currentOverrides().geometry_pair_mode'), 'zip');
run('bindEvents()');
(async () => {
  for (const [id, key] of [['vegetation-spectrum-upload', 'vegetation_spectrum'],
                            ['soil-spectrum-upload', 'soil_spectrum']]) {
    await element(id).handlers.change({ target: { files: [{ name: 'sample.txt', text: async () => '860 .4 .4' }] } });
    assert.equal(run(`spectrumFiles.${key}`), '860 .4 .4');
  }
  for (const [file, identity, cls] of [['index.html', 'TIRT-Static', 'model-static'],
      ['time.html', 'TIRT-Dynamic', 'model-dynamic'], ['angle/index.html', 'TIRT-Static', 'model-static']]) {
    const page = fs.readFileSync(path.join(root, 'gui', file), 'utf8');
    assert(page.includes('aria-label="STREAM TIRT"'));
    assert(page.includes('<span class="logo-stream">STREAM</span>'));
    assert(page.includes('<span class="logo-tirt">TIRT</span>'));
    assert(page.includes(`<h1>${identity}`));
    assert(page.includes(cls));
    assert(page.includes('model_identity.css'));
  }
  assert(html.includes('href="index.html"'));
  console.log('Dynamic payload, uploads, navigation and model identities passed.');
})().catch(error => { console.error(error); process.exitCode = 1; });
