// Actual Edge smoke test, using the browser's DevTools protocol and Node built-ins.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const net = require('node:net');
const { spawn } = require('node:child_process');
const [base, folder] = process.argv.slice(2);
assert(base && folder, 'Usage: node verify_tirt_browser.cjs http://127.0.0.1:PORT OUTPUT');
fs.mkdirSync(folder, { recursive: true });
const delay = ms => new Promise(resolve => setTimeout(resolve, ms));
const profile = fs.mkdtempSync(path.join(os.tmpdir(), 'tirt-browser-'));
let browser, ws;
async function main() {
  const port = await new Promise(resolve => {
    const probe = net.createServer().listen(0, '127.0.0.1', () => {
      const number = probe.address().port; probe.close(() => resolve(number));
    });
  });
  browser = spawn('C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe', [
    '--headless=new', '--no-first-run', '--no-default-browser-check', '--enable-unsafe-swiftshader',
    `--remote-debugging-port=${port}`, `--user-data-dir=${profile}`, 'about:blank'],
    { windowsHide: true, stdio: 'ignore' });
  let target;
  for (let i = 0; i < 100; ++i) {
    try {
      target = await (await fetch(`http://127.0.0.1:${port}/json/new?${encodeURIComponent(base)}`, { method: 'PUT' })).json();
      break;
    } catch { await delay(200); }
  }
  assert(target?.webSocketDebuggerUrl, 'Edge failed to start');
  ws = new WebSocket(target.webSocketDebuggerUrl);
  await new Promise((resolve, reject) => { ws.onopen = resolve; ws.onerror = reject; });
  let sequence = 0;
  const pending = new Map(), errors = [];
  ws.onmessage = event => {
    const msg = JSON.parse(event.data);
    if (msg.method === 'Runtime.exceptionThrown') errors.push(msg.params.exceptionDetails);
    if (pending.has(msg.id)) {
      const { resolve, reject, timeout } = pending.get(msg.id); clearTimeout(timeout); pending.delete(msg.id);
      msg.error ? reject(new Error(JSON.stringify(msg.error))) : resolve(msg.result);
    }
  };
  const send = (method, params = {}) => new Promise((resolve, reject) => {
    const id = ++sequence;
    const timeout = setTimeout(() => reject(new Error(`CDP timeout: ${method}`)), 30000);
    pending.set(id, { resolve, reject, timeout }); ws.send(JSON.stringify({ id, method, params }));
  });
  const evaluate = async expression => {
    const r = await send('Runtime.evaluate', { expression, returnByValue: true, awaitPromise: true });
    assert(!r.exceptionDetails, JSON.stringify(r.exceptionDetails)); return r.result.value;
  };
  const until = async (expression, seconds = 30) => {
    for (let i = 0; i < seconds * 2; ++i) {
      if (await evaluate(expression)) return;
      await delay(500);
    }
    throw new Error(`Browser condition timed out: ${expression}`);
  };
  await send('Runtime.enable'); await send('Page.enable');
  await send('Emulation.setDeviceMetricsOverride', { width: 1680, height: 1050, deviceScaleFactor: 1, mobile: false });
  await until(`document.querySelector('.model-static') && document.querySelector('#run-button')`);
  const staticColor = await evaluate(`getComputedStyle(document.querySelector('.model-static')).backgroundColor`);
  assert.equal(staticColor, 'rgb(15, 114, 199)');
  await until(`document.querySelector('#scene-container canvas')`);
  fs.writeFileSync(path.join(folder, 'static.png'), Buffer.from((await send('Page.captureScreenshot')).data, 'base64'));
  await evaluate(`document.querySelector('a.dynamic-link').click()`);
  await until(`document.querySelector('.model-dynamic') && typeof currentOverrides === 'function'`);
  const dynamicColor = await evaluate(`getComputedStyle(document.querySelector('.model-dynamic')).backgroundColor`);
  assert.equal(dynamicColor, 'rgb(213, 59, 54)');
  // Run the untouched default day before testing edited angles and bands.
  await evaluate(`document.getElementById('run-button').click()`);
  await until(`!document.getElementById('run-button').disabled`, 180);
  const defaultDay = await evaluate(`({steps:window.temporalState.seriesRows.length,
    observations:window.temporalState.observationRows.length,
    toast:document.getElementById('toast').textContent})`);
  assert.equal(defaultDay.steps, 48, JSON.stringify(defaultDay));
  assert.equal(defaultDay.observations, 48);
  console.log('Default dynamic day passed:', JSON.stringify(defaultDay));
  await evaluate(`document.getElementById('time-range').value='214.48-214.54';
    document.getElementById('observation-angles').value='15/77;30/160';
    document.getElementById('wavelengths').value='0.67;0.86;10.5';
    document.getElementById('run-button').click();`);
  await until(`!document.getElementById('run-button').disabled`, 150);
  const result = await evaluate(`({title:document.title, status:document.getElementById('status-pill').textContent,
    toast:document.getElementById('toast').textContent,
    steps:window.temporalState.seriesRows.length, observations:window.temporalState.observationRows.length,
    pairs:[...new Set(window.temporalState.observationRows.map(r=>r.vza+'/'+r.vaa))],
    temperatures:window.temporalState.seriesRows.map(r=>r.Tv_C)})`);
  assert(result.steps > 0, JSON.stringify(result));
  assert.equal(result.observations, result.steps * 2);
  assert.deepEqual(result.pairs, ['15/77', '30/160']);
  assert(result.temperatures.every(Number.isFinite));
  await until(`document.querySelector('#scene-container canvas')`);
  fs.writeFileSync(path.join(folder, 'dynamic.png'), Buffer.from((await send('Page.captureScreenshot')).data, 'base64'));
  await evaluate(`document.querySelector('a.static-link').click()`);
  await until(`document.querySelector('.model-static')`);
  assert.equal(errors.length, 0, JSON.stringify(errors));
  const report = { ok: true, staticColor, dynamicColor, defaultDay, dynamic: result, uncaughtErrors: errors };
  fs.writeFileSync(path.join(folder, 'browser_verification.json'), JSON.stringify(report, null, 2));
  console.log(JSON.stringify(report, null, 2));
  await send('Browser.close').catch(() => {});
}
main().catch(error => { console.error(error); process.exitCode = 1; }).finally(() => {
  ws?.close();
  if (browser && !browser.killed) browser.kill();
  // Keep the temporary profile until Edge's child processes have released it.
});
