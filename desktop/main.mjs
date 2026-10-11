import { app, BrowserWindow, dialog, screen, shell } from 'electron'
import { spawn, spawnSync } from 'node:child_process'
import { existsSync } from 'node:fs'
import { mkdir, readFile, rm } from 'node:fs/promises'
import { dirname, join } from 'node:path'
import { fileURLToPath } from 'node:url'

const sourceRoot = dirname(dirname(fileURLToPath(import.meta.url)))
const smokeTest = process.argv.includes('--smoke-test')
let mainWindow = null
let backendProcess = null
let backendError = ''
let readyFile = null
let isQuitting = false
let backendCleaned = false

if (!app.requestSingleInstanceLock()) {
  app.quit()
} else {
  app.on('second-instance', () => {
    if (!mainWindow) return
    if (mainWindow.isMinimized()) mainWindow.restore()
    mainWindow.show()
    mainWindow.focus()
  })
  app.whenReady().then(createWindow).catch((error) => {
    const message = error?.stack || String(error)
    console.error(message)
    if (!smokeTest && !isQuitting) dialog.showErrorBox('TiRT 启动失败', message)
    requestQuit()
  })
  app.on('before-quit', () => {
    isQuitting = true
    stopBackend()
  })
  app.on('will-quit', stopBackend)
  app.on('window-all-closed', () => app.quit())
}

function resourceRoot() {
  return app.isPackaged ? process.resourcesPath : sourceRoot
}

function backendSpec() {
  const root = resourceRoot()
  const packagedBackend = join(root, 'backend', process.platform === 'win32' ? 'TiRT-service.exe' : 'TiRT-service')
  if (existsSync(packagedBackend)) return { command: packagedBackend, args: [] }

  const localBackend = join(sourceRoot, 'dist', process.platform === 'win32' ? 'TiRT-service.exe' : 'TiRT-service')
  if (existsSync(localBackend)) return { command: localBackend, args: [] }

  const python = process.env.TIRT_PYTHON || (process.platform === 'win32' ? 'python.exe' : 'python3')
  return { command: python, args: [join(sourceRoot, 'gui', 'launcher.py')] }
}

async function startBackend() {
  const tempRoot = join(app.getPath('temp'), 'tirt-desktop')
  await mkdir(tempRoot, { recursive: true })
  readyFile = join(tempRoot, `ready-${process.pid}.json`)
  await rm(readyFile, { force: true })

  const spec = backendSpec()
  backendProcess = spawn(spec.command, [
    ...spec.args,
    '--host', '127.0.0.1',
    '--port', '0',
    '--server-only',
    '--ready-file', readyFile,
  ], {
    cwd: app.isPackaged ? app.getPath('userData') : sourceRoot,
    windowsHide: true,
    stdio: ['ignore', 'ignore', 'pipe'],
    env: { ...process.env, MPLBACKEND: 'Agg' },
  })
  backendProcess.stderr?.on('data', (data) => { backendError = (backendError + String(data)).slice(-8000) })
  backendProcess.once('error', (error) => {
    backendError += `\n${error.message}`
    backendProcess = null
  })
  backendProcess.once('exit', (code, signal) => {
    backendProcess = null
    if (!isQuitting) {
      // Exit code 0 is the in-app Exit button (/api/shutdown), not a failure.
      if (code !== 0 && !smokeTest) {
        const detail = backendError.trim() || `服务已退出（code=${code}, signal=${signal || 'none'}）`
        dialog.showErrorBox('TiRT 服务已停止', detail)
      }
      requestQuit()
    }
  })

  const deadline = Date.now() + 30000
  let lastError = ''
  while (Date.now() < deadline) {
    if (!backendProcess) break
    try {
      const ready = JSON.parse(await readFile(readyFile, 'utf8'))
      const url = `http://${ready.host}:${ready.port}`
      try {
        const response = await fetch(`${url}/api/health`)
        if (response.ok) return url
        lastError = `HTTP ${response.status}`
      } catch (error) {
        lastError = String(error)
      }
    } catch (error) {
      lastError = String(error)
    }
    await new Promise((resolve) => setTimeout(resolve, 150))
  }
  throw new Error(backendError.trim() || `TiRT 本地服务启动超时：${lastError}`)
}

function stopBackend() {
  if (backendCleaned) return
  backendCleaned = true
  const child = backendProcess
  backendProcess = null
  if (readyFile) rm(readyFile, { force: true }).catch(() => {})
  if (!child || child.killed) return
  if (process.platform === 'win32') {
    spawnSync('taskkill.exe', ['/PID', String(child.pid), '/T', '/F'], {
      windowsHide: true,
      stdio: 'ignore',
    })
  } else {
    child.kill('SIGTERM')
  }
}

function requestQuit() {
  isQuitting = true
  stopBackend()
  app.quit()
}

async function createWindow() {
  const url = await startBackend()
  const workArea = screen.getPrimaryDisplay().workArea
  const aspectRatio = 16 / 9
  // Fill the complete usable screen height; the window itself remains 16:9.
  const margin = 0
  const availableWidth = Math.max(1, workArea.width - margin * 2)
  const availableHeight = Math.max(1, workArea.height - margin * 2)
  // Use the largest 16:9 window that fits the usable screen area.
  // Do not cap the width: on a 1920x1080 display the old 1440px cap
  // produced an unnecessarily small 1440x810 window.
  const width = Math.floor(Math.min(availableWidth, availableHeight * aspectRatio))
  const height = Math.floor(width / aspectRatio)
  const minWidth = Math.min(960, width)
  const minHeight = Math.floor(minWidth / aspectRatio)

  mainWindow = new BrowserWindow({
    width,
    height,
    x: workArea.x + Math.round((workArea.width - width) / 2),
    y: workArea.y + Math.round((workArea.height - height) / 2),
    minWidth,
    minHeight,
    resizable: true,
    maximizable: true,
    title: 'TiRT',
    icon: join(resourceRoot(), 'tirt.ico'),
    show: false,
    backgroundColor: '#0b1118',
    autoHideMenuBar: true,
    webPreferences: {
      contextIsolation: true,
      nodeIntegration: false,
      sandbox: true,
    },
  })
  mainWindow.setAspectRatio(aspectRatio)

  mainWindow.webContents.setWindowOpenHandler(({ url: target }) => {
    if (/^https?:/i.test(target)) shell.openExternal(target)
    return { action: 'deny' }
  })
  mainWindow.webContents.on('will-navigate', (event, target) => {
    if (new URL(target).origin !== new URL(url).origin) {
      event.preventDefault()
      if (/^https?:/i.test(target)) shell.openExternal(target)
    }
  })
  mainWindow.once('ready-to-show', () => {
    if (smokeTest) {
      console.log(`TIRT_SMOKE_TEST_OK ${url}`)
      requestQuit()
    } else {
      mainWindow.show()
    }
  })
  mainWindow.on('closed', () => {
    mainWindow = null
    requestQuit()
  })
  await mainWindow.loadURL(url)
}
