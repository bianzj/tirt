// Build a HiSTREAM eRaytracing case from the transverse-dune DEM used by StreamStudio.
// Usage: node build_ridge.mjs --out <runDir> [--res 256] [--samples 32] [--sza 30] [--saa 270] [--vza 0] [--vaa 0]
import { mkdirSync, readFileSync, writeFileSync } from 'node:fs'
import { join, resolve } from 'node:path'
import { pathToFileURL } from 'node:url'

const STREAMSIM = process.env.STREAMSIM_ROOT || 'C:\\work\\streamstudio'
const { createDefaultProject, stringifyProject } = await import(
  pathToFileURL(join(STREAMSIM, 'src', 'renderer', 'src', 'project-schema.js')).href)

const argv = process.argv.slice(2)
const arg = (name, fallback) => { const i = argv.indexOf(`--${name}`); return i >= 0 ? argv[i + 1] : fallback }
const out = resolve(arg('out'))
const resolution = Number(arg('res', 256))
const samples = Number(arg('samples', 32))
const depth = Number(arg('depth', 4))
const SZA = Number(arg('sza', 30)), SAA = Number(arg('saa', 270))
const VZA = Number(arg('vza', 0)), VAA = Number(arg('vaa', 0))
const DEM = arg('dem', 'C:\\work\\streamstudio\\assets\\dem-library\\dune\\transverse_dunes_100m.asc')
const BAND_NM = 10500
const EMIS_TERRAIN = 0.950
const TEMP = [318.15, 303.15]

const text = readFileSync(DEM, 'utf8').split(/\r?\n/)
const header = {}
let firstData = 0
for (let i = 0; i < text.length; i += 1) {
  const line = text[i].trim()
  if (!line) continue
  const m = line.match(/^([A-Za-z_]+)\s+(-?[\d.]+)$/)
  if (m && firstData === 0) { header[m[1].toUpperCase()] = Number(m[2]); continue }
  firstData = i; break
}
const ncols = Math.round(header.NCOLS), nrows = Math.round(header.NROWS), cell = header.CELLSIZE
const z = []
for (let r = 0; r < nrows; r += 1) {
  z.push(text[firstData + r].trim().split(/\s+/).map(Number))
}
if (z.length !== nrows || z[0].length !== ncols) throw new Error(`DEM shape ${z.length}x${z[0]?.length}, expected ${nrows}x${ncols}`)
const L = ncols * cell
const sceneHeight = Math.ceil(z.flat().reduce((a, b) => Math.max(a, b), 0) + 2)

const inputDir = join(out, 'input'), outputDir = join(out, 'output')
mkdirSync(join(inputDir, 'models'), { recursive: true })
mkdirSync(join(inputDir, 'positions'), { recursive: true })
mkdirSync(outputDir, { recursive: true })

const fmt = (v) => Number(v.toFixed(5))
// Engine axes: X north, Y up, Z east. ASC rows run north -> south, columns west -> east.
const xOf = (row) => (nrows - 1 - row) * cell - L / 2
const zOf = (col) => col * cell - L / 2

const lines = ['# transverse dunes heightfield from StreamStudio DEM', 'o heightfield']
let vertices = 0
for (let r = 0; r < nrows - 1; r += 1) {
  for (let c = 0; c < ncols - 1; c += 1) {
    const x0 = xOf(r + 1), x1 = xOf(r), z0 = zOf(c), z1 = zOf(c + 1)
    const y00 = z[r][c], y01 = z[r][c + 1], y11 = z[r + 1][c + 1], y10 = z[r + 1][c]
    // (x0,z0) -> (x0,z1) -> (x1,z1) -> (x1,z0) gives an upward normal.
    const points = [[x0, y00, z0], [x0, y01, z1], [x1, y11, z1], [x1, y10, z0]]
    for (const p of points) lines.push(`v ${fmt(p[0])} ${fmt(p[1])} ${fmt(p[2])}`)
    lines.push(`f ${vertices + 1} ${vertices + 2} ${vertices + 3}`, `f ${vertices + 1} ${vertices + 3} ${vertices + 4}`)
    vertices += 4
  }
}
// Perimeter skirt down to y = 0, so rays leaving the heightfield still hit terrain.
const ring = []
for (let c = 0; c < ncols - 1; c += 1) {
  ring.push([xOf(0), zOf(c), xOf(0), zOf(c + 1)])
  ring.push([xOf(nrows - 1), zOf(c + 1), xOf(nrows - 1), zOf(c)])
}
for (let r = 0; r < nrows - 1; r += 1) {
  ring.push([xOf(r), zOf(0), xOf(r + 1), zOf(0)])
  ring.push([xOf(r + 1), zOf(ncols - 1), xOf(r), zOf(ncols - 1)])
}
for (const [x0, z0, x1, z1] of ring) {
  const y0 = 0, y1 = 0
  const points = [[x0, y0, z0], [x0, -0.05, z0], [x1, -0.05, z1], [x1, y1, z1]]
  for (const p of points) lines.push(`v ${fmt(p[0])} ${fmt(p[1])} ${fmt(p[2])}`)
  lines.push(`f ${vertices + 1} ${vertices + 2} ${vertices + 3}`, `f ${vertices + 1} ${vertices + 3} ${vertices + 4}`)
  vertices += 4
}
const heightfieldPath = join(inputDir, 'models', 'heightfield.obj')
writeFileSync(heightfieldPath, `${lines.join('\n')}\n`, 'utf8')

const positionPath = join(inputDir, 'positions', 'heightfield_position.txt')
writeFileSync(positionPath, `${fmt(L / 2)} ${fmt(L / 2)} 0 1 0\n`, 'utf8')

const spectra = [{ name: 'terrain', label: 'terrain', model: 'custom', reflectance: '0.1', transmittance: '0', refTir: Number((1 - EMIS_TERRAIN).toFixed(6)), tauTir: 0, params: {} }]
const thermals = [{ name: 'terrain_temperature', label: 'terrain_temperature', sunlitTemperature: TEMP[0], shadedTemperature: TEMP[1] }]
const canopies = [{ name: 'rigid_body', label: 'rigid', structureType: 'rigid', lai: 0, density: 0, hc: 2, G: 0.5, LIDFa: 0, LIDFb: 0, hspot: 0, leafwidth: 0 }]
const items = [
  { name: 'heightfield', type: 'Soil', fileName: heightfieldPath, positionFile: positionPath,
    meshes: [{ name: 'heightfield', spectralName: 'terrain', thermalName: 'terrain_temperature', canopyName: 'rigid_body', materialName: 'soilset' }] }
]

const project = createDefaultProject({ mode: 'eRaytracing', name: 'tirt-ridge-dunes', sceneX: L, sceneY: L, sceneHeight, voxelSize: 0.5 })
const c = project.configuration
c.outDir = outputDir
c.scene.background = { spectralName: 'terrain', thermalName: 'terrain_temperature', materialType: 'Soil', materialName: 'soilset', heterogeneityEnabled: false, angularEffectStrength: 0 }
c.scene.terrain = false
c.light = { zenith: SZA, azimuth: SAA, direct: 1, diffuse: 0, skyTemperature: 3 }
Object.assign(c.sensor, { x: resolution, y: resolution, bands: String(BAND_NM), projection: 'parallel', vza: VZA, vaa: VAA,
  continuousBands: false, principalPlane: false, hemisphere: false, cruiseEnabled: false, image: true, process: false,
  radiationProcess: false, energyProcess: false, temperature: false, albedo: false })
Object.assign(c.control, { depth, samples, periodicTraversalCount: 5, heterogeneousVoxel: false, skyboxEnabled: false, radiationSolver: 'traditional' })
c.atmosphere.enabled = false
c.fluid.enabled = false
c.spectra = spectra
c.thermals = thermals
c.canopies = canopies
c.objects = { count: items.length, names: items.map((i) => i.name), items }

writeFileSync(join(inputDir, 'project.json'), stringifyProject(project), 'utf8')
writeFileSync(join(out, 'case.json'), JSON.stringify({ dem: DEM, L, ncols, nrows, cell, sceneHeight, sza: SZA, saa: SAA, vza: VZA, vaa: VAA,
  emissivity: EMIS_TERRAIN, temperatures: TEMP, vertices, triangles: vertices / 2, resolution, samples, depth }, null, 2), 'utf8')
console.log(JSON.stringify({ project: join(inputDir, 'project.json'), L, sceneHeight, vertices }))
