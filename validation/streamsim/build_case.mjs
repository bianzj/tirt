// Builds a HiSTREAM project.json + OBJ/position inputs for one TiRT benchmark case.
// Usage: node build_case.mjs --test T2 --mode eRaytracing --out <runDir> [--samples 16] [--leaf 0.15] [--seed 1]
import { mkdirSync, writeFileSync } from 'node:fs'
import { join, resolve } from 'node:path'
import { pathToFileURL } from 'node:url'

const STREAMSIM = process.env.STREAMSIM_ROOT || 'C:\\work\\streamsim'
const { createDefaultProject, stringifyProject } = await import(
  pathToFileURL(join(STREAMSIM, 'src', 'renderer', 'src', 'project-schema.js')).href)

const argv = process.argv.slice(2)
const arg = (name, fallback) => {
  const index = argv.indexOf(`--${name}`)
  return index >= 0 ? argv[index + 1] : fallback
}
const test = arg('test')
const mode = arg('mode', 'eRaytracing')
const out = resolve(arg('out'))
const samples = Number(arg('samples', 16))
const seed = Number(arg('seed', 1))
const depth = Number(arg('depth', 4))
const resolution = Number(arg('res', 512))
if (!['T1', 'T2', 'T3', 'T4', 'T5az0', 'T5az90', 'TV', 'UV'].includes(test)) throw new Error('--test must be T1..T4, T5az0, T5az90, TV or UV')
if (!['eRaytracing', 'eVoxelRT'].includes(mode)) throw new Error('--mode must be eRaytracing or eVoxelRT')

// TiRT input.csv-matched common parameters.
const SZA = 30, SAA = 0, BAND_NM = 10500
const LAI = 1.5
const HSPOT = 0.15
const EMIS = { soil: 0.955, street: 0.955, leaf: 0.985, roof: 0.950, wall: 0.930, terrain: 0.950 }
const TEMP = {
  soil: [318.15, 303.15], leaf: [306.15, 303.15],
  roof: [318.15, 303.15], wall: [318.15, 303.15], street: [318.15, 303.15]
}
// Decomposition overrides: --blackbody sets every eps=1; --temps "soilSun,soilShade,leafSun,leafShade".
if (argv.includes('--blackbody')) for (const key of Object.keys(EMIS)) EMIS[key] = 1
if (arg('temps')) {
  const [soilSun, soilShade, leafSun, leafShade] = arg('temps').split(',').map(Number)
  TEMP.soil = [soilSun, soilShade]; TEMP.leaf = [leafSun, leafShade]
}
// TiRT terrain+vegetation: every soil surface (ground and cones) has terrain emissivity and soil temperatures.
TEMP.terrain = TEMP.soil
// TiRT urban+vegetation: vegetation is one component at the mean leaf temperature.
if (test === 'UV') TEMP.leaf = [(TEMP.leaf[0] + TEMP.leaf[1]) / 2, (TEMP.leaf[0] + TEMP.leaf[1]) / 2]
// Scattering decomposition: --cold "leaf,..." stops those classes emitting (3 K, B = 0 at 10.5 um);
// --black "terrain,..." stops them reflecting (emissivity 1).
for (const name of (arg('cold') || '').split(',').filter(Boolean)) {
  if (!(name in TEMP)) throw new Error(`--cold: unknown class ${name}`)
  TEMP[name] = [3, 3]
}
for (const name of (arg('black') || '').split(',').filter(Boolean)) {
  if (!(name in EMIS)) throw new Error(`--black: unknown class ${name}`)
  EMIS[name] = 1
}

// Deterministic PRNG (mulberry32) so every rerun builds identical scenes.
function rng(seedValue) {
  let a = seedValue >>> 0
  return () => {
    a = (a + 0x6D2B79F5) >>> 0
    let t = a
    t = Math.imul(t ^ (t >>> 15), t | 1)
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61)
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296
  }
}
const random = rng(seed * 7919 + test.charCodeAt(1))

const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
const norm = (a) => { const l = Math.hypot(...a); return a.map((v) => v / l) }
const fmt = (v) => Number(v.toFixed(5))

function sphericalNormal() {
  const y = random() * 2 - 1
  const phi = random() * 2 * Math.PI
  const r = Math.sqrt(Math.max(0, 1 - y * y))
  return [r * Math.cos(phi), y, r * Math.sin(phi)]
}

// Square leaf of side s with a uniformly random (spherical-LAD) normal.
function squareLeaf(center, s) {
  const n = sphericalNormal()
  const ref = Math.abs(n[1]) < 0.9 ? [0, 1, 0] : [1, 0, 0]
  const u = norm(cross(ref, n))
  const v = norm(cross(n, u))
  const h = s / 2
  return [[-1, -1], [1, -1], [1, 1], [-1, 1]].map(([a, b]) =>
    [center[0] + h * (a * u[0] + b * v[0]), center[1] + h * (a * u[1] + b * v[1]), center[2] + h * (a * u[2] + b * v[2])])
}

// --leaf-shape disk: circular leaf of diameter d (24-gon with the disk's area), spherical-LAD normal.
const LEAF_SHAPE = arg('leaf-shape', 'square')
if (!['square', 'disk'].includes(LEAF_SHAPE)) throw new Error('--leaf-shape must be square or disk')
const DISK_SIDES = 24
function diskLeaf(center, d) {
  const n = sphericalNormal()
  const ref = Math.abs(n[1]) < 0.9 ? [0, 1, 0] : [1, 0, 0]
  const u = norm(cross(ref, n))
  const v = norm(cross(n, u))
  const r = d / 2 * Math.sqrt(2 * Math.PI / (DISK_SIDES * Math.sin(2 * Math.PI / DISK_SIDES)))
  return Array.from({ length: DISK_SIDES }, (_, k) => {
    const a = r * Math.cos(2 * Math.PI * k / DISK_SIDES), b = r * Math.sin(2 * Math.PI * k / DISK_SIDES)
    return [center[0] + a * u[0] + b * v[0], center[1] + a * u[1] + b * v[1], center[2] + a * u[2] + b * v[2]]
  })
}
const makeLeaf = (center, size) => (LEAF_SHAPE === 'disk' ? diskLeaf(center, size) : squareLeaf(center, size))
const leafArea = (size) => (LEAF_SHAPE === 'disk' ? Math.PI * size * size / 4 : size * size)
const leafReach = (size) => (LEAF_SHAPE === 'disk' ? size / 2 : size / Math.SQRT2)

class Obj {
  constructor() { this.lines = ['# TiRT benchmark generated OBJ'], this.vertexCount = 0, this.triangles = 0 }
  object(name) { this.lines.push(`o ${name}`) }
  quad(points) {
    for (const p of points) this.lines.push(`v ${fmt(p[0])} ${fmt(p[1])} ${fmt(p[2])}`)
    const b = this.vertexCount
    this.lines.push(`f ${b + 1} ${b + 2} ${b + 3}`, `f ${b + 1} ${b + 3} ${b + 4}`)
    this.vertexCount += 4
    this.triangles += 2
  }
  polygon(points) {
    for (const p of points) this.lines.push(`v ${fmt(p[0])} ${fmt(p[1])} ${fmt(p[2])}`)
    const b = this.vertexCount
    for (let k = 1; k + 1 < points.length; k += 1) this.lines.push(`f ${b + 1} ${b + k + 1} ${b + k + 2}`)
    this.vertexCount += points.length
    this.triangles += points.length - 2
  }
  save(path) { writeFileSync(path, `${this.lines.join('\n')}\n`, 'utf8'); return path }
}

// Sutherland-Hodgman clip of a planar polygon against axis-aligned half-space p[axis] >= value (sign 1) or <= value (sign -1).
function clipPolygon(points, axis, value, sign) {
  const result = []
  for (let k = 0; k < points.length; k += 1) {
    const a = points[k], b = points[(k + 1) % points.length]
    const da = sign * (a[axis] - value), db = sign * (b[axis] - value)
    if (da >= 0) result.push(a)
    if ((da >= 0) !== (db >= 0)) {
      const t = da / (da - db)
      result.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])])
    }
  }
  return result
}

function polygonArea(points) {
  let area = [0, 0, 0]
  for (let k = 1; k + 1 < points.length; k += 1) {
    const c = cross(sub(points[k], points[0]), sub(points[k + 1], points[0]))
    area = [area[0] + c[0], area[1] + c[1], area[2] + c[2]]
  }
  return Math.hypot(...area) / 2
}

const inputDir = join(out, 'input'), outputDir = join(out, 'output')
mkdirSync(join(inputDir, 'models'), { recursive: true })
mkdirSync(join(inputDir, 'positions'), { recursive: true })
mkdirSync(outputDir, { recursive: true })

function writePositions(name, rows) {
  const path = join(inputDir, 'positions', `${name}_position.txt`)
  writeFileSync(path, `${rows.map((r) => r.map(fmt).join(' ')).join('\n')}\n`, 'utf8')
  return path
}

// Raytracing: explicit 3x3 periodic copies so shadow/scattering rays (which the
// engine does not wrap) see the periodic neighbours. Voxel: only copies that
// overlap the domain; the voxeliser clips everything outside.
function periodicRows(points, L, reach) {
  const rows = []
  for (const p of points) for (const i of [-1, 0, 1]) for (const j of [-1, 0, 1]) {
    const x = p[0] + i * L, z = p[1] + j * L
    if (mode === 'eVoxelRT' && (x < -reach || x > L + reach || z < -reach || z > L + reach)) continue
    rows.push([x, z, p[2], 1, p[3] ?? 0])
  }
  return rows
}

// Ground skirt around the domain (raytracing only) so secondary rays leaving
// the domain still find ground with the background properties.
function skirtObj(L) {
  const obj = new Obj()
  obj.object('skirt')
  const h = L / 2, H = 1.5 * L
  const rects = [[-H, -H, H, -h], [-H, h, H, H], [-H, -h, -h, h], [h, -h, H, h]]
  for (const [x0, z0, x1, z1] of rects) obj.quad([[x0, 0, z0], [x0, 0, z1], [x1, 0, z1], [x1, 0, z0]])
  return obj.save(join(inputDir, 'models', 'skirt.obj'))
}

// Voxel filler: one small triangle at every voxel centre of a shape so the
// voxeliser activates exactly the quantised envelope (density is then set
// from the canopy parameter, independent of leaf triangles).
function voxelFillerObj(name, step, inside, bounds) {
  const obj = new Obj()
  obj.object(name)
  let count = 0
  const t = step * 0.45
  for (let x = bounds[0][0]; x < bounds[0][1] - 1e-9; x += step)
    for (let y = bounds[1][0]; y < bounds[1][1] - 1e-9; y += step)
      for (let z = bounds[2][0]; z < bounds[2][1] - 1e-9; z += step) {
        const c = [x + step / 2, y + step / 2, z + step / 2]
        if (!inside(c)) continue
        obj.quad([[c[0] - t, c[1], c[2] - t], [c[0] - t, c[1], c[2] + t], [c[0] + t, c[1], c[2] + t], [c[0] + t, c[1], c[2] - t]])
        count += 1
      }
  return { path: obj.save(join(inputDir, 'models', `${name}.obj`)), voxels: count }
}

const meta = { test, mode, seed, samples, depth, resolution, sza: SZA, saa: SAA, bandNm: BAND_NM }
const spectra = [], thermals = [], canopies = [], items = []
const addSpectrum = (name, emissivity) => spectra.push({ name, label: name, model: 'custom', reflectance: '0.1', transmittance: '0', refTir: Number((1 - emissivity).toFixed(6)), tauTir: 0, params: {} })
const addThermal = (name, [sunlit, shaded]) => thermals.push({ name, label: name, sunlitTemperature: sunlit, shadedTemperature: shaded })

// Spherical crowns of radius R filled with square leaves (centres [x, z, y, rot], y = crown-centre height),
// split over 4 leaf-layout variants. Leaves lie entirely inside the sphere.
function crownItems(centres, R, leaf, lad, binding) {
  const variants = 4
  const perCrown = Math.round(lad * 4 / 3 * Math.PI * R ** 3 / leafArea(leaf))
  for (let variant = 0; variant < variants; variant += 1) {
    const obj = new Obj(); const name = `crown${variant}`; obj.object(name)
    for (let k = 0; k < perCrown; k += 1) {
      let quad
      do {
        let c
        do c = [(random() * 2 - 1) * R, (random() * 2 - 1) * R, (random() * 2 - 1) * R]
        while (Math.hypot(...c) > R)
        quad = makeLeaf(c, leaf)
      } while (quad.some((p) => Math.hypot(...p) > R))
      obj.polygon(quad)
    }
    const path = obj.save(join(inputDir, 'models', `${name}.obj`))
    const mine = centres.filter((_, index) => index % variants === variant)
    items.push({ name, type: 'Vegetation', fileName: path, positionFile: writePositions(name, periodicRows(mine, L, R)), meshes: [{ name, ...binding }] })
  }
  return { leafSize: leaf, leavesPerCrown: perCrown, crowns: centres.length, lad, leafAreaIndex: centres.length * perCrown * leaf * leaf / (L * L), variants }
}

// Random sequential addition of non-overlapping objects in a periodic L x L tile.
function hardCore(count, overlaps, maxTries = 200000) {
  const placed = []
  for (let tries = 0; placed.length < count && tries < maxTries; tries += 1) {
    const p = [random() * L, random() * L]
    if (!placed.some((q) => overlaps(periodicDelta(p[0] - q[0]), periodicDelta(p[1] - q[1])))) placed.push(p)
  }
  if (placed.length < count) throw new Error(`only ${placed.length} of ${count} objects fit`)
  return placed
}
const periodicDelta = (d) => d - L * Math.round(d / L)

let L, sceneHeight, voxel = Number(arg('voxel', 0.5))
const groundKind = test === 'T4' || test === 'UV' ? 'street' : test === 'TV' ? 'terrain' : 'soil'
addSpectrum(groundKind, EMIS[groundKind]); addThermal(`${groundKind}_temperature`, TEMP[groundKind])

if (test === 'T1') {
  L = Number(arg('L', 20)); sceneHeight = 2
}

if (test === 'T2') {
  L = Number(arg('L', 20)); sceneHeight = 2
  const H = 1.0
  const leaf = Number(arg('leaf', HSPOT * H))  // hspot = leaf size / canopy height
  addSpectrum('leaf', EMIS.leaf); addThermal('leaf_temperature', TEMP.leaf)
  canopies.push({ name: 'leaf_canopy', label: 'leaf', structureType: 'canopy', lai: LAI, density: LAI / H, hc: H, G: 0.5, LIDFa: -0.35, LIDFb: -0.15, hspot: leaf, leafwidth: leaf })
  const binding = { spectralName: 'leaf', thermalName: 'leaf_temperature', canopyName: 'leaf_canopy', materialName: 'leaf_c3' }
  if (mode === 'eRaytracing') {
    const obj = new Obj(); obj.object('leaves')
    const count = Math.round(LAI * L * L / leafArea(leaf))
    for (let k = 0; k < count; k += 1) {
      let quad
      do quad = makeLeaf([(random() - 0.5) * L, random() * H, (random() - 0.5) * L], leaf)
      while (quad.some((p) => p[1] < 0.001))
      obj.polygon(quad)
    }
    const path = obj.save(join(inputDir, 'models', 'leaves.obj'))
    items.push({ name: 'leaves', type: 'Vegetation', fileName: path, positionFile: writePositions('leaves', periodicRows([[L / 2, L / 2, 0]], L, 0)), meshes: [{ name: 'leaves', ...binding }] })
    Object.assign(meta, { leafShape: LEAF_SHAPE, leafSize: leaf, leafCount: count, leafAreaIndex: count * leafArea(leaf) / (L * L), canopyHeight: H })
  } else {
    const filler = voxelFillerObj('slab', voxel, () => true, [[-L / 2, L / 2], [0, H], [-L / 2, L / 2]])
    items.push({ name: 'slab', type: 'Vegetation', fileName: filler.path, positionFile: writePositions('slab', [[L / 2, L / 2, 0, 1, 0]]), meshes: [{ name: 'slab', ...binding }] })
    Object.assign(meta, { voxel, voxels: filler.voxels, lad: LAI / H, hspot: leaf, canopyHeight: H })
  }
}

if (test === 'T3') {
  L = Number(arg('L', 100)); sceneHeight = 10
  const R = 3, CZ = 6, DENSITY = 0.02
  const crownCount = Math.round(DENSITY * L * L)
  const crownVolume = 4 / 3 * Math.PI * R ** 3
  const lad = LAI / (DENSITY * crownVolume)
  const leaf = Number(arg('leaf', HSPOT * 2 * R))  // leaf size / crown depth (2b)
  addSpectrum('leaf', EMIS.leaf); addThermal('leaf_temperature', TEMP.leaf)
  const centres = Array.from({ length: crownCount }, () => [random() * L, random() * L, CZ, random() * 360])
  const binding = { spectralName: 'leaf', thermalName: 'leaf_temperature', canopyName: 'crown_canopy', materialName: 'leaf_c3' }
  if (mode === 'eRaytracing') {
    canopies.push({ name: 'crown_canopy', label: 'crown', structureType: 'canopy', lai: LAI, density: lad, hc: 2 * R, G: 0.5, LIDFa: -0.35, LIDFb: -0.15, hspot: leaf, leafwidth: leaf })
    Object.assign(meta, crownItems(centres, R, leaf, lad, binding))
  } else {
    const half = Math.ceil(R / voxel) * voxel
    const filler = voxelFillerObj('crown', voxel, (c) => Math.hypot(...c) <= R, [[-half, half], [-half, half], [-half, half]])
    // Preserve total leaf area: quantised sphere volume differs from 4/3 pi R^3,
    // and overlapping crowns share voxels (union), so the runner may pass the
    // LAD measured from the engine's inserted-voxel count via --lad.
    const quantisedLad = Number(arg('lad', LAI * L * L / (crownCount * filler.voxels * voxel ** 3)))
    canopies.push({ name: 'crown_canopy', label: 'crown', structureType: 'canopy', lai: LAI, density: quantisedLad, hc: 2 * R, G: 0.5, LIDFa: -0.35, LIDFb: -0.15, hspot: leaf, leafwidth: leaf })
    // Voxel grid placement uses floor(position/step): snap so the quantised envelope is centred.
    const snapped = centres.map(([x, z, y, rot]) => [Math.round(x / voxel) * voxel, Math.round(z / voxel) * voxel, y, 0])
    items.push({ name: 'crown', type: 'Vegetation', fileName: filler.path, positionFile: writePositions('crown', periodicRows(snapped, L, R)), meshes: [{ name: 'crown', ...binding }] })
    Object.assign(meta, { voxel, voxelsPerCrown: filler.voxels, crowns: crownCount, lad, quantisedLad, hspot: leaf })
  }
}

// Boxes of D (x) x W (z) x HB with separate roof/wall meshes, one instance per centre [x, z, y, rot].
function buildingItem(centres, W, D, HB) {
  addSpectrum('roof', EMIS.roof); addThermal('roof_temperature', TEMP.roof)
  addSpectrum('wall', EMIS.wall); addThermal('wall_temperature', TEMP.wall)
  const obj = new Obj()
  const x0 = -D / 2, x1 = D / 2, z0 = -W / 2, z1 = W / 2
  obj.object('roof')
  obj.quad([[x0, HB, z0], [x0, HB, z1], [x1, HB, z1], [x1, HB, z0]])
  obj.object('wall')
  // Outward normals matter: for walls exactly parallel to the sun azimuth the
  // engine's sunlit test offsets the shadow ray along -normal (into the box).
  obj.quad([[x1, 0, z0], [x1, HB, z0], [x1, HB, z1], [x1, 0, z1]])  // north face (+X)
  obj.quad([[x0, 0, z0], [x0, 0, z1], [x0, HB, z1], [x0, HB, z0]])  // south face (-X)
  obj.quad([[x0, 0, z1], [x1, 0, z1], [x1, HB, z1], [x0, HB, z1]])  // east face (+Z)
  obj.quad([[x0, 0, z0], [x0, HB, z0], [x1, HB, z0], [x1, 0, z0]])  // west face (-Z)
  const path = obj.save(join(inputDir, 'models', 'building.obj'))
  items.push({
    name: 'building', type: 'Building', fileName: path,
    positionFile: writePositions('building', periodicRows(centres, L, Math.hypot(W, D) / 2)),
    meshes: [
      { name: 'roof', spectralName: 'roof', thermalName: 'roof_temperature', canopyName: 'rigid_body', materialName: 'building_surface' },
      { name: 'wall', spectralName: 'wall', thermalName: 'wall_temperature', canopyName: 'rigid_body', materialName: 'building_surface' }
    ]
  })
  canopies.push({ name: 'rigid_body', label: 'rigid', structureType: 'rigid', lai: 0, density: 0, hc: 2, G: 0.5, LIDFa: 0, LIDFb: 0, hspot: 0, leafwidth: 0 })
}

if (test === 'T4') {
  L = Number(arg('L', 200)); sceneHeight = 22
  const W = 10, D = 10, HB = 20, DENSITY = 0.0015
  const count = Math.round(DENSITY * L * L)
  const centres = Array.from({ length: count }, () => [random() * L, random() * L, 0, 0])
  buildingItem(centres, W, D, HB)
  Object.assign(meta, { buildings: count, planDensity: count / (L * L), buildingSize: [D, W, HB] })
}

if (test === 'UV' || test === 'TV') {
  if (mode !== 'eRaytracing') throw new Error(`${test} is only built for eRaytracing`)
  addSpectrum('leaf', EMIS.leaf); addThermal('leaf_temperature', TEMP.leaf)
}

// Urban + vegetation (TiRT Urban_Veg): non-overlapping buildings, spherical crowns (R 3 m, centre 6 m) only on streets,
// crown density and LAI per unit street area, vegetation at one temperature.
if (test === 'UV') {
  L = Number(arg('L', 300)); sceneHeight = 22
  const W = 10, D = 10, HB = 20, DENSITY = 0.0015, R = 3, CZ = 6, STAND = 0.02
  const leaf = Number(arg('leaf', 0.15))
  const count = Math.round(DENSITY * L * L)
  const buildings = hardCore(count, (dx, dz) => Math.abs(dx) < D && Math.abs(dz) < W).map(([x, z]) => [x, z, 0, 0])
  buildingItem(buildings, W, D, HB)
  const streetArea = L * L - count * W * D
  const crownCount = Math.round(STAND * streetArea)
  const clear = (x, z) => buildings.every(([bx, bz]) => {
    const dx = Math.max(Math.abs(periodicDelta(x - bx)) - D / 2, 0), dz = Math.max(Math.abs(periodicDelta(z - bz)) - W / 2, 0)
    return Math.hypot(dx, dz) >= R
  })
  const centres = []
  while (centres.length < crownCount) {
    const x = random() * L, z = random() * L
    if (clear(x, z)) centres.push([x, z, CZ, random() * 360])
  }
  const lad = LAI / (STAND * 4 / 3 * Math.PI * R ** 3)
  canopies.push({ name: 'crown_canopy', label: 'crown', structureType: 'canopy', lai: LAI, density: lad, hc: 2 * R, G: 0.5, LIDFa: -0.35, LIDFb: -0.15, hspot: leaf, leafwidth: leaf })
  const binding = { spectralName: 'leaf', thermalName: 'leaf_temperature', canopyName: 'crown_canopy', materialName: 'leaf_c3' }
  Object.assign(meta, crownItems(centres, R, leaf, lad, binding),
    { buildings: count, planDensity: count / (L * L), buildingSize: [D, W, HB], streetArea, crownsPerStreetArea: crownCount / streetArea, leafAreaPerStreetArea: LAI })
}

// Compound terrain + vegetation (TiRT Terrain_Veg): non-overlapping cones on flat ground, Poisson spherical crowns
// per horizontal area everywhere, each standing along the local surface normal (crown centre 6 m out).
if (test === 'TV') {
  L = Number(arg('L', 600)); sceneHeight = 32
  const HC = 20, RC = 30, CONE_DENSITY = 1.5e-4, R = 3, CZ = 6, STAND = 0.02, SEGMENTS = 128
  const leaf = Number(arg('leaf', 0.15))
  const cones = hardCore(Math.round(CONE_DENSITY * L * L), (dx, dz) => Math.hypot(dx, dz) < 2 * RC)
  const cone = new Obj(); cone.object('cone')
  for (let k = 0; k < SEGMENTS; k += 1) {
    const a0 = 2 * Math.PI * k / SEGMENTS, a1 = 2 * Math.PI * (k + 1) / SEGMENTS
    cone.polygon([[RC * Math.cos(a0), 0, RC * Math.sin(a0)], [0, HC, 0], [RC * Math.cos(a1), 0, RC * Math.sin(a1)]])  // outward normal
  }
  items.push({ name: 'cone', type: 'Soil', fileName: cone.save(join(inputDir, 'models', 'cone.obj')),
    positionFile: writePositions('cone', periodicRows(cones.map(([x, z]) => [x, z, 0, 0]), L, RC)),
    meshes: [{ name: 'cone', spectralName: 'terrain', thermalName: 'terrain_temperature', canopyName: 'rigid_body', materialName: 'soilset' }] })
  const slant = Math.hypot(HC, RC), sinB = HC / slant, cosB = RC / slant
  const centres = Array.from({ length: Math.round(STAND * L * L) }, () => {
    const x = random() * L, z = random() * L
    let ground = [0, [0, 1, 0]]
    for (const [cx, cz] of cones) {
      const dx = periodicDelta(x - cx), dz = periodicDelta(z - cz), d = Math.hypot(dx, dz)
      if (d < RC) { ground = [HC * (1 - d / RC), d > 1e-9 ? [sinB * dx / d, cosB, sinB * dz / d] : [0, 1, 0]]; break }
    }
    const [y, n] = ground
    return [x + CZ * n[0], z + CZ * n[2], y + CZ * n[1], random() * 360]
  })
  const lad = LAI / (STAND * 4 / 3 * Math.PI * R ** 3)
  canopies.push({ name: 'crown_canopy', label: 'crown', structureType: 'canopy', lai: LAI, density: lad, hc: 2 * R, G: 0.5, LIDFa: -0.35, LIDFb: -0.15, hspot: leaf, leafwidth: leaf })
  const binding = { spectralName: 'leaf', thermalName: 'leaf_temperature', canopyName: 'crown_canopy', materialName: 'leaf_c3' }
  Object.assign(meta, crownItems(centres, R, leaf, lad, binding),
    { cones: cones.length, coneHeight: HC, coneRadius: RC, coneDensity: cones.length / (L * L), coneCover: cones.length * Math.PI * RC * RC / (L * L) })
}

if (test.startsWith('T5')) {
  // Rows: width 0.5 m, gap 0.5 m (spacing 1 m), height 1 m, in-row LAD 3 => scene LAI 1.5.
  const ROW_W = 0.5, ROW_GAP = 0.5, ROW_H = 1.0, SPACING = ROW_W + ROW_GAP
  const rowAz = test === 'T5az90' ? 90 : 0
  L = Number(arg('L', 20)); sceneHeight = 2
  if (Math.abs(L / SPACING - Math.round(L / SPACING)) > 1e-9) throw new Error('L must be a multiple of the row spacing')
  const rows = Math.round(L / SPACING)
  const inRowLad = LAI * SPACING / (ROW_W * ROW_H)
  const leaf = Number(arg('leaf', HSPOT * ROW_H))
  addSpectrum('leaf', EMIS.leaf); addThermal('leaf_temperature', TEMP.leaf)
  // Row-local (u along row, w across row) -> engine (X north, Z east).
  const toXZ = rowAz === 0 ? (u, y, w) => [u, y, w] : (u, y, w) => [w, y, u]
  const rowCentres = Array.from({ length: rows }, (_, k) => -L / 2 + SPACING / 2 + k * SPACING)
  const binding = { spectralName: 'leaf', thermalName: 'leaf_temperature', canopyName: 'row_canopy', materialName: 'leaf_c3' }
  canopies.push({ name: 'row_canopy', label: 'row', structureType: 'canopy', lai: LAI, density: inRowLad, hc: ROW_H, G: 0.5, LIDFa: -0.35, LIDFb: -0.15, hspot: leaf, leafwidth: leaf })
  if (mode === 'eRaytracing') {
    // Leaf centres are drawn in the envelope grown by the leaf half-diagonal and each
    // square is clipped to the exact row box, so LAD is uniform up to the walls.
    const grow = leafReach(leaf)
    const expandedArea = (ROW_W + 2 * grow) * (ROW_H + 2 * grow)
    const perRow = Math.round(inRowLad * expandedArea * L / leafArea(leaf))
    const obj = new Obj(); obj.object('rows')
    let area = 0, pieces = 0
    for (const wc of rowCentres) for (let k = 0; k < perRow; k += 1) {
      const centre = [(random() - 0.5) * L, -grow + random() * (ROW_H + 2 * grow), wc - ROW_W / 2 - grow + random() * (ROW_W + 2 * grow)]
      let poly = makeLeaf(centre, leaf)  // local (u, y, w)
      poly = clipPolygon(poly, 1, 0, 1); poly = clipPolygon(poly, 1, ROW_H, -1)
      poly = clipPolygon(poly, 2, wc - ROW_W / 2, 1); poly = clipPolygon(poly, 2, wc + ROW_W / 2, -1)
      if (poly.length < 3) continue
      const a = polygonArea(poly)
      if (a < 1e-8) continue
      area += a; pieces += 1
      obj.polygon(poly.map(([u, y, w]) => toXZ(u, y, w)))
    }
    const path = obj.save(join(inputDir, 'models', 'rows.obj'))
    items.push({ name: 'rows', type: 'Vegetation', fileName: path, positionFile: writePositions('rows', periodicRows([[L / 2, L / 2, 0]], L, 0)), meshes: [{ name: 'rows', ...binding }] })
    Object.assign(meta, { rowAzimuth: rowAz, rows, rowWidth: ROW_W, rowGap: ROW_GAP, rowHeight: ROW_H, inRowLad, leafShape: LEAF_SHAPE, leafSize: leaf, leafPieces: pieces, leafAreaIndex: area / (L * L), rowCoverEnvelope: ROW_W / SPACING })
  } else {
    // 0.25 m voxels: row edges (0.25/0.75 m from each 1 m cell origin) fall on voxel faces.
    voxel = Number(arg('voxel', 0.25))
    const inside = (c) => {
      const w = rowAz === 0 ? c[2] : c[0]
      return c[1] <= ROW_H && rowCentres.some((wc) => Math.abs(w - wc) <= ROW_W / 2)
    }
    const filler = voxelFillerObj('rows', voxel, inside, [[-L / 2, L / 2], [0, ROW_H], [-L / 2, L / 2]])
    items.push({ name: 'rows', type: 'Vegetation', fileName: filler.path, positionFile: writePositions('rows', [[L / 2, L / 2, 0, 1, 0]]), meshes: [{ name: 'rows', ...binding }] })
    Object.assign(meta, { rowAzimuth: rowAz, rows, voxel, voxels: filler.voxels, inRowLad, hspot: leaf, leafAreaIndex: filler.voxels * voxel ** 3 * inRowLad / (L * L), rowCoverEnvelope: ROW_W / SPACING })
  }
}

if (mode === 'eRaytracing') {
  items.push({ name: 'skirt', type: 'Soil', fileName: skirtObj(L), positionFile: writePositions('skirt', [[L / 2, L / 2, 0, 1, 0]]),
    meshes: [{ name: 'skirt', spectralName: groundKind, thermalName: `${groundKind}_temperature`, canopyName: 'rigid_body', materialName: 'soilset' }] })
  if (!canopies.some((c) => c.name === 'rigid_body'))
    canopies.push({ name: 'rigid_body', label: 'rigid', structureType: 'rigid', lai: 0, density: 0, hc: 2, G: 0.5, LIDFa: 0, LIDFb: 0, hspot: 0, leafwidth: 0 })
}
if (!canopies.length) canopies.push({ name: 'rigid_body', label: 'rigid', structureType: 'rigid', lai: 0, density: 0, hc: 2, G: 0.5, LIDFa: 0, LIDFb: 0, hspot: 0, leafwidth: 0 })

const project = createDefaultProject({ mode, name: `tirt-benchmark-${test}-${mode}`, sceneX: L, sceneY: L, sceneHeight, voxelSize: voxel })
const c = project.configuration
c.outDir = outputDir
c.scene.background = { spectralName: groundKind, thermalName: `${groundKind}_temperature`, materialType: 'Soil', materialName: 'soilset', heterogeneityEnabled: false, angularEffectStrength: 0 }
c.scene.terrain = false
c.light = { zenith: SZA, azimuth: SAA, direct: 1, diffuse: 0, skyTemperature: 3 }
Object.assign(c.sensor, {
  x: resolution, y: resolution, bands: String(BAND_NM), projection: 'parallel', vza: 0, vaa: 0,
  continuousBands: false, principalPlane: true, hemisphere: false, cruiseEnabled: false,
  image: true, process: false, radiationProcess: false, energyProcess: false, temperature: false, albedo: false
})
Object.assign(c.control, { depth, samples, periodicTraversalCount: 5, heterogeneousVoxel: false, skyboxEnabled: false, radiationSolver: 'traditional' })
c.atmosphere.enabled = false
c.fluid.enabled = false
c.spectra = spectra
c.thermals = thermals
c.canopies = canopies
c.objects = { count: items.length, names: items.map((i) => i.name), items }

writeFileSync(join(inputDir, 'project.json'), stringifyProject(project), 'utf8')
writeFileSync(join(out, 'case.json'), JSON.stringify({ ...meta, L, sceneHeight, emissivity: EMIS, temperature: TEMP }, null, 2), 'utf8')
console.log(JSON.stringify({ project: join(inputDir, 'project.json'), ...meta, L }))
