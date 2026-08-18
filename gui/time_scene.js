import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";

const get = (id) => document.getElementById(id);
const sceneSize = 100;
const halfScene = sceneSize / 2;
let scene;
let camera;
let renderer;
let controls;
let objects;

function number(id, fallback = 0) { const value = Number(get(id)?.value); return Number.isFinite(value) ? value : fallback; }
function material(color, roughness = .8) { return new THREE.MeshStandardMaterial({ color, roughness, metalness: .05 }); }
function addMesh(geometry, mat, position = [0, 0, 0], rotation = [0, 0, 0]) { const mesh = new THREE.Mesh(geometry, mat); mesh.position.set(...position); mesh.rotation.set(...rotation); objects.add(mesh); return mesh; }
function addLine(start, end, color, opacity = 1) { const geometry = new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(...start), new THREE.Vector3(...end)]); objects.add(new THREE.Line(geometry, new THREE.LineBasicMaterial({ color, transparent: opacity < 1, opacity }))); }
function seeded(index) { const value = Math.sin(index * 12.9898 + 78.233) * 43758.5453; return value - Math.floor(value); }

function heightAt(x, z, surface) {
  if (surface === "slope") return -Math.tan(number("terrain-slope", 30) * Math.PI / 180) * (x * Math.sin(number("terrain-aspect") * Math.PI / 180) + z * Math.cos(number("terrain-aspect") * Math.PI / 180));
  if (surface !== "terrain") return 0;
  const radius = Math.max(number("terrain-radius", 100), 5); const height = Math.max(number("terrain-height", 20), 0); const density = Math.max(number("terrain-density", 1), .1); const moundCount = Math.min(18, Math.max(4, Math.round(density * 5)));
  let result = 0;
  for (let index = 0; index < moundCount; index++) { const mx = seeded(index * 3 + 20) * (sceneSize - radius * .3) - halfScene; const mz = seeded(index * 3 + 21) * (sceneSize - radius * .3) - halfScene; const mr = radius * (.15 + seeded(index * 3 + 23) * .25); const distance = Math.hypot(x - mx, z - mz) / mr; if (distance < 1) result = Math.max(result, height * (.35 + seeded(index + 2) * .65) * (1 - distance * distance)); }
  return result;
}

function fitScene() {
  const bounds = new THREE.Box3().setFromObject(objects); if (bounds.isEmpty()) return; const center = bounds.getCenter(new THREE.Vector3()); const size = bounds.getSize(new THREE.Vector3()); const maxDimension = Math.max(size.x, size.y, size.z, 1); const distance = maxDimension / (2 * Math.tan(camera.fov * Math.PI / 360)) * 1.2; const direction = camera.position.clone().sub(controls.target).normalize(); controls.target.copy(center); camera.position.copy(center).add(direction.multiplyScalar(Math.max(distance, 85))); camera.near = .1; camera.far = Math.max(distance * 5, 300); camera.updateProjectionMatrix(); controls.update(); }

function renderScene() {
  if (!objects) return;
  objects.clear();
  const surface = get("surface-model").value; const vegetation = get("vegetation-model").value; const base = surface === "urban" ? 0x3a4642 : surface === "terrain" || surface === "slope" ? 0x5c6552 : 0x384d42; const ground = material(base);
  const segments = surface === "plane" || surface === "urban" ? 1 : 28; const geometry = new THREE.PlaneGeometry(sceneSize, sceneSize, segments, segments); geometry.rotateX(-Math.PI / 2); const positions = geometry.attributes.position; for (let index = 0; index < positions.count; index++) positions.setY(index, heightAt(positions.getX(index), positions.getZ(index), surface)); positions.needsUpdate = true; geometry.computeVertexNormals(); addMesh(geometry, ground);
  if (surface === "plane" || surface === "urban") { const grid = new THREE.GridHelper(sceneSize, 20, 0x45615b, 0x263732); grid.position.y = .02; objects.add(grid); }
  if (surface === "urban") { const building = material(0x79827b); const h = Math.max(number("building-height", 20), 1); const l = Math.max(number("building-length", 10), 1); const w = Math.max(number("building-width", 10), 1); [[-20, -12], [4, -8], [24, 12], [-12, 20]].forEach(([x, z], index) => addMesh(new THREE.BoxGeometry(l * (index % 2 ? .8 : 1), h * (.35 + seeded(index) * .5), w * (index % 2 ? 1.15 : .8)), building, [x, h * (.18 + seeded(index) * .25), z], [0, number("terrain-aspect") * Math.PI / 180, 0])); }
  const green = material(0x5ca67c); const trunk = material(0x725a45);
  if (vegetation === "row") { const width = Math.max(number("row-width", .8), .1); const gap = Math.max(number("row-blank", .8), .05); const rowHeight = Math.max(number("row-height", 1.5), .2); const angle = number("row-azimuth", 90) * Math.PI / 180; const along = [Math.cos(angle), Math.sin(angle)]; for (let cross = -halfScene; cross <= halfScene; cross += width + gap) for (let distance = -halfScene; distance <= halfScene; distance += 2) { const x = along[0] * distance - along[1] * cross; const z = along[1] * distance + along[0] * cross; addMesh(new THREE.ConeGeometry(Math.max(width * .28, .15), rowHeight, 6), green, [x, heightAt(x, z, surface) + rowHeight / 2, z]); } }
  else if (vegetation === "crown") { const radA = Math.max(number("crown-radius", 1.5), .3); const radB = Math.max(number("crown-width", 3), .3); const crownHeight = Math.max(number("crown-height", 3), 1); const density = Math.max(number("crown-density", .1), .01); const count = Math.min(34, Math.max(8, Math.round(density * 100))); for (let index = 0; index < count; index++) { const x = seeded(index * 2) * 90 - 45; const z = seeded(index * 2 + 1) * 90 - 45; const y = heightAt(x, z, surface); addMesh(new THREE.CylinderGeometry(.12, .18, crownHeight, 6), trunk, [x, y + crownHeight / 2, z]); const crown = addMesh(new THREE.SphereGeometry(radA, 12, 8), green, [x, y + crownHeight + radB * .5, z]); crown.scale.set(1, radB / radA, 1); } }
  else if (vegetation === "hom") { const count = 900 + Math.round(Math.max(number("lai", 2), .1) * 70); const facet = new THREE.CylinderGeometry(.2, .2, .15, 8); const instances = new THREE.InstancedMesh(facet, green, count); const dummy = new THREE.Object3D(); for (let index = 0; index < count; index++) { const x = seeded(index * 2 + 1) * 96 - 48; const z = seeded(index * 2 + 2) * 96 - 48; dummy.position.set(x, heightAt(x, z, surface) + .3 + seeded(index * 2 + 3) * Math.max(number("height", 3), 1), z); dummy.rotation.set((seeded(index + 4) - .5) * .6, seeded(index + 5) * Math.PI, (seeded(index + 6) - .5) * .6); dummy.updateMatrix(); instances.setMatrixAt(index, dummy.matrix); } instances.instanceMatrix.needsUpdate = true; instances.frustumCulled = false; objects.add(instances); }
  const observations = window.temporalState?.observationRows || []; const last = observations.length ? observations[observations.length - 1] : {}; const sza = Number(last.sza) || 30; const saa = Number(last.saa) || 0; const originY = heightAt(0, 0, surface) + 1; const distance = 45; const sun = [Math.sin(saa * Math.PI / 180) * Math.cos(sza * Math.PI / 180) * distance, originY + Math.sin(sza * Math.PI / 180) * distance, Math.cos(saa * Math.PI / 180) * Math.cos(sza * Math.PI / 180) * distance]; addMesh(new THREE.SphereGeometry(.45, 16, 12), material(0xedb866), sun); addLine([0, originY, 0], sun, 0xedb866, .6); const vza = number("vza"); const vaa = number("vaa"); const view = [Math.sin(vaa * Math.PI / 180) * Math.cos(vza * Math.PI / 180) * 32, originY + Math.sin(vza * Math.PI / 180) * 32, Math.cos(vaa * Math.PI / 180) * Math.cos(vza * Math.PI / 180) * 32]; addLine([0, originY, 0], view, 0xed8e7d, .7);
  get("scene-title").textContent = `${get("surface-model").selectedOptions[0].textContent} · ${get("vegetation-model").selectedOptions[0].textContent}`; get("scene-readout").textContent = `${get("surface-model").selectedOptions[0].textContent} / ${get("vegetation-model").selectedOptions[0].textContent}`; get("solar-readout").textContent = `SZA ${sza.toFixed(1)}° · SAA ${saa.toFixed(1)}°`; fitScene(); get("scene-loading").classList.add("is-hidden");
}

function initScene() { const container = get("scene-container"); scene = new THREE.Scene(); scene.background = new THREE.Color(0x0c1413); camera = new THREE.PerspectiveCamera(42, container.clientWidth / container.clientHeight, .1, 300); camera.position.set(70, 55, 75); renderer = new THREE.WebGLRenderer({ antialias: true }); renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2)); renderer.setSize(container.clientWidth, container.clientHeight); container.appendChild(renderer.domElement); controls = new OrbitControls(camera, renderer.domElement); controls.enableDamping = true; controls.target.set(0, 2, 0); scene.add(new THREE.HemisphereLight(0xbed5ce, 0x1c2623, 1.8)); const light = new THREE.DirectionalLight(0xffe1b0, 2.2); light.position.set(8, 13, 6); scene.add(light); objects = new THREE.Group(); scene.add(objects); renderScene(); const resize = () => { renderer.setSize(container.clientWidth, container.clientHeight); camera.aspect = container.clientWidth / container.clientHeight; camera.updateProjectionMatrix(); renderScene(); }; window.addEventListener("resize", resize); const animate = () => { requestAnimationFrame(animate); controls.update(); renderer.render(scene, camera); }; animate(); }

window.refreshTemporalScene = renderScene;
initScene();
