export const SCENE_THEME = Object.freeze({
  background: 0xf1f4f7,
  ground: 0xaebecd,
  gridMajor: 0x5f7f9f,
  gridMinor: 0x91a6ba,
  building: 0x686f6b,
  vegetation: 0x327d56,
  trunk: 0x624f3b,
  sun: 0xe7ae3f,
  view: 0xd6534b,
  north: 0x1c5fae,
  skyLight: 0xe8f1fb,
  groundLight: 0xb8c8d8,
  sunlight: 0xfff7e8,
});

export function configureStreamSimRendering(THREE, scene, renderer) {
  scene.background = new THREE.Color(SCENE_THEME.background);
  scene.fog = null;
  renderer.setClearColor(SCENE_THEME.background, 1);
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 0.98;
}
