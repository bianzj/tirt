# StreamSim / HiSTREAM side of the TiRT pilot benchmark

This folder holds the scripts that generate HiSTREAM reference simulations for the TiRT tests (T1–T5, slope S1, terrain TV, urban+vegetation UV). It does not modify StreamSim. The engine (`C:\work\streamsim\models\bin_x64\Release\histream.exe`) and `createDefaultProject()` from `C:\work\streamsim\src\renderer\src\project-schema.js` are only referenced.

| File | Purpose |
|---|---|
| `build_case.mjs` | Builds one case: `project.json` (via `createDefaultProject`), OBJ geometry and instance-position files under `runs/<case>/input/`. |
| `run_benchmark.py` | Builds each case, runs `histream.exe <mode> project.json` the same way `server.mjs` does for non-facet modes, then reduces every view image. The engine runs with cwd = input dir, its own `outDir`, and a timeout that kills the process tree. Slope cases (S1) are rendered as the flat scene in slope-local coordinates with one engine run per view, because the engine takes only a single free view direction. |
| `decompose.py` | Four-component (sunlit/shaded soil and leaf) view-fraction decomposition from blackbody indicator runs. |
| `engine_check.py` | Writes `runs/engine_manifest.csv` (engine sha256 for every case) and `runs/old_vs_new.md` (every result CSV compared with `runs/_old_engine/`). |
| `summarize.py` | Writes `runs/<test>_streamsim.csv` (all main modes) and `runs/summary.md` (noise, mode spread, difference to `runs/tirt/*.csv`). |
| `tirt_side.py` | TiRT side of the benchmark: same angles, band, emissivities and temperatures, written to `runs/tirt/<test>.csv`. |
| `compare.py` | Plots TiRT against every StreamSim mode (principal and perpendicular planes) into `runs/compare/comparison.png`, and writes bias/RMSE/max and angular range vs eRaytracing to `runs/compare/summary.csv`. Variant runs (`_seed2`, `_depth1`, …) appear in the statistics but not in the plot. |
| `overview_plot.py` | Validation figures for GitHub: for every benchmarked scene (forest with both 0.9 m and 0.15 m leaves), HiSTREAM eRaytracing, a 1D homogeneous run (TiRT on flat ground with a horizontally homogeneous layer of the same LAI, i.e. without terrain, buildings or canopy structure) and TiRT. Each figure shows scene sketch, solar principal plane and perpendicular plane on one brightness-temperature axis: `docs/validation/<nn>_<scene>.png`, with statistics in `docs/validation/overview.csv`. |
| `scatter_breakdown.py` | Splits the scattered radiance of TV and UV into four paths (within forest, forest emission reflected by the macro surface, within macro surface, macro emission reflected by the forest) for HiSTREAM and TiRT: `runs/compare/scatter_breakdown.csv`, `runs/compare/scatter/<scene>.png`. HiSTREAM paths come from cold (3 K) emitters and black (ε = 1) reflectors (`build_case.mjs --cold/--black`, reduced with `--keep-zero`); TiRT paths drop the Planck term of switched-off sources from Σ ε_eff·B(T), because TiRT's `planck()` reads temperatures below 100 as °C. |
| `crown_fractions_plot.py` | Forest (T3) four-component fractions: TiRT (hspot = leaf size / 6 m crown depth) against every eRaytracing decomposition, written to `runs/compare/crown_fractions.png` and `runs/tirt/T3_fractions_TiRT.csv`. `--old <csv>` overlays an earlier TiRT fractions file. |

## Engine snapshot

The runs use a frozen copy of the engine in `runs/_engine/`: histream.exe, its DLLs, `shader/` and `defined/`, without PDB/LIB files (about 107 MB). StreamSim's own `models/bin_x64/Release` was rebuilt by someone else during this work.
- Current engine: sha256 `2564ea57…`, built 2026-10-11 10:11 (StreamSim 1.1.2). Every eRaytracing case in `runs/` was rerun with it (see `runs/engine_manifest.csv`; the log is `runs/rerun_2564ea57.log`).
- Its changes are in VoxelRT/VoxelEB only: turbid cells are integrated segment by segment, the Kuusk hotspot uses each canopy's own density and leaf size, and diffuse scattering is done by collision order. The eRaytracing shaders are byte-identical, and every eRaytracing result and fraction file is unchanged to the printed precision (`runs/old_vs_new.md`).
- Previous engine `dac2b149…` (2026-10-11 00:38): snapshot in `runs/_engine_dac2b149/`, results in `runs/_old_engine/`.
- Earlier engine `3b28dfe8…`: results in `runs/_old_engine_3b28dfe8/`. Against it, the `dac2b149` full-depth results differed by ≤ 0.011 K and the fractions by ≤ 0.0008.
- The `_depth1` results of `dac2b149` rose by +0.09 to +0.19 K over `3b28dfe8`. The older engine left about 0.5% of camera samples at 0 at depth 1 (the periodic wrap used up the only path step); later engines no longer do.

To recreate the snapshot:

```powershell
robocopy C:\work\streamsim\models\bin_x64\Release C:\work\tirt\benchmark\streamsim\runs\_engine /E /XF *.pdb *.lib /XD .idea .vscode
```

`run_benchmark.py` uses `runs/_engine/histream.exe` when it exists. Otherwise it falls back to the StreamSim Release engine; the `HISTREAM_EXE` environment variable overrides both.

## Rerun

```powershell
cd C:\work\tirt\benchmark\streamsim
$env:PYTHONDONTWRITEBYTECODE = '1'
$py = 'C:\work\miniconda\envs\python311\python.exe'
& $py tirt_side.py                                       # TiRT side
& $py run_benchmark.py                                   # T1..T4 main cases (32 spp, depth 4)
& $py run_benchmark.py --cases T2:eRaytracing T3:eRaytracing T4:eRaytracing --depth 1 --tag _depth1
& $py run_benchmark.py --cases T2:eRaytracing T3:eRaytracing T4:eRaytracing T2:eVoxelRT T3:eVoxelRT --seed 2 --tag _seed2
& $py run_benchmark.py --cases T2:eRaytracing T3:eRaytracing T4:eRaytracing --samples 128 --tag _spp128
& $py run_benchmark.py --cases T3:eRaytracing T3:eVoxelRT --extra "--leaf 0.3" --tag _leaf0.3
& $py run_benchmark.py --cases T3:eRaytracing --extra "--leaf 0.15" --tag _leaf0.15
& $py run_benchmark.py --cases T5az0:eRaytracing T5az90:eRaytracing T5az0:eVoxelRT T5az90:eVoxelRT
& $py run_benchmark.py --cases T5az0:eRaytracing T5az90:eRaytracing --depth 1 --tag _depth1
& $py run_benchmark.py --cases S1hom:eRaytracing S1crown:eRaytracing UV:eRaytracing TV:eRaytracing
& $py run_benchmark.py --cases TV:eRaytracing UV:eRaytracing --depth 1 --tag _depth1
# scattering paths (TV: macro = terrain; UV: macro = street,wall,roof), each at --depth 2 and 4 (and 1 for _srcF/_srcM)
& $py run_benchmark.py --cases TV:eRaytracing --depth 4 --keep-zero --extra "--cold terrain" --tag _srcF_d4
& $py run_benchmark.py --cases TV:eRaytracing --depth 4 --keep-zero --extra "--cold terrain --black terrain" --tag _srcF_Mblack_d4
& $py run_benchmark.py --cases TV:eRaytracing --depth 4 --keep-zero --extra "--cold leaf" --tag _srcM_d4
& $py run_benchmark.py --cases TV:eRaytracing --depth 4 --keep-zero --extra "--cold leaf --black leaf" --tag _srcM_Fblack_d4
& $py scatter_breakdown.py
& $py summarize.py
& $py engine_check.py                                    # after decompose.py runs (see below)
& $py compare.py
& $py overview_plot.py
```

Each case takes 2–40 s on the RTX PRO 6000. The scripts need Node (≥18), plus numpy and Pillow in the Python env. `HISTREAM_EXE`, `STREAMSIM_ROOT` and `NODE` can be overridden through environment variables.

## Outputs

- `runs/<test>_<mode>[tag].csv` contains `test,mode,vza,vaa,sza,saa,radiance,bt_K`.
  - `radiance` is the mean over all valid pixels of the orthorectified image, in W m-2 sr-1 um-1.
  - `bt_K` is the inverse Planck of that mean at 10.5 µm, with c1 = 11910.439340652e4 and c2 = 14388.291040407 µm K.
  - `vaa`/`saa` are compass azimuths of the sensor and sun (TiRT convention). The engine renders nadir once; that row is repeated for vaa = 0/90/180/270.
- `runs/<case>/diagnostics.csv` has, per view: pixel std, valid-pixel fraction, and the difference between the means of the two image halves.
- `runs/<case>/case.json` records the realised geometry (leaf count, LAI, crown count, LAD, building count…). `runs/<case>/engine.log` is the engine stdout.

## Settings (all cases)

- **Band:** 10500 nm only. `sensor.temperature = false`, so images are radiance.
- **Sun:** SZA 30, SAA 0. Radiation fractions are direct = 1, diffuse = 0 (shortwave only; irrelevant at 10.5 µm).
- **Views:** `principalPlane` (VZA 0–75 step 5 in the solar plane and the perpendicular plane), parallel projection, 512 × 512.
- **Sky and atmosphere:** `skyTemperature` = 3 K, so the reflected sky is ≈ 0 (B(10.5 µm, 3 K) underflows to 0). Atmosphere off, skybox off, background Hapke/heterogeneity off, heterogeneous voxels off.
- **No direct sun at 10.5 µm:** in `offline.rgen` (eRaytracing), solar next-event estimation is applied only for λ ≤ 2.5 µm. In VoxelRT, thermal bands (> 5 µm) also exclude the direct term.
- **Emissivity:** ε = 1 − refTir with tauTir = 0. Soil/street 0.955, leaf 0.985, roof 0.950, wall 0.930.
- **Temperatures (sunlit / shaded, K):** soil, roof, wall, street 318.15 / 303.15; leaf 306.15 / 303.15.
- **Sampling:** `periodicTraversalCount` = 5, `depth` = 4 (multiple scattering between surfaces), 32 samples per pixel.

## Test geometry

| Test | Domain | eRaytracing | eVoxelRT |
|---|---|---|---|
| T1 plane + bare | 20 m | background soil only | same |
| T2 plane + homogeneous | 20 m | 26 667 square leaves, side 0.15 m, LAI 1.50, layer 0–1 m, uniform-sphere normals (G = 0.5) | voxel 0.5 m slab 0–1 m, LAD 1.5, G 0.5, Kuusk hspot 0.15 m |
| T3 plane + crowns | 100 m | 200 spheres, R = 3 m, centres at 6 m, uniform random positions; 93 square leaves of 0.9 m per crown (LAD 0.663, LAI 1.507); 4 crown variants with random rotation | quantised sphere (912 voxels of 0.5 m), LAD rescaled to 0.813 so scene LAI = 1.5 after crown overlap; hspot 0.9 m |
| T5 rows (`T5az0`, `T5az90`) | 20 m | 20 rows, 0.5 m wide, 0.5 m gap, 1 m high, along compass azimuth 0 or 90, exactly aligned with the domain. Square leaves of 0.15 m (hspot 0.15 × row height 1 m) with uniform-sphere normals. Leaf centres are drawn in the row box grown by the leaf half-diagonal, and every square is clipped to the exact 0.5 × 1 m cross-section, so the in-row LAD is uniform (3.0). Realised scene LAI 1.506; row cover 0.5 | 0.25 m voxels; rows are exactly 2 × 4 voxels; LAD 3.0, hspot 0.15 m |
| T4 urban + bare | 200 m | 60 boxes of 10 × 10 × 20 m (0.0015 m⁻²), uniform random positions, azimuth 0; separate `roof`/`wall` meshes with outward normals; street = background | not run |
| S1 slope + homogeneous / crowns (`S1hom`, `S1crown`) | as T2 / T3 | 30° slope facing east (aspect 90, TiRT `slope1` convention: aspect = azimuth the slope faces). TiRT runs the flat canopy model in slope-local angles (trees perpendicular to the slope, density and LAI per slope area), so the reference is the T2 scene and the T3 scene with 0.15 m leaves, rendered with slope-local sun (zenith 41.4°) and view directions, one engine run per view. The 4 views with local zenith > 89° (slope seen from behind) are skipped, as in TiRT | not run |
| TV terrain + crowns | 600 m | 54 non-overlapping cones (random sequential addition), height 20 m, base radius 30 m (slope 33.7°, cover 0.424), 128-segment meshes; 7200 spherical crowns (R 3 m, 0.02 m⁻² of horizontal area, LAI 1.5) of 0.15 m leaves, each with its centre 6 m out along the local surface normal (TiRT evaluates facet crowns in facet-local angles). All soil (ground and cones) has TiRT's terrain emissivity 0.950 and the soil temperatures | not run |
| UV urban + crowns | 300 m | 135 non-overlapping boxes of 10 × 10 × 20 m (0.0015 m⁻²); 1530 spherical crowns (centre 6 m, R 3 m) of 0.15 m leaves placed only on streets, at least 3 m from any building, 0.02 crowns and LAI 1.5 per m² of street; leaves at one temperature, 31.5 °C, the mean leaf temperature TiRT uses | not run |

- **Hotspot / leaf size:** I used the TiRT definition hspot = leaf size / canopy depth with H normalised to 1. T2 therefore uses 0.15 m leaves (H = 1 m). T3 uses 0.15 × crown depth 6 m = 0.9 m. The `_leaf0.3` run shows the sensitivity of T3 to that choice.
- **Periodic boundary:** eRaytracing wraps only camera rays at the domain edge; shadow, scattering and sky rays are not wrapped. All objects are therefore instanced in an explicit 3 × 3 periodic tiling, and a ground "skirt" extends the background to ±1.5 L. The voxel cases include only copies that overlap the domain; the engine clips voxels outside the domain and wraps the rest itself.

## Differences from TiRT that cannot be removed

- **Sunlit/shaded classification.** eRaytracing uses a binary shadow ray from each hit point towards the sun's centre (no penumbra), so surfaces facing away from the sun are shaded. E/W walls are exactly parallel to the sun at SAA 0 (zero irradiance); with outward normals they are classed as shaded. VoxelRT mixes sunlit and shaded radiance per voxel, using a Kuusk bidirectional gap with the canopy `hspot`. The soil under a voxel canopy uses hard-coded `DENSITY`/`HSPOT` constants (`voxelrad_image.comp`), so the voxel hotspot is much weaker (T2 at 30°/0°: VoxelRT 306.4 K, eRaytracing 310.3 K).
- **Multiple scattering.** With depth 4 the TIR reflections between surfaces are included (and sky ≈ 0). The `_depth1` runs show the emission-only part, which is 1.1–1.6 K lower. `scatter_breakdown.py` shows where the scattering in TV and UV comes from (averages over all views, HiSTREAM / TiRT with the forest–background step that `Terrain_Veg` and `Urban_Veg` now include):
  - TV, 1.16 / 1.34 K in total: leaf emission reflected by the ground and cones 0.76 / 0.87 K, within the forest 0.27 / 0.37 K, within the macro surface 0.09 / 0.08 K, macro emission reflected by leaves 0.03 / 0.02 K. Without the forest–background step TiRT had only 0.06 and 0.01 K for the first two, and TV was −1.25 K biased.
  - UV, 1.58 / 1.91 K in total: within the macro surface 0.80 / 1.03 K, leaf emission reflected by streets and walls 0.56 / 0.58 K, within the forest 0.17 / 0.23 K, macro emission reflected by leaves 0.04 / 0.07 K. The remaining UV excess is the building–street exchange, which in TiRT is not blocked by the street trees.
  - Almost all of it is single scattering: depth 2 gives 1.14 K of the 1.16 K (TV) and 1.52 K of the 1.58 K (UV).
- **Voxel quantisation.** VoxelRT uses 0.5 m voxels. Overlapping crowns share voxels (a union), so their LAD is rescaled to keep the scene LAI at 1.5, rather than adding densities where crowns overlap. Spheres are quantised to 912 voxels each (113.1 m³).
- **Scene realisation.** These are single random scenes; the `_seed2` runs give the realisation spread. The Monte-Carlo part is negligible (`_spp128` vs 32 spp: ≤ 0.005 K). Renders are deterministic for a given seed.

## Component-fraction decomposition (T3, T5)

These are component view fractions for T3 from eRaytracing, using the same scene as the main run (seed 1, same crown positions and leaves).
- Every surface is a blackbody (`--blackbody`, refTir = 0) and the sky is 3 K.
- For each fraction there is one indicator run: that component is set to 320 K and the other three to 280 K (`--temps soilSun,soilShade,leafSun,leafShade`).
- The fractions are then f_i = (L_i - B(280)) / (B(320) - B(280)).
- The runs use depth 4. Because ε = 1, every reflected path is zero, so the result is still emission-only. At depth 1 the periodic wrap used up the only path step, and about 0.5% of samples came back as 0 (the fractions summed to 0.986-0.996).

```powershell
cd C:\work\tirt\benchmark\streamsim
$env:PYTHONDONTWRITEBYTECODE = '1'
& C:\work\miniconda\envs\python311\python.exe decompose.py --test T3 --leaf 0.9 0.3 0.15
& C:\work\miniconda\envs\python311\python.exe decompose.py --test T5az0
& C:\work\miniconda\envs\python311\python.exe decompose.py --test T5az90
& C:\work\miniconda\envs\python311\python.exe crown_fractions_plot.py
```

Leaf placement in T3: a leaf is rejected if any corner lies outside the 3 m sphere, so the leaf area is packed into a smaller core and the crowns are more transparent than a turbid sphere of the same LAI. The nadir gap fraction is 0.709 (0.9 m leaves), 0.687 (0.3 m) and 0.679 (0.15 m), approaching the turbid-crown value 0.673 that TiRT uses. With 0.15 m leaves (TiRT hspot 0.025) TiRT − eRaytracing is +0.03 K bias, 0.15 K RMSE; with the default 0.9 m leaves it is −0.38 K, mostly from that extra soil.

Outputs:
- `runs/T3_fractions_eRaytracing.csv` (0.9 m leaves), `runs/T3_fractions_eRaytracing_leaf0.3.csv` and `runs/T3_fractions_eRaytracing_leaf0.15.csv`, with columns `vza,vaa,f_soil_sunlit,f_soil_shaded,f_leaf_sunlit,f_leaf_shaded,sum`.
- `runs/T3_fractions_geometry.json`: realised LAI, crown cover fraction (vertical projection, including periodic copies, on a 2000² grid), and nadir gap fraction (f_soil_sunlit + f_soil_shaded at VZA 0).
- Rows: `runs/T5_row_az0_fractions_eRaytracing.csv` and `runs/T5_row_az90_fractions_eRaytracing.csv`, plus `runs/T5_row_az*_fractions_geometry.json` (realised LAI, row cover 0.5, nadir gap fraction).
- The indicator runs themselves are in `runs/<test>_eRaytracing_frac_*`.
