# HiSTREAM reference runs for the TiRT validation

This folder generates the HiSTREAM eRaytracing references used in [the TiRT test document](../../docs/TiRT_validation.md). It does not modify StreamSim; `createDefaultProject()` from `C:\work\streamsim\src\renderer\src\project-schema.js` is only referenced.

| File | Purpose |
|---|---|
| `build_case.mjs` | Builds one case: `project.json` (via `createDefaultProject`), OBJ geometry and instance-position files under `runs/<case>/input/`. |
| `run_benchmark.py` | Builds each case, runs `histream.exe <mode> project.json`, and reduces every view image to mean radiance and brightness temperature. Slope cases are rendered as the flat scene in slope-local coordinates, one engine run per view. `scripts/validate_current_model.py` also imports its slope-frame helpers. |

## Results

`results/` keeps only the statistics used in the comparison: one CSV per scene with the mean radiance and brightness temperature for every view direction at 10.5 µm. These are the frozen thermal references in `validation/reference/histream_reference.csv` (identical per direction). The raw run directories (view images, OBJ geometry, engine logs) are not kept; rerun the cases below to regenerate them.

| Scene | Result |
|---|---|
| Plane bare | `T1_plane_bare_eRaytracing.csv` |
| Plane homogeneous | `T2_plane_hom_eRaytracing.csv` |
| Plane forest (0.15 m leaves) | `T3_plane_crown_eRaytracing_leaf0.15.csv` |
| Rows, azimuth 0 / 90 | `T5_row_az0_eRaytracing.csv`, `T5_row_az90_eRaytracing.csv` |
| Slope + homogeneous / forest | `S1_slope_hom_eRaytracing.csv`, `S1_slope_crown_eRaytracing.csv` |
| Urban bare | `T4_urban_bare_eRaytracing.csv` |
| Urban + forest | `UV_urban_crown_eRaytracing.csv` |

The compound-terrain reference (100 m hill, 2 m crowns) and the optical (670/860 nm) references were produced outside this repository; their sources and SHA256 are listed in `validation/reference/reference_manifest.json`.

## Engine

The references were rendered with HiSTREAM `histream.exe` sha256 `2564ea57cab7b0af75bd09960ed32401fed2c1dd6760bad62b311e704fa07d55` (StreamSim 1.1.2, built 2026-10-11). `run_benchmark.py` uses `runs/_engine/histream.exe` when present, otherwise the StreamSim Release engine; `HISTREAM_EXE` overrides both. To freeze the engine before rerunning:

```powershell
robocopy C:\work\streamsim\models\bin_x64\Release C:\work\tirt\validation\streamsim\runs\_engine /E /XF *.pdb *.lib /XD .idea .vscode
```

## Rerun

```powershell
cd C:\work\tirt\validation\streamsim
$py = 'C:\work\miniconda\envs\python311\python.exe'
& $py run_benchmark.py --cases T1:eRaytracing T2:eRaytracing T4:eRaytracing
& $py run_benchmark.py --cases T3:eRaytracing --extra "--leaf 0.15" --tag _leaf0.15
& $py run_benchmark.py --cases T5az0:eRaytracing T5az90:eRaytracing
& $py run_benchmark.py --cases S1hom:eRaytracing S1crown:eRaytracing UV:eRaytracing
```

Each case takes 2–40 s on the RTX PRO 6000 and needs Node (≥18), numpy and Pillow. Results are written to `runs/` (not version-controlled); copy the `<test>_<mode>[tag].csv` files to `results/` to update the stored statistics. `HISTREAM_EXE`, `STREAMSIM_ROOT` and `NODE` can be overridden through environment variables.

## Outputs of a rerun

- `runs/<test>_<mode>[tag].csv` (stored copies in `results/`): `test,mode,vza,vaa,sza,saa,radiance,bt_K`. `radiance` is the mean over valid pixels of the orthorectified image (W m⁻² sr⁻¹ µm⁻¹); `bt_K` is its inverse Planck at 10.5 µm. `vaa`/`saa` are compass azimuths; nadir is rendered once and repeated for vaa = 0/90/180/270.
- `runs/<case>/diagnostics.csv`: per-view pixel std, valid-pixel fraction and half-image difference. `runs/<case>/case.json`: realised geometry. `runs/<case>/engine.log`: engine stdout.

## Settings (all cases)

- Band 10500 nm, radiance images, 512 × 512 parallel projection, VZA 0–75 step 5 in the solar and perpendicular planes.
- Sun SZA 30, SAA 0. Sky temperature 3 K (reflected sky ≈ 0), atmosphere off. No direct sun at 10.5 µm.
- Emissivity ε = 1 − refTir, tauTir = 0: soil/street 0.955, leaf 0.985, roof 0.950, wall 0.930.
- Temperatures (sunlit / shaded, K): soil, roof, wall, street 318.15 / 303.15; leaf 306.15 / 303.15 (urban leaves 31.5 °C).
- `periodicTraversalCount` 5, depth 4 (multiple scattering), 32 samples per pixel.

## Test geometry

| Test | Domain | Scene |
|---|---|---|
| T1 plane + bare | 20 m | background soil only |
| T2 plane + homogeneous | 20 m | 26 667 square leaves of 0.15 m, LAI 1.50, layer 0–1 m, uniform-sphere normals |
| T3 plane + crowns | 100 m | 200 spheres, R = 3 m, centres at 6 m, random positions; square leaves of 0.15 m, LAI 1.5 |
| T5 rows | 20 m | 20 rows, 0.5 m wide, 0.5 m gap, 1 m high, azimuth 0 or 90; 0.15 m leaves clipped to the row cross-section, in-row LAD 3.0 |
| T4 urban + bare | 200 m | 60 boxes of 10 × 10 × 20 m (0.0015 m⁻²), azimuth 0; street = background |
| S1 slope | as T2 / T3 | 30° slope facing east, rendered as the flat T2 / T3 scene with slope-local sun (zenith 41.4°) and views; views with local zenith > 89° are skipped |
| UV urban + crowns | 300 m | 135 boxes of 10 × 10 × 20 m; 1530 crowns (centre 6 m, R 3 m) of 0.15 m leaves on streets only |

Periodic boundary: objects are instanced in an explicit 3 × 3 tiling with a ground skirt, because eRaytracing wraps only camera rays.
