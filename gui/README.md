# GUI

`gui/` contains the unified TiRT pages. Directional, energy-balance and multi-angle functions share one local service and one desktop window.

## Start

From the project root:

```bash
python gui/launcher.py --server-only --port 8765
```

For development, open `http://127.0.0.1:8765/` in a browser. For normal use, run `start_tirt.bat` or `npm run desktop`; Electron provides the application window without browser tabs or address bars.

## Layout

- Input settings: scene structure, observation angles, spectra and component temperatures. LAI accepts one value or a semicolon-separated list (for example `1;3;5`); multiple LAIs are overlaid on the principal-plane charts and selected through the LAI control.
- The input panel is grouped as scene structure, sensor settings (observation geometry and wavelengths), and spectral/temperature inputs. Terrain, urban, row and crown controls appear only when their corresponding scene choice is active.
- Spectra and temperatures are separate modules. Spectra contains wavelength bands and component emissivities with file inputs; temperatures contains sunlit and shaded values for each visible component.
- Scene preview: Three.js view of the ground, terrain, buildings, vegetation, sun and observation rays.
- Directional results: the selected band decides the displayed quantity — reflectance (BRF) below 2.5 µm, brightness temperature (°C) at 10.5 µm — plus output records, polar map, parallel solar principal plane and perpendicular solar principal plane.

The form always generates the solar principal plane and the vertical solar principal plane from one VZA step. Additional observation directions can be selected as:

- File: one `vza vaa` direction per line, or `vza vaa sza saa`.
- Manual: enter paired `VZA` and `VAA` lists, for example `0;10;20` and `0;90;180`.

The backend `geometry_mode` still follows the project input convention:

- `0`: use the given `vza` and `vaa` lists.
- `1`: generate relative solar azimuths; `geometry_view` can limit this to `parallel` or `perpendicular`.
- `3`: generate both required principal planes and append file/manual observation directions.
- `2`: read a TXT direction file, or use an uploaded local TXT file.

Three.js and `OrbitControls` are vendored in `gui/vendor/`, so the 3D scene works through the local server without a browser network dependency.

## Energy balance GUI

The energy-balance (time series) page is served by the same TiRT server:

```bash
python gui/server.py --port 8765
```

The TiRT desktop window opens the directional model (`/`). Use the `TIRT-Dynamic` button to open the time-series page (`/time.html`), which uses the energy-balance input template and the local `/api/run_time` endpoint. No second server is required. The page now shows a fixed brightness-temperature time series without a quantity selector; `Save processes.csv` still exports every temperature and flux column.
