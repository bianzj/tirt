# TiRT Directional Three.js GUI

This folder contains the browser interface for the TiRT directional runner.

## Start

From the project root:

```bash
python gui/server.py --port 8765
```

Open `http://127.0.0.1:8765/` in a browser. The server uses the TiRT root `run()` function, so the result panel reflects the current directional implementation and the selected input values.

## Layout

- Input settings: scene structure, observation angles, spectra and component temperatures.
- The input panel is grouped as scene structure, sensor settings (observation geometry and wavelengths), and spectral/temperature inputs. Terrain, urban, row and crown controls appear only when their corresponding scene choice is active.
- Spectra and temperatures are separate modules. Spectra contains wavelength bands and component emissivities with file inputs; temperatures contains sunlit and shaded values for each visible component.
- Scene preview: Three.js view of the ground, terrain, buildings, vegetation, sun and observation rays.
- Directional results: brightness temperature, radiance, output records, polar map, parallel solar principal plane and perpendicular solar principal plane.

The form always generates the solar principal plane and the vertical solar principal plane from one VZA step. Additional observation directions can be selected as:

- File: one `vza vaa` direction per line, or `vza vaa sza saa`.
- Manual: enter paired `VZA` and `VAA` lists, for example `0;10;20` and `0;90;180`.

The backend `geometry_mode` still follows the project input convention:

- `0`: use the given `vza` and `vaa` lists.
- `1`: generate relative solar azimuths; `geometry_view` can limit this to `parallel` or `perpendicular`.
- `3`: generate both required principal planes and append file/manual observation directions.
- `2`: read a TXT direction file, or use an uploaded local TXT file.

Three.js and `OrbitControls` are vendored in `gui/vendor/`, so the 3D scene works through the local server without a browser network dependency.
