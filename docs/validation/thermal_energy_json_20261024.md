# Thermal and energy-balance JSON entry validation — 2026-10-24

Canonical source: `C:\work\tirt`. Environment: Conda `python311`.

## Implemented interface

- `thermal_run.py`: 16 thermal JSON scenes, editable `infile`, `run(..., plots=True)`.
- `tirteb_run.py`: homogeneous, row and crown energy-balance JSON examples with the same calling style.
- Existing CSV solver calls remain supported. The numerical energy-balance solver body is unchanged.
- JSON contains scientific inputs; Python handles output location, CSV exports and optional figures.

## Executed cases

`python examples/run_all.py` and `python examples/verify_outputs.py`:

- 26 supported radiation cases; 2,138 valid rows; 52 PNGs.
- All 52 spectral/angular CSVs are byte-identical to the pre-change reference outputs.
- Six unsupported optical slope/terrain bare/hom/row combinations remain explicitly unsupported.
- Both thermal and energy-balance no-argument entries also executed from an external directory containing spaces.

`python tirteb_run.py --input examples/energy_balance/<scene>.json`:

| Scene | Time steps | Closed | Observation rows | Figures |
|---|---:|---:|---:|---:|
| plane_hom | 48 | 48 | 432 | 2 |
| plane_row | 48 | 48 | 432 | 2 |
| plane_crown | 48 | 48 | 432 | 2 |

The temperature and flux export rows equal their corresponding original output.csv columns.
Source/input hashes, finite temperatures/fluxes and all six energy-balance PNGs were checked.
Representative thermal spectral/angular and energy-balance temperature/flux figures were visually inspected.

Regression checks compare JSON execution against the original CSV solver at absolute tolerance `1e-6`.
`python -m unittest discover -s tests -q` passed all 17 tests (56.844 seconds).
This comparison tolerance does not alter the solver's energy-balance closure criterion.
The tests cover invalid domains, BOM input, CSV compatibility, output defaults, no-plots reruns,
an external working directory and an output path containing spaces and Chinese characters.

Local machine-readable receipts:

- `examples/results/outputs_verified.json`
- `examples/energy_balance/results/<scene>/verification.json`
- `cases/_thermal_energy_json_20261024/verification.json`

These are interface and numerical-regression checks using demonstration inputs, not independent physical validation.
