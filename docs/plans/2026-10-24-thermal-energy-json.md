# Thermal and energy-balance JSON entrypoints

- Keep the existing numerical solvers and CSV calling conventions.
- Give thermal_run.py and tirteb_run.py an editable infile and run(..., plots=True).
- Migrate all 16 thermal examples to JSON, retaining Python shortcuts.
- Add homogeneous, row and crown energy-balance JSON examples with explicit inputs.
- Thermal: spectrum/angles CSV and PNG. Energy balance: temperature/flux time-series CSV and PNG, plus existing solver outputs and closure status.
- Validate JSON before writing outputs; resolve referenced data relative to the JSON.
- Check actual solver results, legacy parity (absolute tolerance 1e-6), no-plot execution, external working directories and bad inputs. Review code and inspect plots.
