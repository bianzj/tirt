# Optical JSON entry

User requirement: select homogeneous, row, discrete-crown and other supported
scenes through JSON files, all run with optical_run.py. Code/model is primary;
no new package delivery or GUI dependency. Continue prior two-plot/LAI behavior.

Implementation: migrate common example simulation/export/plot code into
base/scenario_runner.py; leave examples/common.py as a compatibility facade.
Add schema_version=1 JSON configs beside the ten optical scripts, making JSON
their sole parameter source. optical_run.run supports the new format while
preserving the previous settings JSON and --output CSV contract. New jobs use
--output-dir and optional --no-plots/--show, defaulting to results/<JSON stem>
beside the input JSON. Record input hashes; custom scene names need no .py file.

Checks: all ten JSONs through the real CLI; homogeneous/row/crown direct-call
parity; BOM, Unicode/spaces and independent cwd; malformed domain/version;
legacy input compatibility; existing example and thermal tests. No physics edits.

Completed: all ten JSON files ran through the real optical_run.py CLI from a
different cwd, generating both plots. All 26 scenario cases reran successfully;
all 52 CSV files remain byte-for-byte identical to the previous implementation.
13 targeted tests passed. Independent review additionally blocked imports of
examples and confirmed a custom JSON job still runs. No new distribution bundle.
Evidence: cases/_optical_json_entry_20261024/verification.json.
