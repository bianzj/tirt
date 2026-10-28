# TIRTEB against SCOPE

Both scripts drive TIRTEB (`plane` + `hom`) with the inputs of a SCOPE verification run and compare against SCOPE's shipped outputs; SCOPE itself is not run.

| File | Reference | Comparison |
|---|---|---|
| `compare_scope.py` | SCOPE 1.73 `output/verificationdata` (`--scope <SCOPE_v1.73>`) | 24 h half-hourly series (DOY 169, LAI 3): Rn, LE, H, G and their canopy/soil parts, soil and canopy temperatures, nadir 10.5 µm brightness temperature |
| `compare_scope2.py` | SCOPE 2 `output/verificationdata` (`--scope <SCOPE2 folder with SCOPE.m>`) | 30 Latin-hypercube scenarios; each is one TIRTEB time step at the scenario's sun zenith (equator, equinox), with SCOPE's diffuse fraction, z0/d (`zo_and_d.m`) and soil spectrum. Also feeds SCOPE's temperatures to TiRT to isolate the radiative transfer |

Outputs go to `runs/` (`scope_*` and `scope2_*` csv/png); the figures and summary tables used in the main README are also written to `docs/validation/tirteb_scope173.*` and `docs/validation/tirteb_scope2.*`.

Reading the SCOPE 2 statistics:
- The hypercube samples Ta and ea independently; 17 of 30 scenarios have ea above saturation. The summary reports all scenarios and the 13 unsaturated ones (filled markers in the figure).
- SCOPE ran the scenarios as a daily series, so its soil temperature and G carry over from the previous scenario; TIRTEB starts every scenario from the air temperature. Soil temperature and G are therefore not comparable.
- SCOPE's `fluxes.csv` net radiation is on average 65 W m-2 above its own radiation budget (`radiation.csv`: in - reflected - emitted). `Rn_budget` compares against the latter.

SCOPE 1.73 applies its stable correction only for L > 4000 m while it limits L to 100 m, so its nights are neutral; TIRTEB keeps a Monin–Obukhov stable correction (as SCOPE 2). This changes nighttime turbulent exchange, but the net canopy/soil difference also depends on the radiation partition and the initial soil thermal state.
