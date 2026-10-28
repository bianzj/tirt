# TiRT optical and thermal scenario examples

## Design

User request: one runnable Python file per scene and radiation domain; editable
LAI, angles and multispectral inputs; two plots with wavelength and view angle
on the horizontal axes. No energy balance, as established in this conversation.

Use `examples/optical/` and `examples/thermal/`, with explicit editable parameters
in every script and a shared runner/plotter. Call `base.optical.reflectance` and
the public `run.run` thermal entrypoint. Do not change model physics. Keep all
work in the canonical C:/work/tirt checkout and preserve existing dirty files.

Coverage: 16 surface/vegetation combinations for thermal; 10 allowed by the
public optical interface. Optical slope/terrain bare/hom/row remain explicitly
unsupported, with no silently substituted scene. Bare scenes use LAI=0 only.

Defaults: vegetation LAI 1/3/5; optical 450/550/670/705/860/1610/2200 nm;
thermal 8/9/10/10.5/11/12 um; solar zenith 30 degrees; spectral view at nadir;
angles -60 to 60 degrees in steps of 10, in both solar principal planes;
angular reference bands 860 nm and 10.5 um. Thermal temperatures are prescribed.
Use illustrative material spectra, not sensor-response integration or measured
validation. Preserve backface omissions as missing rows/status rather than zero.

## Implementation plan

- [x] Add meaningful checks for dispatch, bare limits, angle conventions,
  public thermal parity, unsupported modes and output completeness.
- [x] Implement shared simulation/export/plot functions and 26 explicit scripts.
- [x] Add run-all entrypoint and Chinese usage/coverage documentation.
- [x] Execute all default examples, inspect figures, verify output counts and
  finite values, and exercise an extracted copy from a path containing spaces.
- [x] Review new code and report exact delivery/unsupported scope.

No design confirmation: user explicitly opted out of confirmation prompts.
No worktree: established project preference requires this sole canonical checkout.

## Verification record

9 targeted tests passed; 26 default examples generated 2,138 records and 52 PNGs.
Inspected all figures in six contact sheets. Fresh review found the all-backface
thermal sweep case; fixed and covered by a regression test. All 26 examples were
rerun from a clean extracted directory containing spaces, with an unrelated cwd;
all 52 CSVs exactly matched the reference outputs.

User explicitly requested absolute replay tolerance 1e-6 on 2026-10-24.
Updated scripts/validate_current_model.py value/radiance replay checks and its
receipt accordingly. Replay: 1,830 requested, 1,806 valid rows, maximum absolute
prediction difference 7.555229331046576e-9. Passed at the requested tolerance.
