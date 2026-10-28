# TiRT / TIRTEB physics corrections

Goal: correct the confirmed 2026-10-26 audit failures, with coupled non-neutral Monin–Obukhov iteration explicitly required by the user.
Specification: tmp/audit_20261026/审查报告.md and the user's instruction to implement, without additional confirmation.
Constraints: work only in C:/work/tirt; preserve existing dirty changes and frozen validation references; no StreamSim changes or launches; no release packaging or commits.

- [x] Radiation: common reciprocal conservative node operator, phase and bounded hotspot; 10 targeted tests. Coarse-crown accuracy limitations explicitly retained in report.
- [x] Thermal history: physical tests and consistent heat capacities, implicit diffusion, 48-hour actual-interval history and first-step dtime.
- [x] Physiology/spectra: separate leaf temperature/PAR/bulk resistance and broadband energy emissivity independent of output sampling.
- [x] Urban: reciprocal exchange, isothermal and ground-area controls, full street soil plus leaf budgets.
- [x] Terrain: slope projection and statistical compound facets with conservative exposure/exchange; four shared temperatures and area conventions documented.
- [x] Stability: both signs, neutral and warm starts; inverse-L bracket plus bounded temperature/gustiness and joint inverse-L fallback; final common closure checked.
- [x] Integrate: 72 tests+39 subtests, 144/144 time steps, default CSV48/48, all16 parity combinations. Saved HiSTREAM comparison preserves original hashes. Corrected method docs and durable receipt under docs/validation/physics_corrections_20261026.

Review focus: zero LAI/building density; equivalent rotated geometry; variable meteo interval and first-step dtime; stable low-wind nights; conservative/high-reflectance optical limits.
