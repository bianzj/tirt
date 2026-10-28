# Shortwave reciprocity/conservation correction, 2026-10-26

> Reverted on 2026-10-27: `base/reciprocal_optical.py` was removed and the optical first order again uses the thermal four components. This file documents the withdrawn method only; see `docs/TiRT_theory_and_methods.md` sections 4.6-4.7 for the current method.

## Files and method

`base/canopy_closure.py`, new `base/reciprocal_optical.py`, `base/node_radiosity.py`, and `base/optical_terrain_standard.py`.

The old optical first order used the sunlit visible four-component leaf approximation, while higher orders and the energy budget used a different node operator. Replacing the first order by a reciprocal node escape kernel removes that mismatch. No total BRF clipping, direction-swap averaging, or fitted benchmark constants are used.

For node i, p_i(w) is the directional first-hit fraction, A_i its capacity, and q_i(w)=p_i(w) cos(theta) dOmega / integral(p_i cos(theta) dOmega). The separable first escape is c_i p_i(s)p_i(v), with c_i=soil_rho/A_i for soil, and c_i=2(2rho+tau)/(3 A_i) for leaves. Its marginal is exactly the existing node first-escape budget under the common quadrature.

A bounded symmetric hotspot h(s,v)=((1+s dot v)/2)^a uses a=min(512,max(1,1/hspot^2)); hspot<=0 disables it. Center each node's kernel as C_i=h-E_q[h(s,.)]-E_q[h(.,v)]+E_qE_q[h]. The multiplier 1+C_i/2 is nonnegative because h is in [0,1], and has unit marginal. Soil uses it directly. The exponent cap is an explicit angular-resolution floor (hspot below about 0.0442 has the same kernel width); this is a new approximate angular closure, not the former exact joint-gap interpretation of hspot.

Leaf scattering retains the spherical bi-Lambertian phase. For gamma=acos(s dot v), B=2/(3pi)[sin(gamma)+(pi/2-gamma)cos(gamma)], the reflection/transmission kernels are 2[B+cos(gamma)/3] and 2[B-cos(gamma)/3]. This is the analytic spherical-normal integral corresponding to spherical_coefficients/(ks*kv); it avoids the old 64-inclination quadrature. It correctly gives zero first-order backscatter for purely transmitting leaves. Multiply this phase by the bounded hotspot, then symmetrically balance each leaf node: d_i*(K_i@(q_i*d_i))=1. The positive fixed-point iteration must converge below 1e-11 within 200 iterations, otherwise it raises. Off-grid scales are evaluated from the same integral equation, not interpolated. Thus material-dependent phase, nonnegativity, reciprocity and node marginals coexist.

Geometry moments use 32 Gauss-Legendre mu points by 96 full azimuths and are cached separately from the material-dependent balanced phase. Node covered/open soil quadrature now uses the same 32 mu by 24 quadrant directions; row reflection symmetry makes this equivalent to 96 full azimuths. This also changes thermal row multiple scattering slightly (maximum saved-reference model change 0.00116 K). A grazing-row 0/0 transmission-underflow bug was separately reproduced and fixed. Empty crown density now returns bare soil in directional BRF and albedo APIs.

Compound terrain now calls the common local all-order canopy operator for its plane/facet first+multiple terms. Its existing macro single-bounce terrain exchange is retained, so the local canopy proof does not establish global compound-terrain conservation/reciprocity. The legacy view-only Hapke `soil_direct` correction is retained and is explicitly outside the Lambertian guarantee.

## Physical tests

`C:/work/miniconda/python.exe -m pytest tests/test_optical_conservation.py -q --tb=short`

Final retained crown2+phase implementation: 10 tests, 30 subtests passed (39.34 s cold run). The same controls also passed in the rejected crown3 experiment. Tests cover standard/upgrade hom,row,crown reciprocity; direct, mixed and diffuse hemisphere vs independently computed node absorption; nonabsorbing limit; full absorption; shadow; hspot dependence; no first-order transmission backscatter; zero LAI/crown density; grazing rows; zero/vanishing terrain density continuity. Initial tests failed on the audited 0.02423 hom and 0.03275 row energy mismatches, reciprocity failures, grazing NaN, disconnected zero-density conservative singularity, and false pure-transmission backscatter 0.144529.

`optical_controls.py` / `optical_controls.json` record a separate 64 mu by 144 azimuth integral of the retained crown2+phase implementation:

| canopy | exchanged BRF | exchange residual | direct NIR integral minus budget | conservative direct integral minus 1 |
|---|---:|---:|---:|---:|
| hom | 0.3582204054053518 | 0 | +2.34e-10 | +2.80e-10 |
| row | 0.3404862801961447 | 1.11e-16 | -1.40e-5 | -5.51e-6 |
| crown | 0.2121430645886654 | 0 | +2.63e-11 | -1.57e-11 |

Parameters: LAI2,hspot0.1,rho0.45,tau0.40,soil0.25; sun40/130 and view70/280. The old hom swapped pair was 0.3776705316/0.4342328972. Conservative optics are rho=tau0.5,soil1. Diffuse conservative integrals differ from 1 by <=7e-16. The remaining row residual is angular integration error of the piecewise gap geometry, not a percent-level energy excess. Independent test integrals use 48x96, distinct from the 32x96 balancing quadrature.

Hspot still changes backscatter: hom hspot0.05/0.1/0.5 gives 0.44976210/0.44687305/0.40397894; row 0.41780463/0.41518017/0.37812599; crown 0.29762271/0.29550562/0.25967320.

An intermediate whole-project pytest run reported 37 passed and one failure in `test_temperature_and_stability_converge_together[0.0-250.0-15.0-1]` while the root agent was implementing stability. This was reported immediately; the root agent owns the final integrated suite, and no repeated broad runs were launched here.

## Frozen HiSTREAM comparison and limitation

`compare_optical.py` writes `optical_comparison/{all_results.csv,reference_comparison.csv,accuracy_before_after.csv,hotspot_accuracy_before_after.csv,verification.json}`. It evaluates all 1830 requests, with 1806 valid saved HiSTREAM pairs, 10 scenes and 3 bands. All frozen input/reference checksums are verified both before and after; no ray tracing or reference rewriting occurs. The retained phase+crown2 comparison took 13.70 s; the earlier audit replay took 4.26 s, and the simpler phase-free correction took 6.11 s. These are separate cold runs during concurrent development, not a controlled performance benchmark.

Physical constraints improved, but frozen optical accuracy has materially regressed; do not describe this as an accuracy upgrade:

| scene | nm | old RMSE | retained RMSE |
|---|---:|---:|---:|
| plane hom | 670 | .001818 | .006529 |
| plane hom | 860 | .007458 | .014209 |
| plane crown | 670 | .002484 | .010032 |
| plane crown | 860 | .021405 | .052918 |
| row az0 | 860 | .014513 | .027512 |
| row az90 | 860 | .011644 | .014759 |
| slope crown | 860 | .020652 | .063851 |
| terrain crown | 860 | .024082 | .041896 |
| urban crown | 860 | .027558 | .044290 |

Phase restoration improved hom NIR RMSE from .025033 to .014209 and row-az90 from .026419 to .014759 relative to the initial coherent node-only phase closure; it did not remove crown bias. Fixed node escape marginals limit how much angular redistribution alone can recover. Existing crown node transport mixes coarse within-crown and open-ground budgets; replacing that transport requires independent geometric/transport validation. The largest final optical change from the frozen predictions is .131999 BRF. This is an explicit unresolved applicability/accuracy limitation, not a hidden rewritten baseline. Hotspot-region (angular separation <=5deg) scores are included separately in the CSV.

A physically motivated crown3 covered/open soil split was investigated, then rejected as the default: crown860 RMSE improved .052918 to .034348, but crown670 worsened .010032 to .022731, and shared thermal crown RMSE worsened .176419 to .407881 K. The global NODES crown mapping is restored to crown2. Experiment results remain in `optical_comparison_crown3_phase`; intermediate phase-only crown2 and initial phase-free outputs remain in `optical_comparison_crown2_phase` and `optical_comparison_prephase`. No parameter was fitted to the HiSTREAM curves.
