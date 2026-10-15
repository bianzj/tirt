"""Recompute the published ten-scene TIRT/HiSTREAM validation from bundled inputs.

python scripts/validate_current_model.py --plots
No external HiSTREAM executable, private case directory or network is needed.
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / 'docs/validation/current_model_20261015'
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'benchmark/streamsim'))
import run
import run_benchmark
from base.optical import reflectance
from base.radiation import directional_reflectance
from base.physicsF import inv_planck


def digest(path):
    # A Git checkout may use LF or CRLF; the model text is the same.
    return hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def optical(settings, config, wavelength, vza, vaa):
    if config['surface_model'] == 'slope':
        frame = run_benchmark.slope_frame(float(config['slope_angle']), float(config['slope_aspect']))
        view, view_azimuth = run_benchmark.to_local(vza, vaa, frame)
        solar, solar_azimuth = run_benchmark.to_local(30., 0., frame)
        return reflectance(dict(settings, surface_model='plane'), wavelength,
                           view, view_azimuth, solar, solar_azimuth)
    if config['surface_model'] == 'terrain':
        i = settings['wavelengths'].index(wavelength)
        rho, tau, soil = (settings[key][i] for key in ('leaf_reflectance_spectrum',
                         'leaf_transmittance_spectrum', 'soil_reflectance_spectrum'))
        return directional_reflectance(settings, 30., 0., vza, vaa, rho, tau, soil,
                                       use_hapke=False, diffuse_fraction=0.)
    return reflectance(settings, wavelength, vza, vaa, 30., 0.)


def plots(data, scenes, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.sans-serif': ['Microsoft YaHei', 'DejaVu Sans'],
                         'axes.unicode_minus': False})
    folder = output / 'figures'
    folder.mkdir(exist_ok=True)
    axes = []
    for scene in scenes:
        for wave in (10500, 670, 860):
            frame = data[(data.scene == scene['stem']) & (data.wavelength_nm == wave)]
            fig, panels = plt.subplots(1, 2, figsize=(11, 4), sharex=True, sharey=True)
            all_values = []
            for panel, (positive, negative), title in zip(panels, [(0., 180.), (90., 270.)],
                                                          ['太阳主平面', '垂直太阳主平面']):
                for key, label, style in [('value', '当前 TIRT', '-'), ('reference', 'HiSTREAM', 'o-')]:
                    ordered = pd.concat([frame[(frame.vaa == negative) & (frame.vza > 0)].sort_values('vza', ascending=False).assign(sign=-1),
                                         frame[frame.vza == 0].assign(sign=1),
                                         frame[(frame.vaa == positive) & (frame.vza > 0)].sort_values('vza').assign(sign=1)])
                    y = ordered[key] - (273.15 if wave == 10500 else 0.)
                    panel.plot(ordered.vza * ordered.sign, y, style, label=label, lw=1.5, ms=3)
                    all_values.extend(y[np.isfinite(y)])
                panel.set(title=title, xlabel='带符号观测天顶角 (°)', xlim=(-80, 80))
                panel.grid(alpha=.3)
            low, high = min(all_values), max(all_values)
            padding = max(.06 * (high - low), .25 if wave == 10500 else .005)
            panels[0].set_ylim(low-padding if wave == 10500 else max(0., low-padding), high+padding)
            panels[0].set_ylabel('亮温 (°C)' if wave == 10500 else '反射因子 BRF')
            panels[0].legend()
            assert panels[0].get_ylim() == panels[1].get_ylim()
            axes.append(dict(scene=scene['stem'], wavelength_nm=wave, shared_ylim=list(panels[0].get_ylim())))
            fig.suptitle(f"{scene['label']} | {wave / 1000:g} μm | SZA 30°, SAA 0°")
            fig.tight_layout()
            for suffix in ('png', 'pdf'):
                fig.savefig(folder / f"{scene['stem']}_{wave}.{suffix}", dpi=150)
            plt.close(fig)
    return axes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, default=BUNDLE)
    parser.add_argument('--output', type=Path, default=ROOT / 'cases/_current_model_validation')
    parser.add_argument('--plots', action='store_true')
    args = parser.parse_args()
    bundle, output = args.bundle.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((bundle / 'manifest.json').read_text(encoding='utf-8'))
    for name, expected in manifest['bundle_sha256'].items():
        assert digest(bundle / name) == expected, f'Changed bundled input/reference: {name}'
    expected = pd.read_csv(bundle / 'all_results.csv')
    inputs = {scene['stem']: json.loads((bundle / 'inputs' / f"{scene['stem']}.json").read_text(encoding='utf-8'))
              for scene in manifest['scenes']}
    rows = []
    start = time.perf_counter()
    for row in expected.itertuples():
        inp = inputs[row.scene]
        config, settings = inp['thermal_config'], inp['optical_settings']
        geometry = (row.vza, row.vaa, row.sza, row.saa)
        hidden = config['surface_model'] == 'slope' and run._slope_backface(config, geometry)
        result, radiation = np.nan, np.nan
        if not hidden:
            if row.wavelength_nm == 10500:
                radiation, result = run._surface_result(config, config['surface_model'],
                    config['vegetation_model'], 10.5, geometry, inp['emissivity'])
                np.testing.assert_allclose(inv_planck(10.5, radiation), result, atol=1e-9, rtol=0)
            else:
                result = optical(settings, config, row.wavelength_nm, row.vza, row.vaa)
            assert np.isfinite(result) and result >= 0
        rows.append(dict(scene=row.scene, wavelength_nm=row.wavelength_nm, vza=row.vza, vaa=row.vaa,
                         sza=row.sza, saa=row.saa, value=result, radiance=radiation,
                         status='slope_backface' if hidden else 'ok'))
    current = pd.DataFrame(rows)
    assert len(current) == 1830
    assert current.status.equals(expected.status)
    active = current.status == 'ok'
    np.testing.assert_allclose(current.loc[active, 'value'], expected.loc[active, 'value'], rtol=0, atol=1e-10)
    thermal = active & (current.wavelength_nm == 10500)
    np.testing.assert_allclose(current.loc[thermal, 'radiance'], expected.loc[thermal, 'radiance'], rtol=0, atol=1e-10)
    keys = ['scene', 'wavelength_nm', 'vza', 'vaa']
    reference = pd.read_csv(bundle / 'histream_reference.csv')
    paired = current[active].merge(reference, on=keys, how='left', validate='one_to_one')
    assert len(paired) == 1806 and np.isfinite(paired.reference).all()
    paired['difference'] = paired.value - paired.reference
    scores = []
    for (scene, wave), group in paired.groupby(['scene', 'wavelength_nm'], sort=False):
        d = group.difference
        scores.append(dict(scene=scene, wavelength_nm=wave, n=len(d), rmse=float(np.sqrt(np.mean(d*d))),
                           bias=float(d.mean()), max_abs=float(abs(d).max()), unit='K' if wave == 10500 else 'BRF'))
    current.to_csv(output / 'all_results.csv', index=False)
    paired.to_csv(output / 'reference_comparison.csv', index=False)
    pd.DataFrame(scores).to_csv(output / 'accuracy.csv', index=False)
    report = dict(rows=len(current), valid_predictions=int(active.sum()), reference_pairs=len(paired),
                  scene_count=len(inputs), band_count=3, model_recomputed=True, new_raytracing=False,
                  max_abs_prediction_replay_error=float(abs(current.loc[active, 'value']-expected.loc[active, 'value']).max()),
                  prediction_tolerance=1e-10, seconds=time.perf_counter()-start,
                  bundle_checksums_verified=True,
                  source_sha256={name: digest(ROOT / name) for name in manifest['model_source_sha256']})
    if args.plots:
        report['axes'] = plots(current.merge(reference, on=keys, how='left', validate='one_to_one'), manifest['scenes'], output)
    (output / 'verification.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k not in ('source_sha256', 'axes')}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
