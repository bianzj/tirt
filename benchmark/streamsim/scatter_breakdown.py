"""Where the scattered radiance comes from in the terrain+forest (TV) and urban+forest (UV) scenes, HiSTREAM vs TiRT.

Radiance is linear in the emitters, so each scattering path is isolated by switching emitters off (3 K) and
reflectors off (emissivity 1):
  within forest        leaves emit, macro surface cold and black   (leaf emission reflected only by leaves)
  forest -> macro      leaves emit, macro cold  minus the line above  (leaf emission reflected at least once by the macro surface)
  within macro         macro emits, leaves cold and black           (macro emission reflected only by the macro surface)
  macro -> forest      macro emits, leaves cold minus the line above  (macro emission reflected at least once by leaves)
HiSTREAM runs come from run_benchmark.py (tags _srcF/_srcF_Mblack/_srcM/_srcM_Fblack with _d1/_d2/_d4); depth 2 is
single scattering. TiRT is evaluated with the same switches; its "emission only" is TiRT with the scattering terms zeroed.
Each radiance difference is expressed in kelvin with dB/dT at the full scene's brightness temperature of that view.

Writes runs/compare/scatter_breakdown.csv and runs/compare/scatter/<scene>.png.
"""
import csv
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNS = HERE / "runs"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))
warnings.simplefilter("ignore")

import run as R  # noqa: E402
import run_benchmark as rb  # noqa: E402
import tirt_side  # noqa: E402
import base.terrain_veg as terrain_veg_module  # noqa: E402
import base.urban_veg as urban_veg_module  # noqa: E402
from base.inputio import read_input  # noqa: E402
from base.terrain_veg import Terrain_Veg  # noqa: E402
from base.urban_veg import Urban_Veg  # noqa: E402

SCENES = {
    "TV_terrain_crown": {"label": "复合地形 + 森林", "macro_emis": ("terrain_e",)},
    "UV_urban_crown": {"label": "城市 + 森林", "macro_emis": ("street_e", "wall_e", "roof_e")},
}
CATEGORIES = [("forest", "森林内部（叶片之间）"), ("forest_to_macro", "森林发射 → 宏观地表反射"),
              ("macro", "宏观地表内部"), ("macro_to_forest", "宏观地表发射 → 森林反射")]


def histream(test, tag):
    """Mean radiance per view straight from the images, counting zero pixels (cold and black surfaces render as 0)."""
    values = {}
    for row in rb.reduce_images(RUNS / f"{test}_eRaytracing{tag}" / "output", keep_zero=True):
        for vaa in ((0.0, 90.0, 180.0, 270.0) if row["vza"] == 0 else (row["vaa"] % 360,)):
            values[(row["vza"], vaa)] = row["radiance"]
    return values


def tirt(test, views, cold=(), black=(), scatter=True):
    """TiRT radiance per view with some sources switched off and some reflectors made black.

    TiRT's planck() reads temperatures below 100 as Celsius, so sources are switched off by dropping their
    Planck term from sum(effective emissivity x B(T)) rather than by setting them to 3 K."""
    spec = SCENES[test]
    config = dict(read_input(R.DEFAULT_INPUT))
    config.update({key: str(value) for key, value in {**tirt_side.COMMON, **tirt_side.TESTS[test]}.items()})
    spectra = R._spectra(config, R.DEFAULT_INPUT, np.array([10.5]))
    band = {key: values[0] for key, values in spectra.items()}
    for name in black:
        for key in (spec["macro_emis"] if name == "macro" else ("leaf_e",)):
            band[key] = 1.0
    if test.startswith("TV"):
        model_class, surface = Terrain_Veg, "terrain"
        temperatures = np.array(R._vegetation_temperatures(config))
        is_leaf = np.array([False, False, True, True])
    else:
        model_class, surface = Urban_Veg, "urban"
        leaf_sun, leaf_shade = R._vegetation_temperatures(config)[2:]
        temperatures = np.array([R._temperature(config, f"{part}_{state}", 45.0 if state == "sunlit" else 30.0)
                                 for part in ("wall", "street", "roof") for state in ("sunlit", "shaded")] + [0.5 * (leaf_sun + leaf_shade)])
        is_leaf = np.array([False] * 6 + [True])
    emitting = ~is_leaf if "leaf" in cold else is_leaf if "macro" in cold else np.ones(is_leaf.size, bool)
    captured = []
    original_effective = model_class.calculate_effective_component_emissivity

    def capture(self, ifP=0):
        result = original_effective(self, ifP)
        captured.append(np.asarray(result, dtype=float).reshape(-1))
        return result

    originals = Terrain_Veg.calculate_component_scatter_emissivity, Urban_Veg.calculate_scattering_emissivity
    forest_background = terrain_veg_module.emissivity_scatter_analytical, urban_veg_module.emissivity_scatter_analytical
    model_class.calculate_effective_component_emissivity = capture
    if not scatter:
        Terrain_Veg.calculate_component_scatter_emissivity = lambda self, kangle, ifP=0: (0.0, 0.0, 0.0, 0.0)
        Urban_Veg.calculate_scattering_emissivity = lambda self, kangle, ifP=0: (0.0,) * 7
        no_scatter = lambda lai, el, es, vza, sza: (np.zeros((4, np.size(vza))), np.zeros(np.size(vza)))
        terrain_veg_module.emissivity_scatter_analytical = urban_veg_module.emissivity_scatter_analytical = no_scatter
    try:
        radiance = {}
        for view in views:
            R._surface_result(config, surface, "crown", 10.5, (view[0], view[1], 30.0, 0.0), band)
            emissivity = captured[-1]
            planck = np.asarray(R.planck(10.5, temperatures[:emissivity.size]), dtype=float)
            radiance[view] = float(np.sum((emissivity * planck)[emitting[:emissivity.size]]))
        return radiance
    finally:
        model_class.calculate_effective_component_emissivity = original_effective
        Terrain_Veg.calculate_component_scatter_emissivity, Urban_Veg.calculate_scattering_emissivity = originals
        terrain_veg_module.emissivity_scatter_analytical, urban_veg_module.emissivity_scatter_analytical = forest_background


def kelvin(delta, full):
    """Radiance differences per view -> kelvin, with dB/dT at the brightness temperature of the full scene."""
    out = {}
    for view, value in delta.items():
        t = float(R.inv_planck(10.5, full[view]))
        slope = (R.planck(10.5, t + 0.05) - R.planck(10.5, t - 0.05)) / 0.1
        out[view] = value / float(slope)
    return out


def paths(src_f, src_f_mblack, src_m, src_m_fblack, emission_f, emission_m):
    keys = src_f.keys()
    return {
        "forest": {k: src_f_mblack[k] - emission_f[k] for k in keys},
        "forest_to_macro": {k: src_f[k] - src_f_mblack[k] for k in keys},
        "macro": {k: src_m_fblack[k] - emission_m[k] for k in keys},
        "macro_to_forest": {k: src_m[k] - src_m_fblack[k] for k in keys},
    }


def main():
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    out_dir = RUNS / "compare" / "scatter"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for test, spec in SCENES.items():
        full = histream(test, "")
        views = sorted(full)
        h = {depth: paths(*(histream(test, f"{tag}_d{depth}") for tag in ("_srcF", "_srcF_Mblack", "_srcM", "_srcM_Fblack")),
                          histream(test, "_srcF_d1"), histream(test, "_srcM_d1")) for depth in (2, 4)}
        for depth, full_tag in ((1, "_depth1"), (4, "")):
            reference = histream(test, full_tag)
            source_f, source_m = histream(test, f"_srcF_d{depth}"), histream(test, f"_srcM_d{depth}")
            check = max(abs(source_f[k] + source_m[k] - reference[k]) / reference[k] for k in views)
            print(f"{test}: depth {depth}: forest-only + macro-only vs full scene, max relative error {check:.2e}")

        t_full = tirt(test, views)
        t = paths(tirt(test, views, cold=("macro",)), tirt(test, views, cold=("macro",), black=("macro",)),
                  tirt(test, views, cold=("leaf",)), tirt(test, views, cold=("leaf",), black=("leaf",)),
                  tirt(test, views, cold=("macro",), scatter=False), tirt(test, views, cold=("leaf",), scatter=False))

        results = {}
        for key, label in CATEGORIES:
            single = kelvin(h[2][key], full)
            total = kelvin(h[4][key], full)
            model = kelvin(t[key], t_full)
            results[key] = (np.mean(list(single.values())), np.mean(list(total.values())), np.mean(list(model.values())))
            rows.append({"scene": spec["label"], "test": test, "path": label,
                         "histream_single_K": results[key][0], "histream_all_orders_K": results[key][1], "tirt_K": results[key][2],
                         "histream_all_orders_nadir_K": total[(0.0, 0.0)], "tirt_nadir_K": model[(0.0, 0.0)]})
        rows.append({"scene": spec["label"], "test": test, "path": "合计",
                     **{column: sum(row[column] for row in rows[-4:]) for column in
                        ("histream_single_K", "histream_all_orders_K", "tirt_K", "histream_all_orders_nadir_K", "tirt_nadir_K")}})

        fig, ax = plt.subplots(figsize=(9, 4.2))
        x = np.arange(len(CATEGORIES))
        for offset, (index, label, color) in zip((-0.27, 0.0, 0.27), ((0, "HiSTREAM 单次散射", "#f5b971"), (1, "HiSTREAM 全部阶次", "tab:orange"),
                                                                      (2, "TiRT", "tab:blue"))):
            values = [results[key][index] for key, _ in CATEGORIES]
            ax.bar(x + offset, values, width=0.26, color=color, label=label)
            for xi, value in zip(x + offset, values):
                ax.text(xi, value + 0.01, f"{value:.2f}", ha="center", va="bottom", fontsize=8)
        ax.set_xticks(x, [label for _, label in CATEGORIES], fontsize=9)
        ax.set_ylabel("对亮温的贡献 (K，所有观测角平均)")
        ax.set_title(f"{spec['label']}：散射贡献按路径分解（10.5 µm，SZA 30°）", fontsize=11)
        ax.grid(axis="y", alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out_dir / f"{test}.png", dpi=110)
        plt.close(fig)

    with (RUNS / "compare" / "scatter_breakdown.csv").open("w", newline="", encoding="utf-8") as sink:
        writer = csv.DictWriter(sink, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows({k: (f"{v:.3f}" if isinstance(v, float) else v) for k, v in row.items()} for row in rows)
    for row in rows:
        print(f"{row['test']:17s} {row['path']:14s} HiSTREAM single {row['histream_single_K']:+.3f}  all {row['histream_all_orders_K']:+.3f}"
              f"  TiRT {row['tirt_K']:+.3f} K   (nadir {row['histream_all_orders_nadir_K']:+.3f} / {row['tirt_nadir_K']:+.3f})")


if __name__ == "__main__":
    main()
