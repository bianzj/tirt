"""Engine provenance: runs/engine_manifest.csv (engine sha256 per case) and
runs/old_vs_new.md (difference of every result CSV against runs/_old_engine/)."""

import csv
import json
from pathlib import Path

RUNS = Path(__file__).resolve().parent / "runs"
OLD = RUNS / "_old_engine"


def table(path, value_columns):
    with path.open(encoding="utf-8") as source:
        rows = list(csv.DictReader(source))
    columns = [c for c in value_columns if rows and c in rows[0]]
    return {(float(r["vza"]), float(r["vaa"])): [float(r[c]) for c in columns] for r in rows}, columns


def main():
    with (RUNS / "engine_manifest.csv").open("w", newline="", encoding="utf-8") as sink:
        writer = csv.writer(sink)
        writer.writerow(["case", "engine_sha256", "engine_path"])
        for case in sorted(RUNS.glob("*/case.json")):
            engine = json.loads(case.read_text(encoding="utf-8")).get("engine", {})
            writer.writerow([case.parent.name, engine.get("sha256", ""), engine.get("path", "")])
    lines = ["| result CSV | quantity | mean | rms | max abs |", "|---|---|---|---|---|"]
    for old in sorted(OLD.glob("*.csv")):
        new = RUNS / old.name
        if not new.exists() or "_frac_" in old.name:
            continue
        with old.open(encoding="utf-8") as source:
            if "vza" not in (csv.DictReader(source).fieldnames or []):
                continue
        quantity = ["bt_K"] if "fractions" not in old.name else ["f_soil_sunlit", "f_soil_shaded", "f_leaf_sunlit", "f_leaf_shaded"]
        a, columns = table(new, quantity)
        b, _ = table(old, quantity)
        diffs = [x - y for key in a if key in b for x, y in zip(a[key], b[key])]
        if not diffs:
            continue
        rms = (sum(d * d for d in diffs) / len(diffs)) ** 0.5
        unit = "K" if columns == ["bt_K"] else "fraction"
        lines.append(f"| {old.name} | {unit} | {sum(diffs) / len(diffs):+.4f} | {rms:.4f} | {max(abs(d) for d in diffs):.4f} |")
    (RUNS / "old_vs_new.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
