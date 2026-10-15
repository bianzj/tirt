"""Regression checks for replayable inputs and explicit input errors."""

import csv
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from base.inputio import read_geometry, read_input, read_spectrum, read_wavelengths, write_input_snapshot
from run import _geometry, _spectra, run


def radiances(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return [float(row["radiance"]) for row in csv.DictReader(stream)]


class InputSnapshotTests(unittest.TestCase):
    def test_override_snapshot_replays_actual_model_and_parameters(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            output = run(ROOT / "input.csv", work / "original" / "output.csv",
                         model="row", lai=2.2, row_width=0.8, geometry_mode=0,
                         vza="20;40", vaa="30")
            snapshot = output.with_name("input.csv")
            config = read_input(snapshot)
            self.assertEqual(config["vegetation_model"], "row")
            self.assertEqual(config["lai"], "2.2")
            replay = run(snapshot, work / "replay" / "output.csv")
            np.testing.assert_allclose(radiances(output), radiances(replay), rtol=1e-12)

    def test_snapshot_freezes_uploaded_files(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            spectrum = work / "soil.txt"
            spectrum.write_text("8 .7\n10.5 .8\n12 .9\n", encoding="utf-8")
            output = run(ROOT / "input.csv", work / "original" / "output.csv",
                         vegetation_model="bare", geometry_mode=0, vza="30", vaa="0",
                         soil_spectrum_source="file", soil_spectrum_file=str(spectrum))
            spectrum.write_text("8 .1\n10.5 .1\n12 .1\n", encoding="utf-8")
            snapshot = output.with_name("input.csv")
            self.assertFalse(Path(read_input(snapshot)["soil_spectrum_file"]).is_absolute())
            replay = run(snapshot, work / "replay" / "output.csv")
            np.testing.assert_allclose(radiances(output), radiances(replay), rtol=1e-12)

    def test_hemisphere_and_additional_file_are_both_used_and_packaged(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            extra = work / "extra.txt"
            extra.write_text("25 12\n", encoding="utf-8")
            config = read_input(ROOT / "input.csv")
            config.update(geometry_mode="3", include_principal="0", include_hemisphere="1",
                          principal_vza_max="30", additional_geometry_source="file",
                          additional_geometry_file=str(extra))
            expected = _geometry(config, ROOT / "input.csv")
            self.assertIn((25.0, 12.0, 30.0, 0.0), expected)
            self.assertGreater(len(expected), 1)
            snapshot = write_input_snapshot(config, ROOT / "input.csv", work / "bundle" / "input.csv")
            self.assertEqual(expected, _geometry(read_input(snapshot), snapshot))

    def test_nonfinite_and_backside_angles_are_rejected(self):
        config = read_input(ROOT / "input.csv")
        for changes in ({"vza": "nan"}, {"vza": "110"}, {"sza": "100"},
                        {"geometry_mode": "2.5"}, {"geometry_mode": "3", "principal_vza_max": "100"}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                _geometry({**config, **changes}, ROOT / "input.csv")
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "output.csv"
            output.write_text("previous result\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "未找到可见的观测方向"):
                run(ROOT / "input.csv", output, surface_model="slope", slope_angle=60,
                    slope_aspect=0, geometry_mode=0, vza=60, vaa=180)
            self.assertEqual(output.read_text(encoding="utf-8"), "previous result\n")

    def test_scoped_sources_and_overrides_survive_replay(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            (work / "soil.txt").write_text("10.5 .7\n", encoding="utf-8")
            source = work / "source.csv"
            for soil_source, changes in (("file", {}), ("direct", {"soil_spectrum_source": "file"}),
                                         ("direct", {"soil_spectrum_source": " file "})):
                with self.subTest(soil_source=soil_source, changes=changes):
                    source.write_text(
                        "surface,vegetation_model,bare,\n"
                        "spectral,wavelengths,10.5,\n"
                        f"soil,spectrum_source,{soil_source},\n"
                        "soil,spectrum_file,soil.txt,\n"
                        "vegetation,spectrum_source,direct,\n", encoding="utf-8")
                    output = run(source, work / "original" / "output.csv", **changes)
                    snapshot = output.with_name("input.csv")
                    count = len(read_input(snapshot))
                    replay = run(snapshot, work / "replay" / "output.csv")
                    np.testing.assert_allclose(radiances(output), radiances(replay), rtol=1e-12)
                    self.assertEqual(count, len(read_input(replay.with_name("input.csv"))))

    def test_moved_snapshot_defaults_to_its_own_result_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            output = run(ROOT / "input.csv", work / "original" / "output.csv",
                         vegetation_model="bare", geometry_mode=0, vza="20", vaa="0")
            moved = work / "moved"
            shutil.copytree(output.parent, moved)
            replay = run(moved / "input.csv")
            self.assertEqual(replay.resolve(), (moved / "output.csv").resolve())
            np.testing.assert_allclose(radiances(output), radiances(replay), rtol=1e-12)

    def test_nonfinite_wavelength_files_and_invalid_emissivity_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bands.txt"
            path.write_text("10.5\nnan\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                read_wavelengths(path)
        config = read_input(ROOT / "input.csv")
        config["emissivity_soil"] = "0.9;1.1;0.9"
        with self.assertRaises(ValueError):
            _spectra(config, ROOT / "input.csv", np.array([8.0, 10.5, 12.0]))

    def test_invalid_numeric_file_rows_are_not_silently_dropped(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.txt"
            for reader, contents in (
                (read_geometry, "vza vaa\n10 20\n30 wrong\n40 50\n"),
                (read_spectrum, "wavelength emissivity\n8 .8\n10.5 wrong\n12 .9\n"),
            ):
                with self.subTest(reader=reader.__name__):
                    path.write_text(contents, encoding="utf-8")
                    with self.assertRaisesRegex(ValueError, ":3"):
                        reader(path)


if __name__ == "__main__":
    unittest.main()
