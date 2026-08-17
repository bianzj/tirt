# TiRT plots

`plot` reads the unified runner's `output.csv` and provides three view-angle plots:

- `polar`: view hemisphere, azimuth clockwise from north, radius is VZA.
- `parallel`: solar principal plane, using relative azimuth 0 and 180 degrees.
- `perpendicular`: perpendicular solar principal plane, using relative azimuth 90 and 270 degrees.

Examples from the project root:

```bash
python -m plot.run_plot --input cases/file_geometry_spectrum/output.csv --mode polar --wavelength 10.5
python -m plot.run_plot --input cases/file_geometry_spectrum/output.csv --mode parallel --wavelength 10.5
python -m plot.run_plot --input cases/file_geometry_spectrum/output.csv --mode perpendicular --wavelength 10.5
python -m plot.run_plot --input cases/file_geometry_spectrum/output.csv --mode all --wavelength 10.5
```

Images are written to `cases/<case>/plots/` by default. The default plotted quantity is `brightness_temperature_C`; use `--quantity brightness_temperature_K`, `--quantity radiance`, or any numeric output column.
