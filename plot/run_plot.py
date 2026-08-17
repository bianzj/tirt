import argparse
from pathlib import Path

from .plots import load_output, plot_all, plot_parallel, plot_perpendicular, plot_polar


def main():
    parser = argparse.ArgumentParser(description="Plot TiRT output.csv angular results")
    parser.add_argument("--input", type=Path, required=True, help="path to output.csv")
    parser.add_argument("--mode", choices=("polar", "parallel", "perpendicular", "all"), default="all")
    parser.add_argument("--wavelength", type=float, help="thermal infrared wavelength in um")
    parser.add_argument("--quantity", default="brightness_temperature_C", help="output column to plot")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--angle-tolerance", type=float, default=1.0e-6)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    rows, selected, quantity = load_output(args.input, args.wavelength, args.quantity)
    output_dir = args.output_dir or args.input.parent / "plots"
    tag = f"{selected:g}".replace(".", "p")
    if args.mode == "all":
        result = plot_all(args.input, output_dir, selected, quantity, args.angle_tolerance, args.show)
    else:
        functions = {"polar": plot_polar, "parallel": plot_parallel, "perpendicular": plot_perpendicular}
        path = output_dir / f"{args.mode}_{tag}.png"
        result = functions[args.mode](rows, path, quantity, tolerance=args.angle_tolerance, show=args.show) if args.mode != "polar" else functions[args.mode](rows, path, quantity, show=args.show)
    if isinstance(result, dict):
        for name, path in result.items():
            if path is not None:
                print(f"{name}: {path}")
    else:
        print(result)


if __name__ == "__main__":
    main()
