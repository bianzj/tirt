"""Multi-angle / multi-band runner for energy-balance inputs (input_tirteb.csv)."""

import argparse
from pathlib import Path

from run import OUTPUT_FIELDS, ROOT, run

__all__ = ("DEFAULT_INPUT", "OUTPUT_FIELDS", "run_angle")

DEFAULT_INPUT = ROOT / "input_tirteb.csv"


def run_angle(input_path=DEFAULT_INPUT, output_path=None, model=None, **overrides):
    return run(input_path, output_path, model, **overrides)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--model", choices=("bare", "hom", "row", "crown"))
    args = parser.parse_args()
    print(run_angle(args.input, args.output, args.model))


if __name__ == "__main__":
    main()
