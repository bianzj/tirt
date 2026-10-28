"""Render every $$...$$ display formula of the given Markdown files to PNG (matplotlib mathtext).

Run with a Python that has matplotlib. Images go to docs/build/formulas/<sha1>.png; md_to_docx.py looks them up by hash.
"""
import hashlib
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["mathtext.fontset"] = "cm"
OUT = Path(__file__).resolve().parent / "build" / "formulas"
BLOCK = re.compile(r"^\$\$\s*\n(.*?)\n\$\$\s*$", re.S | re.M)


def key(tex):
    return hashlib.sha1(tex.strip().encode("utf-8")).hexdigest()[:16]


def lines_of(tex):
    return [line.strip() for line in tex.strip().splitlines() if line.strip()]


def render(line, path):
    fig = plt.figure(figsize=(0.01, 0.01))
    fig.text(0, 0, f"${line}$", fontsize=13)
    fig.savefig(path, dpi=220, bbox_inches="tight", pad_inches=0.05, facecolor="white")
    plt.close(fig)


def main(files):
    OUT.mkdir(parents=True, exist_ok=True)
    count = errors = 0
    for name in files:
        text = Path(name).read_text(encoding="utf-8-sig")
        for match in BLOCK.finditer(text):
            for line in lines_of(match.group(1)):
                path = OUT / f"{key(line)}.png"
                if path.exists():
                    continue
                try:
                    render(line, path)
                    count += 1
                except Exception as error:  # report every failing formula, keep going
                    errors += 1
                    print(f"{name}: formula failed: {error}\n{line}\n")
    print(f"rendered {count}, failed {errors}")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
