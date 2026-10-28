"""Normalise $$ formulas in the report Markdown to the matplotlib mathtext subset.

Fixes: \\le/\\ge -> \\leq/\\geq, \\iff -> \\Leftrightarrow, \\bmod -> \\operatorname{mod},
\\big| -> |, \\tfrac -> \\frac, and splits multi-line environments (cases/aligned) into
separate display lines so each can be rendered as a standalone formula.
Chinese inside \\text{} is moved out of the math by a table of explicit rewrites.
"""
import re
import sys
from pathlib import Path

BLOCK = re.compile(r"\$\$\s*\n(.*?)\n\$\$", re.S)

REWRITES = {
    r"\left|\,((\varphi_v - \varphi_s + 180)\bmod 360) - 180\,\right|":
        r"\left| (( \varphi_v - \varphi_s + 180 ) \operatorname{mod} 360) - 180 \right|",
    r"\operatorname{overlap}\big((1-c_s)P_r,\ (1-c_v)P_r,\ \tfrac{s_s-s_v}{2},\ P_r\big)":
        r"\operatorname{overlap}\left((1-c_s)P_r,\ (1-c_v)P_r,\ (s_s-s_v)/2,\ P_r\right)",
    r"\mathbf{J} = \left(\mathbf{I} - \mathbf{S}\mathbf{A}^{\mathsf T}\right)^{-1}\mathbf{J}_1":
        r"\mathbf{J} = \left( \mathbf{I} - \mathbf{S} \mathbf{A}^{\mathrm{T}} \right)^{-1} \mathbf{J}_1",
    r"\varepsilon_{eff} = \frac{L\big|_{\text{全部组分 }300\,\mathrm K}}{B(\lambda,300\,\mathrm K)}":
        r"\varepsilon_{eff} = \frac{L_{300}}{B(\lambda, 300\,\mathrm{K})}",
    r"r_{max} = 400\ \text{s m}^{-1}\ (\texttt{max\_aerodynamic\_resistance})":
        r"r_{max} = 400\ \mathrm{s\,m^{-1}}",
    r"\frac{\left|c\right|\cdot10^{-4}}{u_*^{3}}":
        r"\frac{|c| \cdot 10^{-4}}{u_*^{3}}",
    r"\left|x - x_{target}\right| \le 10^{-7} + 10^{-6}\left|x_{target}\right| + \frac{|c| \cdot 10^{-4}}{u_*^{3}}":
        r"\left|x - x_{target}\right| \leq 10^{-7} + 10^{-6}\left|x_{target}\right| + \frac{|c| \cdot 10^{-4}}{u_*^{3}}",
    r"\cos i = \cos\beta\cos\theta + \sin\beta\sin\theta\cos(\varphi - A),\qquad \text{可见/受照} \iff \cos i \ge 0":
        r"\cos i = \cos\beta\cos\theta + \sin\beta\sin\theta\cos(\varphi - A)",
    r"P_{ss} = P_v\exp\left[-n\,(J_s - K)\right],\qquad P_s = P_{ss}\big|_{\theta_v = \theta_s}":
        r"P_{ss} = P_v\exp\left[-n\,(J_s - K)\right],\qquad P_s = \left. P_{ss} \right|_{\theta_v = \theta_s}",
}

# Multi-line environments: each line becomes its own display formula.
CASES = {
    r"\psi_m = \begin{cases}": None,
}


def fix_block(text):
    for old, new in REWRITES.items():
        text = text.replace(old, new)
    text = re.sub(r"\\begin\{aligned\}|\\end\{aligned\}", "", text)
    text = re.sub(r"\\begin\{cases\}|\\end\{cases\}", "", text)
    text = re.sub(r"\\le(?![a-zA-Z])", r"\\leq", text)
    text = re.sub(r"\\ge(?![a-zA-Z])", r"\\geq", text)
    text = text.replace(r"\iff", r"\Leftrightarrow")
    text = text.replace(r"\bmod", r"\operatorname{mod}")
    text = text.replace(r"\tfrac", r"\frac")
    text = text.replace(r"\texttt", r"\mathrm")
    text = re.sub(r"\\big\|", "|", text)
    # unbraced \mathrm X is rejected by mathtext; brace single letters
    text = re.sub(r"\\(mathrm|text|mathbf|mathsf|operatorname)\s+([A-Za-z])(?![A-Za-z])", r"\\\1{\2}", text)
    # end-of-line \\ in split environments is meaningless once lines are separate
    text = re.sub(r"\\\\\s*$", "", text, flags=re.M)
    text = re.sub(r"\\\\\[[0-9.]+pt\]", "", text)
    lines = []
    for line in text.splitlines():
        for segment in line.split(r"\\"):
            stripped = segment.strip()
            stripped = re.sub(r"\\text\{[^}]*[\u4e00-\u9fff][^}]*\}", "", stripped)
            stripped = re.sub(r",?\s*\\(qquad|quad)\s*$", "", stripped)
            stripped = re.sub(r"^\\(qquad|quad)\s*", "", stripped)
            stripped = re.sub(r",?\s*\\(qquad|quad)\s*", r",\\ ", stripped)
            stripped = re.sub(r",?\s*&\s*", r"\\quad ", stripped)
            stripped = re.sub(r",\s*$", "", stripped)
            if stripped:
                lines.append(stripped)
    return "\n".join(lines)


def main(paths):
    for name in paths:
        path = Path(name)
        text = path.read_text(encoding="utf-8")
        out, last = [], 0
        for match in BLOCK.finditer(text):
            out.append(text[last:match.start()])
            out.append("$$\n" + fix_block(match.group(1)) + "\n$$")
            last = match.end()
        out.append(text[last:])
        path.write_text("".join(out), encoding="utf-8")
        print("normalised", path)


if __name__ == "__main__":
    main(sys.argv[1:])
