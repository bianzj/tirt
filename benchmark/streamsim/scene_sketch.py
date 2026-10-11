"""3D scene sketches for the benchmark figures: geometry of each test, the sun direction and the two observation planes.

Axes: x = east, y = north, z = up. The sun is at SZA 30, SAA 0 (north), so the solar principal plane is the
north-south plane and the perpendicular plane the east-west plane.
"""
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

SZA, SAA = 30.0, 0.0
GROUND, LEAF, TRUNK, ROOF, WALL, TERRAIN = "#d9c7a3", "#4f9a3a", "#7a5230", "#b0b0b0", "#d8d8d8", "#c9b48a"


def _plane(ax, half, z=0.0, color=GROUND):
    xs = np.array([[-half, half], [-half, half]])
    ys = np.array([[-half, -half], [half, half]])
    ax.plot_surface(xs, ys, np.full_like(xs, z, dtype=float), color=color, alpha=0.6, linewidth=0, shade=False)


def _sphere(ax, centre, radius, color=LEAF):
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 14), np.linspace(0, np.pi, 8))
    ax.plot_surface(centre[0] + radius * np.cos(u) * np.sin(v), centre[1] + radius * np.sin(u) * np.sin(v),
                    centre[2] + radius * np.cos(v), color=color, linewidth=0, alpha=0.9, shade=True)


def _box(ax, x0, y0, z0, dx, dy, dz, color, alpha=0.9, edge="0.35"):
    x1, y1, z1 = x0 + dx, y0 + dy, z0 + dz
    faces = [[(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0)], [(x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)],
             [(x0, y0, z0), (x1, y0, z0), (x1, y0, z1), (x0, y0, z1)], [(x0, y1, z0), (x1, y1, z0), (x1, y1, z1), (x0, y1, z1)],
             [(x0, y0, z0), (x0, y1, z0), (x0, y1, z1), (x0, y0, z1)], [(x1, y0, z0), (x1, y1, z0), (x1, y1, z1), (x1, y0, z1)]]
    ax.add_collection3d(Poly3DCollection(faces, facecolor=color, edgecolor=edge, linewidths=0.3, alpha=alpha))


def _tree(ax, base, normal, trunk=3.0, radius=3.0):
    base, normal = np.asarray(base, float), np.asarray(normal, float)
    top = base + trunk * normal
    ax.plot(*zip(base, top), color=TRUNK, lw=1.2)
    _sphere(ax, base + (trunk + radius) * normal, radius)


def _cone(ax, cx, cy, height, radius):
    r, a = np.meshgrid(np.linspace(0, radius, 6), np.linspace(0, 2 * np.pi, 28))
    ax.plot_surface(cx + r * np.cos(a), cy + r * np.sin(a), height * (1 - r / radius), color=TERRAIN, linewidth=0, alpha=0.95, shade=True)


def _scatter(rng, count, half, margin, clear=lambda x, y: True):
    points = []
    while len(points) < count:
        x, y = rng.uniform(-half + margin, half - margin, 2)
        if clear(x, y) and all(np.hypot(x - p[0], y - p[1]) > 2 * margin for p in points):
            points.append((x, y))
    return points


def _decorate(ax, half, ztop, zbottom=0.0, ground=lambda x, y: 0.0 * x):
    sun = np.array([np.sin(np.radians(SZA)) * np.sin(np.radians(SAA)), np.sin(np.radians(SZA)) * np.cos(np.radians(SAA)), np.cos(np.radians(SZA))])
    length = 0.55 * (ztop - zbottom)
    target = np.array([0.0, 0.0, zbottom + 0.3 * (ztop - zbottom)])
    tail = target + sun * length
    ax.quiver(*tail, *(-sun * length), color="#e8a200", lw=2.2, arrow_length_ratio=0.18)
    ax.text(*(tail + sun * 0.08 * length), "太阳", color="#b07800", fontsize=8)
    line = np.linspace(-half, half, 20)
    lift = 0.01 * (ztop - zbottom)
    ax.plot(0 * line, line, ground(0 * line, line) + lift, color="tab:red", lw=1.4, ls="--")
    ax.plot(line, 0 * line, ground(line, 0 * line) + lift, color="tab:purple", lw=1.4, ls="--")
    ax.text(0, half * 1.08, float(ground(0.0, half)), "N", color="tab:red", fontsize=8, ha="center")
    ax.text(half * 1.08, 0, float(ground(half, 0.0)), "E", color="tab:purple", fontsize=8, ha="center")
    ax.set_xlim(-half, half)
    ax.set_ylim(-half, half)
    ax.set_zlim(zbottom, ztop)
    ax.set_box_aspect((2 * half, 2 * half, ztop - zbottom), zoom=1.1)
    ax.view_init(elev=24, azim=-58)
    ax.set_axis_off()


def draw(ax, test):
    """Sketch `test` (a benchmark test name) into a 3D axis."""
    rng = np.random.default_rng(3)
    up = (0.0, 0.0, 1.0)
    if test == "T1_plane_bare":
        half = 5.0
        _plane(ax, half)
        _decorate(ax, half, 5.0)
    elif test == "T2_plane_hom":
        half = 4.0
        _plane(ax, half)
        _box(ax, -half, -half, 0, 2 * half, 2 * half, 1.0, LEAF, alpha=0.55, edge=LEAF)
        _decorate(ax, half, 4.5)
    elif test == "T3_plane_crown":
        half = 18.0
        _plane(ax, half)
        for x, y in _scatter(rng, 7, half, 3.2):
            _tree(ax, (x, y, 0), up)
        _decorate(ax, half, 16.0)
    elif test in ("T5_row_az0", "T5_row_az90"):
        half = 3.0
        _plane(ax, half)
        for centre in np.arange(-half + 0.5, half, 1.0):
            if test == "T5_row_az0":
                _box(ax, centre - 0.25, -half, 0, 0.5, 2 * half, 1.0, LEAF, alpha=0.75, edge="#2f6b22")
            else:
                _box(ax, -half, centre - 0.25, 0, 2 * half, 0.5, 1.0, LEAF, alpha=0.75, edge="#2f6b22")
        _decorate(ax, half, 3.6)
    elif test in ("S1_slope_hom", "S1_slope_crown"):
        half = 5.0 if test == "S1_slope_hom" else 16.0
        tilt = np.radians(30.0)
        normal = np.array([np.sin(tilt), 0.0, np.cos(tilt)])  # slope facing east
        xs, ys = np.meshgrid([-half, half], [-half, half])
        zs = -np.tan(tilt) * xs + np.tan(tilt) * half
        ax.plot_surface(xs, ys, zs, color=GROUND, alpha=0.7, linewidth=0, shade=False)
        if test == "S1_slope_hom":
            ax.plot_surface(xs + normal[0], ys, zs + normal[2], color=LEAF, alpha=0.45, linewidth=0, shade=False)
        else:
            for x, y in _scatter(rng, 6, half, 3.4):
                _tree(ax, (x, y, -np.tan(tilt) * x + np.tan(tilt) * half), normal)
        _decorate(ax, half, 2 * np.tan(tilt) * half + (3.0 if test == "S1_slope_hom" else 14.0),
                  ground=lambda x, y: np.tan(tilt) * (half - x) + 0.0 * y)
    elif test == "TV_terrain_crown":
        half = 70.0
        _plane(ax, half, color=TERRAIN)
        cones = [(-30.0, 30.0), (32.0, 15.0), (-12.0, -38.0)]
        for cx, cy in cones:
            _cone(ax, cx, cy, 20.0, 30.0)
        beta = np.arctan(20.0 / 30.0)
        for x, y in _scatter(rng, 26, half, 4.0):
            base, normal = (x, y, 0.0), up
            for cx, cy in cones:
                d = np.hypot(x - cx, y - cy)
                if d < 30.0:
                    base = (x, y, 20.0 * (1 - d / 30.0))
                    normal = (np.sin(beta) * (x - cx) / max(d, 1e-9), np.sin(beta) * (y - cy) / max(d, 1e-9), np.cos(beta))
                    break
            _tree(ax, base, normal)
        _decorate(ax, half, 45.0)
    elif test in ("T4_urban_bare", "UV_urban_crown"):
        half = 30.0
        _plane(ax, half, color="#c8c8c8")
        buildings = [(-18.0, 12.0), (8.0, 18.0), (16.0, -10.0), (-10.0, -16.0), (-20.0, -4.0)]
        for x, y in buildings:
            _box(ax, x - 5, y - 5, 0, 10, 10, 20, WALL, alpha=0.95)
            _box(ax, x - 5, y - 5, 20, 10, 10, 0.01, ROOF, alpha=1.0)
        if test == "UV_urban_crown":
            clear = lambda x, y: all(max(abs(x - bx) - 5, 0) ** 2 + max(abs(y - by) - 5, 0) ** 2 >= 9 for bx, by in buildings)
            for x, y in _scatter(rng, 8, half, 3.2, clear):
                _tree(ax, (x, y, 0), up)
        _decorate(ax, half, 30.0)
    else:
        raise ValueError(f"no sketch for {test}")


PLANES = {"太阳主平面（红虚线，+：太阳一侧）": (0.0, 180.0), "垂直太阳主平面（紫虚线，+：东）": (90.0, 270.0)}
FIGSIZE, DPI = (9.8, 4.2), 150  # 1470 x 630 px; no bbox_inches="tight", which would change it


def scene_figure(test, title, series, out, note=""):
    """Sketch + solar principal plane + perpendicular plane, both planes on the same brightness-temperature axis.

    series: list of (values {(vza, vaa): BT_K}, matplotlib format, plot kwargs, label)."""
    import matplotlib.pyplot as plt
    import compare

    labels = [label for *_, label in series]
    ncol = len(series) if sum(len(label) for label in labels) <= 60 else 1
    legend_rows = -(-len(series) // ncol)
    bottom = 0.15 + 0.05 * legend_rows

    fig = plt.figure(figsize=FIGSIZE, dpi=DPI)
    grid = fig.add_gridspec(1, 3, width_ratios=[0.9, 1.0, 1.0], wspace=0.08, left=0.0, right=0.985,
                            top=0.80 if note else 0.85, bottom=bottom)
    sketch = fig.add_subplot(grid[0], projection="3d")
    draw(sketch, test)
    sketch.set_title("场景示意", fontsize=9)
    principal = fig.add_subplot(grid[1])
    perpendicular = fig.add_subplot(grid[2], sharex=principal, sharey=principal)
    panels = [principal, perpendicular]
    values = []
    for ax, (plane, (pos, neg)) in zip(panels, PLANES.items()):
        for data, fmt, kwargs, label in series:
            x, y = compare.signed_curve(data, pos, neg)
            ax.plot(x, y - 273.15, fmt, label=label, **kwargs)
            values.extend(v for v in y - 273.15 if np.isfinite(v))
        ax.set_title(plane, fontsize=9)
        ax.set_xlabel("带符号观测天顶角 (°)", fontsize=8.5)
        ax.tick_params(labelsize=8)
        ax.grid(alpha=0.3)
    principal.set_ylabel("亮温 (°C)", fontsize=8.5)
    perpendicular.tick_params(labelleft=False)
    principal.set_xlim(-80, 80)
    low, high = min(values), max(values)
    pad = max(0.06 * (high - low), 0.25)
    principal.set_ylim(low - pad, high + pad)
    fig.legend(*principal.get_legend_handles_labels(), loc="lower center", ncol=ncol, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, 0.0))
    fig.suptitle(f"{title}（10.5 µm，SZA 30°，SAA 0°）" + (f"\n{note}" if note else ""), fontsize=10)
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
