"""Create short MP4 introductions from TIRTEB screenshots."""

from __future__ import annotations

from pathlib import Path

import cv2
import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageFilter

ROOT = Path(__file__).resolve().parents[1]
SCREEN = ROOT / "docs" / "screenshots"
OUT = ROOT / "videos"
W, H, FPS, SECONDS = 1280, 720, 24, 3.0


def get_font(size: int, bold: bool = False):
    names = [
        Path(r"C:\Windows\Fonts\msyhbd.ttc") if bold else Path(r"C:\Windows\Fonts\msyh.ttc"),
        Path(r"C:\Windows\Fonts\arial.ttf"),
        Path("/mnt/c/Windows/Fonts/msyh.ttc"),
    ]
    for path in names:
        if path.exists():
            return ImageFont.truetype(str(path), size)
    return ImageFont.load_default()


TITLE, BODY, SMALL = get_font(44, True), get_font(25), get_font(18)


def background(shot: str | None) -> Image.Image:
    if shot and (SCREEN / shot).exists():
        image = Image.open(SCREEN / shot).convert("RGB")
        image.thumbnail((W, H))
        canvas = Image.new("RGB", (W, H), (12, 20, 19))
        canvas.paste(image, ((W - image.width) // 2, (H - image.height) // 2))
        canvas = canvas.filter(ImageFilter.GaussianBlur(1.0))
        return Image.blend(canvas, Image.new("RGB", (W, H), (6, 12, 12)), 0.48)
    return Image.new("RGB", (W, H), (12, 20, 19))


def wrap(draw: ImageDraw.ImageDraw, text: str, font, width: int) -> list[str]:
    lines = []
    for raw in text.split("\n"):
        line = ""
        for char in raw:
            if draw.textbbox((0, 0), line + char, font=font)[2] <= width:
                line += char
            else:
                if line:
                    lines.append(line)
                line = char
        if line:
            lines.append(line)
    return lines


def slide(title: str, body: list[str], shot: str | None = None, accent=(87, 211, 198)) -> Image.Image:
    image = background(shot).convert("RGBA")
    image = Image.alpha_composite(image, Image.new("RGBA", (W, H), (5, 12, 12, 165 if shot else 0))).convert("RGB")
    draw = ImageDraw.Draw(image)
    draw.rectangle((0, 0, W, 8), fill=accent)
    draw.text((54, 42), "TIRTEB", font=SMALL, fill=accent)
    draw.text((54, 82), title, font=TITLE, fill=(238, 245, 242))
    y = 175
    for item in body:
        draw.ellipse((60, y + 11, 70, y + 21), fill=accent)
        lines = wrap(draw, item, BODY, 625)
        for line in lines:
            draw.text((86, y), line, font=BODY, fill=(218, 231, 226))
            y += 35
        y += 16
    if shot and (SCREEN / shot).exists():
        source = Image.open(SCREEN / shot).convert("RGB")
        source.thumbnail((430, 520))
        card = Image.new("RGB", (450, 540), (18, 29, 27))
        card.paste(source, ((450 - source.width) // 2, (540 - source.height) // 2))
        image.paste(card, (760, 108))
        draw.rounded_rectangle((760, 108, 1210, 648), radius=14, outline=(82, 122, 112), width=2)
    return image


def write_video(name: str, slides: list[Image.Image]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    writer = cv2.VideoWriter(str(path), cv2.VideoWriter_fourcc(*"mp4v"), FPS, (W, H))
    for image in slides:
        frame = cv2.cvtColor(np.asarray(image), cv2.COLOR_RGB2BGR)
        for _ in range(int(FPS * SECONDS)):
            writer.write(frame)
    writer.release()
    print(path)


def main() -> None:
    videos = {
        "01_theory_overview.mp4": [
            slide("Theory and method", ["Time-step energy balance updates soil, leaf, roof, wall and street temperatures.", "Radiative transfer combines solar, sky and directional component proportions.", "The result contains fluxes, temperatures, reflectance and thermal brightness temperature."], "02_default_results.png"),
            slide("Model workflow", ["input_tirteb.csv -> defined and meteo files -> solar geometry -> resistance and physiology -> energy balance -> observations.csv", "State is carried from one meteorological time step to the next."], "03_urban_crown_settings.png", (237, 184, 102)),
        ],
        "02_gui_quick_start.mp4": [
            slide("Open the local GUI", ["Run TIRTEB.exe and use the local address ending in /time.html.", "The packaged executable contains gui, data and input_tirteb.csv."], "01_home.png"),
            slide("Run a time series", ["Set the time range and click Run time series.", "The charts show brightness temperature, component temperature and energy flux.", "Save processes.csv and observations.csv after completion."], "02_default_results.png"),
        ],
        "03_scenes_and_files.mp4": [
            slide("Compare scenes", ["Plane, terrain and urban surface modes are supported.", "Bare soil, homogeneous canopy, crop rows and forest crowns select the active structure."], "03_urban_crown_settings.png"),
            slide("Use external files", ["Upload meteorological, direction and wavelength files from the input panel.", "Measured leaf, soil and urban spectra are interpolated to the internal spectral grid."], "05_manual_geometry_results.png"),
        ],
        "04_results_and_distribution.mp4": [
            slide("Read and export results", ["Use the time-series chart for fluxes and temperatures.", "Use observations.csv for directional multi-band analysis and analysis/plot_observations.py for polar figures."], "04_urban_crown_results.png"),
            slide("Distribute the application", ["Copy or zip dist/TIRTEB.exe.", "The recipient does not need Python, NumPy, SciPy or the source data for the packaged GUI.", "Keep input_tirteb.csv and data/ only when running from source or changing default files."], "06_mobile_view.png", (237, 142, 125)),
        ],
    }
    for name, slides in videos.items():
        write_video(name, slides)


if __name__ == "__main__":
    main()
