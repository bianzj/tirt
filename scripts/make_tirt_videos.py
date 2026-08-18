from __future__ import annotations

from pathlib import Path

import cv2
import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageFilter

ROOT = Path(__file__).resolve().parent.parent
SCREEN = ROOT / "docs" / "screenshots"
OUT = ROOT / "videos"
W, H = 1280, 720
FPS = 24
SECONDS_PER_SLIDE = 3.2

FONT_CANDIDATES = [
    Path(r"C:\Windows\Fonts\msyh.ttc"),
    Path(r"C:\Windows\Fonts\simhei.ttf"),
    Path(r"C:\Windows\Fonts\arial.ttf"),
]


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [Path(r"C:\Windows\Fonts\msyhbd.ttc")] if bold else []
    candidates += FONT_CANDIDATES
    for path in candidates:
        if path.exists():
            return ImageFont.truetype(str(path), size=size)
    return ImageFont.load_default()


F_TITLE = font(48, True)
F_SUB = font(28)
F_BODY = font(25)
F_SMALL = font(19)
F_LABEL = font(18, True)


def cover_image(path: Path | None, size: tuple[int, int], dim: float = 0.35) -> Image.Image:
    if path and path.exists():
        img = Image.open(path).convert("RGB")
        img.thumbnail((int(size[0] * 1.08), int(size[1] * 1.08)))
        canvas = Image.new("RGB", size, (16, 23, 22))
        x = (size[0] - img.width) // 2
        y = (size[1] - img.height) // 2
        canvas.paste(img, (x, y))
        canvas = canvas.filter(ImageFilter.GaussianBlur(1.2))
        overlay = Image.new("RGB", size, (8, 12, 12))
        return Image.blend(canvas, overlay, dim)
    return Image.new("RGB", size, (16, 23, 22))


def wrap(draw: ImageDraw.ImageDraw, text: str, fnt, max_width: int) -> list[str]:
    lines: list[str] = []
    for raw in text.split("\n"):
        line = ""
        for char in raw:
            candidate = line + char
            if draw.textbbox((0, 0), candidate, font=fnt)[2] <= max_width:
                line = candidate
            else:
                if line:
                    lines.append(line)
                line = char
        if line:
            lines.append(line)
    return lines


def panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], fill=(22, 32, 30), outline=(67, 98, 91)) -> None:
    draw.rounded_rectangle(box, radius=18, fill=fill, outline=outline, width=2)


def add_text(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, fnt, fill=(235, 244, 241), max_width: int | None = None, line_gap: int = 8) -> int:
    x, y = xy
    lines = wrap(draw, text, fnt, max_width) if max_width else text.split("\n")
    for line in lines:
        draw.text((x, y), line, font=fnt, fill=fill)
        y += draw.textbbox((0, 0), line, font=fnt)[3] + line_gap
    return y


def screenshot_card(base: Image.Image, path: Path, box: tuple[int, int, int, int]) -> None:
    if not path.exists():
        return
    img = Image.open(path).convert("RGB")
    target_w = box[2] - box[0]
    target_h = box[3] - box[1]
    img.thumbnail((target_w, target_h))
    card = Image.new("RGB", (target_w, target_h), (18, 27, 25))
    card.paste(img, ((target_w - img.width) // 2, (target_h - img.height) // 2))
    mask = Image.new("L", (target_w, target_h), 0)
    ImageDraw.Draw(mask).rounded_rectangle((0, 0, target_w - 1, target_h - 1), radius=16, fill=255)
    base.paste(card, (box[0], box[1]), mask)
    d = ImageDraw.Draw(base)
    d.rounded_rectangle(box, radius=16, outline=(82, 122, 112), width=2)


def slide(title: str, subtitle: str = "", bullets: list[str] | None = None, shot: str | None = None, accent=(87, 211, 198)) -> Image.Image:
    bg = cover_image(SCREEN / shot if shot else None, (W, H), 0.55 if shot else 0.0)
    overlay = Image.new("RGBA", (W, H), (7, 13, 13, 170 if shot else 0))
    img = Image.alpha_composite(bg.convert("RGBA"), overlay).convert("RGB")
    draw = ImageDraw.Draw(img)
    draw.rectangle((0, 0, W, 8), fill=accent)
    draw.text((54, 44), "TiRT", font=F_LABEL, fill=accent)
    add_text(draw, (54, 82), title, F_TITLE, max_width=720)
    if subtitle:
        add_text(draw, (56, 158), subtitle, F_SUB, fill=(202, 218, 213), max_width=700)
    if bullets:
        y = 250
        for item in bullets:
            draw.ellipse((58, y + 10, 68, y + 20), fill=accent)
            y = add_text(draw, (84, y), item, F_BODY, max_width=610, line_gap=6) + 8
    if shot:
        screenshot_card(img, SCREEN / shot, (760, 96, 1218, 620))
    return img


def write_video(name: str, frames: list[Image.Image]) -> None:
    OUT.mkdir(exist_ok=True)
    path = OUT / name
    writer = cv2.VideoWriter(str(path), cv2.VideoWriter_fourcc(*"mp4v"), FPS, (W, H))
    for frame in frames:
        arr = cv2.cvtColor(np.asarray(frame), cv2.COLOR_RGB2BGR)
        for _ in range(int(FPS * SECONDS_PER_SLIDE)):
            writer.write(arr)
    writer.release()
    print(path)


def make_all() -> None:
    videos = {
        "01_overview_theory.mp4": [
            slide("TiRT 理论概览", "热红外方向性辐射传输模拟", ["输入地表结构、光谱、温度和观测几何", "输出方向辐亮度、亮温和组分发射率", "适用于植被、坡面、地形和城市场景"], "02_default_results.png"),
            slide("从组分到亮温", "核心计算流程", ["普朗克函数把组分温度转换为辐亮度", "几何光学估计可见比例和日照/阴影比例", "反普朗克函数把总辐亮度转换为亮温"], None, (237, 184, 102)),
            slide("支持的场景组合", "Surface × Vegetation", ["Plane / Slope / Terrain / Urban", "Bare soil / Turbid veg / Crop rows / Forest crowns", "统一由 run(input.csv) 生成输出表"], "03_urban_crown_settings.png"),
            slide("方向性结果", "为什么同一地表不同角度亮温不同", ["观测角度改变可见组分比例", "太阳方位改变日照与阴影分布", "结构越复杂，方向差异越明显"], "04_urban_crown_results.png"),
        ],
        "02_gui_quick_start.mp4": [
            slide("快速上手", "打开 TiRT.exe 后直接进入本地 GUI", ["无需额外复制 data 或 input.csv", "默认数据已经内置到 exe", "浏览器访问 127.0.0.1 本地服务"], "01_home.png"),
            slide("默认运行", "点击 Run simulation", ["默认场景为平面 + 均质植被", "默认三波段：8、10.5、12 um", "运行完成后状态显示 Complete"], "02_default_results.png"),
            slide("查看指标", "结果区自动更新", ["Mean temperature：当前波段平均亮温", "Directional span：方向性变化范围", "Observations：当前波段方向数量"], "02_default_results.png"),
            slide("保存结果", "导出 observations.csv", ["保存所有方向和波段结果", "可切换摄氏温度、K 温度或辐亮度", "适合后续绘图和统计分析"], "02_default_results.png"),
        ],
        "03_scene_and_input_settings.mp4": [
            slide("场景结构", "选择地表和植被类型", ["Surface mode 控制平面、坡面、地形、城市", "Vegetation type 控制裸土、均质植被、垄行和树冠", "界面会自动显示对应参数"], "03_urban_crown_settings.png"),
            slide("角度设置", "主平面、半球和手动角度", ["Solar zenith / Solar azimuth 设置太阳位置", "Principal plane 生成主平面方向", "Manual 可输入 VZA/VAA 方向对"], "05_manual_geometry_results.png"),
            slide("光谱设置", "直接输入或上传文件", ["波段和发射率用分号分隔", "文件模式会插值到目标波段", "发射率数量必须与波段数量一致"], "01_home.png"),
            slide("温度设置", "日照与阴影组分温度", ["叶片、土壤、地形、屋顶、墙体、街道分别设置", "温度单位为摄氏度", "结构和温差共同决定方向性"], "04_urban_crown_results.png"),
        ],
        "04_results_export_distribution.mp4": [
            slide("结果图解读", "Polar / Parallel / Vertical", ["Polar 图展示半球方向冷热分布", "Parallel 曲线展示太阳主平面响应", "Vertical 曲线展示垂直主平面响应"], "02_default_results.png"),
            slide("切换输出量", "Temperature 或 Radiance", ["Temperature (C)：常用亮温结果", "Temperature (K)：绝对温度单位", "Radiance：热辐射亮度"], "04_urban_crown_results.png"),
            slide("导出与复现", "保存 observations.csv 和 input.csv", ["observations.csv 保存模拟结果", "input.csv 保存当前参数", "两者可用于论文图表和复现实验"], "05_manual_geometry_results.png"),
            slide("分发给别人", "压缩 dist 目录", ["当前 dist 里只有 TiRT.exe", "默认 gui/data/input.csv 已内置", "对方双击即可打开本地 GUI"], None, (237, 142, 125)),
        ],
    }
    for name, frames in videos.items():
        write_video(name, frames)


if __name__ == "__main__":
    make_all()
