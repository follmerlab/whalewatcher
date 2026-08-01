"""Rasterise the whalewatcher icon into PNGs, a multi-resolution .ico, and a
legibility contact sheet.

The master SVG is rendered once at 1024 px through headless Chromium, then
Lanczos downsampled. Rendering the vector natively at 16 px throws away the
antialiasing that makes a small icon readable; downsampling from a high-res
render is what actually holds up.

cairosvg is unusable here (no libcairo on Windows), so Playwright rasterises
and Pillow does everything downstream.

    python make_icon.py                 # the shipped mark
    python make_icon.py binoculars eyes # plus superseded variants, to compare
"""
import io
import sys
import pathlib

from PIL import Image, ImageDraw
from playwright.sync_api import sync_playwright

HERE = pathlib.Path(__file__).parent
MASTER_PX = 1024
PNG_SIZES = [512, 256, 128, 64, 48, 32, 16]

# The shipped mark. Rendered by default; feeds every .ico slot.
ARTWORK = {"icon": "whalewatcher_icon.svg"}

# Superseded explorations, kept as sources so the design history is inspectable.
# Rendered only when named explicitly: `python make_icon.py binoculars`
VARIANTS = {
    "binoculars": "variants/whalewatcher_icon_binoculars.svg",
    "eyes": "variants/whalewatcher_icon_eyes.svg",
    "flat": "variants/whalewatcher_icon_small.svg",
}

# .ico slot -> artwork.
# The scope mark carries every size on its own: its dark barrel ring gives a
# crisp boundary at 16 px that a flat disc lacks, and dropping the in-scene
# binoculars freed enough field for the whale to read at 32 px. Earlier passes
# needed separate reduced artwork in the 16/32 slots; that is no longer
# necessary, though the per-slot mechanism is left in place.
ICO_PLAN = {s: "icon" for s in (16, 32, 48, 64, 128, 256)}

DARK_BG = "#1a1a2e"   # the app's 3D canvas colour
LIGHT_BG = "#ffffff"


def render_master(svg_path: pathlib.Path, px: int) -> Image.Image:
    svg = svg_path.read_text(encoding="utf-8")
    html = (
        "<!doctype html><html><head><meta charset='utf-8'><style>"
        "html,body{margin:0;padding:0;background:transparent;}"
        f"svg{{display:block;width:{px}px;height:{px}px;}}"
        "</style></head><body>" + svg + "</body></html>"
    )
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page(viewport={"width": px, "height": px},
                                device_scale_factor=1)
        page.set_content(html, wait_until="load")
        buf = page.screenshot(omit_background=True)
        browser.close()
    return Image.open(io.BytesIO(buf)).convert("RGBA")


def downsample(master: Image.Image, size: int) -> Image.Image:
    return master.resize((size, size), Image.LANCZOS)


def flatten(img: Image.Image, bg_hex: str) -> Image.Image:
    bg = Image.new("RGBA", img.size, bg_hex)
    return Image.alpha_composite(bg, img).convert("RGB")


def contact_sheet(masters: dict, out: pathlib.Path) -> None:
    """One row per (artwork, background). Each cell shows the icon at true
    scale with a 4x nearest-neighbour blow-up beneath, so the pixel grid is
    actually judgeable rather than asserted."""
    sizes = [16, 32, 64, 128]
    pad, gap, label_h = 30, 30, 20
    cell = 132

    rows = [(name, bg, lbl)
            for name in masters
            for bg, lbl in ((LIGHT_BG, "#ffffff"), (DARK_BG, DARK_BG))]

    row_h = label_h + cell + 6 + cell + 18
    W = pad * 2 + len(sizes) * cell + (len(sizes) - 1) * gap
    H = pad * 2 + len(rows) * row_h + 30

    sheet = Image.new("RGB", (W, H), "#efefef")
    draw = ImageDraw.Draw(sheet)
    draw.text((pad, 10), "whalewatcher icon - true scale over 4x nearest-neighbour",
              fill="#222222")

    for r, (name, bg_hex, bglbl) in enumerate(rows):
        y0 = pad + 26 + r * row_h
        draw.text((pad, y0), f"{name}  on {bglbl}", fill="#222222")
        for c, s in enumerate(sizes):
            x0 = pad + c * (cell + gap)
            flat = flatten(downsample(masters[name], s), bg_hex)

            sw = Image.new("RGB", (cell, cell), bg_hex)
            sw.paste(flat, ((cell - s) // 2, (cell - s) // 2))
            sheet.paste(sw, (x0, y0 + label_h))
            draw.rectangle([x0, y0 + label_h, x0 + cell, y0 + label_h + cell],
                           outline="#999999")
            draw.text((x0 + 3, y0 + label_h + 2), f"{s}px", fill="#666666")

            # blow-up factor chosen so it always fits the cell (no cropping,
            # which silently broke the 64/128 cells in the first pass)
            factor = max(1, cell // s)
            zoom = flat.resize((s * factor, s * factor), Image.NEAREST)
            zc = Image.new("RGB", (cell, cell), bg_hex)
            zc.paste(zoom, ((cell - zoom.width) // 2, (cell - zoom.height) // 2))
            zy = y0 + label_h + cell + 6
            sheet.paste(zc, (x0, zy))
            draw.rectangle([x0, zy, x0 + cell, zy + cell], outline="#999999")
            draw.text((x0 + 3, zy + 2), f"{factor}x", fill="#666666")

    sheet.save(out)
    print(f"contact sheet -> {out.name}  ({sheet.width}x{sheet.height})")


def main() -> None:
    requested = sys.argv[1:]
    plan = dict(ARTWORK)
    for name in requested:
        if name in VARIANTS:
            plan[name] = VARIANTS[name]
        elif name not in ARTWORK:
            print(f"unknown artwork {name!r}; "
                  f"known variants: {list(VARIANTS)}")

    masters = {}
    for stem, svg_name in plan.items():
        svg = HERE / svg_name
        if not svg.exists():
            print(f"SKIP {svg_name} (not found)")
            continue
        print(f"rendering {svg_name} at {MASTER_PX}px...")
        masters[stem] = render_master(svg, MASTER_PX)
        masters[stem].save(HERE / f"whalewatcher_{stem}_1024.png")
        for s in PNG_SIZES:
            downsample(masters[stem], s).save(
                HERE / f"whalewatcher_{stem}_{s}.png")
        print(f"  png -> whalewatcher_{stem}_{{{','.join(map(str, PNG_SIZES))}}}.png")

    if set(ICO_PLAN.values()) <= set(masters):
        frames = []
        for size in sorted(ICO_PLAN):
            frames.append(downsample(masters[ICO_PLAN[size]], size))
        ico = HERE / "whalewatcher.ico"
        # append_images carries the remaining frames; each keeps its own artwork
        frames[-1].save(ico, format="ICO",
                        sizes=[(f.width, f.height) for f in frames],
                        append_images=frames[:-1])
        summary = ", ".join(f"{s}:{ICO_PLAN[s]}" for s in sorted(ICO_PLAN))
        print(f"ico -> {ico.name}  ({summary})")

    if masters:
        contact_sheet(masters, HERE / "whalewatcher_contactsheet.png")
    print("done")


if __name__ == "__main__":
    main()
