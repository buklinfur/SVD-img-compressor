#!/usr/bin/env python3
"""Regenerate the README figures and tables from a real svdc run.

usage: python3 scripts/make_figures.py [--svdc build/svdc] [--image examples/test2.png]

Writes into docs/:
  grid.jpg              original next to reconstructions at k = 5, 20, 50
  rate_distortion.png   PSNR vs bits per pixel: svdc (8 and 16-bit factors) and JPEG
  spectrum.png          singular values and captured energy of each channel
  bench.csv             the raw svdc bench output (8-bit factors)
and prints a Markdown table for the README.

Needs numpy, Pillow and matplotlib (docs only, not needed to build svdc).
"""

import argparse
import csv
import io
import os
import subprocess
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image, ImageDraw, ImageFont  # noqa: E402

KS = [1, 2, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300]
GRID_KS = [5, 20, 50]

# validated categorical palette (dataviz reference palette, light mode)
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID = "#1f1f1e", "#6b6a63", "#e4e3de"


def psnr(a, b):
    mse = np.mean((a.astype(np.float64) - b.astype(np.float64)) ** 2)
    return float("inf") if mse == 0 else 10 * np.log10(255.0 ** 2 / mse)


def bench(svdc, image, bits, workdir, save=None):
    out = os.path.join(workdir, "bench_%d.csv" % bits)
    cmd = [svdc, "bench", image, "--ks", ",".join(map(str, KS)), "--bits", str(bits), "--csv", out]
    if save:
        cmd += ["--save", save]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    with open(out) as f:
        return list(csv.DictReader(f))


def jpeg_curve(pixels):
    points = []
    for q in [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95]:
        buf = io.BytesIO()
        Image.fromarray(pixels).save(buf, "JPEG", quality=q)
        decoded = np.asarray(Image.open(io.BytesIO(buf.getvalue())).convert("RGB"))
        points.append((len(buf.getvalue()), psnr(pixels, decoded), q))
    return points


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def make_grid(original, prefix, path):
    tiles = [("original", original)]
    for k in GRID_KS:
        tiles.append(("k = %d" % k, Image.open("%s_k%d.png" % (prefix, k)).convert("RGB")))
    h = 360
    scaled = [(name, im.resize((round(im.width * h / im.height), h), Image.LANCZOS)) for name, im in tiles]
    gap, label = 8, 30
    width = sum(im.width for _, im in scaled) + gap * (len(scaled) - 1)
    canvas = Image.new("RGB", (width, h + label), "white")
    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.load_default(size=18)  # Pillow >= 10.1
    except TypeError:
        font = ImageFont.load_default()
    x = 0
    for name, im in scaled:
        canvas.paste(im, (x, label))
        tw = draw.textlength(name, font=font)
        draw.text((x + (im.width - tw) / 2, 4), name, fill=INK, font=font)
        x += im.width + gap
    canvas.save(path, quality=92)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--svdc", default="build/svdc")
    parser.add_argument("--image", default="examples/test2.png")
    parser.add_argument("--out", default="docs")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    original = Image.open(args.image).convert("RGB")
    pixels = np.asarray(original)
    npx = pixels.shape[0] * pixels.shape[1]

    with tempfile.TemporaryDirectory() as work:
        prefix = os.path.join(work, "rec")
        rows8 = bench(args.svdc, args.image, 8, work, save=prefix)
        rows16 = bench(args.svdc, args.image, 16, work)
        make_grid(original, prefix, os.path.join(args.out, "grid.jpg"))
        with open(os.path.join(work, "bench_8.csv")) as src, \
                open(os.path.join(args.out, "bench.csv"), "w") as dst:
            dst.write(src.read())

    # rate-distortion: same x axis (bits per pixel) and y axis (PSNR) for every codec
    jpeg = jpeg_curve(pixels)
    fig, ax = plt.subplots(figsize=(7.2, 4.2), dpi=150)
    style(ax)
    series = [
        ("svdc, 8-bit factors", rows8, BLUE, "o"),
        ("svdc, 16-bit factors", rows16, AQUA, "s"),
    ]
    for name, rows, color, marker in series:
        x = [8 * int(r["bytes"]) / npx for r in rows]
        y = [float(r["psnr_db"]) for r in rows]
        ax.plot(x, y, color=color, linewidth=2, marker=marker, markersize=5, label=name)
        for r, xi, yi in zip(rows, x, y):
            if name.startswith("svdc, 8") and int(r["k"]) in (5, 20, 50, 100):
                ax.annotate("k=%s" % r["k"], (xi, yi), textcoords="offset points",
                            xytext=(-6, 6), ha="right", fontsize=8, color=MUTED)
    ax.plot([8 * b / npx for b, _, _ in jpeg], [p for _, p, _ in jpeg], color=ORANGE,
            linewidth=2, marker="^", markersize=5, label="JPEG (Pillow, q 5-95)")
    ax.set_xscale("log")
    ax.set_xlabel("bits per pixel (log scale)", color=INK)
    ax.set_ylabel("PSNR, dB", color=INK)
    ax.set_title("Quality vs size, %s (%dx%d)" % (os.path.basename(args.image),
                 pixels.shape[1], pixels.shape[0]), color=INK, fontsize=11, loc="left")
    ax.legend(frameon=False, fontsize=9, labelcolor=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "rate_distortion.png"))
    plt.close(fig)

    # singular value spectrum of each channel (computed with numpy, docs only)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.2), dpi=150)
    for ax in (a1, a2):
        style(ax)
    for c, (name, color) in enumerate([("red", ORANGE), ("green", AQUA), ("blue", BLUE)]):
        s = np.linalg.svd(pixels[..., c].astype(np.float64), compute_uv=False)
        a1.semilogy(np.arange(1, len(s) + 1), s, color=color, linewidth=1.5, label=name)
        a2.plot(np.arange(1, len(s) + 1), np.cumsum(s ** 2) / np.sum(s ** 2), color=color,
                linewidth=1.5, label=name)
    a1.set_title("singular values", color=INK, fontsize=10, loc="left")
    a1.set_xlabel("index i", color=INK)
    a2.set_title("energy kept by the first k", color=INK, fontsize=10, loc="left")
    a2.set_xlabel("k", color=INK)
    a2.set_xscale("log")
    a2.legend(frameon=False, fontsize=8, labelcolor=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "spectrum.png"))
    plt.close(fig)

    # markdown table for the README
    raw = pixels.size
    print("| k | file size | vs raw | bits/pixel | PSNR, dB | compress, s |")
    print("|--:|--:|--:|--:|--:|--:|")
    for r in rows8:
        b = int(r["bytes"])
        print("| %s | %.1f KiB | %.1fx | %.2f | %.2f | %.3f |" % (
            r["k"], b / 1024, raw / b, 8 * b / npx, float(r["psnr_db"]), float(r["compress_s"])))
    print("\nJPEG reference:")
    for b, p, q in jpeg:
        print("  q=%d: %.1f KiB, %.2f bpp, %.2f dB" % (q, b / 1024, 8 * b / npx, p))


if __name__ == "__main__":
    main()
