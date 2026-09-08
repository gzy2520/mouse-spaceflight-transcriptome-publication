#!/usr/bin/env python3
"""Create a local figure index and compact review sheet for the selected release."""
from html import escape
from pathlib import Path
import sys

from PIL import Image, ImageDraw


def build(root: Path) -> None:
    figures = sorted(root.glob("Main/*.png")) + sorted(root.glob("Suppl/*.png"))
    figures += sorted(root.glob("Suppl/Fig_S6_per_tissue/*.png"))
    cards = []
    for path in figures:
        relative = path.relative_to(root).as_posix()
        cards.append(
            f'<article><h2>{escape(path.stem)}</h2>'
            f'<a href="{escape(relative)}"><img loading="lazy" src="{escape(relative)}" '
            f'alt="{escape(path.stem)}"></a>'
            f'<p><a href="{escape(relative)}">PNG</a> · '
            f'<a href="{escape(relative[:-4] + ".pdf")}">PDF</a></p></article>'
        )
    (root / "index.html").write_text(
        '<!doctype html><html lang="zh-CN"><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        '<title>最终图表查验</title><style>body{font-family:system-ui;margin:32px;'
        'background:#f6f7f8;color:#203040}main{display:grid;grid-template-columns:'
        'repeat(auto-fit,minmax(440px,1fr));gap:24px}article{background:white;'
        'padding:16px;border:1px solid #ddd}img{width:100%;height:320px;object-fit:contain}'
        'h2{font-size:18px}a{color:#245b8b}@media(max-width:520px){main{display:block}}'
        '</style><h1>最终图表查验</h1><p>16 个主图面板、9 个补充图面板、26 个组织相关性图。'
        '点击图像查看原分辨率；PDF 可单独打开。Fig. 1a 不在此包中。</p>'
        '<p><a href="README.md">方法与目录说明</a> · '
        '<a href="provenance/output_sha256.csv">文件校验清单</a></p><main>'
        + ''.join(cards) + '</main></html>', encoding='utf-8'
    )
    main = sorted(root.glob("Main/*.png"))
    sheet = Image.new("RGB", (2000, 1200), "white")
    draw = ImageDraw.Draw(sheet)
    for i, path in enumerate(main):
        x, y = (i % 4) * 500, (i // 4) * 300
        with Image.open(path) as original:
            thumb = original.convert("RGB")
            thumb.thumbnail((490, 270))
            sheet.paste(thumb, (x + (500 - thumb.width) // 2, y + 25))
        draw.text((x + 12, y + 8), path.stem, fill="#203040")
    sheet.save(root / "provenance/main_figure_contact_sheet.jpg", quality=90)


if __name__ == "__main__":
    build(Path(sys.argv[1]).resolve(strict=True))
