#!/usr/bin/env python
"""Create a simple local HTML gallery from generated PNG figures."""
from __future__ import annotations
from pathlib import Path
import html

def main():
    root=Path('.')
    fig=root/'figures'
    items=sorted(fig.glob('*.png'))
    blocks=[]
    for p in items:
        title=p.stem.replace('_',' ')
        blocks.append(f'<figure><img src="figures/{html.escape(p.name)}" loading="lazy"><figcaption>{html.escape(title)}</figcaption></figure>')
    page='''<!doctype html><html><head><meta charset="utf-8"><title>Geometric Cusp Validation Gallery</title>
<style>body{font-family:system-ui,sans-serif;margin:2rem}main{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:1rem}figure{margin:0}img{width:100%;height:auto;border:1px solid #bbb}figcaption{padding:.4rem 0 .8rem;font-size:.9rem}</style></head><body><h1>Geometric Cusp Validation Gallery</h1><main>'''+''.join(blocks)+'''</main></body></html>'''
    (root/'figures'/'index.html').write_text(page)
    print(root/'figures'/'index.html')
if __name__=='__main__': main()
