#!/usr/bin/env python3
"""Parse the abundant-phylum network gml into a colors CSV (Id,Label,type,Color) that feeds
get_network_legend.R. Colors come straight from the gml graphics.fill (auto, always matches the figure)."""
import re, csv

GML = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/network_99_revised/whole_network2_abundant_phylum.gml"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/network_99_revised/network_colors.csv"

rows = []
cur = {}
state = None          # None | "node" | "edge"
in_graphics = False
for ln in open(GML):
    s = ln.strip()
    if s == "node [":
        state = "node"; cur = {}; in_graphics = False
    elif s == "edge [":
        state = "edge"
    elif s == "graphics [":
        in_graphics = True
    elif s == "]":
        if in_graphics:
            in_graphics = False
        elif state == "node":
            rows.append(cur); state = None
        elif state == "edge":
            state = None
    elif state == "node":
        m = re.match(r'(\w+)\s+"?([^"]*)"?$', s)
        if not m:
            continue
        k, v = m.group(1), m.group(2)
        if in_graphics and k == "fill":
            cur["Color"] = v
        elif k in ("id", "label", "type", "color"):
            cur[{"id": "Id", "label": "Label", "type": "type", "color": "color"}[k]] = v

# prefer graphics fill; fall back to flat color
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["Id", "Label", "type", "Color"])
    for r in rows:
        w.writerow([r.get("Id", ""), r.get("Label", ""), r.get("type", ""),
                    r.get("Color") or r.get("color", "")])
print(f"wrote {OUT} ({len(rows)} nodes)")
from collections import Counter
cats = Counter(r.get("type") for r in rows)
print("categories:", dict(cats))
