#!/usr/bin/env python3
"""
Reshape the per-site source-data CSV (from plot_persite_methylation.py) into the input format
of the mod_site_vis React/D3 web app, and print a ready-to-paste CONTIGS snippet.

mod_site_vis expects one CSV per (sample, contig, motif) at
  public/data/{sample}/{csvPrefix}_{motif}.csv
with header:  name,contig,type,start,stop,strand,score,legend
and treats a site as MODIFIED iff score == 1.0.

We have one sample, one motif (ACC), six contigs (tracks), each a 20 kb window. Positions are
shifted to window-relative coordinates (0..20 kb) so each track fills its axis like the static
PDF; the absolute genomic coordinate is preserved inside `name`.

Run:
  /home/shuaiw/miniconda3/envs/modifi/bin/python benchmark/borg/persite_to_modsitevis.py
"""
import os, csv, argparse
from collections import OrderedDict

SRC = "/home/shuaiw/MODIFI/tmp/rev_figs/borg/persite_ACC_1_SR-VP_07_25_2022_A1_100cm.sourcedata.csv"
APP_DATA_DIR = "/home/shuaiw/borg/revision/borg/mod_site_vis/public/data"
SAMPLE = "soil_100"
MOTIF = "ACC"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=SRC)
    ap.add_argument("--out_dir", default=APP_DATA_DIR)
    ap.add_argument("--sample", default=SAMPLE)
    ap.add_argument("--motif", default=MOTIF)
    args = ap.parse_args()

    out_sample_dir = os.path.join(args.out_dir, args.sample)
    os.makedirs(out_sample_dir, exist_ok=True)

    # group rows by track, preserving first-seen order
    tracks = OrderedDict()
    with open(args.src) as f:
        for row in csv.DictReader(f):
            tracks.setdefault(row["track"], []).append(row)

    contigs_cfg = []
    for track, rows in tracks.items():
        win_start = int(rows[0]["window_start"])
        win_end = int(rows[0]["window_end"])
        width = win_end - win_start
        out_csv = os.path.join(out_sample_dir, f"{track}_{args.motif}.csv")
        nmod = 0
        with open(out_csv, "w", newline="") as fo:
            w = csv.writer(fo)
            w.writerow(["name", "contig", "type", "start", "stop", "strand", "score", "legend"])
            for r in rows:
                pos = int(r["position_1based"])
                rel = pos - win_start                     # window-relative coordinate
                strand = r["strand"]
                modified = r["status"] == "modified"
                nmod += int(modified)
                name = f"{track}:{pos}{strand}_{args.motif}"   # absolute coord kept here
                score = "1.0" if modified else "0.5"
                legend = "yes" if modified else "no"
                w.writerow([name, track, args.motif, rel, rel, strand, score, legend])
        print(f"wrote {out_csv}  ({len(rows)} sites, {nmod} modified, window {win_start}-{win_end})")
        contigs_cfg.append((track, width))

    # print a CONTIGS array ready to paste into src/types.ts
    print("\n// ---- paste into src/types.ts (replace CONTIGS) ----")
    print("export const CONTIGS: ContigConfig[] = [")
    for track, width in contigs_cfg:
        print(f'  {{ id: "{track}", label: "{track}", csvPrefix: "{track}", '
              f'length: {width}, defaultViewEnd: {width} }},')
    print("];")
    print(f"\n// SAMPLES = [\"{args.sample}\"];  MOTIFS = [\"{args.motif}\"];")


if __name__ == "__main__":
    main()
