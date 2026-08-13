#!/usr/bin/env python3
"""Shared plotting helpers for the simu_meta figures.

Bacterial species/genus names must render in italics (scientific nomenclature), while
strain/isolate IDs, GTDB placeholder epithets (`sp002270295`), and ranks above genus stay
upright. Use these when a figure puts species names on axes/legends/annotations.
"""
import re

# GTDB placeholder epithet like "sp002270295" / "sp." — genus italic, epithet upright
_GTDB_SP = re.compile(r"^(sp\d+|sp\.?)$")


def mathtext_species(name):
    r"""Return a matplotlib mathtext string that italicizes a bacterial name.

    'Escherichia coli'      -> r'$\it{Escherichia\ coli}$'
    'E. coli'               -> r'$\it{E.\ coli}$'
    'Klebsiella sp002270295'-> r'$\it{Klebsiella}$ sp002270295'   (genus only)
    'Enterobacteriaceae'    -> 'Enterobacteriaceae'               (single word/rank: upright)
    Non-species text is returned unchanged.
    """
    if not isinstance(name, str) or not name.strip():
        return name
    parts = name.split()
    if len(parts) == 1:
        return name                                   # bare genus/rank: leave upright (safe default)
    genus, epithet = parts[0], parts[1]
    esc = lambda w: w.replace(" ", r"\ ")
    if _GTDB_SP.match(epithet):
        return rf"$\it{{{esc(genus)}}}$ {' '.join(parts[1:])}"   # italic genus + upright placeholder
    binom = esc(f"{genus} {epithet}")
    tail = (" " + " ".join(parts[2:])) if len(parts) > 2 else ""
    return rf"$\it{{{binom}}}$" + tail


def italicize_ticklabels(ax, axis="x"):
    """Set existing tick labels on `axis` ('x'/'y') to italic species mathtext."""
    getter = ax.get_xticklabels if axis == "x" else ax.get_yticklabels
    setter = ax.set_xticklabels if axis == "x" else ax.set_yticklabels
    setter([mathtext_species(t.get_text()) for t in getter()])
