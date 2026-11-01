#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Erzeugt die Markdown-Tafel der Netzreihe AUS den Ergebnisdateien.

Warum ein eigenes Programm fuer eine Tabelle: Abgeschriebene Zahlen sind eine
Fehlerquelle, die sich nicht von selbst meldet. Steht die Tafel im Methodenblatt und
stammt aus der Hand, so kann sie von den Ergebnisdateien abweichen, ohne dass es jemand
merkt. Hier wird sie aus genau den Dateien gebildet, die die Laeufe geschrieben haben.

    python3 tafel_netzreihe.py                nur anzeigen
    python3 tafel_netzreihe.py --eintragen    in METHODE_UND_NACHWEIS.md einsetzen
                                              (ersetzt die Zeile <!-- TABELLE-NETZREIHE -->
                                               bzw. die zuvor eingesetzte Tafel)
"""
FASSUNG = "1.0"
import glob
import json
import math
import os
import sys

import wirbelstrasse as W

HIER = os.path.dirname(os.path.abspath(__file__))
BLATT = os.path.join(os.path.dirname(HIER), "METHODE_UND_NACHWEIS.md")
MARKE_AUF = "<!-- TABELLE-NETZREIHE -->"
MARKE_ZU = "<!-- ENDE-TABELLE-NETZREIHE -->"
MARKE_T_AUF = "<!-- TABELLE-TREPPE -->"
MARKE_T_ZU = "<!-- ENDE-TABELLE-TREPPE -->"


def _komma(x, n=4):
    """Deutsche Schreibweise: Dezimalkomma und echtes Minuszeichen, kein Bindestrich."""
    return f"{x:.{n}f}".replace(".", ",").replace("-", "\u2212")


def _prozent(x):
    return f"{x:+.1f} %".replace(".", ",").replace("-", "\u2212")


def stufen(muster, schluessel):
    aus = []
    for d in sorted(glob.glob(os.path.join(HIER, muster))):
        if "zweiwege" in d:
            continue
        e = json.load(open(d))
        if schluessel not in e:
            continue
        if any(not math.isfinite(float(e[k])) for k in (schluessel,)):
            continue
        aus.append((e["nx"] / W.LX * W.DZ, e))
    aus.sort(key=lambda s: s[0])
    return aus


def tafel():
    reihe = stufen("ergebnis_*x*.json", "cD")
    dps = stufen("ergebnis_dp_*x*.json", "dp")
    z = []
    z.append("| Zellen über $D$ | Netz | $c_D^{\\max}$ | $c_L^{\\max}$ | $\\mathrm{St}$ | $\\Delta p$ |")
    z.append("|--:|:--|--:|--:|--:|--:|")
    dp_nach_zpd = {round(k, 3): e["dp"] for k, e in dps}
    for zpd, e in reihe:
        dp = dp_nach_zpd.get(round(zpd, 3))
        z.append(f"| {zpd:.0f} | {e['nx']} × {e['ny']} | {_komma(e['cD'])} | "
                 f"{_komma(e['cL'])} | {_komma(e['St'])} | "
                 f"{_komma(dp) if dp is not None else '—'} |")
    # Die Referenzzeile traegt in der Netz-Spalte nur einen Strich: „Schaefer-Turek 2D-2"
    # bricht im Satz um und stoert die Tafel. Woher die Werte stammen, steht im Text.
    ref = [("Referenz", "—", "3,22 … 3,24", "0,99 … 1,01",
            "0,295 … 0,305", "2,46 … 2,50")]
    for r in ref:
        z.append("| **" + "** | **".join(r) + "** |")
    z.append("")
    z.append("Abweichung von der Mitte des jeweiligen Referenzintervalls:")
    z.append("")
    z.append("| Zellen über $D$ | $c_D$ | $c_L$ | $\\mathrm{St}$ | $\\Delta p$ |")
    z.append("|--:|--:|--:|--:|--:|")
    mitten = {"cD": 3.23, "cL": 1.00, "St": 0.30, "dp": 2.48}
    for zpd, e in reihe:
        dp = dp_nach_zpd.get(round(zpd, 3))
        felder = []
        for k in ("cD", "cL", "St"):
            felder.append(_prozent(100 * (e[k] / mitten[k] - 1)))
        felder.append(_prozent(100 * (dp / mitten["dp"] - 1))
                      if dp is not None else "—")
        z.append(f"| {zpd:.0f} | " + " | ".join(felder) + " |")
    return "\n".join(z)


def tafel_treppe():
    """Die Streuprobe: gleiche Geometrie, gleiche Aufloesung, andere Lage der Treppe."""
    import numpy as np
    from gitter import Gitter
    zeilen = []
    for d in sorted(glob.glob(os.path.join(HIER, "streuung_*x*.json"))):
        zeilen.append(json.load(open(d)))
    # der Lauf aus der Netzreihe gehoert dazu: er ist die vierte Lage
    grund = os.path.join(HIER, "ergebnis_440x82.json")
    if os.path.isfile(grund):
        e = json.load(open(grund))
        g = Gitter(W.LX, W.LY, e["nx"], e["ny"])
        e = dict(e, dx=g.dx, dy=g.dy,
                 lage_x=(W.XC / g.dx) % 1.0, lage_y=(W.YC / g.dy) % 1.0)
        zeilen.insert(0, e)
    if len(zeilen) < 2:
        return "*(noch keine Streuprobe gerechnet)*"
    z = ["| Netz | $\\Delta y$ | Mittelpunkt im Raster | $c_D^{\\max}$ | $c_L^{\\max}$ | $\\mathrm{St}$ |",
         "|:--|--:|:--|--:|--:|--:|"]
    for e in zeilen:
        z.append(f"| {e['nx']} × {e['ny']} | {_komma(e['dy'], 6)} | "
                 f"x {_komma(e['lage_x'], 2)} Δx, y {_komma(e['lage_y'], 2)} Δy | "
                 f"{_komma(e['cD'])} | {_komma(e['cL'])} | {_komma(e['St'])} |")
    import statistics as st
    z.append("")
    for schl, name in (("cD", "c_D"), ("cL", "c_L"), ("St", "St")):
        w = [e[schl] for e in zeilen]
        spanne = max(w) - min(w)
        z.append(f"- **${name}$**: {_komma(min(w))} … {_komma(max(w))}, "
                 f"Spanne {_komma(spanne)} "
                 f"({_prozent(100*spanne/st.fmean(w)).lstrip('+')} vom Mittelwert)")
    return "\n".join(z)


def _einsetzen(blatt, auf, zu, inhalt):
    s = open(blatt, encoding="utf-8").read()
    if zu in s:
        a = s.index(auf); b = s.index(zu) + len(zu)
        s = s[:a] + auf + "\n\n" + inhalt + "\n\n" + zu + s[b:]
    elif auf in s:
        s = s.replace(auf, auf + "\n\n" + inhalt + "\n\n" + zu, 1)
    else:
        raise SystemExit("FEHLER: Marke nicht gefunden: " + auf)
    open(blatt, "w", encoding="utf-8").write(s)


if __name__ == "__main__":
    t = tafel()
    print(t)
    tt = tafel_treppe()
    print()
    print(tt)
    if "--eintragen" in sys.argv:
        _einsetzen(BLATT, MARKE_AUF, MARKE_ZU, t)
        _einsetzen(BLATT, MARKE_T_AUF, MARKE_T_ZU, tt)
        print(f"\n  beide Tafeln in {os.path.basename(BLATT)} eingetragen")
