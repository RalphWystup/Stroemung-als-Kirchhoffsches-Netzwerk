#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Wieviel streuen die Beiwerte allein durch die Lage der Treppe?

Anlass: Die Richardson-Extrapolation aus drei Netzstufen laeuft ueber die Referenz hinaus,
und die beobachteten Ordnungen sind untereinander unvertraeglich (0,85 / 1,70 / 0,90). Bei
einem Verfahren mit EINEM beherrschenden Fehleranteil muessten sie aehnlich sein. Die
Vermutung ist, dass die Treppendarstellung des runden Zylinders den Fehler traegt und dass
der nicht glatt mit dem Zellabstand faellt, weil die abgedeckten Zellen SPRINGEN.

Eine Vermutung ist aber keine Begruendung. Gemessen wird deshalb so:

Die Geometrie bleibt UNVERAENDERT -- Kanal, Zylinder, Stoffwerte, alles wie im Benchmark.
Nur die Zellzahl wird um eins und zwei erhoeht. Der Zellabstand aendert sich dadurch um
weniger als zweieinhalb Prozent, die Aufloesung ist also praktisch dieselbe. Was sich stark
aendert, ist die LAGE des Zylindermittelpunkts im Zellraster und damit das Treppenmuster.

Streuen die Beiwerte zwischen diesen drei Laeufen deutlich, so ist die Treppe die Ursache
und nicht die Aufloesung -- und dann ist auch klar, warum eine Extrapolation, die einen
glatten Fehlerverlauf voraussetzt, hier nicht traegt.

Geschrieben wird nach streuung_<nx>x<ny>.json, NICHT nach ergebnis_*, damit die Netzreihe
davon unberuehrt bleibt.
"""
FASSUNG = "1.0"
import json
import sys

import numpy as np

import wirbelstrasse as W
from gitter import Gitter


def lage_im_raster(nx, ny):
    """Wo liegt der Zylindermittelpunkt zwischen den Zellmitten? 0 = auf einer Mitte."""
    g = Gitter(W.LX, W.LY, nx, ny)
    return (W.XC / g.dx) % 1.0, (W.YC / g.dy) % 1.0, g.dx, g.dy


if __name__ == "__main__":
    nx = int(sys.argv[1]); ny = int(sys.argv[2])
    t_ende = float(sys.argv[3]) if len(sys.argv) > 3 else 12.0
    fx, fy, dx, dy = lage_im_raster(nx, ny)
    _, _, mp = W.koerper_masken(Gitter(W.LX, W.LY, nx, ny))
    print("=" * 92)
    print(f" TREPPENSTREUUNG -- Fassung {FASSUNG}   Netz {nx}x{ny}")
    print(f" dx = {dx:.6f}  dy = {dy:.6f}   Mittelpunkt im Raster: "
          f"x {fx:.3f} dx, y {fy:.3f} dy   {int(mp.sum())} Zellen ganz im Koerper")
    print("=" * 92)
    g, u, v, p, verlauf, dt = W.lauf(nx, ny, t_ende=t_ende,
                                     melde=max(1, int(2.0 / 1e-3)))
    cD, cL, St, T = W.auswerten(verlauf, dt)
    print(f"\n  c_D max {cD:.4f}   c_L max {cL:.4f}   St {St:.4f}   T {T*1000:.1f} ms")
    json.dump({"nx": nx, "ny": ny, "dx": dx, "dy": dy, "lage_x": fx, "lage_y": fy,
               "zellen_im_koerper": int(mp.sum()),
               "cD": cD, "cL": cL, "St": St, "T": T, "dt": dt},
              open(f"streuung_{nx}x{ny}.json", "w"), indent=1)
