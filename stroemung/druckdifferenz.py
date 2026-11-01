#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Die vierte Referenzgroesse: Druckdifferenz ueber dem Zylinder.

Schaefer und Turek geben neben c_D, c_L und St noch

    dp = p(0,15 ; 0,20) - p(0,25 ; 0,20)   zum Zeitpunkt t_0 + 1/(2f)       2,46 ... 2,50

mit t_0 als dem Zeitpunkt des groessten Auftriebsbeiwerts. Diese Groesse ist von den
Kraftbeiwerten UNABHAENGIG: Sie liest den Druck an zwei Punkten ab, waehrend c_D und c_L
aus dem Impulsstrom gebildet werden. Trifft die Rechnung beide, so ist das kein Zufall
eines Auswerteweges.

Die beiden Punkte liegen GENAU auf der Zylinderoberflaeche. Auf dem kartesischen Gitter
liegt dort keine Zellmitte, und die naechstgelegene Zelle kann bereits im Koerper liegen --
dort ist der Druck bedeutungslos (Einheitszeile). Deshalb wird aus den beiden naechsten
FLUIDZELLEN auf der Achse linear auf die Oberflaeche extrapoliert. Das ist ein
Extrapolationsfehler erster Ordnung im Zellabstand; er wird mit dem Netz kleiner und ist
in der Netzreihe sichtbar. Geschaetzt wird nichts.
"""
FASSUNG = "1.0"
import json
import sys

import numpy as np

import wirbelstrasse as W


def _p_auf_achse(g, p, mp):
    """Druck entlang y = YC als Funktion von x, und wo er gueltig ist.

    In y wird zwischen den beiden Zellreihen interpoliert, die y = YC einschliessen.
    Eine Spalte gilt nur, wenn BEIDE beteiligten Zellen ausserhalb des Koerpers liegen.
    """
    j = np.searchsorted(g.yp, W.YC) - 1
    j = min(max(j, 0), g.ny - 2)
    y0, y1 = g.yp[j], g.yp[j + 1]
    w = (W.YC - y0) / (y1 - y0)
    pa = (1 - w) * p[:, j] + w * p[:, j + 1]
    gueltig = (~mp[:, j]) & (~mp[:, j + 1])
    return pa, gueltig


def _auf_rand(g, pa, gueltig, x_ziel, von_links):
    """Aus den beiden naechsten Fluidzellen linear auf x_ziel extrapolieren."""
    if von_links:
        kand = np.where(gueltig & (g.xp < x_ziel))[0]
        if len(kand) < 2:
            return float("nan")
        i2, i1 = kand[-1], kand[-2]          # i1 weiter weg, i2 naeher
    else:
        kand = np.where(gueltig & (g.xp > x_ziel))[0]
        if len(kand) < 2:
            return float("nan")
        i2, i1 = kand[0], kand[1]
    x1, x2 = g.xp[i1], g.xp[i2]
    return float(pa[i2] + (pa[i2] - pa[i1]) * (x_ziel - x2) / (x2 - x1))


def lauf(nx, ny, t_ende=12.0):
    vorn, hinten = W.XC - W.DZ / 2.0, W.XC + W.DZ / 2.0
    dp = []

    def schauen(k, t, g, u, v, p):
        _, _, mp = schauen.masken
        pa, gu = _p_auf_achse(g, p, mp)
        dp.append(_auf_rand(g, pa, gu, vorn, True) - _auf_rand(g, pa, gu, hinten, False))

    from gitter import Gitter
    schauen.masken = W.koerper_masken(Gitter(W.LX, W.LY, nx, ny))
    g, u, v, p, verlauf, dt = W.lauf(nx, ny, t_ende=t_ende, beobachter=schauen)
    return np.array(dp), verlauf, dt


def auswerten(dp, verlauf, dt, ab=0.5):
    """dp zum Zeitpunkt t_0 + T/2, mit t_0 = groesstes c_L im eingeschwungenen Teil."""
    n = len(verlauf)
    k0 = int(ab * n)
    t, cL = verlauf[k0:, 0], verlauf[k0:, 2]
    kreuz = np.where((cL[:-1] < 0) & (cL[1:] >= 0))[0]
    if len(kreuz) < 2:
        return float("nan"), float("nan"), float("nan")
    T = float(np.median(np.diff(t[kreuz])))
    # letzte volle Schwingung, darin das groesste c_L
    a, b = kreuz[-2], kreuz[-1]
    k_max = k0 + a + int(np.argmax(cL[a:b]))
    k_ziel = k_max + int(round(0.5 * T / dt))
    if k_ziel >= len(dp):
        k_ziel = k_max - int(round(0.5 * T / dt))       # eine Periode zurueck
    return float(dp[k_ziel]), T, float(verlauf[k_max, 0])


if __name__ == "__main__":
    nx = int(sys.argv[1]) if len(sys.argv) > 1 else 220
    ny = int(sys.argv[2]) if len(sys.argv) > 2 else 41
    t_ende = float(sys.argv[3]) if len(sys.argv) > 3 else 12.0
    print("=" * 92)
    print(f" DRUCKDIFFERENZ UEBER DEM ZYLINDER -- Fassung {FASSUNG}   Netz {nx}x{ny}")
    print(" dp = p(0,15;0,20) - p(0,25;0,20)  bei  t_0 + T/2,  t_0 = groesstes c_L")
    print("=" * 92)
    dp, verlauf, dt = lauf(nx, ny, t_ende)
    wert, T, t0 = auswerten(dp, verlauf, dt)
    lo, hi = 2.46, 2.50
    ok = lo <= wert <= hi
    print(f"  {len(dp)} Zeitschritte, dt = {dt:.3e} s")
    print(f"  Periode T = {T*1000:.1f} ms,  t_0 = {t0:.3f} s")
    print(f"  dp = {wert:.4f}      Referenz {lo} … {hi}      "
          f"{'im Intervall' if ok else 'daneben'}")
    print(f"  (Spanne von dp im eingeschwungenen Teil: "
          f"{np.nanmin(dp[len(dp)//2:]):.3f} … {np.nanmax(dp[len(dp)//2:]):.3f})")
    json.dump({"nx": nx, "ny": ny, "dp": wert, "T": T, "dt": dt},
              open(f"ergebnis_dp_{nx}x{ny}.json", "w"), indent=1)
