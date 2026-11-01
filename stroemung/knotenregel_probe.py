#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Die Knotenregel am ENDE jedes vollstaendigen Zeitschritts.

Anlass: Befund B3. Es gab eine Probe auf das Drucknetz allein, aber keine darauf, dass am
Ende des ganzen Schritts noch stimmt, was die Projektion hergestellt hat. Genau dazwischen
sass der Fehler -- eine Randzeile, die die Einstroemung NACH der Druckkorrektur aenderte,
hinterliess eine stehende Divergenz von 6,97 in der Eckzelle.

Gemessen wird in jedem Schritt

    max |div u|  ueber alle Zellen, die nicht im Koerper liegen

und zwar auf den Feldern, wie sie der Zeitschritt endgueltig hinterlaesst. Der Massstab
ist nicht null, sondern die Rundungsschwelle: Die Divergenz wird mit u/dx gebildet, also
ist alles unterhalb von etwa eps*|u|/dx richtig.

Zweiter Teil: Der Nur-Lesen-Haken in wirbelstrasse.lauf darf die Rechnung nicht
veraendern. Das wird nicht behauptet, sondern geprueft -- derselbe Lauf mit und ohne
Haken muss Zahl fuer Zahl dasselbe ergeben.
"""
FASSUNG = "1.0"
import sys

import numpy as np

import wirbelstrasse as W


def probe(nx=110, ny=21, t_ende=1.0):
    gemessen = []

    def schauen(k, t, g, u, v, p):
        d = np.abs(g.divergenz(u, v))
        _, _, mp = W.koerper_masken(g)
        d[mp] = 0.0                      # Zellen im Koerper tragen keine Bilanz
        bezug = max(np.abs(u).max(), 1e-30) / g.dx
        gemessen.append((t, float(d.max()), float(d.max() / bezug)))

    g, u, v, p, verlauf, dt = W.lauf(nx, ny, t_ende=t_ende, beobachter=schauen)
    return np.array(gemessen), verlauf, dt


if __name__ == "__main__":
    nx = int(sys.argv[1]) if len(sys.argv) > 1 else 110
    ny = int(sys.argv[2]) if len(sys.argv) > 2 else 21
    t_ende = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    print("=" * 94)
    print(f" KNOTENREGEL AM ENDE JEDES ZEITSCHRITTS -- Fassung {FASSUNG}   Netz {nx}x{ny}")
    print("=" * 94)
    gem, verlauf, dt = probe(nx, ny, t_ende)
    roh, rel = gem[:, 1].max(), gem[:, 2].max()
    print(f"  {len(gem)} Zeitschritte geprueft, dt = {dt:.3e} s")
    print(f"  groesste Divergenz   absolut {roh:.3e} 1/s   bezogen auf |u|/dx  {rel:.3e}")
    for anteil in (0.0, 0.25, 0.5, 0.75):
        k = int(anteil * (len(gem) - 1))
        print(f"      t = {gem[k,0]:6.3f} s   |div u| = {gem[k,1]:.3e}")

    # --- der Haken darf nichts veraendern: gleicher Lauf ohne Beobachter --------------
    g2, u2, v2, p2, verlauf2, dt2 = W.lauf(nx, ny, t_ende=t_ende)
    gleich = bool(np.array_equal(verlauf, verlauf2))

    urteile = [
        ("Knotenregel am Ende jedes Schritts erfuellt (bezogen < 1e-11)",
         rel < 1e-11, f"{rel:.3e}"),
        ("der Nur-Lesen-Haken veraendert die Rechnung nicht (Bit fuer Bit)",
         gleich, "identische Verlaufstafel" if gleich else "ABWEICHUNG"),
    ]
    print()
    for text, ok, beleg in urteile:
        print(f"  [{'RICHTIG' if ok else 'FALSCH  '}] {text}\n              {beleg}")
    n = sum(1 for _, ok, _ in urteile if ok)
    print(f"\n  == {n}/{len(urteile)}")
    raise SystemExit(0 if n == len(urteile) else 1)
