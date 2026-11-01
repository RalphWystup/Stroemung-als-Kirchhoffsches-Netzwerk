#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Sprachprobe des Drucknetzes MIT Stoerkoerper und Ausstromrand.

Die Probe in ns_loeser.py prueft den nackten Laplace-Operator. Der Fall der Wirbelstrasse
verlangt aber zwei Dinge mehr, und genau an beiden sass je ein Fehler (Befund B2/B3):

  * aufgetrennte Zweige dort, wo die Flaeche im Koerper liegt
  * ein Zweig ans Bezugspotential in der letzten Zellspalte (Umgebungsdruck am Ausstrom)

Geprueft wird, dass die beiden Aufbauwege denselben Operator ergeben:

  Weg A  drucknetz_matrix        -- Zweig fuer Zweig in eine duennbesetzte Matrix gestempelt
  Weg B  drucknetz_wirkung_voll  -- Flaechenfluesse als ganze Felder, ohne jede Matrix

Dazu drei Eigenschaften, die ohne Rechnung gelten muessen.
"""
FASSUNG = "1.0"
import numpy as np

from gitter import Gitter
from ns_loeser import drucknetz_matrix, drucknetz_wirkung_voll
from wirbelstrasse import koerper_masken, LX, LY


if __name__ == "__main__":
    print("=" * 96)
    print(f" DRUCKNETZ MIT STOERKOERPER UND AUSSTROM -- Fassung {FASSUNG}")
    print(" Weg A (Zweigstempel -> Matrix)  gegen  Weg B (Feldfluesse, matrixfrei)")
    print("=" * 96)
    rng = np.random.default_rng(5)
    urteile = []
    for nx, ny in ((110, 21), (220, 41), (440, 82)):
        g = Gitter(LX, LY, nx, ny)
        mu, mv, mp = koerper_masken(g)
        A, leer = drucknetz_matrix(g, ausstrom_rechts=True, u_maske=mu, v_maske=mv)
        p = rng.standard_normal((nx, ny))
        rA = (A @ p.ravel()).reshape(nx, ny)
        rB = drucknetz_wirkung_voll(g, p, leer=leer, ausstrom_rechts=True,
                                    u_maske=mu, v_maske=mv)
        rel = np.abs(rA - rB).max() / np.abs(rA).max()
        urteile.append((f"{nx}x{ny}: beide Wege, {int(leer.sum())} Zellen im Koerper, "
                        f"{int(mu.sum()+mv.sum())} Zweige getrennt", rel < 1e-12, f"{rel:.3e}"))
        if nx == 220:
            # 1 Konstanter Druck: ohne Ausstromzweig stromlos, MIT Ausstromzweig fliesst es
            A0, leer0 = drucknetz_matrix(g, ausstrom_rechts=False, u_maske=mu, v_maske=mv)
            eins = np.ones(nx * ny); eins[leer0] = 0.0
            r0 = A0 @ eins
            # Nicht mitgezaehlt werden die Knoten, die KEINE Knotenregel tragen, sondern
            # eine Bedingung: die Zellen ganz im Koerper (Einheitszeile) und -- im
            # geschlossenen Fall -- der eine festgehaltene Bezugsknoten 0. Ohne diese
            # Ausnahme meldet die Probe die Bezugszeile als Strom (gemessen: genau 1,0).
            r0[leer0] = 0.0
            r0[0] = 0.0
            urteile.append(("konstanter Druck treibt im geschlossenen Netz keinen Strom "
                            "(ohne Bezugs- und Koerperknoten)",
                            np.abs(r0).max() < 1e-10, f"{np.abs(r0).max():.3e}"))
            eins = np.ones(nx * ny); eins[leer] = 0.0
            r1 = (A @ eins).reshape(nx, ny)
            soll = g.dy / (g.dx / 2.0)
            innen = np.abs(r1[:-1, :]).max()
            urteile.append(("am Ausstrom treibt er genau den Zweig ans Bezugspotential",
                            innen < 1e-10 and abs(r1[-1, 0] - soll) < 1e-10,
                            f"innen {innen:.3e}, letzte Spalte {r1[-1,0]:.6g} "
                            f"(erwartet {soll:.6g})"))
            # 2 Symmetrie: ein Widerstandsnetz hat keine Vorzugsrichtung
            uns = abs(A - A.T).max() / abs(A).max()
            urteile.append(("die Drucknetzmatrix ist symmetrisch", uns < 1e-14, f"{uns:.3e}"))
    for text, ok, beleg in urteile:
        print(f"  [{'RICHTIG' if ok else 'FALSCH  '}] {text}\n              {beleg}")
    n = sum(1 for _, ok, _ in urteile if ok)
    print(f"\n  == {n}/{len(urteile)}")
    raise SystemExit(0 if n == len(urteile) else 1)
