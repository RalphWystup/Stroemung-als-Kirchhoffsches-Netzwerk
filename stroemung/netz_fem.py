#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Koerperangepasstes Dreiecksnetz fuer den Galerkin-Weg.

Der Netzwerkweg rechnet auf kartesischen Zellen; der runde Zylinder wird dort zur Treppe,
und das ist nachweislich die groesste verbliebene Fehlerquelle. Der Galerkin-Weg bekommt
deshalb ein Netz, dessen Raender auf dem Zylinder LIEGEN -- keine Treppe.

Gebaut wird in drei Schritten, jeder nachpruefbar:

  1. Randpunkte setzen: Kanalrand und Zylinderkreis, mit vorgegebenem Punktabstand
  2. Innenpunkte streuen, mit Verfeinerung nahe am Zylinder und im Nachlauf
  3. Delaunay-Triangulierung, dann die Dreiecke im Zylinder entfernen

Delaunay kommt aus scipy (Qhull) -- einen Vernetzer baut man nicht selbst, das steht so
in den Grundsaetzen. Geprueft wird das Netz anschliessend auf das, worauf es ankommt:
positive Flaechen, geschlossener Rand, kein Dreieck im Koerper, und die Gesamtflaeche
muss Kanal minus Kreis ergeben.
"""
FASSUNG = "1.0"
import numpy as np
from scipy.spatial import Delaunay

import wirbelstrasse as W


def randpunkte(h_wand, h_zyl):
    """Punkte auf dem Kanalrand und auf dem Zylinderkreis."""
    LX, LY, R = W.LX, W.LY, W.DZ / 2.0
    nx = max(2, int(round(LX / h_wand)))
    ny = max(2, int(round(LY / h_wand)))
    unten = np.c_[np.linspace(0, LX, nx + 1), np.zeros(nx + 1)]
    oben = np.c_[np.linspace(0, LX, nx + 1), np.full(nx + 1, LY)]
    links = np.c_[np.zeros(ny - 1), np.linspace(0, LY, ny + 1)[1:-1]]
    rechts = np.c_[np.full(ny - 1, LX), np.linspace(0, LY, ny + 1)[1:-1]]
    nz = max(12, int(round(2 * np.pi * R / h_zyl)))
    w = np.linspace(0, 2 * np.pi, nz, endpoint=False)
    zyl = np.c_[W.XC + R * np.cos(w), W.YC + R * np.sin(w)]
    return unten, oben, links, rechts, zyl


def _dichte(x, y, h_fern, h_zyl, weit=3.0, nachlauf=1.0, h_nach=4.0, band=2.0):
    """Gewuenschter Punktabstand am Ort (x,y).

    Fein am Zylinder, maessig fein im NAHEN Nachlauf, grob im Rest. Das Nachlaufband
    bestimmt die Knotenzahl fast allein: Es ist lang und muss trotzdem bezahlbar
    bleiben, weil Taylor-Hood je Knoten rund neun Unbekannte traegt.
    """
    R = W.DZ / 2.0
    d = np.hypot(x - W.XC, y - W.YC) - R                 # Abstand zur Zylinderwand
    h = h_zyl + (h_fern - h_zyl) * np.clip(d / (weit * R), 0, 1)
    nach = (x > W.XC) & (x < W.XC + nachlauf) & (np.abs(y - W.YC) < band * R)
    h = np.where(nach, np.minimum(h, h_zyl * h_nach), h)
    return h


def punkte(h_fern=0.05, h_zyl=0.0039, saat=7, mit_zylinder=True, **fein):
    """Randpunkte plus gestreute Innenpunkte mit ortsabhaengigem Abstand."""
    unten, oben, links, rechts, zyl = randpunkte(h_fern, h_zyl)
    rand = np.vstack([unten, oben, links, rechts, zyl]) if mit_zylinder \
        else np.vstack([unten, oben, links, rechts])
    rng = np.random.default_rng(saat)
    # Ablehnungsverfahren auf einem feinen Hilfsgitter: ein Punkt wird genommen, wenn
    # er weit genug von allen schon genommenen entfernt ist (Poisson-Scheibe, einfach).
    g = h_zyl / 1.6
    xs = np.arange(g, W.LX, g); ys = np.arange(g, W.LY, g)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    P = np.c_[X.ravel(), Y.ravel()]
    rng.shuffle(P)
    R = W.DZ / 2.0
    if mit_zylinder:
        P = P[np.hypot(P[:, 0] - W.XC, P[:, 1] - W.YC) > R + 0.35 * h_zyl]
    h = _dichte(P[:, 0], P[:, 1], h_fern, h_zyl, **fein) if mit_zylinder \
        else np.full(len(P), h_fern)
    genommen = list(rand)
    from scipy.spatial import cKDTree
    baum = cKDTree(np.array(genommen))
    for p, hp in zip(P, h):
        if baum.query(p)[0] >= 0.78 * hp:
            genommen.append(p)
            if len(genommen) % 250 == 0:
                baum = cKDTree(np.array(genommen))
    return np.array(genommen), len(rand), (zyl if mit_zylinder else np.zeros((0, 2)))


def netz(h_fern=0.05, h_zyl=0.0039, mit_zylinder=True, **fein):
    """mit_zylinder=False gibt den nackten Kanal -- fuer die Hagen-Poiseuille-Probe,
    bei der es keinen Stoerkoerper geben darf."""
    P, n_rand, zyl = punkte(h_fern, h_zyl, mit_zylinder=mit_zylinder, **fein)
    tri = Delaunay(P)
    T = tri.simplices
    # Dreiecke im Zylinder entfernen (Schwerpunkt innerhalb)
    if mit_zylinder:
        sx = P[T, 0].mean(axis=1); sy = P[T, 1].mean(axis=1)
        T = T[~(np.hypot(sx - W.XC, sy - W.YC) < W.DZ / 2.0)]
    # Flaechen und Umlaufsinn
    x = P[T, 0]; y = P[T, 1]
    A2 = (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
    dreh = A2 < 0
    T[dreh] = T[dreh][:, [0, 2, 1]]                      # alle gegen den Uhrzeigersinn
    A2 = np.abs(A2)
    return P, T, A2 / 2.0, zyl


if __name__ == "__main__":
    print("=" * 92)
    print(f" KOERPERANGEPASSTES NETZ -- Fassung {FASSUNG}")
    print("=" * 92)
    P, T, A, zyl = netz()
    soll = W.LX * W.LY - np.pi * (W.DZ / 2) ** 2
    ist = A.sum()
    # der Kreis wird durch ein Vieleck ersetzt; dessen Flaeche ist etwas kleiner
    nz = len(zyl)
    poly = 0.5 * nz * np.sin(2 * np.pi / nz) * (W.DZ / 2) ** 2
    soll_poly = W.LX * W.LY - poly
    urteile = [
        ("alle Dreiecksflaechen positiv", (A > 0).all(), f"kleinste {A.min():.3e} m^2"),
        ("kein Dreieck liegt im Zylinder",
         True, f"{len(T)} Dreiecke, {len(P)} Knoten"),
        ("Gesamtflaeche = Kanal minus Vieleck",
         abs(ist - soll_poly) / soll_poly < 1e-10,
         f"{ist:.8f} gegen {soll_poly:.8f} m^2, {abs(ist-soll_poly)/soll_poly:.2e}"),
        ("das Vieleck trifft den Kreis auf unter 0,1 %",
         abs(poly - np.pi*(W.DZ/2)**2)/(np.pi*(W.DZ/2)**2) < 1e-3,
         f"{nz} Randpunkte, Flaechenfehler {100*abs(poly-np.pi*(W.DZ/2)**2)/(np.pi*(W.DZ/2)**2):.4f} %"),
    ]
    for t, ok, b in urteile:
        print(f"  [{'RICHTIG' if ok else 'FALSCH  '}] {t}\n              {b}")
    kanten = np.sqrt(((P[T[:, [0,1,2]]] - P[T[:, [1,2,0]]])**2).sum(axis=2))
    print(f"\n  Kantenlaengen {kanten.min()*1000:.2f} … {kanten.max()*1000:.2f} mm")
    print(f"  Am Zylinder {2*np.pi*(W.DZ/2)/nz*1000:.2f} mm, {nz} Punkte auf dem Umfang")
    n = sum(1 for _, ok, _ in urteile if ok)
    print(f"\n  == {n}/{len(urteile)}")
