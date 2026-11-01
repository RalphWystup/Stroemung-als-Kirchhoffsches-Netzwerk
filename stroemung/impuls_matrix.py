#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Weg A: der Impulsoperator als MATRIX aus Zweigstempeln -- Knotenanalyse.

Das Gegenstueck zu konvektion.py (Weg B, Feldoperationen). Beide diskretisieren dasselbe,
werden aber voellig verschieden aufgebaut:

    Weg B   Fluesse als Feldoperationen, keine Matrix, geschoben mit numpy-Schnitten
    Weg A   jeder Zweig wird EINZELN gestempelt: zwei Leitwerte, die sich um den
            Volumenstrom unterscheiden, in eine duennbesetzte Matrix eingetragen

Das ist die Sprachprobe aus dem Feld-Labor: Wenn zwei getrennt gebaute Wege dieselbe
Wirkung ergeben, ist der Operator belegt und nicht nur plausibel.

In der Netzwerksprache ist Weg A das Eigentliche: eine Knotenadmittanzmatrix, deren
UNSYMMETRIE genau die gesteuerten Quellen traegt. Fuer u_f -> 0 wird sie symmetrisch und
das Netz ist das passive Widerstandsnetz des Feld-Labors.

Randbehandlung: Ein Geisterwert, der eine Kopie (freier Ausstrom) oder eine Spiegelung
(Haftbedingung, Faktor -1) eines INNEREN Unbekannten ist, wird als Matrixeintrag auf
diesen Knoten gelegt -- nicht als Konstante in b. Nur das vorgegebene Einstroemprofil ist
eine echte Konstante. Damit ist M die vollstaendige Operatormatrix und nicht bloss eine
am aktuellen Zustand abgelesene Naeherung.
"""
FASSUNG = "1.1"
import numpy as np
import scipy.sparse as sp

from gitter import Gitter
from konvektion import _rand_u, _rand_v, sg_leitwerte


class _Stempelbuch:
    """Sammelt Matrixeintraege und Konstanten; ein Zweig nach dem anderen."""

    def __init__(self, n):
        self.zi, self.zj, self.zw = [], [], []
        self.b = np.zeros(n)
        self.n = n

    def knoten(self, a, j, w):
        self.zi.append(a); self.zj.append(j); self.zw.append(w)

    def fest(self, a, w):
        self.b[a] += w

    def matrix(self):
        return sp.coo_matrix((self.zw, (self.zi, self.zj)),
                             shape=(self.n, self.n)).tocsr()


def impuls_matrix_u(g, u, v, nu, einstroem=None):
    """M, b fuer die u-Zweige, so dass M @ u.ravel() + b = konv_zaeh_u(...).ravel().

    u hat die Form (nx+1, ny); Knotennummer a = i*ny + j.
    """
    nx, ny, dx, dy = g.nx, g.ny, g.dx, g.dy
    U = _rand_u(g, u, einstroem)
    V = _rand_v(g, v)
    uf_x = 0.5 * (U[:-1, 1:-1] + U[1:, 1:-1])      # (nx+2, ny)  x-Flaechen
    vf_y = 0.5 * (V[:-1, 1:-1] + V[1:, 1:-1])      # (nx+1, ny+1) y-Flaechen
    ghx, grx = sg_leitwerte(uf_x, nu, dx, dy)
    ghy, gry = sg_leitwerte(vf_y, nu, dy, dx)
    iVz = 1.0 / (dx * dy)

    B = _Stempelbuch((nx + 1) * ny)
    nr = lambda i, j: i * ny + j

    for i in range(nx + 1):
        for j in range(ny):
            a = nr(i, j)
            # ---- x-Zweige -------------------------------------------------------
            # Zufluss ueber die linke Flaeche i:  + gh*U[i,j+1] - gr*u[i,j]
            gh, gr = ghx[i, j] * iVz, grx[i, j] * iVz
            B.knoten(a, a, -gr)
            if i >= 1:
                B.knoten(a, nr(i - 1, j), gh)                 # U[i,j+1] = u[i-1,j]
            elif einstroem is None:
                B.knoten(a, nr(0, j), gh)                     # Geist = Kopie von u[0,j]
            else:
                B.fest(a, gh * float(einstroem[j]))           # vorgegeben: echte Konstante
            # Abfluss ueber die rechte Flaeche i+1: - gh*u[i,j] + gr*U[i+2,j+1]
            gh, gr = ghx[i + 1, j] * iVz, grx[i + 1, j] * iVz
            B.knoten(a, a, -gh)
            if i <= nx - 1:
                B.knoten(a, nr(i + 1, j), gr)                 # U[i+2,j+1] = u[i+1,j]
            else:
                B.knoten(a, nr(nx, j), gr)                    # freier Ausstrom: Kopie
            # ---- y-Zweige -------------------------------------------------------
            # Zufluss ueber die untere Flaeche j: + gh*U[i+1,j] - gr*u[i,j]
            gh, gr = ghy[i, j] * iVz, gry[i, j] * iVz
            B.knoten(a, a, -gr)
            if j >= 1:
                B.knoten(a, nr(i, j - 1), gh)
            else:
                B.knoten(a, nr(i, 0), -gh)                    # Haftbedingung: Spiegel, -1
            # Abfluss ueber die obere Flaeche j+1: - gh*u[i,j] + gr*U[i+1,j+2]
            gh, gr = ghy[i, j + 1] * iVz, gry[i, j + 1] * iVz
            B.knoten(a, a, -gh)
            if j <= ny - 2:
                B.knoten(a, nr(i, j + 1), gr)
            else:
                B.knoten(a, nr(i, ny - 1), -gr)               # Haftbedingung oben
    return B.matrix(), B.b


def impuls_matrix_v(g, u, v, nu):
    """M, b fuer die v-Zweige, so dass M @ v.ravel() + b = konv_zaeh_v(...).ravel().

    v hat die Form (nx, ny+1); Knotennummer a = i*(ny+1) + j.
    """
    nx, ny, dx, dy = g.nx, g.ny, g.dx, g.dy
    U = _rand_u(g, u)
    V = _rand_v(g, v)
    uf_x = 0.5 * (U[1:-1, :-1] + U[1:-1, 1:])      # (nx+1, ny+1)
    vf_y = 0.5 * (V[1:-1, :-1] + V[1:-1, 1:])      # (nx, ny+2)
    ghx, grx = sg_leitwerte(uf_x, nu, dx, dy)
    ghy, gry = sg_leitwerte(vf_y, nu, dy, dx)
    iVz = 1.0 / (dx * dy)

    m = ny + 1
    B = _Stempelbuch(nx * m)
    nr = lambda i, j: i * m + j

    for i in range(nx):
        for j in range(m):
            a = nr(i, j)
            # ---- x-Zweige: Flaeche i (links) und i+1 (rechts) --------------------
            gh, gr = ghx[i, j] * iVz, grx[i, j] * iVz
            B.knoten(a, a, -gr)                               # V[i+1,j+1] = v[i,j]
            if i >= 1:
                B.knoten(a, nr(i - 1, j), gh)                 # V[i,j+1] = v[i-1,j]
            else:
                B.knoten(a, nr(0, j), gh)                     # V[0,*] = V[1,*]: Kopie
            gh, gr = ghx[i + 1, j] * iVz, grx[i + 1, j] * iVz
            B.knoten(a, a, -gh)
            if i <= nx - 2:
                B.knoten(a, nr(i + 1, j), gr)                 # V[i+2,j+1] = v[i+1,j]
            else:
                B.knoten(a, nr(nx - 1, j), gr)                # V[nx+1,*] = V[nx,*]
            # ---- y-Zweige: Flaeche j (unten) und j+1 (oben) ----------------------
            gh, gr = ghy[i, j] * iVz, gry[i, j] * iVz
            B.knoten(a, a, -gr)                               # V[i+1,j+1] = v[i,j]
            if j >= 1:
                B.knoten(a, nr(i, j - 1), gh)                 # V[i+1,j] = v[i,j-1]
            else:
                B.knoten(a, nr(i, 0), -gh)                    # V[i+1,0] = -v[i,0]
            gh, gr = ghy[i, j + 1] * iVz, gry[i, j + 1] * iVz
            B.knoten(a, a, -gh)
            if j <= m - 2:
                B.knoten(a, nr(i, j + 1), gr)                 # V[i+1,j+2] = v[i,j+1]
            else:
                B.knoten(a, nr(i, m - 1), -gr)                # V[i+1,-1] = -v[i,-1]
    return B.matrix(), B.b


if __name__ == "__main__":
    from konvektion import konv_zaeh_u, konv_zaeh_v
    print("=" * 94)
    print(" SPRACHPROBE IMPULSOPERATOR -- Fassung " + FASSUNG)
    print(" Weg A (Matrix, Zweig fuer Zweig gestempelt) gegen Weg B (Feldoperationen)")
    print("=" * 94)
    rng = np.random.default_rng(23)
    urteile = []
    for (nx, ny, nu) in ((12, 8, 1e-1), (25, 13, 1e-2), (40, 20, 1e-3)):
        g = Gitter(2.2, 0.41, nx, ny)
        u = rng.standard_normal(g.leer_u().shape) * 0.7
        v = rng.standard_normal(g.leer_v().shape) * 0.7
        for name, ein in (("u-Zweige, Einstroemprofil vorgegeben", rng.standard_normal(ny) * 0.3),
                          ("u-Zweige, freier Einstrom", None)):
            rB = konv_zaeh_u(g, u, v, nu, einstroem=ein)
            M, b = impuls_matrix_u(g, u, v, nu, einstroem=ein)
            rA = (M @ u.ravel() + b).reshape(nx + 1, ny)
            rel = np.abs(rA - rB).max() / max(np.abs(rB).max(), 1e-30)
            urteile.append((f"{nx:3d}x{ny:<3d} nu={nu:<6g} {name}", rel < 1e-11, rel))
        rB = konv_zaeh_v(g, u, v, nu)
        M, b = impuls_matrix_v(g, u, v, nu)
        rA = (M @ v.ravel() + b).reshape(nx, ny + 1)
        rel = np.abs(rA - rB).max() / max(np.abs(rB).max(), 1e-30)
        urteile.append((f"{nx:3d}x{ny:<3d} nu={nu:<6g} v-Zweige", rel < 1e-11, rel))
        # Die Unsymmetrie ist die gesteuerte Quelle: bei verschwindender Stroemung
        # muss die Matrix symmetrisch werden (passives Widerstandsnetz).
        M0, _ = impuls_matrix_u(g, u * 1e-9, v * 1e-9, nu, einstroem=None)
        uns = abs(M0 - M0.T).max() / abs(M0).max()
        urteile.append((f"{nx:3d}x{ny:<3d} nu={nu:<6g} u->0 macht die Matrix symmetrisch",
                        uns < 1e-7, uns))
    for text, ok, z in urteile:
        print(f"  [{'RICHTIG' if ok else 'FALSCH  '}] {text}   -> {z:.3e}")
    n = sum(1 for _, ok, _ in urteile if ok)
    print(f"\n  == {n}/{len(urteile)} " +
          ("beide Wege ergeben denselben Impulsoperator" if n == len(urteile) else "ABWEICHUNG"))
    raise SystemExit(0 if n == len(urteile) else 1)
