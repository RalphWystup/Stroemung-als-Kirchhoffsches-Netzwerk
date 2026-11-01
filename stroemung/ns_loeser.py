#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Inkompressible Navier-Stokes auf dem versetzten Gitter -- ZWEI getrennte Aufbauwege.

Sein Wort 31.10.2026: „zwei getrennt aufgebaute Loeser auf dasselbe Problem ... genau so!"

Beide Wege diskretisieren DASSELBE, aber auf voellig verschiedene Art:

    Weg A  NETZWERK   Jeder Zweig wird einzeln gestempelt (Leitwert plus gesteuerte Quelle),
                      daraus entsteht eine duennbesetzte Matrix. Kirchhoff, Knotenanalyse.
    Weg B  FINITE VOLUMEN
                      Die Fluesse werden unmittelbar als Feldoperationen gebildet, ohne
                      Matrix. Klassischer Aufbau.

Gemeinsam sind nur Geometrie und Stoffwerte. Stimmen beide ueberein, ist das ein Nachweis,
kein Zufall -- die Sprachprobe aus dem Feld-Labor.

Zeitverfahren: Projektionsverfahren (Chorin).
    1. Vorhersage   u* aus Konvektion und Zaehigkeit, ohne Druck
    2. Drucknetz    Laplace p = div(u*)/dt          <- reines Widerstandsnetz
    3. Korrektur    u = u* - dt grad p              <- erfuellt die Knotenregel

Schritt 2 ist genau der Operator, den das Feld-Labor schon rechnet.
"""
FASSUNG = "1.2"
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

from gitter import Gitter
from sg_stempel import bernoulli


# =====================================================================================
#  Drucknetz -- bei beiden Wegen dasselbe Widerstandsnetz, aber getrennt aufgebaut
# =====================================================================================

def drucknetz_matrix(g, periodisch_x=False, ausstrom_rechts=False,
                     u_maske=None, v_maske=None):
    """Weg A: Laplace-Operator aus ZWEIGLEITWERTEN zusammengestempelt.

    Jeder innere Zweig zwischen zwei Zellen traegt den Leitwert A_f/d_f. Die Knotenregel
    an jeder Zelle liefert eine Zeile. Das ist buchstaeblich Knotenanalyse.

    Fassung 1.1, Befund 31.10.2026 -- zwei Dinge muss das Drucknetz WISSEN, sonst kann die
    Projektion die Knotenregel nicht herstellen:

      ``ausstrom_rechts``  Am Ausstrom liegt Umgebungsdruck. Im Netzwerk heisst das: die
                           letzte Zellspalte haengt ueber einen Zweig am BEZUGSPOTENTIAL
                           (Leitwert dy/(dx/2), weil der halbe Zellabstand zaehlt). Ohne
                           diesen Zweig ist das Netz ringsum stromlos aufgehaengt, der
                           Ausstrom laesst sich nicht korrigieren, und die Divergenz bleibt
                           stehen -- gemessen 6,97 statt 1e-13.
      ``u_maske/v_maske``  Zweige, deren Flaeche im Koerper liegt, werden AUFGETRENNT
                           (Leitwert null). Sonst rechnet das Drucknetz durch den Koerper
                           hindurch, und das nachtraegliche Nullsetzen der Zweige zerstoert
                           die gerade hergestellte Knotenregel -- gemessen 17,0.

    Zellen, die dadurch ganz ohne Zweig dastehen (vollstaendig im Koerper), bekommen eine
    Einheitszeile; ihr Druck ist bedeutungslos und darf das System nicht singulaer machen.
    """
    nx, ny, dx, dy = g.nx, g.ny, g.dx, g.dy
    N = nx * ny
    nr = lambda i, j: i * ny + j
    zi, zj, zw = [], [], []
    grad = np.zeros(N)            # Zahl der angeschlossenen Zweige je Knoten

    def zweig(i1, j1, i2, j2, G):
        a, b = nr(i1, j1), nr(i2, j2)
        zi.extend([a, a, b, b]); zj.extend([a, b, b, a]); zw.extend([G, -G, G, -G])
        grad[a] += 1; grad[b] += 1

    for i in range(nx):
        for j in range(ny):
            if i + 1 < nx:
                # der Zweig zwischen (i,j) und (i+1,j) liegt auf der u-Flaeche Nr. i+1
                if u_maske is None or not u_maske[i + 1, j]:
                    zweig(i, j, i + 1, j, dy / dx)
            elif periodisch_x:
                zweig(i, j, 0, j, dy / dx)
            if j + 1 < ny:
                if v_maske is None or not v_maske[i, j + 1]:
                    zweig(i, j, i, j + 1, dx / dy)

    A = sp.coo_matrix((zw, (zi, zj)), shape=(N, N)).tolil()
    if ausstrom_rechts:
        # letzte Spalte an das Bezugspotential: nur die Diagonale, die Gegenseite ist 0
        for j in range(ny):
            a = nr(nx - 1, j)
            A[a, a] += dy / (dx / 2.0)
            grad[a] += 1
    # Knoten ohne jeden Zweig (ganz im Koerper): Einheitszeile
    for n in np.where(grad == 0)[0]:
        A[n, :] = 0.0; A[n, n] = 1.0
    if not ausstrom_rechts:
        # reines Neumann-Problem: ein Knoten wird festgehalten
        A[0, :] = 0.0; A[0, 0] = 1.0
    return A.tocsr(), (grad == 0)


def drucknetz_wirkung(g, p, periodisch_x=False):
    """Weg B: dieselbe Wirkung, aber als Feldoperation ohne Matrix."""
    nx, ny, dx, dy = g.nx, g.ny, g.dx, g.dy
    r = np.zeros_like(p)
    # Fluesse in x
    fx = np.zeros((nx + 1, ny))
    fx[1:nx, :] = (p[1:, :] - p[:-1, :]) * (dy / dx)
    if periodisch_x:
        rand = (p[0, :] - p[-1, :]) * (dy / dx)
        fx[0, :] = rand; fx[nx, :] = rand
    fy = np.zeros((nx, ny + 1))
    fy[:, 1:ny] = (p[:, 1:] - p[:, :-1]) * (dx / dy)
    r = -((fx[1:, :] - fx[:-1, :]) + (fy[:, 1:] - fy[:, :-1]))
    return r


def drucknetz_wirkung_voll(g, p, leer=None, ausstrom_rechts=False,
                           u_maske=None, v_maske=None):
    """Weg B fuer den vollen Fall: dieselbe Wirkung wie drucknetz_matrix, OHNE Matrix.

    Fassung 1.2, 31.10.2026. drucknetz_wirkung (1.1) kennt weder Ausstromrand noch
    Stoerkoerper und taugt damit nicht fuer die Wirbelstrasse. Die geprueften Fassungen
    bleiben unangetastet; dies ist die Erweiterung daneben (Hausregel „Geprueftes kapseln").

    Aufgebaut ist das hier bewusst anders als in drucknetz_matrix: dort wird Zweig fuer
    Zweig in eine Matrix gestempelt, hier werden die Flaechenfluesse als ganze Felder
    gebildet und bilanziert. Stimmen beide ueberein, ist der Operator belegt.

        r[i] = Summe ueber die Zweige am Knoten i von  G * (p_i - p_Nachbar)

    Das ist +A*p mit dem A aus drucknetz_matrix (dort steht +G auf der Diagonalen).
    """
    nx, ny, dx, dy = g.nx, g.ny, g.dx, g.dy
    # --- Leitwerte der inneren Zweige als Felder ------------------------------------
    Gx = np.zeros((nx + 1, ny)); Gx[1:nx, :] = dy / dx
    Gy = np.zeros((nx, ny + 1)); Gy[:, 1:ny] = dx / dy
    if u_maske is not None:
        Gx[u_maske] = 0.0                      # Zweig im Koerper: aufgetrennt
    if v_maske is not None:
        Gy[v_maske] = 0.0
    # --- Fluesse ueber die Flaechen --------------------------------------------------
    Fx = np.zeros((nx + 1, ny)); Fy = np.zeros((nx, ny + 1))
    Fx[1:nx, :] = Gx[1:nx, :] * (p[1:, :] - p[:-1, :])
    Fy[:, 1:ny] = Gy[:, 1:ny] * (p[:, 1:] - p[:, :-1])
    r = (Fx[:-1, :] - Fx[1:, :]) + (Fy[:, :-1] - Fy[:, 1:])
    if ausstrom_rechts:
        r[-1, :] += (dy / (dx / 2.0)) * p[-1, :]     # Zweig ans Bezugspotential
    if leer is not None:
        m = leer.reshape(nx, ny)
        r[m] = p[m]                                   # Einheitszeile
    return r


# =====================================================================================
#  Zaehigkeit und Konvektion
# =====================================================================================

def zaeh_u(g, u, nu, oben_unten="haft"):
    """Zaehigkeitsterm nu*Laplace(u) auf den u-Zweigen. Feldoperation (Weg B)."""
    nx, ny, dx, dy = g.nx, g.ny, g.dx, g.dy
    L = np.zeros_like(u)
    L[1:-1, :] += (u[2:, :] - 2 * u[1:-1, :] + u[:-2, :]) / dx**2
    # in y: Haftbedingung an Wand -> Spiegelwert -u am Rand
    uy = np.empty((u.shape[0], ny + 2))
    uy[:, 1:-1] = u
    if oben_unten == "haft":
        uy[:, 0] = -u[:, 0]          # u_Wand = 0 zwischen Zelle und Spiegel
        uy[:, -1] = -u[:, -1]
    else:                            # frei
        uy[:, 0] = u[:, 0]; uy[:, -1] = u[:, -1]
    L += (uy[:, 2:] - 2 * uy[:, 1:-1] + uy[:, :-2]) / dy**2
    return nu * L


def zaeh_v(g, v, nu):
    nx, ny, dx, dy = g.nx, g.ny, g.dx, g.dy
    L = np.zeros_like(v)
    L[:, 1:-1] += (v[:, 2:] - 2 * v[:, 1:-1] + v[:, :-2]) / dy**2
    vx = np.empty((nx + 2, v.shape[1]))
    vx[1:-1, :] = v
    vx[0, :] = v[-1, :]; vx[-1, :] = v[0, :]        # periodisch in x
    L += (vx[2:, :] - 2 * vx[1:-1, :] + vx[:-2, :]) / dx**2
    return nu * L


if __name__ == "__main__":
    print("=" * 84)
    print(" SPRACHPROBE DRUCKNETZ: Weg A (Matrix aus Zweigstempeln) gegen Weg B (Feldoperation)")
    print("=" * 84)
    rng = np.random.default_rng(7)
    schlecht = 0
    for (nx, ny) in ((8, 6), (17, 11), (40, 24)):
        g = Gitter(1.3, 0.7, nx, ny)
        A, _ = drucknetz_matrix(g)
        p = rng.standard_normal((nx, ny))
        rA = (A @ p.ravel()).reshape(nx, ny)
        rB = drucknetz_wirkung(g, p)
        # Zeile 0 ist bei Weg A der festgehaltene Knoten -- dort wird nicht verglichen
        d = np.abs(rA - rB); d[0, 0] = 0.0
        rel = d.max() / max(np.abs(rB).max(), 1e-30)
        ok = rel < 1e-12
        schlecht += 0 if ok else 1
        print(f"  [{'RICHTIG' if ok else 'FALSCH  '}] {nx:3d}x{ny:<3d}  "
              f"groesste Abweichung {d.max():.3e}  relativ {rel:.3e}")
    print(f"\n  == {'beide Wege bauen denselben Operator' if not schlecht else 'ABWEICHUNG'}")
    raise SystemExit(0 if not schlecht else 1)
