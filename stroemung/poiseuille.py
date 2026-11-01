#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Stufe 2: der Loeser gegen Hagen-Poiseuille -- eine EXAKTE Loesung der NS-Gleichungen.

Aufbau: ebener Kanal der Hoehe h, in x periodisch, Haftbedingung oben und unten,
angetrieben durch eine Volumenkraft G = -dp/dx. Im ausgebildeten Zustand verschwindet der
konvektive Term identisch; geprueft werden also Zaehigkeit, Haftbedingung, Zeitintegration
und die Druckkorrektur -- ohne dass die Konvektion das Ergebnis verdecken koennte.

Exakte Loesung:

    u(y) = (G / 2 nu) * y * (h - y)        u_max = G h^2 / (8 nu)      u_quer = G h^2 / (12 nu)

Das ist keine Naeherung und kein Benchmark, sondern eine geschlossene Loesung der
vollstaendigen Gleichungen. Trifft der Loeser sie nicht, ist er falsch.

Geprueft wird zweierlei:
  * das PROFIL punktweise (Form)
  * der MITTELWERT (Integral)  -- beides muss stimmen, nicht nur eines
Dazu die Ordnung ueber eine Netzreihe: der Fehler muss mit h^2 fallen.
"""
FASSUNG = "1.0"
import numpy as np
import scipy.sparse.linalg as spl

from gitter import Gitter
from ns_loeser import drucknetz_matrix, zaeh_u, zaeh_v


def rechne(ny, nu=1.0, h=1.0, G=1.0, nx=4, schritte=200000, genau=1e-12, melde=False):
    """Bis zum stationaeren Zustand rechnen und u(y) zurueckgeben."""
    g = Gitter(h * nx / ny, h, nx, ny)          # x-Laenge beliebig (periodisch)
    u, v = g.leer_u(), g.leer_v()
    A, _ = drucknetz_matrix(g, periodisch_x=True)
    loes = spl.factorized(A.tocsc())
    dt = 0.2 * min(g.dx, g.dy) ** 2 / nu        # Stabilitaetsgrenze der Diffusion
    vor = None
    for k in range(schritte):
        # 1 Vorhersage: Zaehigkeit und Antrieb (Konvektion ist hier identisch null)
        us = u + dt * (zaeh_u(g, u, nu) + G)
        vs = v + dt * zaeh_v(g, v, nu)
        us[0, :] = us[-1, :] = 0.5 * (us[0, :] + us[-1, :])     # periodisch
        vs[:, 0] = vs[:, -1] = 0.0                               # Wand undurchlaessig
        # 2 Drucknetz: Knotenregel herstellen
        d = g.divergenz(us, vs)
        # Vorzeichen und Zellvolumen, Befund 31.10.2026:
        # drucknetz_matrix baut A = -Laplace * V (Leitwerte A_f/d_f enthalten die Flaeche
        # schon). Verlangt ist Laplace p = div(u*)/dt, also A p = -div * V / dt.
        # Fassung 1.0 schrieb b = div/dt -- falsches Vorzeichen UND der Faktor V = dx*dy
        # (hier 3,9e-4) fehlte. Der Druck war damit 2560-fach zu gross und wirkte in die
        # falsche Richtung; der Loeser explodierte im ersten Zeitschritt.
        #
        # Warum Stufe 2 das nicht gezeigt hat: In der ausgebildeten Kanalstroemung ist die
        # Divergenz identisch null, der Druckterm also wirkungslos. Hagen-Poiseuille prueft
        # Zaehigkeit und Haftbedingung -- die Druckkopplung prueft erst dieser Fall.
        b = -d.ravel() * (g.dx * g.dy) / dt; b[0] = 0.0
        p = loes(b).reshape(g.nx, g.ny)
        # 3 Korrektur
        u = us.copy(); vneu = vs.copy()
        u[1:-1, :] -= dt * (p[1:, :] - p[:-1, :]) / g.dx
        rand = dt * (p[0, :] - p[-1, :]) / g.dx
        u[0, :] -= rand; u[-1, :] -= rand
        vneu[:, 1:-1] -= dt * (p[:, 1:] - p[:, :-1]) / g.dy
        v = vneu
        if k % 200 == 0:
            m = float(np.abs(u).max())
            if vor is not None and abs(m - vor) < genau * max(m, 1e-30):
                if melde:
                    print(f"      stationaer nach {k} Schritten")
                break
            vor = m
    return g, u, v


def exakt(y, G, nu, h):
    return G / (2.0 * nu) * y * (h - y)


if __name__ == "__main__":
    print("=" * 88)
    print(" STUFE 2: HAGEN-POISEUILLE -- exakte Loesung der Navier-Stokes-Gleichungen")
    print(" ebener Kanal, Haftbedingung, Antrieb G;  u(y) = G/(2nu) y (h-y)")
    print("=" * 88)
    G, nu, h = 1.0, 1.0, 1.0
    u_max_ex = G * h * h / (8 * nu)
    u_quer_ex = G * h * h / (12 * nu)
    print(f"\n exakt:  u_max = {u_max_ex:.10f}   u_quer = {u_quer_ex:.10f}\n")
    print(f"{'ny':>5s} {'u_max':>13s} {'Fehler max':>12s} {'u_quer':>13s} "
          f"{'Fehler quer':>12s} {'Profil L2':>11s} {'Ordnung':>8s}")
    vorher = None
    zeilen = []
    for ny in (8, 16, 32, 64):
        g, u, v = rechne(ny, nu=nu, h=h, G=G)
        uy = u.mean(axis=0)                        # ueber x gemittelt (ist konstant in x)
        ue = exakt(g.yu, G, nu, h)
        f_prof = float(np.sqrt(np.mean((uy - ue) ** 2)))
        u_max = float(uy.max()); u_quer = float(uy.mean())
        ordn = "" if vorher is None else f"{np.log2(vorher / f_prof):8.2f}"
        print(f"{ny:5d} {u_max:13.9f} {abs(u_max-u_max_ex):12.3e} {u_quer:13.9f} "
              f"{abs(u_quer-u_quer_ex):12.3e} {f_prof:11.3e} {ordn:>8s}")
        vorher = f_prof
        zeilen.append((ny, f_prof, abs(u_quer - u_quer_ex)))

    # Urteil: zweite Ordnung und feinstes Netz unter 1e-4
    ordnungen = [np.log2(zeilen[i][1] / zeilen[i + 1][1]) for i in range(len(zeilen) - 1)]
    o_ok = all(o > 1.8 for o in ordnungen)
    f_ok = zeilen[-1][1] < 1e-4
    print()
    print(f"  [{'RICHTIG' if o_ok else 'FALSCH  '}] Fehler faellt mit h^2 "
          f"(Ordnungen {', '.join(f'{o:.2f}' for o in ordnungen)}, alle > 1,8)")
    print(f"  [{'RICHTIG' if f_ok else 'FALSCH  '}] feinstes Netz unter 1e-4 "
          f"({zeilen[-1][1]:.3e})")
    raise SystemExit(0 if (o_ok and f_ok) else 1)
