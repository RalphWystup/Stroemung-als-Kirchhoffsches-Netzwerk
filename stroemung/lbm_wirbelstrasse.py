#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Die Wirbelstrasse mit dem UNABHAENGIGEN Loeser -- Lattice-Boltzmann.

Derselbe Fall wie wirbelstrasse.py (Schaefer-Turek 2D-2, Re = 100), aber mit einem
Verfahren, das mit dem Netzwerkloeser nichts teilt. Alle Benchmarkgroessen sind
dimensionslos und werden in GITTEREINHEITEN gebildet; eine Umrechnung findet nicht statt
und kann deshalb auch nicht schiefgehen.

    H        Zellen ueber der Kanalhoehe; daraus folgt alles andere
    u_quer   mittlere Geschwindigkeit in Gittereinheiten (Mach-Zahl klein halten)
    D_l      = D/Ly * H                 Zylinderdurchmesser in Zellen
    nu_l     = u_quer * D_l / Re        damit ist Re getroffen
    tau      = 0,5 + 3 nu_l

Rand: Einstrom nach Zou-He mit dem parabolischen Profil, Ausstrom mit Nullgradient,
Waende und Zylinder mit Rueckprall auf halbem Weg. Die Kraft kommt aus dem
Impulsaustausch an den Verbindungen.

Aufruf:  python3 lbm_wirbelstrasse.py [H] [t_ende in s] [u_quer]
"""
FASSUNG = "1.0"
#
# STAND 31.10.2026 -- NICHT WEITERVERFOLGT, und zwar aus zwei Gruenden:
#
#   1. Lattice-Boltzmann gibt die Navier-Stokes-Gleichungen nur im Grenzfall kleiner
#      Mach-Zahl wieder; der Fehler ist von der Ordnung Ma^2, hier rund 0,75 %. Als
#      Gegenprobe fuer eine Groesse, die auf +-0,3 % angegeben ist, bringt das eine
#      ZUSAETZLICHE Naeherungsschicht mit, statt nur anders zu rechnen.
#   2. Die Unabhaengigkeit zweier Wege sitzt in der HERLEITUNG des Ansatzes, nicht im
#      Gleichungsloeser. Die Paarung, die zu diesem Projekt gehoert, ist Flussbilanz am
#      Kontrollvolumen (Kirchhoff) gegen Minimierung ueber Ansatzfunktionen (Ritz /
#      Galerkin) -- so wie im Feld-Labor. Die bauen wir, nicht ein kinetisches Verfahren.
#
# Hinzu kam, dass der Lauf bei tau = 0,515 instabil wurde (BGK ist dicht an tau = 0,5
# anfaellig). Das waere mit MRT oder regularisiertem Stoss zu beheben -- es lohnt aber
# nicht, solange Grund 1 und 2 stehen.
#
# lbm.py selbst bleibt: Seine beiden analytischen Proben (Taylor-Green, Hagen-Poiseuille)
# gehen auf, und das Modul ist damit ein geprueftes Werkzeug fuer spaeter.
import json
import sys
import time

import numpy as np

from lbm import E, W, GEGEN, CS2, gleichgewicht, stoss_stroemen, kraft_impulsaustausch
import wirbelstrasse as W_FALL          # LX, LY, DZ, XC, YC, RE, REFERENZ


def aufbau(H, u_quer):
    dx = W_FALL.LY / H
    ny = H + 2                                   # je eine Festreihe oben und unten
    nx = int(round(W_FALL.LX / dx)) + 1
    D_l = W_FALL.DZ / dx
    nu_l = u_quer * D_l / W_FALL.RE
    tau = 0.5 + nu_l / CS2
    j = np.arange(ny)
    y_phys = (j - 0.5) * dx                      # Wand auf halbem Weg: y=0 bei j=0,5
    x_phys = np.arange(nx) * dx
    fest = np.zeros((nx, ny), bool)
    fest[:, 0] = True; fest[:, -1] = True
    X, Y = np.meshgrid(x_phys, y_phys, indexing="ij")
    fest |= ((X - W_FALL.XC) ** 2 + (Y - W_FALL.YC) ** 2) < (W_FALL.DZ / 2) ** 2
    fest[0, 1:-1] = False; fest[-1, 1:-1] = False
    y_l = j - 0.5                                # 0 an der unteren Wand, H an der oberen
    u_ein = 4.0 * 1.5 * u_quer * y_l * (H - y_l) / (H * H)
    u_ein[0] = 0.0; u_ein[-1] = 0.0
    return dict(dx=dx, nx=nx, ny=ny, D_l=D_l, nu_l=nu_l, tau=tau, fest=fest,
                u_ein=u_ein, dt=dx * u_quer / W_FALL.U_QUER, u_quer=u_quer, H=H)


def zou_he_einstrom(f, u_ein):
    """Geschwindigkeitsrand nach Zou und He am Westrand: die drei fehlenden
    Verteilungen folgen aus Dichte und Impuls, nicht aus einer Extrapolation."""
    ux = u_ein; uy = 0.0
    f0, f2, f3, f4, f6, f7 = f[0, 0], f[2, 0], f[3, 0], f[4, 0], f[6, 0], f[7, 0]
    rho = (f0 + f2 + f4 + 2.0 * (f3 + f6 + f7)) / (1.0 - ux)
    f[1, 0] = f3 + (2.0 / 3.0) * rho * ux
    f[5, 0] = f7 - 0.5 * (f2 - f4) + (1.0 / 6.0) * rho * ux + 0.5 * rho * uy
    f[8, 0] = f6 + 0.5 * (f2 - f4) + (1.0 / 6.0) * rho * ux - 0.5 * rho * uy
    return f


def lauf(H=41, t_ende=12.0, u_quer=0.05, melde=None, speichern=None):
    A = aufbau(H, u_quer)
    nx, ny, fest, tau = A["nx"], A["ny"], A["fest"], A["tau"]
    n = int(round(t_ende / A["dt"]))
    rho = np.ones((nx, ny))
    ux = np.tile(A["u_ein"], (nx, 1)); uy = np.zeros((nx, ny))
    ux[fest] = 0.0
    f = gleichgewicht(rho, ux, uy)
    nenner = 0.5 * u_quer * u_quer * A["D_l"]
    verlauf = []
    t0 = time.time()
    for k in range(n):
        f = zou_he_einstrom(f, A["u_ein"])
        f[:, -1] = f[:, -2]                                   # Ausstrom: Nullgradient
        f, fs, _, _, _ = stoss_stroemen(f, tau, fest)
        Fx, Fy = kraft_impulsaustausch(fs, fest)
        verlauf.append((k * A["dt"], Fx / nenner, Fy / nenner))
        if melde and k % melde == 0:
            print(f"    t={k*A['dt']:6.3f}  cD={Fx/nenner:7.3f}  cL={Fy/nenner:7.3f}  "
                  f"({time.time()-t0:5.0f} s)", flush=True)
    verlauf = np.array(verlauf)
    if speichern:
        rho, uxx, uyy = (f.sum(axis=0), None, None)
        np.savez_compressed(speichern, f=f.astype("f4"), fest=fest, verlauf=verlauf,
                            dx=A["dx"], dt=A["dt"], H=H)
    return A, f, verlauf


if __name__ == "__main__":
    H = int(sys.argv[1]) if len(sys.argv) > 1 else 41
    t_ende = float(sys.argv[2]) if len(sys.argv) > 2 else 12.0
    u_quer = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05
    A = aufbau(H, u_quer)
    print("=" * 92)
    print(f" WIRBELSTRASSE MIT LATTICE-BOLTZMANN -- Fassung {FASSUNG}")
    print(f" H = {H} Zellen ueber dem Kanal, Gitter {A['nx']}x{A['ny']}, "
          f"D = {A['D_l']:.1f} Zellen")
    print(f" u_quer = {u_quer} (Mach {u_quer/np.sqrt(CS2):.3f}), nu = {A['nu_l']:.6f}, "
          f"tau = {A['tau']:.4f}, dt = {A['dt']:.3e} s")
    print("=" * 92)
    A, f, verlauf = lauf(H, t_ende, u_quer, melde=max(1, int(0.5 / A["dt"])),
                         speichern=f"lbm_feld_H{H}.npz")
    cD, cL, St, T = W_FALL.auswerten(verlauf, A["dt"])
    print(f"\n  {len(verlauf)} Zeitschritte")
    print(f"\n{'Groesse':>10s} {'gerechnet':>12s} {'Referenz':>16s} {'Urteil':>14s}")
    for name, wert, (lo, hi) in (("c_D max", cD, W_FALL.REFERENZ["cD"]),
                                 ("c_L max", cL, W_FALL.REFERENZ["cL"]),
                                 ("St", St, W_FALL.REFERENZ["St"])):
        print(f"{name:>10s} {wert:12.4f} {f'{lo} … {hi}':>16s} "
              f"{'im Intervall' if lo <= wert <= hi else 'daneben':>14s}")
    json.dump({"H": H, "nx": A["nx"], "ny": A["ny"], "D_l": A["D_l"], "tau": A["tau"],
               "u_quer": u_quer, "cD": cD, "cL": cL, "St": St, "T": T, "dt": A["dt"]},
              open(f"ergebnis_lbm_H{H}.json", "w"), indent=1)
