#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Erzeugt die Vergleichsdaten, gegen die der Browserkern geprueft wird.

Der Rechenkern des Labors ist in JavaScript neu geschrieben. Damit ist er ein ZWEITES
Programm -- und ein zweites Programm ist so lange eine zweite Fehlerquelle, bis es gegen
das erste gehalten wurde. Dieses Programm rechnet mit den geprueften Python-Modulen eine
Reihe von Faellen und schreibt Eingang und Ausgang nach vergleichsdaten.json.
pruefe_kern.js rechnet dieselben Faelle im Kern und vergleicht Zahl fuer Zahl.

Die Zufallsfelder werden HIER erzeugt und mitgeschrieben, damit beide Seiten wirklich
dasselbe rechnen. Ein eigener Zufall auf jeder Seite waere kein Vergleich.
"""
FASSUNG = "1.0"
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "stroemung"))

from sg_stempel import bernoulli, loese, genau            # noqa: E402
from gitter import Gitter                                  # noqa: E402
from konvektion import konv_zaeh_u, konv_zaeh_v            # noqa: E402
from impuls_matrix import impuls_matrix_u, impuls_matrix_v  # noqa: E402
from ns_loeser import drucknetz_matrix, drucknetz_wirkung_voll  # noqa: E402
import wirbelstrasse as W                                   # noqa: E402

aus = {"fassung": FASSUNG, "faelle": []}

# --- 1 Bernoulli ueber den ganzen Zahlenbereich -------------------------------------
xs = [0.0, 1e-12, -1e-12, 1e-6, 0.5, -0.5, 1.0, -1.0, 37.0, -37.0, 499.0, -499.0,
      501.0, -501.0, 1e6, -1e6]
aus["faelle"].append({"art": "bernoulli", "x": xs,
                      "y": [float(bernoulli(np.array([x]))[0]) for x in xs]})

# --- 2 Stationaere 1D-Konvektion-Diffusion, drei Stempel ----------------------------
for Pe in (1.0, 5.0, 20.0, 100.0, 1e3, 1e6):
    n, L, D = 20, 1.0, 1.0
    u = Pe * D / L
    for st in ("sg", "zentral", "aufwind"):
        aus["faelle"].append({"art": "loese1d", "n": n, "L": L, "u": u, "D": D,
                              "stempel": st, "c": loese(n, L, u, D, st).tolist()})
    aus["faelle"].append({"art": "genau1d", "L": L, "Pe": Pe,
                          "x": [i * L / n for i in range(n + 1)],
                          "c": [float(genau(i * L / n, L, Pe)) for i in range(n + 1)]})

# --- 3 Feldoperationen und Matrizen auf einem festen Zufallsfeld --------------------
rng = np.random.default_rng(4711)
for (nx, ny, nu) in ((12, 8, 1e-1), (25, 13, 1e-2)):
    g = Gitter(W.LX, W.LY, nx, ny)
    u = rng.standard_normal(g.leer_u().shape) * 0.7
    v = rng.standard_normal(g.leer_v().shape) * 0.7
    p = rng.standard_normal((nx, ny))
    ein = W.einstroemprofil(g.yu)
    mu, mv, mp = W.koerper_masken(g)
    A_p, leer = drucknetz_matrix(g, ausstrom_rechts=True, u_maske=mu, v_maske=mv)
    Mu, bu = impuls_matrix_u(g, u, v, nu, einstroem=ein)
    Mv, bv = impuls_matrix_v(g, u, v, nu)
    # Auch der Fall OHNE Ausstrom und OHNE Koerper -- die Seite bietet ihn an,
    # also muss er geprueft sein (Befund B11).
    A_nackt, leer_nackt = drucknetz_matrix(g)
    aus["faelle"].append({
        "art": "drucknetz_nackt", "nx": nx, "ny": ny, "Lx": W.LX, "Ly": W.LY,
        "p": p.ravel().tolist(),
        "A": (A_nackt @ p.ravel()).tolist(),
        "B": drucknetz_wirkung_voll(g, p, leer=leer_nackt, ausstrom_rechts=False).ravel().tolist(),
    })
    aus["faelle"].append({
        "art": "felder", "nx": nx, "ny": ny, "nu": nu,
        "Lx": W.LX, "Ly": W.LY,
        "u": u.ravel().tolist(), "v": v.ravel().tolist(), "p": p.ravel().tolist(),
        "ein": ein.tolist(),
        "mu": mu.ravel().astype(int).tolist(), "mv": mv.ravel().astype(int).tolist(),
        "leer": leer.astype(int).tolist(),
        "konv_u": konv_zaeh_u(g, u, v, nu, einstroem=ein).ravel().tolist(),
        "konv_v": konv_zaeh_v(g, u, v, nu).ravel().tolist(),
        "impuls_u": (Mu @ u.ravel() + bu).tolist(),
        "impuls_v": (Mv @ v.ravel() + bv).tolist(),
        "druck_A": (A_p @ p.ravel()).tolist(),
        "druck_B": drucknetz_wirkung_voll(g, p, leer=leer, ausstrom_rechts=True,
                                          u_maske=mu, v_maske=mv).ravel().tolist(),
        "divergenz": g.divergenz(u, v).ravel().tolist(),
    })

ziel = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vergleichsdaten.json")
with open(ziel, "w") as f:
    json.dump(aus, f)
print(f"  {len(aus['faelle'])} Faelle nach {os.path.basename(ziel)} "
      f"({os.path.getsize(ziel)/1024:.0f} KB)")
