#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Wo bleibt die Zeit in einem Zeitschritt?

Anlass: Im Pruefplan steht die Aussage, der Jetson braeuchte hier nicht hinzugezogen zu
werden. Eine solche Aussage ist nur so viel wert wie ihre Messung. Gemessen wird deshalb,
welcher Anteil eines Zeitschritts auf welchen Teil entfaellt:

    Konvektion und Zaehigkeit   Feldoperationen, gut vektorisiert, parallelisierbar
    Drucknetz loesen            Vorwaerts-/Rueckwaertseinsetzen der LU-Zerlegung
    Rest                        Divergenz, Korrektur, Kraftbilanz

Nur was auf das Loesen entfaellt, koennte eine GPU ueberhaupt angreifen -- und auch das nur
dann, wenn man vom direkten Loeser auf ein iteratives Verfahren wechselt. Das Einsetzen in
eine Dreiecksmatrix ist von Natur aus sequentiell.
"""
FASSUNG = "1.0"
import sys
import time

import numpy as np
import scipy.sparse.linalg as spl

import wirbelstrasse as W
from gitter import Gitter
from ns_loeser import drucknetz_matrix
from konvektion import konv_zaeh_u, konv_zaeh_v


def messen(nx, ny, schritte=40):
    g = Gitter(W.LX, W.LY, nx, ny)
    mu_m, mv_m, mp_m = W.koerper_masken(g)
    ein = W.einstroemprofil(g.yu)
    u, v = g.leer_u(), g.leer_v()
    u[:, :] = ein[None, :]
    u[mu_m] = 0.0; v[mv_m] = 0.0
    t0 = time.perf_counter()
    A, leer = drucknetz_matrix(g, ausstrom_rechts=True, u_maske=mu_m, v_maske=mv_m)
    t_bau = time.perf_counter() - t0
    t0 = time.perf_counter()
    loes = spl.factorized(A.tocsc())
    t_zerlegung = time.perf_counter() - t0
    dt = min(0.4 * min(g.dx, g.dy) / (2.0 * W.UM),
             0.2 * min(g.dx, g.dy) ** 2 / W.NU)
    t_konv = t_loes = t_rest = 0.0
    for k in range(schritte):
        t0 = time.perf_counter()
        us = u + dt * konv_zaeh_u(g, u, v, W.NU, einstroem=ein)
        vs = v + dt * konv_zaeh_v(g, u, v, W.NU)
        t_konv += time.perf_counter() - t0
        t0 = time.perf_counter()
        vs[:, 0] = vs[:, -1] = 0.0; vs[0, :] = 0.0
        us[0, :] = ein; us[-1, :] = us[-2, :]; vs[-1, :] = vs[-2, :]
        us[mu_m] = 0.0; vs[mv_m] = 0.0
        b = -g.divergenz(us, vs).ravel() * (g.dx * g.dy) / dt
        b[leer] = 0.0
        t_rest += time.perf_counter() - t0
        t0 = time.perf_counter()
        p = loes(b).reshape(nx, ny)
        t_loes += time.perf_counter() - t0
        t0 = time.perf_counter()
        u = us.copy(); v = vs.copy()
        u[1:-1, :] -= dt * (p[1:, :] - p[:-1, :]) / g.dx
        v[:, 1:-1] -= dt * (p[:, 1:] - p[:, :-1]) / g.dy
        u[-1, :] -= dt * (0.0 - p[-1, :]) / (g.dx / 2.0)
        u[mu_m] = 0.0; v[mv_m] = 0.0
        u[0, :] = ein
        t_rest += time.perf_counter() - t0
    return t_bau, t_zerlegung, t_konv / schritte, t_loes / schritte, t_rest / schritte


if __name__ == "__main__":
    nx = int(sys.argv[1]) if len(sys.argv) > 1 else 440
    ny = int(sys.argv[2]) if len(sys.argv) > 2 else 82
    print("=" * 92)
    print(f" ZEITAUFTEILUNG EINES ZEITSCHRITTS -- Fassung {FASSUNG}   Netz {nx}x{ny}")
    print("=" * 92)
    bau, zerl, konv, loes, rest = messen(nx, ny)
    ges = konv + loes + rest
    print(f"  einmalig:  Netz stempeln {bau*1000:7.1f} ms,  LU-Zerlegung {zerl*1000:7.1f} ms")
    print(f"\n  je Zeitschritt (Mittel aus 40):")
    for name, w in (("Konvektion und Zaehigkeit", konv), ("Drucknetz loesen", loes),
                    ("Rest (Divergenz, Korrektur, Kraft)", rest)):
        print(f"      {name:36s} {w*1000:7.2f} ms   {100*w/ges:5.1f} %")
    print(f"      {'zusammen':36s} {ges*1000:7.2f} ms")
    print(f"\n  Nur der Loeseanteil ({100*loes/ges:.1f} %) kaeme fuer eine GPU ueberhaupt in Frage.")
    print(f"  Selbst wenn er auf null faellt, wird der Schritt hoechstens "
          f"{1/(1-loes/ges):.2f}-mal schneller.")
