#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Konvektion auf dem versetzten Gitter -- in beiden Aufbauwegen.

Die Konvektion ist der Teil, den ein passives Widerstandsnetz NICHT kann: Sie ist gerichtet
und erzeugt eine unsymmetrische Matrix. Im Netzwerk ist das die gesteuerte Quelle.

Geschrieben wird der Term in Divergenzform, weil nur sie die Erhaltung sichert:

    (u grad)u  =  div(u u)  -  u div(u)        und div(u) = 0

Auf dem versetzten Gitter heisst das: Fuer den u-Zweig wird ueber die Flaechen der
u-Kontrollzelle bilanziert, die zwischen zwei Druckzellen liegt.

Der Zweigstempel ist Scharfetter-Gummel aus sg_stempel.py: exakt bei jeder Peclet-Zahl.
Fuer die Impulsgleichung ist nu die Diffusivitaet, Pe_h = u*h/nu.
"""
FASSUNG = "1.1"
import numpy as np

from gitter import Gitter
from sg_stempel import bernoulli


def _rand_u(g, u, einstroem=None):
    """u mit einer Geisterschicht ringsum (Haft oben/unten, vorgegeben links, frei rechts)."""
    nx, ny = g.nx, g.ny
    U = np.empty((nx + 3, ny + 2))
    U[1:-1, 1:-1] = u
    U[1:-1, 0] = -u[:, 0]                 # Haftbedingung unten: u_Wand = 0
    U[1:-1, -1] = -u[:, -1]               # oben
    U[0, :] = U[1, :] if einstroem is None else np.r_[0, einstroem, 0]
    U[-1, :] = U[-2, :]                   # freier Ausstrom
    return U


def _rand_v(g, v):
    nx, ny = g.nx, g.ny
    V = np.empty((nx + 2, ny + 3))
    V[1:-1, 1:-1] = v
    V[1:-1, 0] = -v[:, 0]
    V[1:-1, -1] = -v[:, -1]
    V[0, :] = V[1, :]
    V[-1, :] = V[-2, :]
    return V


def sg_leitwerte(u_f, nu, d, A):
    """Zweigleitwerte nach Scharfetter-Gummel aus der FLAECHENGESCHWINDIGKEIT u_f.

    Fassung 1.1, Befund 31.10.2026: Fassung 1.0 nahm hier den Volumenstrom u*A statt der
    Geschwindigkeit. Die Peclet-Zahl war damit um den Faktor A zu klein (bei dy = 0,0195
    also fuenfzigfach), der Stempel entartete zu reiner Diffusion und der konvektive Anteil
    war um denselben Faktor zu schwach. Aufgefallen ist es erst an der Eigenschaftsprobe
    g_hin - g_rueck = u_f*A, die in Fassung 1.0 nur deshalb bestand, weil dort A = 1
    gesetzt war -- eine Pruefung, die ihren eigenen Fehler nicht sehen konnte.

        Pe = u_f * d / nu            (dimensionslos, wie es sein muss)
        g  = (nu*A/d) * B(-+Pe)
        g_hin - g_rueck = (nu*A/d)*Pe = u_f*A = Volumenstrom      <- die gesteuerte Quelle
    """
    Pe = u_f * d / nu
    g_hin = (nu * A / d) * bernoulli(-Pe)
    g_rueck = (nu * A / d) * bernoulli(Pe)
    return g_hin, g_rueck


def konv_zaeh_u(g, u, v, nu, einstroem=None):
    """Konvektion und Zaehigkeit auf den u-Zweigen, als Feldoperation (Weg B).

    Rueckgabe: d(u)/dt-Beitrag, also (Zufluss - Abfluss)/Volumen.
    """
    nx, ny, dx, dy = g.nx, g.ny, g.dx, g.dy
    U = _rand_u(g, u, einstroem)
    V = _rand_v(g, v)
    # Massenstroeme durch die Flaechen der u-Kontrollzelle
    # x-Flaechen liegen in den Druckzellmitten: u dort = Mittel der beiden Nachbarzweige
    uf_x = 0.5 * (U[:-1, 1:-1] + U[1:, 1:-1])               # (nx+2, ny)
    # y-Flaechen: v an der Ecke = Mittel in x
    vf_y = 0.5 * (V[:-1, 1:-1] + V[1:, 1:-1])               # (nx+1, ny+1)

    gh, gr = sg_leitwerte(uf_x, nu, dx, dy)
    # Fluss durch die x-Flaechen: J = gh*u_links - gr*u_rechts
    Jx = gh * U[:-1, 1:-1] - gr * U[1:, 1:-1]               # (nx+2, ny)
    gh, gr = sg_leitwerte(vf_y, nu, dy, dx)
    Jy = gh * U[1:-1, :-1] - gr * U[1:-1, 1:]               # (nx+1, ny+1)

    V_zelle = dx * dy
    return (Jx[:-1, :] - Jx[1:, :] + Jy[:, :-1] - Jy[:, 1:]) / V_zelle


def konv_zaeh_v(g, u, v, nu):
    nx, ny, dx, dy = g.nx, g.ny, g.dx, g.dy
    U = _rand_u(g, u)
    V = _rand_v(g, v)
    uf_x = 0.5 * (U[1:-1, :-1] + U[1:-1, 1:])               # (nx+1, ny+1)
    vf_y = 0.5 * (V[1:-1, :-1] + V[1:-1, 1:])               # (nx, ny+2)
    gh, gr = sg_leitwerte(uf_x, nu, dx, dy)
    Jx = gh * V[:-1, 1:-1] - gr * V[1:, 1:-1]               # (nx+1, ny+1)
    gh, gr = sg_leitwerte(vf_y, nu, dy, dx)
    Jy = gh * V[1:-1, :-1] - gr * V[1:-1, 1:]               # (nx, ny+2)
    return (Jx[:-1, :] - Jx[1:, :] + Jy[:, :-1] - Jy[:, 1:]) / (dx * dy)


if __name__ == "__main__":
    print("=" * 86)
    print(" KONVEKTIONSSTEMPEL -- Eigenschaften, die ohne Rechnung gelten muessen")
    print("=" * 86)
    g = Gitter(2.0, 1.0, 24, 12)
    rng = np.random.default_rng(11)
    urteile = []

    # 1 Grenzfall nu -> gross: der Stempel muss in reine Diffusion uebergehen
    from ns_loeser import zaeh_u
    u = rng.standard_normal(g.leer_u().shape) * 1e-6
    v = np.zeros(g.leer_v().shape)
    a = konv_zaeh_u(g, u, v, nu=1.0)
    b = zaeh_u(g, u, nu=1.0)
    # Raender weglassen: dort unterscheiden sich die Behandlungen
    d = np.abs(a[2:-2, 1:-1] - b[2:-2, 1:-1]).max() / max(np.abs(b[2:-2, 1:-1]).max(), 1e-30)
    urteile.append(("verschwindende Geschwindigkeit -> reine Diffusion", d < 1e-6, f"{d:.3e}"))

    # 2 Gleichfoermige Stroemung: konstantes u darf sich nicht aendern
    u = np.ones(g.leer_u().shape); v = np.zeros(g.leer_v().shape)
    r = konv_zaeh_u(g, u, v, nu=0.01, einstroem=np.ones(g.ny))
    d = np.abs(r[2:-2, 1:-1]).max()
    urteile.append(("gleichfoermige Stroemung bleibt unveraendert", d < 1e-10, f"{d:.3e}"))

    # 3 Die Unsymmetrie IST der Massenstrom: g_hin - g_rueck = ms
    # A ABSICHTLICH ungleich 1: mit A = 1 bestand diese Probe auch mit dem Fehler
    # aus Fassung 1.0 (Befund 31.10.2026).
    A_f, d_f = 0.0195, 0.02
    for uf in (1e-8, 1e-3, 1.0, 1e3, 1e6):
        gh, gr = sg_leitwerte(np.array([uf]), 1e-3, d_f, A_f)
        soll = uf * A_f
        urteile.append((f"g_hin - g_rueck = Volumenstrom bei u_f={uf:g}, A={A_f}",
                        abs(float(gh[0] - gr[0]) - soll) < 1e-9 * max(soll, 1.0),
                        f"{float(gh[0]-gr[0]):.6g} gegen {soll:.6g}"))

    for text, ok, beleg in urteile:
        print(f"  [{'RICHTIG' if ok else 'FALSCH  '}] {text}")
        print(f"              {beleg}")
    n = sum(1 for _, ok, _ in urteile if ok)
    print(f"\n  == {n}/{len(urteile)}")
    raise SystemExit(0 if n == len(urteile) else 1)
