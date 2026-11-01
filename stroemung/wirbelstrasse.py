#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Stufe 3: Kaermaensche Wirbelstrasse am Zylinder -- Schaefer-Turek-Benchmark 2D-2.

Aufbau nach Schaefer und Turek (DFG 1996), Fall 2D-2 (instationaer):

    Kanal          2,2 x 0,41
    Zylinder       Mittelpunkt (0,20 | 0,20), Durchmesser D = 0,1
                   leicht aussermittig (0,20 gegen H/2 = 0,205), damit die Abloesung anlaeuft
    Einstroemung   u(y) = 4 U_m y (H-y) / H^2 ,  U_m = 1,5  ->  u_quer = 1,0
    Stoffwerte     rho = 1 , nu = 1e-3        ->  Re = u_quer * D / nu = 100

Veroeffentlichte Referenzwerte (von mehreren Gruppen unabhaengig bestaetigt):

    c_D max   3,22 ... 3,24
    c_L max   0,99 ... 1,01
    St        0,295 ... 0,305

Der Zylinder wird als eingetauchter Koerper gesetzt (Geschwindigkeit null in den
abgedeckten Zweigen). Das ist auf kartesischem Gitter nur naeherungsweise rund; die
Strouhalzahl ist davon wenig betroffen, c_D und c_L dagegen empfindlich. Deshalb wird eine
NETZREIHE gefahren und der Trend berichtet -- nicht ein einzelner Lauf, der zufaellig im
Intervall liegt.
"""
FASSUNG = "1.4"
import json
import sys
import time

import numpy as np
import scipy.sparse.linalg as spl

from gitter import Gitter
from ns_loeser import drucknetz_matrix
from konvektion import konv_zaeh_u, konv_zaeh_v

# --- Benchmark-Festlegungen (Schaefer-Turek 2D-2) ------------------------------------
LX, LY = 2.2, 0.41
XC, YC, DZ = 0.2, 0.2, 0.1
UM = 1.5
NU = 1e-3
RHO = 1.0
U_QUER = 2.0 / 3.0 * UM          # = 1,0
RE = U_QUER * DZ / NU            # = 100

REFERENZ = {"cD": (3.22, 3.24), "cL": (0.99, 1.01), "St": (0.295, 0.305)}


def einstroemprofil(y):
    return 4.0 * UM * y * (LY - y) / LY**2


def koerper_masken(g):
    """Welche Zweige und Zellen liegen im Zylinder?"""
    XU, YU = np.meshgrid(g.xu, g.yu, indexing="ij")
    XV, YV = np.meshgrid(g.xv, g.yv, indexing="ij")
    XP, YP = np.meshgrid(g.xp, g.yp, indexing="ij")
    r = DZ / 2.0
    mu = (XU - XC) ** 2 + (YU - YC) ** 2 < r * r
    mv = (XV - XC) ** 2 + (YV - YC) ** 2 < r * r
    mp = (XP - XC) ** 2 + (YP - YC) ** 2 < r * r
    return mu, mv, mp


def kraefte(g, u, v, p, mu_m, mv_m, mp_m):
    """Kraft auf den Koerper aus dem Impulsverlust -- die saubere Art beim eingetauchten Koerper.

    Statt Spannungen auf einer Treppenkontur zu integrieren (dort ist die Normale
    unbestimmt), wird die Kraft als der Impuls gemessen, den das Verfahren je Zeitschritt
    aus den abgedeckten Zweigen nimmt. Das ist exakt die Reaktionskraft des Koerpers.
    Dieser Weg ist unter dem Namen „momentum forcing" ueblich.
    """
    return None   # wird im Lauf direkt gebildet (siehe unten)


def lauf(nx, ny, t_ende=8.0, cfl=0.4, melde=None, speichern=None, beobachter=None):
    g = Gitter(LX, LY, nx, ny)
    mu_m, mv_m, mp_m = koerper_masken(g)
    u, v = g.leer_u(), g.leer_v()
    ein = einstroemprofil(g.yu)
    u[:, :] = ein[None, :]
    u[mu_m] = 0.0; v[mv_m] = 0.0
    A, leer = drucknetz_matrix(g, ausstrom_rechts=True,
                               u_maske=mu_m, v_maske=mv_m)
    loes = spl.factorized(A.tocsc())
    dt = cfl * min(g.dx, g.dy) / (2.0 * UM)
    dt = min(dt, 0.2 * min(g.dx, g.dy) ** 2 / NU)
    n = int(t_ende / dt)
    verlauf = []
    t0 = time.time()
    for k in range(n):
        t = k * dt
        # 1 Vorhersage
        us = u + dt * konv_zaeh_u(g, u, v, NU, einstroem=ein)
        vs = v + dt * konv_zaeh_v(g, u, v, NU)
        # Fassung 1.2, Befund 31.10.2026 -- zwei Fehler an derselben Stelle:
        #
        # 1. `us[:, 0] = us[:, -1] = 0` war falsch. u sitzt bei y = dy/2 und y = Ly - dy/2,
        #    also INNEN, nicht an der Wand. Die Haftbedingung gehoert in den Spiegelwert des
        #    Zaehigkeitsterms (_rand_u in konvektion.py macht das bereits). Hier hat sie die
        #    Geschwindigkeit in der ersten und letzten Zellreihe erschlagen.
        # 2. Die Zeile stand NACH `us[0,:] = ein` und setzte damit die Einstroemung in der
        #    Eckzelle wieder auf null. Nach der Projektion wurde sie am Ende erneut auf ein
        #    gesetzt -- der Zufluss aenderte sich also NACH der Druckkorrektur, und genau
        #    diese Aenderung blieb als Divergenz stehen: gemessen 6,97 in Zelle (0,0).
        #
        # v dagegen sitzt bei y = 0 und y = Ly genau auf der Wand; dort ist null richtig.
        vs[:, 0] = vs[:, -1] = 0.0                      # undurchlaessige Kanalwand
        vs[0, :] = 0.0                                  # Einstroemung ohne Querkomponente
        us[0, :] = ein                                  # Einstroemung
        us[-1, :] = us[-2, :]                           # freier Ausstrom
        vs[-1, :] = vs[-2, :]
        # 2 Koerper: Impuls herausnehmen und als Kraft buchen
        #
        # Fassung 1.3, Befund 31.10.2026: Hier wurde nur der Impuls gezaehlt, den
        # Konvektion und Zaehigkeit in den Koerper tragen. Der DRUCKANTEIL fehlte --
        # und bei einem umstroemten Zylinder bei Re = 100 ist gerade er der groessere.
        # Gemessen ergab das c_D = 0,60 statt 3,23, also Faktor 5,4 zu klein.
        #
        # Richtig ist: Der Koerper nimmt in JEDEM Teilschritt Impuls aus der Stroemung.
        # Gezaehlt wird deshalb beides -- was nach Konvektion/Zaehigkeit in den
        # abgedeckten Zweigen steht UND was die Druckkorrektur dort hineinschiebt.
        fx1 = float(np.sum(us[mu_m])); fy1 = float(np.sum(vs[mv_m]))
        us[mu_m] = 0.0; vs[mv_m] = 0.0
        # Fassung 1.1: Die Massenbilanz wird NICHT mehr von Hand erzwungen. Der Ausstrom
        # haengt jetzt am Bezugspotential (drucknetz_matrix(ausstrom_rechts=True)); die
        # Druckkorrektur stellt die Bilanz damit von selbst her. Das Zurechtskalieren der
        # letzten Spalte hatte die gerade hergestellte Knotenregel wieder zerstoert.
        # 3 Drucknetz
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
        b = -d.ravel() * (g.dx * g.dy) / dt
        b[leer] = 0.0                      # Zellen ganz im Koerper: Einheitszeile
        p = loes(b).reshape(g.nx, g.ny)
        # 4 Korrektur
        u = us.copy(); v = vs.copy()
        u[1:-1, :] -= dt * (p[1:, :] - p[:-1, :]) / g.dx
        v[:, 1:-1] -= dt * (p[:, 1:] - p[:, :-1]) / g.dy
        # Ausstrom: Gegenseite ist das Bezugspotential 0, Abstand dx/2
        u[-1, :] -= dt * (0.0 - p[-1, :]) / (g.dx / 2.0)
        fx2 = float(np.sum(u[mu_m])); fy2 = float(np.sum(v[mv_m]))
        u[mu_m] = 0.0; v[mv_m] = 0.0
        u[0, :] = ein
        # Beiwerte: der GESAMTE je Zeitschritt entnommene Impuls ist die Reaktionskraft
        fx = (fx1 + fx2) * g.dx * g.dy * RHO / dt
        fy = (fy1 + fy2) * g.dx * g.dy * RHO / dt
        nenner = 0.5 * RHO * U_QUER**2 * DZ
        verlauf.append((t, fx / nenner, fy / nenner))
        # Fassung 1.4: NUR-LESEN-Haken. Er rechnet nichts und aendert nichts; er reicht
        # die Felder am ENDE des vollstaendigen Zeitschritts heraus. Noetig geworden durch
        # Befund B3: zwischen „die Projektion rechnet richtig" und „am Ende des Schritts
        # ist die Knotenregel erfuellt" liegt alles, was danach noch an den Feldern gedreht
        # wird -- und genau dort sass der Fehler. Ohne Haken kein laufender Nachweis.
        if beobachter is not None:
            beobachter(k, t, g, u, v, p)
        if melde and k % melde == 0:
            print(f"    t={t:6.3f}  cD={fx/nenner:7.3f}  cL={fy/nenner:7.3f}  "
                  f"({time.time()-t0:5.0f} s)", flush=True)
    if speichern:
        np.savez_compressed(speichern, u=u, v=v, p=p, mp=mp_m,
                            xp=g.xp, yp=g.yp, xu=g.xu, yu=g.yu, xv=g.xv, yv=g.yv,
                            verlauf=np.array(verlauf))
    return g, u, v, p, np.array(verlauf), dt


def auswerten(verlauf, dt, ab=0.5):
    """c_D max, c_L max und Strouhalzahl aus dem eingeschwungenen Teil."""
    n = len(verlauf)
    i0 = int(n * ab)
    t, cD, cL = verlauf[i0:, 0], verlauf[i0:, 1], verlauf[i0:, 2]
    # Strouhal aus dem Auftriebsbeiwert: Nulldurchgaenge aufsteigend
    s = np.signbit(cL)
    kreuz = np.where(~s[1:] & s[:-1])[0]
    if len(kreuz) >= 3:
        perioden = np.diff(t[kreuz])
        T = float(np.median(perioden))
        St = DZ / (U_QUER * T)
    else:
        T, St = float("nan"), float("nan")
    return float(np.max(cD)), float(np.max(cL)), St, T


if __name__ == "__main__":
    nx = int(sys.argv[1]) if len(sys.argv) > 1 else 220
    ny = int(sys.argv[2]) if len(sys.argv) > 2 else 41
    t_ende = float(sys.argv[3]) if len(sys.argv) > 3 else 8.0
    print("=" * 88)
    print(f" STUFE 3: KAERMAENSCHE WIRBELSTRASSE -- Schaefer-Turek 2D-2, Re = {RE:.0f}")
    print(f" Gitter {nx} x {ny}, bis t = {t_ende} s")
    print("=" * 88)
    g, u, v, p, verlauf, dt = lauf(nx, ny, t_ende=t_ende, melde=max(1, int(0.5 / 1e-3)),
                                   speichern=f"feld_{nx}x{ny}.npz")
    cD, cL, St, T = auswerten(verlauf, dt)
    print(f"\n  Zeitschritt dt = {dt:.3e} s, {len(verlauf)} Schritte")
    print(f"\n{'Groesse':>10s} {'gerechnet':>12s} {'Referenz':>16s} {'Urteil':>10s}")
    for name, wert, (lo, hi) in (("c_D max", cD, REFERENZ["cD"]),
                                 ("c_L max", cL, REFERENZ["cL"]),
                                 ("St", St, REFERENZ["St"])):
        ok = lo <= wert <= hi
        print(f"{name:>10s} {wert:12.4f} {f'{lo} … {hi}':>16s} "
              f"{'im Intervall' if ok else 'daneben':>14s}")
    json.dump({"nx": nx, "ny": ny, "cD": cD, "cL": cL, "St": St, "T": T, "dt": dt},
              open(f"ergebnis_{nx}x{ny}.json", "w"), indent=1)
