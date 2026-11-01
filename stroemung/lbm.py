#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Der UNABHAENGIGE zweite Loeser: Lattice-Boltzmann, D2Q9 mit BGK-Stoss.

Warum dieser und nicht noch ein Aufbau desselben Verfahrens. Zwei Wege sind nur dann
unabhaengig, wenn sie VERSCHIEDEN scheitern koennen. wirbelstrasse.py (Weg A und Weg B)
teilt sich Zeitschritt, Randbehandlung, Koerperbehandlung, Kraftbilanz und den
Scharfetter-Gummel-Zweigstempel -- waere einer davon begrifflich falsch, waeren beide
einträchtig falsch. Dieser Loeser teilt mit ihm NICHTS:

    Zustand        Verteilungsfunktionen f_i statt u, v, p
    Zeitschritt    Stoss und Stroemen statt Vorhersage und Projektion
    Druck          gibt es nicht -- p = c_s^2 rho folgt aus der Dichte
    Inkompressibel nur im Grenzfall kleiner Mach-Zahl, nicht durch eine Zwangsbedingung
    Wand           Rueckprall an den Verbindungen statt Spiegelwerte
    Kraft          Impulsaustausch an den Verbindungen statt Impulsbilanz der Zweige

Er ist damit die Gegenprobe, die der Netzwerkloeser braucht: Trifft er dieselben
Beiwerte, so liegt es nicht an einer gemeinsamen Annahme.

Geprueft wird er selbst in drei Stufen, jede gegen eine geschlossene Loesung:
    L01  Zaehigkeit: Taylor-Green-Wirbel klingt mit exp(-2 nu k^2 t) ab
    L02  Hagen-Poiseuille: parabolisches Profil im Kanal
    L03  Wirbelstrasse gegen Schaefer-Turek

Alle Benchmarkgroessen sind dimensionslos und werden in GITTEREINHEITEN gebildet --
damit entfaellt jede Umrechnung und mit ihr eine Fehlerquelle.
"""
FASSUNG = "1.0"
import numpy as np

# --- D2Q9 ---------------------------------------------------------------------------
E = np.array([[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1],
              [1, 1], [-1, 1], [-1, -1], [1, -1]], dtype=int)
W = np.array([4/9] + [1/9]*4 + [1/36]*4)
GEGEN = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])      # entgegengesetzte Richtung
CS2 = 1.0 / 3.0                                     # Schallgeschwindigkeit im Quadrat


def gleichgewicht(rho, ux, uy):
    """f_i^eq = w_i rho [1 + 3(e.u) + 4,5(e.u)^2 - 1,5 u^2]."""
    uu = ux * ux + uy * uy
    f = np.empty((9,) + rho.shape)
    for i in range(9):
        eu = E[i, 0] * ux + E[i, 1] * uy
        f[i] = W[i] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uu)
    return f


def groessen(f):
    """Dichte und Geschwindigkeit aus den Verteilungen -- die Momente."""
    rho = f.sum(axis=0)
    ux = (E[:, 0][:, None, None] * f).sum(axis=0) / rho
    uy = (E[:, 1][:, None, None] * f).sum(axis=0) / rho
    return rho, ux, uy


def kraft_guo(tau, ux, uy, Fx, Fy):
    """Volumenkraft nach Guo: der Zusatz zu f_i, der die Kraft exakt einbringt.

    Ohne diesen Zusatz traegt eine Volumenkraft einen Fehler zweiter Ordnung in die
    Impulsbilanz; fuer die Hagen-Poiseuille-Probe waere das genau der Teil, den die
    Probe nachweisen soll.
    """
    q = np.empty((9,) + ux.shape)
    for i in range(9):
        eu = E[i, 0] * ux + E[i, 1] * uy
        q[i] = W[i] * (1.0 - 1.0 / (2.0 * tau)) * (
            3.0 * ((E[i, 0] - ux) * Fx + (E[i, 1] - uy) * Fy)
            + 9.0 * eu * (E[i, 0] * Fx + E[i, 1] * Fy))
    return q


def stoss_stroemen(f, tau, fest=None, kraft=None):
    """Ein Zeitschritt: BGK-Stoss, dann Stroemen, dann Rueckprall an festen Knoten.

    Rueckprall auf halbem Weg: Die Wand liegt zwischen Fluid- und Festknoten. Kommt eine
    Verteilung von einem Festknoten, so wird stattdessen die entgegengesetzte des eigenen
    Knotens von vor dem Stroemen eingesetzt.
    """
    rho = f.sum(axis=0)
    ux = (E[:, 0][:, None, None] * f).sum(axis=0) / rho
    uy = (E[:, 1][:, None, None] * f).sum(axis=0) / rho
    if kraft is not None:                   # halbe Kraft in die Geschwindigkeit (Guo)
        ux = ux + 0.5 * kraft[0] / rho
        uy = uy + 0.5 * kraft[1] / rho
    fs = f - (f - gleichgewicht(rho, ux, uy)) / tau            # Stoss
    if kraft is not None:
        fs = fs + kraft_guo(tau, ux, uy, kraft[0], kraft[1])
    fn = np.empty_like(fs)
    for i in range(9):
        fn[i] = np.roll(np.roll(fs[i], E[i, 0], axis=0), E[i, 1], axis=1)
    if fest is not None:
        for i in range(9):
            quelle_fest = np.roll(np.roll(fest, E[i, 0], axis=0), E[i, 1], axis=1)
            fn[i] = np.where(quelle_fest, fs[GEGEN[i]], fn[i])
    return fn, fs, rho, ux, uy


def kraft_impulsaustausch(fs, fest):
    """Kraft auf den Koerper: 2 e_i f_i* ueber alle Verbindungen, die ihn treffen.

    Das ist der Impuls, den die zurueckgeworfene Verteilung dem Koerper uebergibt --
    ein voellig anderer Weg als die Impulsbilanz der abgedeckten Zweige im
    Netzwerkloeser, und deshalb eine echte Gegenprobe auch der Kraftermittlung.
    """
    Fx = Fy = 0.0
    fluid = ~fest
    for i in range(1, 9):
        ziel_fest = np.roll(np.roll(fest, -E[i, 0], axis=0), -E[i, 1], axis=1)
        m = fluid & ziel_fest
        if not m.any():
            continue
        s = 2.0 * fs[i][m].sum()
        Fx += E[i, 0] * s
        Fy += E[i, 1] * s
    return Fx, Fy


# =====================================================================================
#  Proben gegen geschlossene Loesungen
# =====================================================================================

def probe_taylor_green(n=64, tau=0.8, schritte=400, U=0.02):
    """L01  Taylor-Green-Wirbel im periodischen Kasten klingt mit exp(-2 nu k^2 t) ab.

    Das prueft den Stoss und die Zaehigkeit nu = c_s^2 (tau - 1/2) -- eine Beziehung,
    die aus der Chapman-Enskog-Entwicklung folgt und nicht eingestellt werden kann.
    """
    nu = CS2 * (tau - 0.5)
    k = 2.0 * np.pi / n
    x = np.arange(n)[:, None] + 0.5
    y = np.arange(n)[None, :] + 0.5
    ux0 = -U * np.cos(k * x) * np.sin(k * y)
    uy0 = U * np.sin(k * x) * np.cos(k * y)
    f = gleichgewicht(np.ones((n, n)), ux0, uy0)
    for _ in range(schritte):
        f, _, _, _, _ = stoss_stroemen(f, tau)
    _, ux, uy = groessen(f)
    gemessen = np.sqrt((ux**2 + uy**2).max())
    erwartet = U * np.exp(-2.0 * nu * k * k * schritte)
    return gemessen, erwartet, nu


def probe_poiseuille(H=40, tau=0.9, schritte=20000, g=1e-6):
    """L02  Ebener Kanal unter Volumenkraft: u(y) = g/(2 nu) y (H - y).

    Das prueft zugleich, WO die Wand liegt: Beim Rueckprall auf halbem Weg sitzt sie
    zwischen Fest- und Fluidknoten. Eine um eine halbe Zelle verschobene Wand faellt
    hier sofort auf, weil das Profil dann nicht mehr passt.
    """
    nx, ny = 8, H + 2
    fest = np.zeros((nx, ny), bool)
    fest[:, 0] = True; fest[:, -1] = True
    nu = CS2 * (tau - 0.5)
    f = gleichgewicht(np.ones((nx, ny)), np.zeros((nx, ny)), np.zeros((nx, ny)))
    Fx = np.full((nx, ny), g); Fy = np.zeros((nx, ny))
    for _ in range(schritte):
        f, _, _, _, _ = stoss_stroemen(f, tau, fest, kraft=(Fx, Fy))
    rho = f.sum(axis=0)
    ux = ((E[:, 0][:, None, None] * f).sum(axis=0) + 0.5 * Fx) / rho
    y = np.arange(ny) - 0.5                      # Wand auf halbem Weg bei y = 0 und y = H
    genau = g / (2.0 * nu) * y * (H - y)
    innen = slice(1, ny - 1)
    gem = ux[nx // 2, innen]
    return gem, genau[innen], nu


if __name__ == "__main__":
    print("=" * 92)
    print(f" LATTICE-BOLTZMANN -- Fassung {FASSUNG}, Proben gegen geschlossene Loesungen")
    print("=" * 92)
    urteile = []
    g_, e_, nu_ = probe_taylor_green()
    rel = abs(g_ - e_) / e_
    urteile.append((f"L01 Taylor-Green klingt mit exp(-2 nu k^2 t) ab (nu = {nu_:.6f})",
                    rel < 0.02, f"gemessen {g_:.6e}, erwartet {e_:.6e}, {100*rel:.3f} % daneben"))
    gem, gen, nu_ = probe_poiseuille()
    rel = np.abs(gem - gen).max() / gen.max()
    urteile.append((f"L02 Hagen-Poiseuille, Profil und Wandlage (nu = {nu_:.6f})",
                    rel < 1e-3, f"groesste Abweichung {100*rel:.4f} % der Spitze"))
    for text, ok, beleg in urteile:
        print(f"  [{'RICHTIG' if ok else 'FALSCH  '}] {text}\n              {beleg}")
    n = sum(1 for _, ok, _ in urteile if ok)
    print(f"\n  == {n}/{len(urteile)}")
    raise SystemExit(0 if n == len(urteile) else 1)
