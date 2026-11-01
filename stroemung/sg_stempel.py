#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Stufe 1: der Scharfetter-Gummel-Zweigstempel -- und der Nachweis, dass er exakt ist.

Sein Wort 31.10.2026: „ohne Fake, rein analytisch und nummerisch bewiesen".

Worum es geht
-------------
Konvektion-Diffusion auf einem Zweig zwischen zwei Knoten. Der Fluss ist

    J = u*c - D*dc/dx

Die naheliegenden Stempel versagen beide, sobald die Gitter-Peclet-Zahl Pe_h = u*h/D gross
wird -- und bei Tinte in Wasser ist sie rund 10^6:

    zentrale Differenz   schwingt, sobald Pe_h > 2
    Aufwind              stabil, aber mit numerischer Diffusion, die das physikalische D
                         um Groessenordnungen uebertreffen kann

Scharfetter-Gummel (1969, aus der Halbleitersimulation) loest die 1D-Gleichung auf dem
Zweig EXAKT und stempelt diese Loesung als Leitwert. Damit ist er auf gleichmaessigem
Gitter KNOTENWEISE EXAKT, bei jedem Pe -- und genau das wird hier nachgewiesen.

    J_{i+1/2} = (D/h) * [ B(-Pe_h)*c_i - B(Pe_h)*c_{i+1} ]      mit B(x) = x/(e^x - 1)

In der Netzwerksprache: B(-Pe_h)*D/h und B(Pe_h)*D/h sind die beiden (unterschiedlichen!)
Leitwerte des Zweigs in Hin- und Rueckrichtung. Ein passiver Widerstand hat nur einen --
die Unsymmetrie IST die gesteuerte Quelle, nach der er gefragt hat.

Nachweis
--------
Stationaer, 1D, c(0)=0, c(L)=1, konstantes u. Geschlossene Loesung:

    c(x) = (exp(Pe*x/L) - 1) / (exp(Pe) - 1)        Pe = u*L/D

Geprueft ueber sieben Groessenordnungen in Pe. SG muss auf Maschinengenauigkeit treffen;
die beiden anderen Stempel werden mitgerechnet, damit der Unterschied sichtbar ist und
nicht behauptet werden muss.
"""
FASSUNG = "1.0"
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl


def bernoulli(x):
    """B(x) = x/(e^x - 1), numerisch sauber ueber den ganzen Bereich.

    Die direkte Auswertung bricht an drei Stellen: bei x -> 0 (0/0), bei grossem positivem x
    (e^x laeuft ueber) und bei grossem negativem x (Ausloeschung). Die Fallunterscheidung
    ist nicht Kosmetik -- ohne sie ist der Stempel genau dort falsch, wo er gebraucht wird.
    """
    x = np.asarray(x, float)
    y = np.empty_like(x)
    klein = np.abs(x) < 1e-10
    gross = x > 500.0
    tief = x < -500.0
    rest = ~(klein | gross | tief)
    y[klein] = 1.0 - x[klein] / 2.0                  # Reihe: B(x) = 1 - x/2 + x^2/12 - ...
    y[gross] = 0.0                                   # x/e^x -> 0
    y[tief] = -x[tief]                               # x/(0-1) = -x
    y[rest] = x[rest] / np.expm1(x[rest])            # expm1 statt exp(x)-1: keine Ausloeschung
    return y


def loese(n, L, u, D, stempel):
    """Stationaere 1D-Konvektion-Diffusion auf n Zellen, c(0)=0, c(L)=1.

    Knoten i = 0..n liegen bei x_i = i*h. Bilanziert wird um die inneren Knoten:
    was durch den linken Zweig hereinkommt, muss durch den rechten hinaus -- Kirchhoffs
    Knotenregel. Der Stempel bestimmt nur, WIE der Zweigfluss aus den beiden Knotenwerten
    gebildet wird.
    """
    h = L / n
    N = n + 1
    Pe_h = u * h / D
    A = sp.lil_matrix((N, N))
    b = np.zeros(N)

    if stempel == "sg":
        g_hin = D / h * bernoulli(np.array([-Pe_h]))[0]    # Leitwert i -> i+1
        g_rueck = D / h * bernoulli(np.array([Pe_h]))[0]   # Leitwert i+1 -> i
    elif stempel == "zentral":
        g_hin = D / h + u / 2.0
        g_rueck = D / h - u / 2.0
    elif stempel == "aufwind":
        g_hin = D / h + max(u, 0.0)
        g_rueck = D / h + max(-u, 0.0)
    else:
        raise ValueError(stempel)

    for i in range(1, n):
        # Fluss durch den rechten Zweig minus Fluss durch den linken = 0
        A[i, i - 1] += -g_hin
        A[i, i] += g_hin + g_rueck
        A[i, i + 1] += -g_rueck
    A[0, 0] = 1.0; b[0] = 0.0
    A[n, n] = 1.0; b[n] = 1.0
    return np.asarray(spl.spsolve(A.tocsr(), b)).ravel()


def genau(x, L, Pe):
    """c(x) = (e^{Pe x/L} - 1)/(e^{Pe} - 1), ueberlaufsicher."""
    if abs(Pe) < 1e-12:
        return x / L
    if Pe > 500:                      # Grenzfall: alles sitzt am rechten Rand
        return np.exp(Pe * (x / L - 1.0))
    if Pe < -500:
        return 1.0 - np.exp(-Pe * (x / L))*0 - (1.0 - np.exp(Pe * x / L))
    return np.expm1(Pe * x / L) / np.expm1(Pe)


if __name__ == "__main__":
    L, D, n = 1.0, 1.0, 20
    print("=" * 86)
    print(" SCHARFETTER-GUMMEL-STEMPEL gegen die geschlossene Loesung")
    print(" stationaer, 1D, c(0)=0, c(L)=1, 20 Zellen -- groesster Knotenfehler")
    print("=" * 86)
    print(f"\n{'Pe':>12s} {'Pe_h':>10s} {'SG':>14s} {'zentral':>14s} {'Aufwind':>14s}")
    urteile = []
    for Pe in (0.0, 1.0, 2.0, 5.0, 20.0, 100.0, 1e3, 1e4, 1e6):
        u = Pe * D / L
        x = np.linspace(0.0, L, n + 1)
        c_ex = np.array([genau(xx, L, Pe) for xx in x])
        zeile = []
        for st in ("sg", "zentral", "aufwind"):
            c = loese(n, L, u, D, st)
            zeile.append(float(np.max(np.abs(c - c_ex))))
        print(f"{Pe:12.3g} {Pe*L/n/L:10.1f} {zeile[0]:14.3e} {zeile[1]:14.3e} {zeile[2]:14.3e}")
        urteile.append((Pe, zeile[0]))

    schlecht = [p for p, e in urteile if e > 1e-12]
    print()
    print(f"  [{'RICHTIG' if not schlecht else 'FALSCH  '}] SG ist bei JEDEM Pe knotenweise "
          f"exakt (Fehler < 1e-12)")
    if schlecht:
        print(f"           nicht erreicht bei Pe = {schlecht}")
    print(f"           groesster SG-Fehler ueber alle Faelle: "
          f"{max(e for _, e in urteile):.3e}")
    raise SystemExit(0 if not schlecht else 1)
