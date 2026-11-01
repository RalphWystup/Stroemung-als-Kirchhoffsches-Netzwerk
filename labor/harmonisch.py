#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Das Zeitgesetz als harmonische Anpassung -- und der Vergleich mit DMD.

Befund 31.10.2026. Die DMD-Anpassung auf eingeschwungenen Daten liefert die richtige
Struktur (Grundfrequenz 2,524 Hz mit ganzzahligen Oberwellen, Eigenwerte auf dem
Einheitskreis), aber die Frequenz bleibt rund 0,17 % daneben. Ueber 6,3 Perioden laeuft
daraus ein Phasenfehler von 3,84 Grad auf, und weil das erste Modenpaar 96 % der
Schwankung traegt, sind das 7 % Feldfehler. Der Fehler WAECHST nicht, er schwingt --
das ist die Unterschrift eines Phasenfehlers, nicht eines Modellfehlers.

Dagegen hilft nicht ein laengerer Lauf, sondern die bekannte Struktur. Ein Grenzzyklus
ist periodisch; seine Beiwerte sind

    a(t) = a_0 + Summe_k [ c_k cos(2 pi k f t) + s_k sin(2 pi k f t) ]

mit EINER Grundfrequenz f und ganzzahligen Oberwellen. Das ist keine Annahme, die man
den Daten ueberstuelpt: Die paarweisen Singulaerwerte und die ganzzahligen
DMD-Frequenzen zeigen es. Angepasst werden c_k, s_k durch Ausgleichsrechnung ueber das
GANZE Fenster, und f selbst durch eine eindimensionale Suche -- damit kann sich keine
Phase davonstehlen, denn das Modell ist exakt periodisch.

Beide Modelle werden gemessen und beide Zahlen ausgewiesen. Welches besser ist, soll
die Messung sagen, nicht die Vorliebe.
"""
FASSUNG = "1.0"
import json
import os
import sys

import numpy as np

HIER = os.path.dirname(os.path.abspath(__file__))


def entwurf(t, f, K):
    """Entwurfsmatrix [1, cos(2pi f t), sin(2pi f t), cos(4pi f t), ...]."""
    S = [np.ones_like(t)]
    for k in range(1, K + 1):
        S.append(np.cos(2 * np.pi * k * f * t))
        S.append(np.sin(2 * np.pi * k * f * t))
    return np.array(S).T


def anpassen(a, t, f, K):
    """Kleinste Quadrate fuer alle Beiwerte auf einmal. a hat die Form (r, Zeiten)."""
    S = entwurf(t, f, K)
    B, *_ = np.linalg.lstsq(S, a.T, rcond=None)
    return B, S @ B


def suche_f(a, t, f0, K, spanne=0.02, stufen=241):
    """Eindimensionale Suche nach der Grundfrequenz: welches f erklaert die Daten am besten?"""
    bestes, bf = np.inf, f0
    for f in np.linspace(f0 * (1 - spanne), f0 * (1 + spanne), stufen):
        _, vor = anpassen(a, t, f, K)
        r = float(np.sum((vor.T - a) ** 2))
        if r < bestes:
            bestes, bf = r, f
    # verfeinern
    for schritt in (1e-4, 1e-5, 1e-6):
        for f in np.linspace(bf * (1 - 10 * schritt), bf * (1 + 10 * schritt), 41):
            _, vor = anpassen(a, t, f, K)
            r = float(np.sum((vor.T - a) ** 2))
            if r < bestes:
                bestes, bf = r, f
    return bf


if __name__ == "__main__":
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    d = np.load(os.path.join(HIER, "surrogat_daten.npz"))
    XA, XB, mittel, U, aA, aB = d["XA"], d["XB"], d["mittel"], d["U"], d["aA"], d["aB"]
    zeiten = d["zeiten"]; t = zeiten - zeiten[0]
    n = XA.shape[1]; spanne = np.abs(XA - mittel).max()
    print("=" * 94)
    print(f" ZEITGESETZ ALS HARMONISCHE ANPASSUNG -- Fassung {FASSUNG}, {K} Oberwellen")
    print("=" * 94)

    # Startwert fuer f aus dem fuehrenden Beiwert
    kr = np.where((aA[0, :-1] < 0) & (aA[0, 1:] >= 0))[0]
    dts = float(t[1] - t[0])
    f0 = 1.0 / (float(np.median(np.diff(kr))) * dts)
    f = suche_f(aA, t, f0, K)
    print(f"  Grundfrequenz: Startwert {f0:.5f} Hz -> angepasst {f:.6f} Hz "
          f"(Periode {1/f:.6f} s)")

    erg = {}
    for name, a, X in (("A", aA, XA), ("B", aB, XB)):
        B, vor = anpassen(a, t, f, K)
        F = U @ vor.T + mittel
        fehl = np.abs(F - X).max() / spanne
        erg[name] = {"B": B, "fehler": float(fehl)}
        print(f"  Weg {name}: Fehler ueber alle {n} Aufnahmen "
              f"{100*fehl:.3f} % der Schwankungsbreite")
    FA = U @ anpassen(aA, t, f, K)[1].T + mittel
    FB = U @ anpassen(aB, t, f, K)[1].T + mittel
    d_sur = np.abs(FA - FB).max() / spanne
    d_echt = np.abs(XA - XB).max() / spanne
    print(f"  Unterschied der beiden harmonischen Surrogate: {d_sur:.3e}  "
          f"(echte Rechnung: {d_echt:.3e})")

    # Vergleich mit DMD aus demselben Datensatz
    a = aA[:, 0].copy(); vor = []
    for k in range(n):
        vor.append(a.copy()); a = d["OpA"] @ a
    f_dmd = np.abs(U @ np.array(vor).T + mittel - XA).max() / spanne
    print(f"\n  {'Modell':34s} {'Fehler':>10s}")
    print(f"  {'DMD (linearer Zeitoperator)':34s} {100*f_dmd:9.3f} %")
    print(f"  {'harmonische Anpassung':34s} {100*erg['A']['fehler']:9.3f} %")
    besser = erg["A"]["fehler"] < f_dmd
    print(f"\n  == die harmonische Anpassung ist {'BESSER' if besser else 'nicht besser'} "
          f"(Faktor {f_dmd/max(erg['A']['fehler'],1e-30):.1f})")

    np.savez_compressed(os.path.join(HIER, "surrogat_harmonisch.npz"),
                        f=np.array([f]), K=np.array([K]),
                        BA=erg["A"]["B"], BB=erg["B"]["B"],
                        mittel=mittel, U=U, zeiten=zeiten)
    json.dump({"f": float(f), "K": K, "fehler_harmonisch": erg["A"]["fehler"],
               "fehler_dmd": float(f_dmd), "unterschied_surrogat": float(d_sur),
               "unterschied_echt": float(d_echt)},
              open(os.path.join(HIER, "harmonisch_bericht.json"), "w"), indent=1)
