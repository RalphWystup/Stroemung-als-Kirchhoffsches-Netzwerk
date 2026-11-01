#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Vergleichsdaten fuer den Galerkin-Kern im Browser.

Der Browserkern kern_fem.js muss dasselbe rechnen wie fem_galerkin.py und
fem_wirbelstrasse.py -- nicht ungefaehr, sondern auf Maschinengenauigkeit. Dieses
Programm erzeugt die Zahlen, gegen die geprueft wird: vergleichsdaten_fem.json.

Dasselbe Vorgehen wie bei vergleichsdaten.py fuer den Netzwerkkern. Der Unterschied:
Hier wird AUCH DAS NETZ mitgegeben. Der Vernetzer im Browser streut seine Innenpunkte
mit einem eigenen Zufallsgenerator und kann deshalb nicht dasselbe Netz erzeugen wie
numpy -- das ist richtig so und keine Schwaeche: Das Netz ist Eingabe, nicht Verfahren.
Geprueft wird der RECHENKERN, und der bekommt dafuer genau das Netz, das hier entstand.
Der Vernetzer wird getrennt geprueft, an seinen eigenen Eigenschaften.

Aufruf:  python3 vergleichsdaten_fem.py
"""
FASSUNG = "1.0"
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "stroemung"))
import netz_fem
import fem_galerkin as G
import fem_wirbelstrasse as FW
import wirbelstrasse as W
import scipy.sparse as sp
import scipy.sparse.linalg as spl

HIER = os.path.dirname(os.path.abspath(__file__))


def liste(a):
    return np.asarray(a, float).ravel().tolist()


def probe(a, n=12):
    """Ein paar Werte und die Kennzahlen -- die ganze Liste waere unnoetig gross."""
    a = np.asarray(a, float).ravel()
    schritt = max(1, len(a) // n)
    return {"n": int(len(a)), "summe": float(a.sum()), "quadratsumme": float((a * a).sum()),
            "min": float(a.min()), "max": float(a.max()),
            "stellen": [int(i) for i in range(0, len(a), schritt)][:n],
            "werte": [float(a[i]) for i in range(0, len(a), schritt)][:n]}


def main():
    d = {"fassung": FASSUNG}

    # ---------------------------------------------------------------- 1. Das Netz
    P, T, A, zyl = netz_fem.netz(0.09, 0.03)
    d["netz"] = {"P": liste(P), "T": [int(i) for i in np.asarray(T).ravel()],
                 "nP": int(len(P)), "nT": int(len(T)), "nzyl": int(len(zyl)),
                 "flaeche": float(A.sum())}

    # ---------------------------------------------------------------- 2. Der P2-Raum
    XY, E6, nv, ne = G.p2_raum(P, T)
    d["p2raum"] = {"n2": int(len(XY)), "nv": int(nv), "nkanten": int(ne),
                   "E6_erste": [int(i) for i in E6[:5].ravel()],
                   "XY_probe": probe(XY)}

    # ---------------------------------------------------------------- 3. Formfunktionen
    punkte = [[1/3, 1/3, 1/3], [0.6, 0.3, 0.1], [1.0, 0.0, 0.0], [0.0, 0.5, 0.5]]
    d["formfunktionen"] = {"punkte": punkte,
                           "werte": [liste(G.p2_werte(np.array(L))) for L in punkte]}
    g, Ae = G.lambda_gradienten(P, T)
    d["lambda"] = {"g_erste": liste(g[:3]), "A_probe": probe(Ae)}
    d["ableitungen"] = {"element": 0,
                        "werte": [liste(G.p2_ableitungen(np.array(L), g[0])) for L in punkte]}

    # ---------------------------------------------------------------- 4. Die Matrizen
    M, K, Bx, By, g2, A2 = G.zusammenstellen(P, T, XY, E6, nv)
    for name, S in (("M", M), ("K", K), ("Bx", Bx), ("By", By)):
        S = S.tocsr()
        d.setdefault("matrizen", {})[name] = {
            "zeilen": int(S.shape[0]), "spalten": int(S.shape[1]),
            "besetzt": int(S.nnz), "summe": float(S.sum()),
            "quadratsumme": float(S.multiply(S).sum()),
            "zeilensumme_probe": probe(np.asarray(S.sum(axis=1)).ravel())}
    # Gegenprobe von Hand: Die Massenmatrix muss in der Summe die Flaeche geben,
    # die Steifigkeit muss konstante Felder vernichten.
    eins = np.ones(len(XY))
    d["matrizen"]["probe_masse_gleich_flaeche"] = float(M.sum())
    d["matrizen"]["probe_steifigkeit_mal_eins"] = float(np.abs(K @ eins).max())
    d["matrizen"]["probe_B_mal_eins"] = [float(np.abs(Bx @ eins).max()),
                                         float(np.abs(By @ eins).max())]

    # ---------------------------------------------------------------- 5. Die Raender
    R = G.raender(XY, nv, mit_zylinder=True)
    d["raender"] = {k: int(np.count_nonzero(R[k])) for k in ("ein", "wand", "zyl", "aus", "fest")}
    d["raender"]["u_probe"] = probe(R["u"])

    # ---------------------------------------------------------------- 6. Die Konvektion
    # ein festes, reproduzierbares Feld -- keine Zufallszahlen in einer Pruefdatei
    x, y = XY[:, 0], XY[:, 1]
    ux = np.sin(3 * x) * np.cos(5 * y)
    uy = np.cos(2 * x) * np.sin(4 * y)
    fx, fy = G.konvektion(P, T, E6, g2, A2, ux, uy)
    d["konvektion"] = {"feld": "ux=sin(3x)cos(5y), uy=cos(2x)sin(4y)",
                       "fx": probe(fx), "fy": probe(fy)}

    # ---------------------------------------------------------------- 7. Dirichlet
    dt = 2e-3
    Sys = G.sattelsystem(K, Bx, By, W.NU, M=M, dt=dt)
    N = 2 * len(XY) + nv
    feste = np.zeros(N, bool); werte = np.zeros(N)
    feste[:len(XY)] = R["fest"]; werte[:len(XY)] = R["u"]
    feste[len(XY):2*len(XY)] = R["fest"]; werte[len(XY):2*len(XY)] = R["v"]
    # Im Druckteil NULL -- so sieht die rechte Seite im echten Zeitschritt aus, und nur
    # dann laesst sie sich mit dem Browserkern vergleichen: Der schreibt die dritte
    # Blockzeile mit -B statt +B (symmetrisch, fuer MINRES). Beide Systeme haben
    # dieselbe Loesung, aber ihre rechten Seiten sind im Druckteil Gegenzahlen. Mit
    # einem von null verschiedenen Eintrag waeren sie es NICHT, und der Vergleich
    # haette eine Abweichung gemeldet, die keine ist.
    rhs0 = np.zeros(N)
    rhs0[:2*len(XY)] = np.arange(2*len(XY), dtype=float) * 1e-3
    Sfest, rhsfest = G.festsetzen(Sys, rhs0, feste, werte)
    d["dirichlet"] = {"dt": dt, "N": int(N), "feste": int(feste.sum()),
                      "S_besetzt": int(Sfest.nnz), "S_summe": float(Sfest.sum()),
                      "S_quadratsumme": float(Sfest.multiply(Sfest).sum()),
                      "rhs_probe": probe(rhsfest)}

    # ---------------------------------------------------------------- 8. Zeitschritte
    S_ = dict(P=P, T=T, XY=XY, E6=E6, nv=nv, ne=ne, M=M, K=K, Bx=Bx, By=By,
              g=g2, A=A2, R=R, n2=len(XY), n1=nv, nzyl=len(zyl))
    ux, uy, pr, verlauf, _ = FW.lauf(S_, t_ende=40 * dt, dt=dt)
    d["schritte"] = {"dt": dt, "anzahl": int(len(verlauf)),
                     "cD": [float(v[1]) for v in verlauf],
                     "cL": [float(v[2]) for v in verlauf],
                     "ux_probe": probe(ux), "uy_probe": probe(uy), "p_probe": probe(pr)}

    # ---------------------------------------------------------------- 9. Hagen-Poiseuille
    # Der nackte Kanal: Die Loesung ist geschlossen bekannt, und ein P2-Ansatz gibt ein
    # Parabelprofil EXAKT wieder. Weicht hier etwas ab, ist der Kern falsch -- nicht
    # ungenau. Das ist die schaerfste Probe, die es fuer diesen Weg gibt.
    Pk, Tk, Ak, _ = netz_fem.netz(0.09, 0.03, mit_zylinder=False)
    XYk, E6k, nvk, nek = G.p2_raum(Pk, Tk)
    Mk, Kk, Bxk, Byk, gk, Aek = G.zusammenstellen(Pk, Tk, XYk, E6k, nvk)
    Rk = G.raender(XYk, nvk, mit_zylinder=False)
    Nk = 2 * len(XYk) + nvk
    Sk = G.sattelsystem(Kk, Bxk, Byk, W.NU)
    fk = np.zeros(Nk, bool); wk = np.zeros(Nk)
    fk[:len(XYk)] = Rk["fest"]; wk[:len(XYk)] = Rk["u"]
    fk[len(XYk):2*len(XYk)] = Rk["fest"]; wk[len(XYk):2*len(XYk)] = Rk["v"]
    Skf, rhsk = G.festsetzen(Sk, np.zeros(Nk), fk, wk)
    loes = spl.spsolve(Skf.tocsc(), rhsk)
    uxk = loes[:len(XYk)]
    soll = W.einstroemprofil(XYk[:, 1])
    fehler = float(np.abs(uxk - soll).max() / soll.max())
    d["poiseuille"] = {"nP": int(len(Pk)), "N": int(Nk), "fehler": fehler,
                       "ux_probe": probe(uxk),
                       # das nackte Netz muss mit, sonst kann der Browser nicht dieselbe
                       # Rechnung machen -- das Netz ist Eingabe, nicht Verfahren
                       "P": liste(Pk), "T": [int(i) for i in np.asarray(Tk).ravel()],
                       "nT": int(len(Tk))}

    # ---------------------------------------------------------------- 10. Die Kennzahlen
    d["fall"] = {"LX": W.LX, "LY": W.LY, "XC": W.XC, "YC": W.YC, "DZ": W.DZ,
                 "UM": W.UM, "NU": W.NU, "RHO": W.RHO, "U_QUER": W.U_QUER,
                 "REFERENZ": {k: list(v) for k, v in W.REFERENZ.items()}}

    # --------------------------------------------------------------- 11. Die Anzeige
    # Die Wirbelstaerke auf das Anzeigegitter zu bringen ist eine ANZEIGE-Operation,
    # keine Rechenoperation -- aber sie entscheidet, was der Betrachter sieht, und
    # Befund B16 hat gezeigt, wie teuer ein Fehler genau dort wird. Also wird auch sie
    # geprueft. Ausgewertet wird auf den GITTERKNOTEN (i*dx, j*dy), dort, wo der
    # Netzwerkweg seine Wirbelstaerke bildet.
    nxA, nyA = 110, 20
    om = FW.wirbel_auf_gitter(S_, ux, uy, nxA, nyA)
    # dieselben Punkte wie der Netzwerkweg: (i+1)*LX/(nxA+1), (j+1)*LY/(nyA+1)
    import matplotlib.tri as mtri
    P2, T2, E62, g2b, A2b = S_["P"], S_["T"], S_["E6"], S_["g"], S_["A"]
    Lc = np.array([1/3, 1/3, 1/3])
    dN = np.empty((len(T2), 6, 2))
    dN[:, 0] = (4*Lc[0]-1) * g2b[:, 0]; dN[:, 1] = (4*Lc[1]-1) * g2b[:, 1]
    dN[:, 2] = (4*Lc[2]-1) * g2b[:, 2]
    dN[:, 3] = 4 * (Lc[2]*g2b[:, 1] + Lc[1]*g2b[:, 2])
    dN[:, 4] = 4 * (Lc[0]*g2b[:, 2] + Lc[2]*g2b[:, 0])
    dN[:, 5] = 4 * (Lc[1]*g2b[:, 0] + Lc[0]*g2b[:, 1])
    om_e = (np.einsum("ei,ei->e", uy[E62], dN[:, :, 0])
            - np.einsum("ei,ei->e", ux[E62], dN[:, :, 1]))
    zahl = np.zeros(len(P2)); summe_ = np.zeros(len(P2))
    for k in range(3):
        np.add.at(summe_, T2[:, k], om_e * A2b)
        np.add.at(zahl, T2[:, k], A2b)
    om_v = summe_ / np.maximum(zahl, 1e-30)
    tri = mtri.Triangulation(P2[:, 0], P2[:, 1], T2)
    ip = mtri.LinearTriInterpolator(tri, om_v)
    xs = (np.arange(nxA) + 1.0) * W.LX / (nxA + 1)
    ys = (np.arange(nyA) + 1.0) * W.LY / (nyA + 1)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    Zk = ip(X, Y)
    Zk = np.where(np.isfinite(Zk), Zk, 0.0)
    d["anzeige"] = {"nx": nxA, "ny": nyA,
                    "mitte": probe(om.ravel()),
                    "knoten": probe(Zk.ravel())}

    ziel = os.path.join(HIER, "vergleichsdaten_fem.json")
    json.dump(d, open(ziel, "w"), separators=(",", ":"))
    print(f"  {os.path.basename(ziel)}: {os.path.getsize(ziel)/1024:.1f} kB")
    print(f"  Netz {len(P)} Ecken, {len(T)} Dreiecke, {len(XY)} P2-Knoten, {N} Unbekannte")
    print(f"  Masse = Flaeche:        {M.sum():.12f} gegen {A.sum():.12f}")
    print(f"  K mal eins:             {np.abs(K @ eins).max():.3e}")
    print(f"  Hagen-Poiseuille:       {fehler:.3e} relativ")
    print(f"  {len(verlauf)} Zeitschritte, cD zuletzt {verlauf[-1][1]:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
