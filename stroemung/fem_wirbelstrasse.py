#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Die Wirbelstrasse auf dem UNABHAENGIGEN Weg: Galerkin, Taylor-Hood, koerperangepasst.

Derselbe Fall wie wirbelstrasse.py, aber aus einer anderen Herleitung:

    Netzwerkweg   Flussbilanz am Kontrollvolumen; Unbekannte sind Zellmittel und
                  Flaechenfluesse; der Zylinder ist eine Treppe
    dieser Weg    schwache Form; Unbekannte sind Knotenwerte eines Polynomansatzes
                  (P2 fuer u, P1 fuer p); der Zylinderrand LIEGT auf Knoten

Zeitschritt halbimplizit: Zaehigkeit und Druck implizit, Konvektion explizit. Damit ist
die Systemmatrix konstant und wird nur EINMAL zerlegt; je Zeitschritt bleibt ein
Einsetzen. Die Konvektion wird als Vektor ausgewertet, nicht als Matrix.

Die Kraft kommt aus der schwachen Form selbst: Testet man die Impulsgleichung mit einer
Funktion, die auf dem Zylinder gleich eins und auf allen anderen Raendern null ist, so
steht auf der einen Seite das Randintegral der Spannungen -- also die Kraft. Dieser Weg
heisst „konsistente" Kraftermittlung und ist deutlich genauer als das Aufsummieren von
Wandspannungen. Er hat mit der Impulsbilanz der abgedeckten Zweige im Netzwerkweg
nichts gemein.

Aufruf:  python3 fem_wirbelstrasse.py [h_zyl] [t_ende] [dt]
"""
FASSUNG = "1.0"
import json
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

import netz_fem
import fem_galerkin as G
import wirbelstrasse as W


def aufbau(h_zyl=0.0039, h_fern=0.05, **fein):
    P, T, A, zyl = netz_fem.netz(h_fern, h_zyl, **fein)
    XY, E6, nv, ne = G.p2_raum(P, T)
    M, K, Bx, By, g, Ae = G.zusammenstellen(P, T, XY, E6, nv)
    R = G.raender(XY, nv, mit_zylinder=True)
    return dict(P=P, T=T, XY=XY, E6=E6, nv=nv, ne=ne, M=M, K=K, Bx=Bx, By=By,
                g=g, A=Ae, R=R, n2=len(XY), n1=nv, nzyl=len(zyl))


def wirbel_auf_gitter(S_, ux, uy, nx, ny):
    """Wirbelstaerke des Galerkin-Feldes auf ein kartesisches Anzeigegitter bringen.

    omega = dv/dx - du/dy ist bei P2 je Element LINEAR und ueber die Kanten unstetig.
    Fuer die Anzeige wird es flaechengewichtet auf die Ecken gemittelt und dann linear
    interpoliert -- das ist eine Anzeigeoperation, keine Rechenoperation, und sie steht
    deshalb hier und nicht im Loeser.
    """
    import matplotlib.tri as mtri
    P, T, E6, g, A = S_["P"], S_["T"], S_["E6"], S_["g"], S_["A"]
    ne = len(T)
    # Gradienten im Schwerpunkt (baryzentrisch 1/3,1/3,1/3)
    L = np.array([1/3, 1/3, 1/3])
    dN = np.empty((ne, 6, 2))
    dN[:, 0] = (4*L[0]-1) * g[:, 0]; dN[:, 1] = (4*L[1]-1) * g[:, 1]
    dN[:, 2] = (4*L[2]-1) * g[:, 2]
    dN[:, 3] = 4 * (L[2]*g[:, 1] + L[1]*g[:, 2])
    dN[:, 4] = 4 * (L[0]*g[:, 2] + L[2]*g[:, 0])
    dN[:, 5] = 4 * (L[1]*g[:, 0] + L[0]*g[:, 1])
    om_e = (np.einsum("ei,ei->e", uy[E6], dN[:, :, 0])
            - np.einsum("ei,ei->e", ux[E6], dN[:, :, 1]))
    zahl = np.zeros(len(P)); summe = np.zeros(len(P))
    for k in range(3):
        np.add.at(summe, T[:, k], om_e * A)
        np.add.at(zahl, T[:, k], A)
    om_v = summe / np.maximum(zahl, 1e-30)
    tri = mtri.Triangulation(P[:, 0], P[:, 1], T)
    ip = mtri.LinearTriInterpolator(tri, om_v)
    xs = (np.arange(nx) + 0.5) * W.LX / nx
    ys = (np.arange(ny) + 0.5) * W.LY / ny
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    Z = ip(X, Y)
    return np.where(np.isfinite(Z), Z, 0.0)


def lauf(S_, t_ende=12.0, dt=5e-4, melde=None, speichern=None, schnapp=None):
    n2, n1 = S_["n2"], S_["n1"]
    M, K, Bx, By, R = S_["M"], S_["K"], S_["Bx"], S_["By"], S_["R"]
    N = 2 * n2 + n1
    Sys = G.sattelsystem(K, Bx, By, W.NU, M=M, dt=dt)
    feste = np.zeros(N, bool); werte = np.zeros(N)
    feste[:n2] = R["fest"]; werte[:n2] = R["u"]
    feste[n2:2*n2] = R["fest"]; werte[n2:2*n2] = R["v"]
    Sfest, _ = G.festsetzen(Sys, np.zeros(N), feste, werte)
    t0 = time.time()
    loeser = spl.factorized(Sfest.tocsc())
    t_zerl = time.time() - t0
    # Testfunktion fuer die Kraft: eins auf dem Zylinder, null sonst
    vz = np.zeros(n2); vz[R["zyl"]] = 1.0
    ux = np.where(R["fest"], R["u"], 0.0)
    uy = np.zeros(n2)
    nenner = 0.5 * W.RHO * W.U_QUER**2 * W.DZ
    verlauf = []
    t1 = time.time()
    n = int(round(t_ende / dt))
    for k in range(n):
        fx, fy = G.konvektion(S_["P"], S_["T"], S_["E6"], S_["g"], S_["A"], ux, uy)
        rhs = np.zeros(N)
        rhs[:n2] = M @ ux / dt - fx
        rhs[n2:2*n2] = M @ uy / dt - fy
        _, rhs = G.festsetzen(Sys, rhs, feste, werte)
        loes = loeser(rhs)
        uxn, uyn, pr = loes[:n2], loes[n2:2*n2], loes[2*n2:]
        # Konsistente Kraft: Residuum der Impulsgleichung, getestet mit vz
        Fx = -(M @ ((uxn - ux) / dt) + W.NU * (K @ uxn) + fx - Bx.T @ pr) @ vz
        Fy = -(M @ ((uyn - uy) / dt) + W.NU * (K @ uyn) + fy - By.T @ pr) @ vz
        ux, uy = uxn, uyn
        verlauf.append((k * dt, Fx / nenner, Fy / nenner))
        if schnapp is not None and k * dt >= schnapp["ab"] and k % schnapp["jeder"] == 0:
            schnapp["bilder"].append(
                wirbel_auf_gitter(S_, ux, uy, schnapp["nx"], schnapp["ny"]).ravel().astype("f4"))
            schnapp["zeiten"].append(k * dt)
            schnapp["cl"].append(Fy / nenner)
        if melde and k % melde == 0:
            print(f"    t={k*dt:6.3f}  cD={Fx/nenner:7.3f}  cL={Fy/nenner:7.3f}  "
                  f"({time.time()-t1:5.0f} s)", flush=True)
    verlauf = np.array(verlauf)
    if speichern:
        np.savez_compressed(speichern, XY=S_["XY"], T=S_["T"], E6=S_["E6"],
                            ux=ux, uy=uy, p=pr, verlauf=verlauf, dt=dt)
    return ux, uy, pr, verlauf, t_zerl


if __name__ == "__main__":
    h_zyl = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0039
    t_ende = float(sys.argv[2]) if len(sys.argv) > 2 else 12.0
    dt = float(sys.argv[3]) if len(sys.argv) > 3 else 5e-4
    S_ = aufbau(h_zyl)
    print("=" * 94)
    print(f" WIRBELSTRASSE NACH GALERKIN -- Fassung {FASSUNG}")
    print(f" {len(S_['P'])} Ecken, {S_['ne']} Kanten, {len(S_['T'])} Dreiecke, "
          f"{S_['nzyl']} Knoten auf dem Zylinder")
    print(f" Unbekannte: {2*S_['n2']} Geschwindigkeit (P2) + {S_['n1']} Druck (P1) "
          f"= {2*S_['n2']+S_['n1']}   dt = {dt:.1e} s")
    print("=" * 94)
    ux, uy, pr, verlauf, tz = lauf(S_, t_ende, dt, melde=max(1, int(0.5/dt)),
                                   speichern=f"fem_feld_{h_zyl}.npz")
    print(f"\n  Zerlegung der Systemmatrix: {tz:.1f} s (einmalig)")
    cD, cL, St, T = W.auswerten(verlauf, dt)
    print(f"\n{'Groesse':>10s} {'gerechnet':>12s} {'Referenz':>16s} {'Urteil':>14s}")
    for name, wert, (lo, hi) in (("c_D max", cD, W.REFERENZ["cD"]),
                                 ("c_L max", cL, W.REFERENZ["cL"]),
                                 ("St", St, W.REFERENZ["St"])):
        print(f"{name:>10s} {wert:12.4f} {f'{lo} … {hi}':>16s} "
              f"{'im Intervall' if lo <= wert <= hi else 'daneben':>14s}")
    json.dump({"h_zyl": h_zyl, "knoten": len(S_["P"]), "dreiecke": len(S_["T"]),
               "nzyl": S_["nzyl"], "unbekannte": 2*S_["n2"]+S_["n1"],
               "cD": cD, "cL": cL, "St": St, "T": T, "dt": dt},
              open(f"ergebnis_fem_{h_zyl}.json", "w"), indent=1)
