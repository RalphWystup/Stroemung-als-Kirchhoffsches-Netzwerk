#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Dieselbe Wirbelstrasse, zweimal gerechnet -- auf zwei getrennt gebauten Wegen.

Das ist die Sprachprobe aus dem Feld-Labor, jetzt auf den VOLLSTAENDIGEN Loeser
angewandt, nicht nur auf einzelne Operatoren.

    Weg A -- NETZWERK NACH KIRCHHOFF
             Jeder Zweig wird einzeln gestempelt. Impuls und Druck stehen als
             duennbesetzte Knotenadmittanzmatrizen da; der Zeitschritt ist ein
             Matrix-Vektor-Produkt, das Drucknetz wird DIREKT geloest (LU).
             Die Unsymmetrie der Impulsmatrix ist die gesteuerte Quelle.

    Weg B -- KLASSISCHE FINITE VOLUMEN
             Keine Matrix wird je aufgestellt. Alle Operatoren sind Feldoperationen
             auf dem versetzten Gitter; das Drucknetz wird MATRIXFREI mit dem
             Verfahren der konjugierten Gradienten geloest.

Gemeinsam sind nur die Aufgabenstellung, das Gitter und der Zeitschritt. Stimmen die
Bahnen von c_D und c_L ueber Hunderte von Schritten ueberein, ist das Ergebnis nicht
die Eigenart eines Programms.

Aufruf:  python3 wirbelstrasse_zweiwege.py [nx] [ny] [t_ende]
"""
FASSUNG = "1.2"
import sys
import time
import json

import numpy as np
import scipy.sparse.linalg as spl
from scipy.sparse.linalg import LinearOperator

from gitter import Gitter
from ns_loeser import drucknetz_matrix, drucknetz_wirkung_voll
from konvektion import konv_zaeh_u, konv_zaeh_v
from impuls_matrix import impuls_matrix_u, impuls_matrix_v
from wirbelstrasse import (LX, LY, DZ, XC, YC, UM, NU, RHO, U_QUER,
                           einstroemprofil, koerper_masken, auswerten)

try:                                   # scipy >= 1.12 nennt es rtol
    _CG_ARG = "rtol"
    spl.cg(np.eye(2), np.ones(2), rtol=1e-10)
except TypeError:
    _CG_ARG = "tol"


def _raender(g, us, vs, ein):
    """Dieselben Randbedingungen fuer beide Wege -- sie gehoeren zur Aufgabe, nicht zum Weg."""
    vs[:, 0] = vs[:, -1] = 0.0
    vs[0, :] = 0.0
    us[0, :] = ein
    us[-1, :] = us[-2, :]
    vs[-1, :] = vs[-2, :]
    return us, vs


def lauf_zweiwege(nx, ny, t_ende=2.0, cfl=0.4, melde=50, beobachter=None):
    g = Gitter(LX, LY, nx, ny)
    mu_m, mv_m, mp_m = koerper_masken(g)
    ein = einstroemprofil(g.yu)
    dt = min(cfl * min(g.dx, g.dy) / (2.0 * UM),
             0.2 * min(g.dx, g.dy) ** 2 / NU)
    n = int(t_ende / dt)
    nenner = 0.5 * RHO * U_QUER**2 * DZ

    # ---- Weg A: alles als Matrix ----------------------------------------------------
    A_p, leer = drucknetz_matrix(g, ausstrom_rechts=True, u_maske=mu_m, v_maske=mv_m)
    loes_direkt = spl.factorized(A_p.tocsc())

    # ---- Weg B: nichts als Matrix ---------------------------------------------------
    N = nx * ny

    def wirkung(x):
        return drucknetz_wirkung_voll(g, x.reshape(nx, ny), leer=leer,
                                      ausstrom_rechts=True,
                                      u_maske=mu_m, v_maske=mv_m).ravel()
    Op = LinearOperator((N, N), matvec=wirkung, dtype=float)

    def start():
        u, v = g.leer_u(), g.leer_v()
        u[:, :] = ein[None, :]
        u[mu_m] = 0.0; v[mv_m] = 0.0
        return u, v

    uA, vA = start()
    uB, vB = start()
    pB_alt = np.zeros(N)
    bahn_A, bahn_B, abweichung = [], [], []
    t0 = time.time()

    for k in range(n):
        t = k * dt
        # ================= Weg A: Knotenanalyse, Matrizen, direkter Loeser =============
        Mu, bu = impuls_matrix_u(g, uA, vA, NU, einstroem=ein)
        Mv, bv = impuls_matrix_v(g, uA, vA, NU)
        us = uA + dt * (Mu @ uA.ravel() + bu).reshape(uA.shape)
        vs = vA + dt * (Mv @ vA.ravel() + bv).reshape(vA.shape)
        us, vs = _raender(g, us, vs, ein)
        fx1 = float(np.sum(us[mu_m])); fy1 = float(np.sum(vs[mv_m]))
        us[mu_m] = 0.0; vs[mv_m] = 0.0
        b = -g.divergenz(us, vs).ravel() * (g.dx * g.dy) / dt
        b[leer] = 0.0
        pA = loes_direkt(b).reshape(nx, ny)
        u2 = us.copy(); v2 = vs.copy()
        u2[1:-1, :] -= dt * (pA[1:, :] - pA[:-1, :]) / g.dx
        v2[:, 1:-1] -= dt * (pA[:, 1:] - pA[:, :-1]) / g.dy
        u2[-1, :] -= dt * (0.0 - pA[-1, :]) / (g.dx / 2.0)
        fx2 = float(np.sum(u2[mu_m])); fy2 = float(np.sum(v2[mv_m]))
        u2[mu_m] = 0.0; v2[mv_m] = 0.0
        u2[0, :] = ein
        uA, vA = u2, v2
        faktor = g.dx * g.dy * RHO / dt / nenner
        bahn_A.append((t, (fx1 + fx2) * faktor, (fy1 + fy2) * faktor))

        # ================= Weg B: Feldoperationen, matrixfreier CG ====================
        us = uB + dt * konv_zaeh_u(g, uB, vB, NU, einstroem=ein)
        vs = vB + dt * konv_zaeh_v(g, uB, vB, NU)
        us, vs = _raender(g, us, vs, ein)
        gx1 = float(np.sum(us[mu_m])); gy1 = float(np.sum(vs[mv_m]))
        us[mu_m] = 0.0; vs[mv_m] = 0.0
        b = -g.divergenz(us, vs).ravel() * (g.dx * g.dy) / dt
        b[leer] = 0.0
        pB, info = spl.cg(Op, b, x0=pB_alt, maxiter=2000, **{_CG_ARG: 1e-13})
        if info != 0:
            raise RuntimeError(f"CG in Weg B nicht konvergiert, info={info}, Schritt {k}")
        pB_alt = pB
        pB = pB.reshape(nx, ny)
        u2 = us.copy(); v2 = vs.copy()
        u2[1:-1, :] -= dt * (pB[1:, :] - pB[:-1, :]) / g.dx
        v2[:, 1:-1] -= dt * (pB[:, 1:] - pB[:, :-1]) / g.dy
        u2[-1, :] -= dt * (0.0 - pB[-1, :]) / (g.dx / 2.0)
        gx2 = float(np.sum(u2[mu_m])); gy2 = float(np.sum(v2[mv_m]))
        u2[mu_m] = 0.0; v2[mv_m] = 0.0
        u2[0, :] = ein
        uB, vB = u2, v2
        bahn_B.append((t, (gx1 + gx2) * faktor, (gy1 + gy2) * faktor))

        # Fassung 1.2: Nur-Lesen-Haken. Er rechnet nichts und aendert nichts; er reicht
        # BEIDE Felder am Ende des vollstaendigen Zeitschritts heraus. Gebraucht fuer das
        # Surrogat, das aus Momentaufnahmen beider Wege angelernt wird.
        if beobachter is not None:
            beobachter(k, t, g, uA, vA, pA, uB, vB, pB)
        bezug = max(np.abs(uA).max(), 1e-30)
        abweichung.append((float(np.abs(uA - uB).max() / bezug),
                           float(np.abs(vA - vB).max() / bezug),
                           float(abs(bahn_A[-1][1] - bahn_B[-1][1])),
                           float(abs(bahn_A[-1][2] - bahn_B[-1][2]))))
        if melde and k % melde == 0:
            du, dv, dcd, dcl = abweichung[-1]
            print(f"    t={t:6.3f}  A: cD={bahn_A[-1][1]:7.3f} cL={bahn_A[-1][2]:7.3f}   "
                  f"B: cD={bahn_B[-1][1]:7.3f} cL={bahn_B[-1][2]:7.3f}   "
                  f"|A-B|: u {du:.2e}  cD {dcd:.2e}   ({time.time()-t0:4.0f} s)", flush=True)
    return (g, np.array(bahn_A), np.array(bahn_B), np.array(abweichung), dt)


if __name__ == "__main__":
    nx = int(sys.argv[1]) if len(sys.argv) > 1 else 110
    ny = int(sys.argv[2]) if len(sys.argv) > 2 else 21
    t_ende = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
    print("=" * 104)
    print(f" WIRBELSTRASSE AUF ZWEI WEGEN -- Fassung {FASSUNG}   Netz {nx}x{ny}, "
          f"t bis {t_ende} s")
    print(" Weg A: Zweigstempel -> Matrix -> direkter Loeser   |   "
          "Weg B: Feldoperationen -> matrixfrei -> CG")
    print("=" * 104)
    g, bA, bB, abw, dt = lauf_zweiwege(nx, ny, t_ende=t_ende,
                                       melde=max(1, int(0.1 / 1e-3)))
    du = abw[:, 0].max(); dv = abw[:, 1].max()
    dcd = abw[:, 2].max(); dcl = abw[:, 3].max()
    print(f"\n  {len(bA)} Zeitschritte zu dt = {dt:.3e} s, beide Wege vom selben Anfang")
    # Fassung 1.1, Befund 31.10.2026: Hier stand bA[:,1].max() -- das ROHE Maximum ueber
    # den ganzen Lauf, also der Anlaufstoss des ersten Zeitschritts (gemessen 222,0). Als
    # A-gegen-B-Vergleich war das gueltig, als Kennzahl war es irrefuehrend. Jetzt wird
    # dieselbe Auswertung benutzt wie beim Produktionsloeser: eingeschwungener Teil,
    # Maximum darin, Strouhalzahl aus den Nulldurchgaengen. Damit ist zusaetzlich
    # nachgewiesen, dass beide Wege dieselben BEIWERTE liefern, nicht nur dieselben Felder.
    from wirbelstrasse import auswerten as _aus
    cDA, cLA, StA, TA = _aus(bA, dt)
    cDB, cLB, StB, TB = _aus(bB, dt)
    # Beiwerte nur melden, wenn die Abloesung wirklich ausgebildet ist. In einem zu kurzen
    # Lauf schwingt c_L noch gar nicht; die Nulldurchgaenge stammen dann aus Rundung, und
    # die daraus gebildete Strouhalzahl ist eine Zufallszahl. Sie dann hinzuschreiben
    # waere genau das, was hier nicht sein soll.
    ausgebildet = cLA > 0.05 and np.isfinite(StA) and TA * 3 < len(bA) * dt
    if ausgebildet:
        print(f"\n  Beiwerte im eingeschwungenen Teil, beide Wege getrennt ausgewertet:")
        print(f"      Weg A   c_D {cDA:.6f}   c_L {cLA:.6f}   St {StA:.6f}")
        print(f"      Weg B   c_D {cDB:.6f}   c_L {cLB:.6f}   St {StB:.6f}")
        print(f"      Unterschied  {abs(cDA-cDB):.3e}   {abs(cLA-cLB):.3e}   {abs(StA-StB):.3e}")
    else:
        print(f"\n  Beiwerte nicht ausgewertet: in {len(bA)*dt:.2f} s ist die Abloesung noch "
              f"nicht ausgebildet (groesstes c_L = {cLA:.4f}).")
        print(f"  Dieser Lauf prueft die UEBEREINSTIMMUNG beider Wege, nicht die Beiwerte.")
        cDA = cLA = StA = cDB = cLB = StB = float("nan")
    np.savez_compressed(f"bahnen_zweiwege_{nx}x{ny}.npz", A=bA, B=bB, abw=abw, dt=dt)

    urteile = [
        ("Geschwindigkeitsfeld u stimmt ueber den ganzen Lauf", du < 1e-9, f"{du:.3e} relativ"),
        ("Geschwindigkeitsfeld v stimmt ueber den ganzen Lauf", dv < 1e-9, f"{dv:.3e} relativ"),
        ("Widerstandsbeiwert c_D stimmt",                       dcd < 1e-6, f"{dcd:.3e} absolut"),
        ("Auftriebsbeiwert c_L stimmt",                         dcl < 1e-6, f"{dcl:.3e} absolut"),
    ]
    for text, ok, beleg in urteile:
        print(f"  [{'RICHTIG' if ok else 'FALSCH  '}] {text}\n              groesste Abweichung {beleg}")
    n_ok = sum(1 for _, ok, _ in urteile if ok)
    print(f"\n  == {n_ok}/{len(urteile)} -- " +
          ("dieselbe Stroemung auf beiden Wegen" if n_ok == len(urteile) else "ABWEICHUNG"))
    json.dump({"nx": nx, "ny": ny, "dt": dt, "schritte": len(bA),
               "du": du, "dv": dv, "dcD": dcd, "dcL": dcl,
               "cD_A": cDA, "cL_A": cLA, "St_A": StA,
               "cD_B": cDB, "cL_B": cLB, "St_B": StB},
              open(f"ergebnis_zweiwege_{nx}x{ny}.json", "w"), indent=1)
    raise SystemExit(0 if n_ok == len(urteile) else 1)
