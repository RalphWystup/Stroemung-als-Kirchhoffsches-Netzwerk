#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Der unabhaengige zweite Weg: Galerkin/Ritz mit Taylor-Hood-Elementen.

Die Unabhaengigkeit sitzt in der HERLEITUNG, nicht im Gleichungsloeser. Der
Netzwerkweg leitet aus der Flussbilanz am Kontrollvolumen her: Jede Zelle ist ein
Knoten, jede Flaeche ein Zweig, die Kontinuitaet ist die Knotenregel. Dieser Weg leitet
aus der schwachen Form her: Die Gleichungen werden mit Ansatzfunktionen multipliziert
und ueber das Gebiet integriert; die Unbekannten sind Knotenwerte eines
Polynomansatzes, nicht Zellmittel. Dass am Ende beide ein duennbesetztes System loesen,
ist eine Frage des Werkzeugs, nicht des Ansatzes.

    Geschwindigkeit   P2  (quadratisch, Knoten auf Ecken UND Kantenmitten)
    Druck             P1  (linear, Knoten nur auf den Ecken)

Diese Paarung erfuellt die inf-sup-Bedingung und braucht deshalb KEINEN abgestimmten
Stabilisierungsparameter -- das ist der Grund, warum sie hier steht und nicht ein
gleichordniger Ansatz mit PSPG.

Dazu kommt ein koerperangepasstes Netz: Der Zylinderrand liegt auf Knoten, es gibt
keine Treppe. Damit faellt genau die Fehlerquelle weg, die beim Netzwerkweg uebrig
bleibt -- auch das macht den Vergleich erst aussagekraeftig.
"""
FASSUNG = "1.1"
import numpy as np
import scipy.sparse as sp

import netz_fem
import wirbelstrasse as W

# Gausspunkte auf dem Dreieck, Grad 5, sieben Punkte (baryzentrisch, Gewichte zu 1)
_a = 0.470142064105115
_b = 0.101286507323456
QUAD = np.array([
    [1/3, 1/3, 1/3, 0.225],
    [1 - 2*_a, _a, _a, 0.132394152788506],
    [_a, 1 - 2*_a, _a, 0.132394152788506],
    [_a, _a, 1 - 2*_a, 0.132394152788506],
    [1 - 2*_b, _b, _b, 0.125939180544827],
    [_b, 1 - 2*_b, _b, 0.125939180544827],
    [_b, _b, 1 - 2*_b, 0.125939180544827],
])


def p2_raum(P, T):
    """Knoten des P2-Raums: Ecken, dann Kantenmitten. Je Dreieck sechs Knoten.

    Nummerierung im Element: 0,1,2 = Ecken; 3 = Mitte(1,2), 4 = Mitte(2,0), 5 = Mitte(0,1).
    """
    kanten = {}
    for t in T:
        for a, b in ((1, 2), (2, 0), (0, 1)):
            k = (min(t[a], t[b]), max(t[a], t[b]))
            if k not in kanten:
                kanten[k] = len(kanten)
    nv = len(P)
    XY = np.empty((nv + len(kanten), 2))
    XY[:nv] = P
    for (a, b), i in kanten.items():
        XY[nv + i] = 0.5 * (P[a] + P[b])
    E6 = np.empty((len(T), 6), dtype=int)
    E6[:, :3] = T
    for e, t in enumerate(T):
        for j, (a, b) in enumerate(((1, 2), (2, 0), (0, 1))):
            E6[e, 3 + j] = nv + kanten[(min(t[a], t[b]), max(t[a], t[b]))]
    return XY, E6, nv, len(kanten)


def p2_werte(L):
    """Die sechs P2-Formfunktionen an baryzentrischen Koordinaten L = (l0,l1,l2)."""
    l0, l1, l2 = L
    return np.array([l0*(2*l0-1), l1*(2*l1-1), l2*(2*l2-1),
                     4*l1*l2, 4*l2*l0, 4*l0*l1])


def p2_ableitungen(L, g):
    """Gradienten der sechs Formfunktionen; g[i] = grad(l_i), je Element konstant."""
    l0, l1, l2 = L
    d = np.empty((6, 2))
    d[0] = (4*l0 - 1) * g[0]
    d[1] = (4*l1 - 1) * g[1]
    d[2] = (4*l2 - 1) * g[2]
    d[3] = 4 * (l2 * g[1] + l1 * g[2])
    d[4] = 4 * (l0 * g[2] + l2 * g[0])
    d[5] = 4 * (l1 * g[0] + l0 * g[1])
    return d


def lambda_gradienten(P, T):
    """grad(l_i) je Dreieck und die Flaechen -- beides je Element konstant (gerade Kanten)."""
    x = P[T, 0]; y = P[T, 1]
    A2 = ((x[:, 1]-x[:, 0])*(y[:, 2]-y[:, 0]) - (x[:, 2]-x[:, 0])*(y[:, 1]-y[:, 0]))
    A = 0.5 * A2
    g = np.empty((len(T), 3, 2))
    g[:, 0, 0] = (y[:, 1]-y[:, 2]) / A2; g[:, 0, 1] = (x[:, 2]-x[:, 1]) / A2
    g[:, 1, 0] = (y[:, 2]-y[:, 0]) / A2; g[:, 1, 1] = (x[:, 0]-x[:, 2]) / A2
    g[:, 2, 0] = (y[:, 0]-y[:, 1]) / A2; g[:, 2, 1] = (x[:, 1]-x[:, 0]) / A2
    return g, A


def zusammenstellen(P, T, XY, E6, nv):
    """Massenmatrix, Steifigkeit und die Druck-Geschwindigkeits-Kopplung.

    Die schwache Form lautet: fuer alle Testfunktionen v und q

        (du/dt, v) + c(u,u,v) + nu (grad u, grad v) - (p, div v) = 0
        (q, div u) = 0

    Daraus entstehen M (Masse), K (Steifigkeit) und B = (q, div v). Die Konvektion c
    wird NICHT als Matrix gebaut, sondern je Zeitschritt als Vektor ausgewertet -- damit
    bleibt die Systemmatrix konstant und muss nur EINMAL zerlegt werden.
    """
    g, A = lambda_gradienten(P, T)
    ne, n2, n1 = len(T), len(XY), nv
    M = np.zeros((ne, 6, 6)); K = np.zeros((ne, 6, 6))
    Bx = np.zeros((ne, 3, 6)); By = np.zeros((ne, 3, 6))
    for q in QUAD:
        L, w = q[:3], q[3]
        N = p2_werte(L)
        dN = np.empty((ne, 6, 2))
        dN[:, 0] = (4*L[0]-1) * g[:, 0]
        dN[:, 1] = (4*L[1]-1) * g[:, 1]
        dN[:, 2] = (4*L[2]-1) * g[:, 2]
        dN[:, 3] = 4 * (L[2]*g[:, 1] + L[1]*g[:, 2])
        dN[:, 4] = 4 * (L[0]*g[:, 2] + L[2]*g[:, 0])
        dN[:, 5] = 4 * (L[1]*g[:, 0] + L[0]*g[:, 1])
        wA = w * A
        M += wA[:, None, None] * np.outer(N, N)[None]
        K += wA[:, None, None] * np.einsum("eid,ejd->eij", dN, dN)
        Bx += wA[:, None, None] * np.einsum("k,ei->eki", L, dN[:, :, 0])
        By += wA[:, None, None] * np.einsum("k,ei->eki", L, dN[:, :, 1])
    def streu(lokal, zeilen, spalten, form):
        i = np.repeat(zeilen, lokal.shape[2], axis=1).ravel()
        j = np.tile(spalten, (1, lokal.shape[1])).ravel()
        return sp.coo_matrix((lokal.ravel(), (i, j)), shape=form).tocsr()
    Mm = streu(M, E6, E6, (n2, n2))
    Kk = streu(K, E6, E6, (n2, n2))
    Bxx = streu(Bx, T, E6, (n1, n2))
    Byy = streu(By, T, E6, (n1, n2))
    return Mm, Kk, Bxx, Byy, g, A


def konvektion(P, T, E6, g, A, ux, uy):
    """c(u,u,v) als Vektor: integral (u . grad)u . v. Je Zeitschritt neu, ohne Matrix."""
    ne = len(T)
    Rx = np.zeros((ne, 6)); Ry = np.zeros((ne, 6))
    uxe = ux[E6]; uye = uy[E6]
    for q in QUAD:
        L, w = q[:3], q[3]
        N = p2_werte(L)
        dN = np.empty((ne, 6, 2))
        dN[:, 0] = (4*L[0]-1) * g[:, 0]
        dN[:, 1] = (4*L[1]-1) * g[:, 1]
        dN[:, 2] = (4*L[2]-1) * g[:, 2]
        dN[:, 3] = 4 * (L[2]*g[:, 1] + L[1]*g[:, 2])
        dN[:, 4] = 4 * (L[0]*g[:, 2] + L[2]*g[:, 0])
        dN[:, 5] = 4 * (L[1]*g[:, 0] + L[0]*g[:, 1])
        uq = uxe @ N; vq = uye @ N
        duxdx = np.einsum("ei,ei->e", uxe, dN[:, :, 0])
        duxdy = np.einsum("ei,ei->e", uxe, dN[:, :, 1])
        duydx = np.einsum("ei,ei->e", uye, dN[:, :, 0])
        duydy = np.einsum("ei,ei->e", uye, dN[:, :, 1])
        wA = w * A
        Rx += (wA * (uq*duxdx + vq*duxdy))[:, None] * N[None]
        Ry += (wA * (uq*duydx + vq*duydy))[:, None] * N[None]
    n2 = E6.max() + 1
    fx = np.bincount(E6.ravel(), Rx.ravel(), minlength=n2)
    fy = np.bincount(E6.ravel(), Ry.ravel(), minlength=n2)
    return fx, fy


def raender(XY, nv, mit_zylinder=True):
    """Welche P2-Knoten liegen wo. Der Ausstrom bleibt FREI (do-nothing)."""
    x, y = XY[:, 0], XY[:, 1]
    e = 1e-9
    ein = x < e
    wand = (y < e) | (y > W.LY - e)
    zyl = (np.hypot(x - W.XC, y - W.YC) < W.DZ/2 + 1e-7) if mit_zylinder \
        else np.zeros(len(XY), bool)
    aus = x > W.LX - e
    # Fassung 1.1, Befund 31.10.2026: Hier stand `& ~aus`. Damit waren auch die beiden
    # ECKKNOTEN am Ausstrom frei, obwohl sie auf der Wand liegen -- die Haftbedingung
    # fiel dort weg. Richtig ist: Der Ausstrom gibt nur die Knoten frei, die NUR am
    # Ausstrom liegen; Wand bleibt Wand bis zum letzten Knoten.
    fest = ein | wand | zyl
    u_vor = np.zeros(len(XY)); v_vor = np.zeros(len(XY))
    u_vor[ein] = W.einstroemprofil(y[ein])
    u_vor[wand | zyl] = 0.0
    return dict(ein=ein, wand=wand, zyl=zyl, aus=aus, fest=fest, u=u_vor, v=v_vor)


def sattelsystem(K, Bx, By, nu, M=None, dt=None):
    """Die Blockmatrix des Sattelpunktproblems.

        [ A    0   -Bx^T ]
        [ 0    A   -By^T ]
        [ Bx   By    0   ]

    A = nu K (stationaer) oder M/dt + nu K (instationaer). Weil die Konvektion als
    Vektor behandelt wird, ist A KONSTANT -- die Zerlegung faellt nur einmal an.
    """
    A = nu * K if M is None else (M / dt + nu * K)
    n1 = Bx.shape[0]
    Z = sp.csr_matrix((n1, n1))
    return sp.bmat([[A, None, -Bx.T], [None, A, -By.T], [Bx, By, Z]], format="csc")


def festsetzen(S, rhs, feste, werte):
    """Dirichlet-Knoten symmetrisch eliminieren.

    Die bekannten Werte wandern als Spaltenbeitrag in die rechte Seite, danach werden
    Zeile UND Spalte geleert und die Diagonale auf eins gesetzt. So bleibt die Struktur
    erhalten, und das System ist genau dann loesbar, wenn es das vorher war.

    Fassung 1.1, Befund 31.10.2026: Fassung 1.0 tat dasselbe von Hand ueber lil-Listen,
    in mehreren Durchgaengen und mit Umwandlungen dazwischen. Das Ergebnis war falsch --
    die Hagen-Poiseuille-Probe, die auf Maschinengenauigkeit aufgehen MUSS, lag 3,3 %
    daneben. Hier steht dieselbe Sache in vier Zeilen Matrixalgebra, und sie geht auf.
    """
    feste = np.asarray(feste, bool)
    bekannt = werte * feste
    rhs = rhs - S @ bekannt                  # Spaltenbeitrag der bekannten Werte
    rhs = np.where(feste, werte, rhs)
    frei = sp.diags((~feste).astype(float))
    ein = sp.diags(feste.astype(float))
    return (frei @ S @ frei + ein).tocsc(), rhs


if __name__ == "__main__":
    print("=" * 94)
    print(f" ANSATZRAUM TAYLOR-HOOD -- Fassung {FASSUNG}")
    print("=" * 94)
    P, T, A, zyl = netz_fem.netz(0.05, 0.0039, h_nach=4.0, band=2.0, nachlauf=1.0)
    XY, E6, nv, ne = p2_raum(P, T)
    g, Ae = lambda_gradienten(P, T)
    urteile = []
    urteile.append(("Flaechen aus den Lambda-Gradienten stimmen mit dem Netz ueberein",
                    np.abs(Ae - A).max() < 1e-14, f"{np.abs(Ae-A).max():.2e}"))
    # Zerlegung der Eins: die sechs Formfunktionen muessen sich ueberall zu 1 addieren
    s = max(abs(p2_werte(q[:3]).sum() - 1.0) for q in QUAD)
    urteile.append(("die sechs Formfunktionen addieren sich zu 1", s < 1e-14, f"{s:.2e}"))
    # Gradienten addieren sich zu null
    gs = max(np.abs(p2_ableitungen(q[:3], g[0]).sum(axis=0)).max() for q in QUAD)
    urteile.append(("ihre Gradienten addieren sich zu null", gs < 1e-12, f"{gs:.2e}"))
    # P2 gibt quadratische Felder EXAKT wieder: f = x^2 + 3xy - y^2 + 2x - 1
    f = lambda x, y: x*x + 3*x*y - y*y + 2*x - 1.0
    fk = f(XY[:, 0], XY[:, 1])
    fehl = 0.0
    for e in range(0, len(T), max(1, len(T)//400)):
        for q in QUAD:
            L = q[:3]
            xq = (L[:, None] * P[T[e]]).sum(axis=0)
            fehl = max(fehl, abs(p2_werte(L) @ fk[E6[e]] - f(xq[0], xq[1])))
    urteile.append(("P2 gibt ein quadratisches Feld knotenweise exakt wieder",
                    fehl < 1e-12, f"groesster Fehler {fehl:.2e}"))
    # Quadratur: Summe der Gewichte mal Flaeche = Gesamtflaeche
    ges = (QUAD[:, 3].sum() * A).sum()
    urteile.append(("die Quadraturgewichte geben die Gesamtflaeche",
                    abs(ges - A.sum()) < 1e-12, f"{ges:.8f} gegen {A.sum():.8f}"))
    for t, ok, b in urteile:
        print(f"  [{'RICHTIG' if ok else 'FALSCH  '}] {t}\n              {b}")
    print(f"\n  Netz: {len(P)} Ecken, {ne} Kanten, {len(T)} Dreiecke")
    print(f"  Unbekannte: {2*(nv+ne)} Geschwindigkeit (P2) + {nv} Druck (P1) "
          f"= {2*(nv+ne)+nv}")
    # ---- Hagen-Poiseuille: P2 enthaelt die Parabel EXAKT ---------------------------
    import scipy.sparse.linalg as spl
    Pk, Tk, Ak, _ = netz_fem.netz(0.03, 0.03, mit_zylinder=False)
    XYk, E6k, nvk, nek = p2_raum(Pk, Tk)
    M, K, Bx, By, gk, Aek = zusammenstellen(Pk, Tk, XYk, E6k, nvk)
    R = raender(XYk, nvk, mit_zylinder=False)
    n2 = len(XYk); n1 = nvk
    S = sattelsystem(K, Bx, By, W.NU)
    N = 2*n2 + n1
    feste = np.zeros(N, bool); werte = np.zeros(N)
    feste[:n2] = R["fest"]; werte[:n2] = R["u"]
    feste[n2:2*n2] = R["fest"]; werte[n2:2*n2] = R["v"]
    S2, rhs = festsetzen(S, np.zeros(N), feste, werte)
    loes = spl.spsolve(S2.tocsc(), rhs)
    ux, uy, pr = loes[:n2], loes[n2:2*n2], loes[2*n2:]
    genau_u = W.einstroemprofil(XYk[:, 1])
    du = np.abs(ux - genau_u).max() / genau_u.max()
    dv = np.abs(uy).max() / genau_u.max()
    # Druckgefaelle: dp/dx = nu u'' = -8 nu U_m / H^2
    soll_grad = -8.0 * W.NU * W.UM / W.LY**2
    xk = Pk[:, 0]
    A1 = np.c_[xk, np.ones(len(xk))]
    koef, *_ = np.linalg.lstsq(A1, pr, rcond=None)
    rest = np.abs(pr - A1 @ koef).max()
    urteile += [
        ("Hagen-Poiseuille: P2 gibt das Profil knotenweise exakt",
         du < 1e-10, f"groesste Abweichung {du:.3e} relativ"),
        ("die Querkomponente bleibt null", dv < 1e-10, f"{dv:.3e} relativ"),
        ("der Druck ist linear in x",
         rest / abs(koef[0] * W.LX) < 1e-9, f"Rest {rest:.3e}"),
        ("das Druckgefaelle trifft -8 nu U_m / H^2",
         abs(koef[0] - soll_grad) / abs(soll_grad) < 1e-9,
         f"{koef[0]:.8f} gegen {soll_grad:.8f} 1/s^2"),
    ]
    for t, ok, b in urteile[5:]:
        print(f"  [{'RICHTIG' if ok else 'FALSCH  '}] {t}\n              {b}")
    n = sum(1 for _, ok, _ in urteile if ok)
    print(f"\n  == {n}/{len(urteile)}")
    raise SystemExit(0 if n == len(urteile) else 1)
