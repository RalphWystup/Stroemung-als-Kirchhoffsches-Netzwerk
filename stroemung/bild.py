#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Das Feldbild: Stoerkoerper und Wirbel.

Gezeigt wird die WIRBELSTAERKE omega = dv/dx - du/dy. Sie ist die richtige Groesse fuer
dieses Bild: Geschwindigkeitsbetrag oder Druck zeigen den Nachlauf nur verwaschen, waehrend
die Wirbelstaerke genau das sichtbar macht, was die Wirbelstrasse ausmacht -- abwechselnd
links- und rechtsdrehende Ballen.

Die Farbskala ist symmetrisch um null, damit Drehrichtung an der Farbe ablesbar ist und
nicht an der Helligkeit geraten werden muss.
"""
FASSUNG = "1.2"
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import wirbelstrasse as W


def wirbelstaerke(u, v, dx, dy):
    r"""omega = dv/dx - du/dy auf den INNEREN Gitterecken, plus die Ausdehnung fuer imshow.

    Fassung 1.1, Befund 31.10.2026. Fassung 1.0 bildete

        dvdx = (v[1:,:] - v[:-1,:])/dx          Form (nx-1, ny+1)
        dudy = (u[:,1:] - u[:,:-1])/dy          Form (nx+1, ny-1)

    und schnitt beide auf die kleinere Form zurecht. Das ist falsch: dvdx[0,0] sitzt auf
    der Ecke (1*dx, 0), dudy[0,0] dagegen auf (0, 1*dy). Voneinander abgezogen wurden also
    Werte von VERSCHIEDENEN Orten -- um je eine ganze Zelle versetzt. Der Fehler ist von
    der Ordnung h und zerstoert das Bild nicht, er verwischt es nur; genau deshalb ist er
    beim Ansehen nicht aufgefallen. Zusaetzlich wurde das Ergebnis mit der Ausdehnung
    [0, LX, 0, LY] gezeichnet, obwohl es gar nicht bis an den Rand reicht.

    Richtig liegen beide Anteile auf derselben inneren Ecke (i*dx, j*dy):

        omega[i,j] = (v[i,j] - v[i-1,j])/dx - (u[i,j] - u[i,j-1])/dy,
        i = 1 ... nx-1 ,  j = 1 ... ny-1

    Auf dem versetzten Gitter ist das die natuerliche Lage der Wirbelstaerke: Sie ist die
    Zirkulation um die Zelle, deren Ecken die vier Geschwindigkeitszweige sind -- im
    Netzwerkbild die Maschenregel, so wie die Divergenz die Knotenregel ist.
    """
    om = (v[1:, 1:-1] - v[:-1, 1:-1]) / dx - (u[1:-1, 1:] - u[1:-1, :-1]) / dy
    nx1, ny1 = om.shape                      # = (nx-1, ny-1)
    ausdehnung = [0.5 * dx, (nx1 + 0.5) * dx, 0.5 * dy, (ny1 + 0.5) * dy]
    return om, ausdehnung


def zeichne(datei, ziel, titel=None):
    d = np.load(datei)
    u, v, p, mp = d["u"], d["v"], d["p"], d["mp"]
    xp, yp = d["xp"], d["yp"]
    dx = float(xp[1] - xp[0]); dy = float(yp[1] - yp[0])
    om, ausdehnung = wirbelstaerke(u, v, dx, dy)
    gr = float(np.percentile(np.abs(om), 99))          # Ausreisser nicht die Skala bestimmen lassen

    fig, ax = plt.subplots(figsize=(13.0, 2.9), dpi=150)
    bild = ax.imshow(om.T, origin="lower", cmap="RdBu_r", vmin=-gr, vmax=gr,
                     extent=ausdehnung, aspect="equal", interpolation="bilinear")
    # Zylinder als echter Kreis darueber, nicht als Treppe
    ax.add_patch(plt.Circle((W.XC, W.YC), W.DZ / 2, facecolor="0.25",
                            edgecolor="black", linewidth=0.8, zorder=5))
    # Der Rahmen zeigt den GANZEN Kanal, auch wenn das Wirbelstaerkefeld die aeusserste
    # halbe Zelle nicht abdeckt: omega sitzt auf den inneren Gitterecken. Eine Wandformel
    # ueber den Spiegelwert waere moeglich, ist aber nur von erster Ordnung genau -- sie
    # haette in diesem Programm nichts zu suchen, das sonst exakt rechnet. Lieber ein
    # schmaler ungemalter Saum als eine Zahl, die anders zustande kommt als der Rest.
    ax.set_xlim(0, W.LX); ax.set_ylim(0, W.LY)
    ax.set_xlabel("x in m"); ax.set_ylabel("y in m")
    ax.set_title(titel or f"Kármánsche Wirbelstraße, Re = {W.RE:.0f} — Wirbelstärke",
                 fontsize=10)
    cb = fig.colorbar(bild, ax=ax, fraction=0.025, pad=0.012)
    cb.set_label("ω in 1/s", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(ziel, bbox_inches="tight")
    plt.close(fig)
    return om, gr


def probe():
    """Die Wirbelstaerke gegen Felder mit bekanntem omega -- analytisch, nicht geschaetzt."""
    from gitter import Gitter
    g = Gitter(W.LX, W.LY, 37, 19)
    XU, YU = np.meshgrid(g.xu, g.yp, indexing="ij")
    XV, YV = np.meshgrid(g.xp, g.yv, indexing="ij")
    faelle = []
    # 1 Starrkoerperdrehung u = -a(y-y0), v = a(x-x0)  ->  omega = 2a, ueberall gleich
    a = 3.7
    faelle.append(("Starrkoerperdrehung a=3,7", -a * (YU - 0.2), a * (XV - 0.2), 2 * a))
    # 2 Reine Scherung u = c*y, v = 0  ->  omega = -c
    c = -1.9
    faelle.append(("reine Scherung c=-1,9", c * YU, np.zeros_like(XV), -c))
    # 3 Dehnstroemung u = k*x, v = -k*y  ->  omega = 0 (drehungsfrei)
    k = 2.3
    faelle.append(("drehungsfreie Dehnstroemung", k * XU, -k * YV, 0.0))
    # 4+5 ORTSABHAENGIGE Wirbelstaerke -- nur sie deckt einen Versatz auf. Bei linearen
    # Feldern ist die Ableitung ueberall gleich; ein um eine Zelle verschobener Ablesepunkt
    # liefert denselben Wert, und der Fehler bleibt unsichtbar. Genau deshalb ist der
    # Versatz in Fassung 1.0 nicht aufgefallen.
    # v = x*y, u = 0  ->  omega = dv/dx = y
    XV2, YV2 = np.meshgrid(g.xp, g.yv, indexing="ij")
    XC2, YC2 = np.meshgrid(g.xu[1:-1], g.yv[1:-1], indexing="ij")
    faelle.append(("v = x*y, u = 0  ->  omega = y (ortsabhaengig)",
                   np.zeros_like(XU), XV2 * YV2, YC2))
    # u = x*y, v = 0  ->  omega = -du/dy = -x
    faelle.append(("u = x*y, v = 0  ->  omega = -x (ortsabhaengig)",
                   XU * YU, np.zeros_like(XV), -XC2))
    urteile = []
    for name, u, v, soll in faelle:
        om, _ = wirbelstaerke(u, v, g.dx, g.dy)
        d = float(np.abs(om - soll).max())
        bez = (f"omega = {soll:g}" if np.isscalar(soll)
               else f"omega ortsabhaengig, {np.min(soll):.4f} … {np.max(soll):.4f}")
        urteile.append((name, d < 1e-11, f"groesste Abweichung {d:.3e} von {bez}"))
    return urteile


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "probe":
        print("=" * 84)
        print(f" WIRBELSTAERKE GEGEN BEKANNTE FELDER -- Fassung {FASSUNG}")
        print("=" * 84)
        urteile = probe()
        for text, ok, beleg in urteile:
            print(f"  [{'RICHTIG' if ok else 'FALSCH  '}] {text}\n              {beleg}")
        n = sum(1 for _, ok, _ in urteile if ok)
        print(f"\n  == {n}/{len(urteile)}")
        raise SystemExit(0 if n == len(urteile) else 1)
    datei = sys.argv[1] if len(sys.argv) > 1 else "feld_220x41.npz"
    ziel = sys.argv[2] if len(sys.argv) > 2 else "wirbelstrasse.png"
    om, gr = zeichne(datei, ziel)
    print(f"  {ziel} geschrieben")
    print(f"  Wirbelstaerke: {om.min():.1f} bis {om.max():.1f} 1/s, Skala +-{gr:.1f}")
