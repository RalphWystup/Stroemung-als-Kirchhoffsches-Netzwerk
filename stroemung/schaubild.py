#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Das Schaubild: Wirbelstrasse, Ablösetakt und Netzkonvergenz auf einem Blatt.

bild.py (Fassung 1.0) bleibt unangetastet -- es ist geprueft und zeigt das nackte Feld.
Dies hier ist das Schaustueck daneben: Es zeigt nicht nur, DASS sich Wirbel ablösen,
sondern belegt im selben Blatt, dass die Rechnung auf die gemessene Wirklichkeit zulaeuft.

    oben       Wirbelstaerke omega = dv/dx - du/dy. Abwechselnd links- und
               rechtsdrehende Ballen; die Farbe ist die Drehrichtung.
    unten links  Auftriebsbeiwert ueber der Zeit. Eine volle Schwingung ist EIN Wirbelpaar;
               daraus kommt die Strouhalzahl.
    unten rechts Netzreihe: wie c_D, c_L und St mit feinerem Netz in die gemessenen
               Referenzbaender von Schaefer-Turek laufen.

Aufruf: python3 schaubild.py feld_880x164.npz wirbelstrasse_schaubild.png
"""
FASSUNG = "1.2"
import glob
import json
import math
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import wirbelstrasse as W
from bild import wirbelstaerke

# Farbwelt: kuehles Blaugrau als Grund, die Drehrichtung traegt die einzige kraeftige Farbe
GRUND = "#F7F8FA"
LINIE = "#2B3440"
MATT = "#8A94A6"
BAND = "#C8D4E8"


def _netzreihe():
    """Alle gerechneten Netzstufen einsammeln, nach Zellen je Durchmesser sortiert."""
    reihe = []
    for d in sorted(glob.glob("ergebnis_*x*.json")):
        try:
            e = json.load(open(d))
        except Exception:
            continue
        if "cD" not in e:
            continue
        if any(not math.isfinite(float(e[k])) for k in ("cD", "cL", "St")):
            continue          # Lauf zu kurz fuer eine volle Ablösung: keine Aussage
        zpd = e["nx"] / W.LX * W.DZ          # Zellen ueber dem Zylinderdurchmesser
        reihe.append((zpd, e))
    reihe.sort(key=lambda s: s[0])
    return reihe


def zeichne(datei, ziel):
    d = np.load(datei)
    u, v = d["u"], d["v"]
    xp, yp = d["xp"], d["yp"]
    dx = float(xp[1] - xp[0]); dy = float(yp[1] - yp[0])
    om, ausdehnung = wirbelstaerke(u, v, dx, dy)
    gr = float(np.percentile(np.abs(om), 99))
    verlauf = d["verlauf"] if "verlauf" in d else None

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.edgecolor": MATT, "axes.labelcolor": LINIE,
        "xtick.color": MATT, "ytick.color": MATT, "text.color": LINIE,
        "axes.linewidth": 0.8,
    })
    fig = plt.figure(figsize=(13.6, 6.4), dpi=170, facecolor=GRUND)
    gs = GridSpec(2, 2, height_ratios=[1.0, 0.72], width_ratios=[1.45, 1.0],
                  hspace=0.42, wspace=0.20, left=0.052, right=0.965,
                  top=0.875, bottom=0.085)

    # ---------- oben: das Feld ------------------------------------------------------
    ax = fig.add_subplot(gs[0, :])
    ax.set_facecolor(GRUND)
    bild = ax.imshow(om.T, origin="lower", cmap="RdBu_r", vmin=-gr, vmax=gr,
                     extent=ausdehnung, aspect="equal",
                     interpolation="bilinear")
    ax.add_patch(plt.Circle((W.XC, W.YC), W.DZ / 2, facecolor="#333A45",
                            edgecolor="#11151B", linewidth=1.0, zorder=5))
    # Der Rahmen zeigt den GANZEN Kanal, auch wenn das Wirbelstaerkefeld die aeusserste
    # halbe Zelle nicht abdeckt: omega sitzt auf den inneren Gitterecken. Eine Wandformel
    # ueber den Spiegelwert waere moeglich, ist aber nur von erster Ordnung genau -- sie
    # haette in diesem Programm nichts zu suchen, das sonst exakt rechnet. Lieber ein
    # schmaler ungemalter Saum als eine Zahl, die anders zustande kommt als der Rest.
    ax.set_xlim(0, W.LX); ax.set_ylim(0, W.LY)
    ax.set_xlabel("x in m", labelpad=2)
    ax.set_ylabel("y in m", labelpad=2)
    for s in ax.spines.values():
        s.set_color(MATT)
    cb = fig.colorbar(bild, ax=ax, fraction=0.0165, pad=0.008)
    cb.set_label("Wirbelstärke ω in 1/s", fontsize=8)
    cb.ax.tick_params(labelsize=7); cb.outline.set_edgecolor(MATT)
    ax.text(0.012, 1.10, "Kármánsche Wirbelstraße am Kreiszylinder",
            transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")
    ny_ = om.shape[1]
    ax.text(0.012, 1.015,
            f"Re = {W.RE:.0f}  ·  Kanal {W.LX} × {W.LY} m  ·  Zylinder D = {W.DZ} m  ·  "
            f"Netz {om.shape[0]+1} × {ny_+1}  ·  Scharfetter-Gummel-Zweigstempel auf "
            f"versetztem Gitter",
            transform=ax.transAxes, fontsize=8.5, color=MATT, va="bottom")

    # ---------- unten links: der Ablösetakt -----------------------------------------
    axl = fig.add_subplot(gs[1, 0]); axl.set_facecolor(GRUND)
    if verlauf is not None and len(verlauf) > 10:
        t, cL = verlauf[:, 0], verlauf[:, 2]
        maske = t >= t.max() - min(2.0, 0.45 * t.max())
        axl.axhline(0, color=MATT, lw=0.6)
        axl.plot(t[maske], cL[maske], color="#B23A3A", lw=1.3)
        # aufsteigende Nulldurchgaenge markieren: ein Punkt je Wirbelpaar
        tt, cc = t[maske], cL[maske]
        kr = np.where((cc[:-1] < 0) & (cc[1:] >= 0))[0]
        if len(kr) >= 2:
            T = float(np.median(np.diff(tt[kr])))
            for kx in tt[kr]:
                axl.axvline(kx, color=MATT, lw=0.6, ls=(0, (3, 3)))
            a, b = tt[kr[0]], tt[kr[1]]
            y = 1.14 * np.max(np.abs(cc))        # ueber der Kurve, nicht darauf
            axl.set_ylim(-1.30 * np.max(np.abs(cc)), 1.42 * np.max(np.abs(cc)))
            axl.annotate("", xy=(b, y), xytext=(a, y),
                         arrowprops=dict(arrowstyle="<->", color=LINIE, lw=1.0))
            axl.text(0.5 * (a + b), y * 1.10,
                     f"eine Ablösung  T = {T*1000:.0f} ms", ha="center",
                     va="bottom", fontsize=8.5, color=LINIE)
        axl.set_xlabel("Zeit t in s", labelpad=2)
        axl.set_ylabel("Auftriebsbeiwert $c_L$", labelpad=2)
    axl.set_title("Der Takt der Ablösung", fontsize=10.5, fontweight="bold",
                  loc="left", pad=6)
    for s in ("top", "right"):
        axl.spines[s].set_visible(False)

    # ---------- unten rechts: die Netzreihe ------------------------------------------
    axr = fig.add_subplot(gs[1, 1]); axr.set_facecolor(GRUND)
    reihe = _netzreihe()
    groessen = [("cD", "$c_D$", W.REFERENZ["cD"]),
                ("cL", "$c_L$", W.REFERENZ["cL"]),
                ("St", "$St$", W.REFERENZ["St"])]
    if reihe:
        # Aufgetragen wird die Abweichung von der Referenzmitte in Prozent. Die drei
        # Referenzintervalle sind VERSCHIEDEN breit (c_D ±0,31 %, c_L ±1,0 %, St ±1,7 %);
        # ein einziges graues Band waere daher falsch. Jede Groesse bekommt ihr eigenes.
        farbe = {"cD": "#B23A3A", "cL": "#2E6F9E", "St": "#4C8C56"}
        for schl, name, (lo, hi) in groessen:
            mitte = 0.5 * (lo + hi)
            axr.axhspan(100 * (lo / mitte - 1), 100 * (hi / mitte - 1),
                        color=farbe[schl], alpha=0.22, lw=0, zorder=0)
        axr.axhline(0.0, color=MATT, lw=0.8, zorder=1)
        for schl, name, (lo, hi) in groessen:
            mitte = 0.5 * (lo + hi)
            x = [s_[0] for s_ in reihe]
            y = [100 * (s_[1][schl] / mitte - 1) for s_ in reihe]
            # Der Endwert steht in der LEGENDE, nicht am Punkt: auf dem feinsten Netz
            # liegen die drei Kurven so dicht beieinander, dass sich Beschriftungen am
            # Punkt gegenseitig verdecken.
            ende = f"{y[-1]:+.1f} %".replace(".", ",").replace("-", "\u2212")
            axr.plot(x, y, "o-", lw=1.4, ms=5, label=f"{name}  {ende}",
                     color=farbe[schl], zorder=3)
        axr.set_xscale("log", base=2)
        axr.set_xticks([s_[0] for s_ in reihe])
        axr.set_xticklabels([f"{s_[0]:.0f}" for s_ in reihe])
        axr.set_xlim(min(x) / 1.25, max(x) * 1.18)
        axr.set_xlabel("Zellen über dem Durchmesser", labelpad=2)
        axr.set_ylabel("Abweichung in %", labelpad=2)
        axr.legend(frameon=False, fontsize=8.5, ncols=1, loc="lower right",
                   handlelength=1.6, labelspacing=0.35)
        axr.text(0.0, -0.30, "farbige Bänder = die gemessenen Referenzintervalle "
                 "(Schäfer-Turek 2D-2); sie sind verschieden breit",
                 transform=axr.transAxes, fontsize=7.5, color=MATT, va="top")
    axr.set_title("Was feineres Netz bringt", fontsize=10.5, fontweight="bold",
                  loc="left", pad=10)
    for s in ("top", "right"):
        axr.spines[s].set_visible(False)

    fig.savefig(ziel, facecolor=GRUND, bbox_inches="tight")
    plt.close(fig)
    return om, gr, reihe


if __name__ == "__main__":
    datei = sys.argv[1] if len(sys.argv) > 1 else "feld_440x82.npz"
    ziel = sys.argv[2] if len(sys.argv) > 2 else "wirbelstrasse_schaubild.png"
    om, gr, reihe = zeichne(datei, ziel)
    print(f"  {ziel} geschrieben")
    print(f"  Wirbelstaerke {om.min():.1f} … {om.max():.1f} 1/s, Skala ±{gr:.1f}")
    for zpd, e in reihe:
        print(f"  Netzstufe {e['nx']:4d}x{e['ny']:<4d} = {zpd:4.1f} Zellen/D   "
              f"cD {e['cD']:.4f}  cL {e['cL']:.4f}  St {e['St']:.4f}")
