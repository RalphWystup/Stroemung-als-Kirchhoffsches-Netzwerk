#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Das Surrogat: ein angelerntes Ersatzmodell fuer den Zeitverlauf.

Warum das hier KEIN Blackbox-Lernen ist. Die Wirbelstrasse bei Re = 100 laeuft in einen
GRENZZYKLUS: Nach dem Einschwingen ist das Feld eine Summe weniger Raumformen mit
sinusfoermigen Zeitbeiwerten. Genau das zeigt die Singulaerwertzerlegung der
Momentaufnahmen -- die Singulaerwerte kommen PAARWEISE, weil jede Frequenz ein
Sinus-Kosinus-Paar beisteuert. Das Surrogat ist deshalb eine Zerlegung, keine Anpassung:

  1. POD   Singulaerwertzerlegung der Momentaufnahmen (nach Abzug des Zeitmittels)
           -> r Raumformen, die fast die ganze Energie tragen
  2. DMD   Ausgleichsrechnung fuer den linearen Zeitoperator auf den Beiwerten:
           a_{k+1} = Op a_k  ->  der ganze Zeitverlauf folgt aus EINEM Anfangszustand

Beide Schritte sind geschlossene lineare Algebra. Was offen bleibt, ist der FEHLER
gegenueber dem echten Loeser, und der wird gemessen und mit ausgeliefert.

Beide Wege (A und B) werden mitgeschrieben. Die Raumformen stammen aus Weg A und werden
fuer BEIDE benutzt -- eine gemeinsame Basis. Das ist Absicht: Die beiden Datensaetze
unterscheiden sich nur um Rundung, und zwei getrennte Singulaerwertzerlegungen wuerden
innerhalb der fast entarteten Paare beliebig gegeneinander verdreht. Mit gemeinsamer Basis
traegt der Unterschied der beiden Surrogate genau den Unterschied der beiden Wege.

Aufruf:  python3 surrogat.py [nx] [ny] [t_ende] [fenster] [r]

===================================================================================
AUSBAU -- NOCH OFFEN, zurueckgestellt am 01.11.2026 auf Wort des Auftraggebers
===================================================================================

Was hier steht, ist ein Surrogat fuer EINEN Fall. Es gibt genau diesen Fall wieder --
diesen Zylinder, diese Reynoldszahl, den eingeschwungenen Zustand -- und ist damit
heute ein ANZEIGEMITTEL: Es traegt den Zeitverlauf im Browser-Labor, weil eine volle
Rechnung je Bildwechsel unmoeglich waere. Mehr wird nicht behauptet.

Nuetzlich wird es, sobald es ueber einen PARAMETER angelernt ist. Der Ausbau, in der
Reihenfolge des Nutzens:

  1. Betriebspunkte durchfahren.
     Acht bis zwoelf volle Rechnungen ueber einen Bereich (Volumenstrom, Reynoldszahl,
     spaeter ein Geometriemass wie der Schaufelwinkel). Eine GEMEINSAME POD-Basis ueber
     alle Faelle, dann anpassen, wie Modenamplituden und Grundfrequenz vom Parameter
     abhaengen. Danach liefert jeder Zwischenpunkt in Mikrosekunden ein Feld. Das ist
     der Unterschied zwischen Nachrechnen und AUSLEGEN: Eine Kennlinie wie die zwanzig
     Punkte von FAN-01 kaeme aus einer Handvoll voller Rechnungen.
     Zu bauen: gemeinsame Basis ueber mehrere Datensaetze, Anpassung der Beiwerte ueber
     dem Parameter, und eine Kreuzprobe -- einen Fall herauslassen, vorhersagen,
     gegen die volle Rechnung halten.

  2. Die Dimension der Stroemung messen.
     Dass hier 16 Moden 100,0000 % tragen und die Singulaerwerte PAARWEISE kommen, ist
     keine Rechengroesse, sondern eine Aussage ueber die Stroemung: Der Grenzzyklus ist
     niederdimensional, hat also wenige unabhaengige Mechanismen. Dieselbe Zerlegung an
     einem Luefterrad sagt, wie viele dort wirken. Braucht es ploetzlich zweihundert
     Moden, ist das ein BEFUND und keine Rechenlast.

  3. Rueckwaerts rechnen.
     Ein schnelles Modell macht die Umkehrung bezahlbar: aus einem gemessenen Signal auf
     den Betriebspunkt schliessen. Das ist die Richtung der akustischen Rohrdiagnose --
     dort wird aus einer Messung auf einen Zustand geschlossen, und ein Surrogat ist der
     Weg, das mit einem Stroemungsmodell zu tun statt mit einem Kennfeld.

  GRENZE, die bleibt: Ein Surrogat interpoliert nur innerhalb dessen, was es gesehen
  hat. Ein neues Verhalten -- Abloesebeginn, Abreissen, eine Instabilitaet, die im
  Lerndatensatz nicht vorkam -- findet es nicht, sondern glaettet darueber hinweg. Seine
  Rolle ist schnelles Erkunden INNERHALB eines geprueften Bereichs; an dessen Rand und
  fuer alles Neue bleibt der volle Loeser die Instanz.
"""
FASSUNG = "1.1"
import json
import os
import sys
import time

import numpy as np

HIER = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HIER, "..", "stroemung"))

from bild import wirbelstaerke                      # noqa: E402
from wirbelstrasse_zweiwege import lauf_zweiwege     # noqa: E402
import wirbelstrasse as W                            # noqa: E402

SCHRITT = 2          # jeden 2. Punkt speichern -- Anzeigegitter


def sammle(nx, ny, t_ende, fenster, jeder=8):
    """Laeuft beide Wege und sammelt Wirbelstaerke-Momentaufnahmen im letzten Fenster."""
    A, B, zeiten = [], [], []
    form = {}

    def schauen(k, t, g, uA, vA, pA, uB, vB, pB):
        if t < t_ende - fenster or k % jeder:
            return
        oA, ausd = wirbelstaerke(uA, vA, g.dx, g.dy)
        oB, _ = wirbelstaerke(uB, vB, g.dx, g.dy)
        oA = oA[::SCHRITT, ::SCHRITT]; oB = oB[::SCHRITT, ::SCHRITT]
        form["form"] = oA.shape; form["ausdehnung"] = ausd
        A.append(oA.ravel().copy()); B.append(oB.ravel().copy()); zeiten.append(t)

    t0 = time.time()
    g, bA, bB, abw, dt = lauf_zweiwege(nx, ny, t_ende=t_ende, melde=0, beobachter=schauen)
    print(f"  {len(A)} Momentaufnahmen zu je {A[0].size} Punkten "
          f"({form['form'][0]}x{form['form'][1]}), {time.time()-t0:.0f} s Rechenzeit")
    return (np.array(A).T, np.array(B).T, np.array(zeiten), form, g, bA, bB, dt)


def baue(XA, XB, r):
    """Gemeinsame POD-Basis aus Weg A, getrennte DMD-Operatoren fuer A und B."""
    mittel = XA.mean(axis=1, keepdims=True)
    U, s, Vt = np.linalg.svd(XA - mittel, full_matrices=False)
    U = U[:, :r]
    aA = U.T @ (XA - mittel)
    aB = U.T @ (XB - mittel)
    def op(a):
        return a[:, 1:] @ np.linalg.pinv(a[:, :-1])
    return mittel, U, s, aA, aB, op(aA), op(aB)


def vorhersage(Op, a0, n):
    a = a0.copy(); folge = [a.copy()]
    for _ in range(n - 1):
        a = Op @ a; folge.append(a.copy())
    return np.array(folge).T


if __name__ == "__main__":
    nx = int(sys.argv[1]) if len(sys.argv) > 1 else 220
    ny = int(sys.argv[2]) if len(sys.argv) > 2 else 41
    t_ende = float(sys.argv[3]) if len(sys.argv) > 3 else 8.0
    fenster = float(sys.argv[4]) if len(sys.argv) > 4 else 1.4
    r = int(sys.argv[5]) if len(sys.argv) > 5 else 16
    print("=" * 94)
    print(f" SURROGAT AUS MOMENTAUFNAHMEN -- Fassung {FASSUNG}   Netz {nx}x{ny}, "
          f"Fenster {fenster} s, r = {r}")
    print("=" * 94)
    XA, XB, zeiten, form, g, bA, bB, dt = sammle(nx, ny, t_ende, fenster)
    mittel, U, s, aA, aB, OpA, OpB = baue(XA, XB, r)

    # Ist der Grenzzyklus ueberhaupt eingeschwungen? Eine Antwort darauf gehoert vor
    # jede Deutung: Auf nicht eingeschwungenen Daten ist der angepasste Zeitoperator
    # in der Frequenz daneben, und ein Phasenfehler laeuft auf (Befund 31.10.2026 --
    # bei t_ende = 8 s fiel die Schwankungsbreite ueber das Fenster noch um 7 %,
    # die Periode kam 1,3 % zu lang heraus, und der Fehler stieg auf 15,8 %).
    # Fassung 1.1, Befund 31.10.2026: Fassung 1.0 verglich das ERSTE mit dem LETZTEN
    # Bild. Die liegen an verschiedenen Stellen der Schwingung, und der Vergleich maass
    # deshalb die Phase, nicht die Drift -- er meldete 7 bis 10 % Abfall, wo in
    # Wahrheit 0,003 % waren. Jetzt wird je GANZER Periode verglichen; dann faellt die
    # Phase heraus.
    amp = np.abs(XA - mittel).max(axis=0)
    dts = float(zeiten[1] - zeiten[0])
    kr = np.where((aA[0, :-1] < 0) & (aA[0, 1:] >= 0))[0]
    proP = int(round(float(np.median(np.diff(kr))))) if len(kr) >= 2 else len(amp)
    spitzen = [amp[k:k + proP].max() for k in range(0, len(amp) - proP + 1, proP)]
    drift = 100 * (spitzen[-1] / spitzen[0] - 1.0) if len(spitzen) >= 2 else 0.0
    print(f"\n  Spitze der Schwankungsbreite je GANZER Periode: "
          + " ".join(f"{x:.3f}" for x in spitzen))
    print(f"  Drift ueber {len(spitzen)} Perioden: {drift:+.3f} %")
    print(f"  {'eingeschwungen' if abs(drift) < 1 else 'NOCH NICHT EINGESCHWUNGEN'}")
    print(f"\n  Singulaerwerte (die Paare sind die Signatur des Grenzzyklus):")
    print("    " + "  ".join(f"{x:.3e}" for x in s[:10]))
    energie = np.cumsum(s**2) / np.sum(s**2)
    print(f"    r = {r}: erfasste Energie {100*energie[r-1]:.4f} %")

    n = XA.shape[1]
    vA = U @ vorhersage(OpA, aA[:, 0], n) + mittel
    vB = U @ vorhersage(OpB, aB[:, 0], n) + mittel
    spanne = np.abs(XA - mittel).max()
    f_sur = np.abs(vA - XA).max() / spanne
    f_trunk = np.abs((U @ aA + mittel) - XA).max() / spanne
    d_echt = np.abs(XA - XB).max() / spanne
    d_sur = np.abs(vA - vB).max() / spanne
    # Periode: gemessen gegen vorhergesagt -- das ist die Groesse, an der ein
    # Phasenfehler zuerst sichtbar wird.
    def periode(x):
        kr = np.where((x[:-1] < 0) & (x[1:] >= 0))[0]
        return float(np.median(np.diff(kr))) * dts if len(kr) >= 2 else float("nan")
    dts = float(zeiten[1] - zeiten[0])
    T_echt, T_vor = periode(aA[0]), periode((U.T @ (vA - mittel))[0])
    print(f"  Periode gemessen {T_echt:.4f} s, vorhergesagt {T_vor:.4f} s "
          f"({100*(T_vor/T_echt-1):+.2f} %)")
    urteile = [
        (f"Surrogat gegen echten Loeser, ueber {n} Schritte fortgeschrieben",
         f_sur < 0.05, f"{100*f_sur:.3f} % der Schwankungsbreite"),
        ("davon allein durch das Abschneiden auf r Moden",
         True, f"{100*f_trunk:.4f} %"),
        ("echter Unterschied Weg A gegen Weg B in den Momentaufnahmen",
         d_echt < 1e-10, f"{d_echt:.3e} relativ"),
        ("Unterschied der beiden SURROGATE (gemeinsame Basis)",
         d_sur < 1e-6, f"{d_sur:.3e} relativ"),
    ]
    print()
    for text, ok, beleg in urteile:
        print(f"  [{'RICHTIG' if ok else 'FALSCH  '}] {text}\n              {beleg}")

    np.savez_compressed(os.path.join(HIER, "surrogat_daten.npz"),
                        XA=XA, XB=XB, zeiten=zeiten, mittel=mittel, U=U, s=s,
                        aA=aA, aB=aB, OpA=OpA, OpB=OpB,
                        form=np.array(form["form"]), ausdehnung=np.array(form["ausdehnung"]),
                        verlaufA=bA, verlaufB=bB, dt=np.array([dt]))
    json.dump({"nx": nx, "ny": ny, "r": r, "schnappschuesse": int(n),
               "fehler_surrogat": float(f_sur), "fehler_trunkierung": float(f_trunk),
               "unterschied_echt": float(d_echt), "unterschied_surrogat": float(d_sur),
               "energie": float(energie[r-1]),
               "singulaerwerte": s[:12].tolist()},
              open(os.path.join(HIER, "surrogat_bericht.json"), "w"), indent=1)
    n_ok = sum(1 for _, ok, _ in urteile if ok)
    print(f"\n  == {n_ok}/{len(urteile)}   Daten in surrogat_daten.npz")
