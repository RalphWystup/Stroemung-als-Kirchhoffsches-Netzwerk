#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Das Pruefprotokoll des Stroemungslösers -- ein Aufruf, alle Stufen.

Jede Stufe ist ein eigenstaendiges Programm, das sich SELBST beurteilt und mit
Rueckgabewert 0 oder 1 endet. Dieses Programm ruft sie der Reihe nach auf und sammelt die
Urteile. Es rechnet selbst nichts nach -- sonst waere es eine zweite Quelle der Wahrheit.

    schnell   die kurzen Stufen (Sekunden bis Minuten)
    ganz      zusaetzlich die langen Laeufe (Zweiwege-Vergleich, Netzreihe)

Aufruf:  python3 pruefprotokoll.py [schnell|ganz]
"""
FASSUNG = "1.1"
import glob
import json
import re
import math
import os
import subprocess
import sys
import time

HIER = os.path.dirname(os.path.abspath(__file__))

# (Kuerzel, Was geprueft wird, Befehl, Zeitgrenze in s, nur bei "ganz")
STUFEN = [
    ("S01", "Zweigstempel nach Scharfetter-Gummel gegen die geschlossene Loesung",
     [sys.executable, "sg_stempel.py"], 300, False),
    ("S02", "Knotenregel auf dem versetzten Gitter (divergenzfreies Feld)",
     [sys.executable, "gitter.py"], 120, False),
    ("S03", "Drucknetz: Zweigstempel gegen Feldoperation (Sprachprobe)",
     [sys.executable, "ns_loeser.py"], 300, False),
    ("S04", "Drucknetz MIT Stoerkoerper und Ausstrom, beide Wege",
     [sys.executable, "drucknetz_koerper_probe.py"], 600, False),
    ("S05", "Ebener Kanal gegen Hagen-Poiseuille, Konvergenzordnung h^2",
     [sys.executable, "poiseuille.py"], 600, False),
    ("S06", "Konvektionsstempel: Eigenschaften, die ohne Rechnung gelten muessen",
     [sys.executable, "konvektion.py"], 300, False),
    ("S07", "Impulsoperator: Matrix gegen Feldoperation, Symmetrie im Grenzfall",
     [sys.executable, "impuls_matrix.py"], 600, False),
    ("S08", "Knotenregel am ENDE jedes vollstaendigen Zeitschritts",
     [sys.executable, "knotenregel_probe.py", "110", "21", "1.0"], 900, False),
    ("S09", "Wirbelstaerke im Bild gegen Felder mit bekanntem omega",
     [sys.executable, "bild.py", "probe"], 300, False),
    ("S10", "Vollstaendiger Loeser auf zwei getrennt gebauten Wegen (kurz)",
     [sys.executable, "wirbelstrasse_zweiwege.py", "110", "21", "0.5"], 900, False),
    ("S11", "Vollstaendiger Loeser auf zwei Wegen, bis in die Abloesung",
     [sys.executable, "wirbelstrasse_zweiwege.py", "220", "41", "8.0"], 7200, True),
    ("S12", "Lattice-Boltzmann gegen Taylor-Green und Hagen-Poiseuille",
     [sys.executable, "lbm.py"], 1200, False),
    ("S13", "Koerperangepasstes Netz: Flaechen, Rand, Kreistreue",
     [sys.executable, "netz_fem.py"], 900, False),
    ("S14", "Galerkin-Ansatzraum und Hagen-Poiseuille EXAKT (P2 enthaelt die Parabel)",
     [sys.executable, "fem_galerkin.py"], 1800, False),
]


def netzreihe_urteil():
    """S15: Laeuft die Netzreihe auf die Referenz zu? Geprueft wird der TREND, nicht ein Lauf."""
    import wirbelstrasse as W
    stufen = []
    for d in sorted(glob.glob(os.path.join(HIER, "ergebnis_*x*.json"))):
        if "zweiwege" in d:
            continue
        e = json.load(open(d))
        if "cD" not in e or any(not math.isfinite(float(e[k])) for k in ("cD", "cL", "St")):
            continue
        stufen.append((e["nx"] / W.LX * W.DZ, e))
    stufen.sort(key=lambda s: s[0])
    if len(stufen) < 2:
        return None, "weniger als zwei auswertbare Netzstufen vorhanden", []
    zeilen, alle_besser = [], True
    for schl, (lo, hi) in (("cD", W.REFERENZ["cD"]), ("cL", W.REFERENZ["cL"]),
                           ("St", W.REFERENZ["St"])):
        mitte = 0.5 * (lo + hi)
        abw = [abs(e[schl] / mitte - 1.0) for _, e in stufen]
        besser = all(abw[i + 1] < abw[i] for i in range(len(abw) - 1))
        alle_besser &= besser
        zeilen.append((schl, [s[0] for s in stufen], [e[schl] for _, e in stufen],
                       [100 * a for a in abw], besser))
    text = ("jede Groesse kommt mit jeder Verfeinerung naeher an die Referenz"
            if alle_besser else "mindestens eine Groesse entfernt sich")
    return alle_besser, text, zeilen


def richardson(stufen):
    """Beobachtete Konvergenzordnung und Grenzwert aus DREI Netzstufen.

    Mit zwei Stufen laesst sich nur extrapolieren, WENN man die Ordnung vorgibt -- und
    genau das waere Raten. Mit drei Stufen gleichen Verfeinerungsverhaeltnisses r folgt
    die Ordnung aus den Daten selbst:

        p = ln( |f_4h - f_2h| / |f_2h - f_h| ) / ln r
        f_exakt = f_h + (f_h - f_2h) / (r^p - 1)

    Vorausgesetzt ist, dass die drei Werte monoton und im asymptotischen Bereich liegen;
    das zeigt sich daran, dass der Quotient groesser als eins ist (die Aenderungen werden
    kleiner). Ist er es nicht, wird kein Wert zurueckgegeben statt eines erfundenen.
    """
    import wirbelstrasse as W
    if len(stufen) < 3:
        return None
    (z4, e4), (z2, e2), (z1, e1) = stufen[-3], stufen[-2], stufen[-1]
    r = z1 / z2
    if abs(r - z2 / z4) > 1e-6:
        return None                       # ungleiche Verfeinerungsschritte
    aus = []
    for schl, (lo, hi) in (("cD", W.REFERENZ["cD"]), ("cL", W.REFERENZ["cL"]),
                           ("St", W.REFERENZ["St"])):
        d2, d1 = e2[schl] - e4[schl], e1[schl] - e2[schl]
        if d1 == 0 or (d2 / d1) <= 1.0:
            aus.append((schl, None, None, (lo, hi), "nicht im asymptotischen Bereich"))
            continue
        pord = math.log(abs(d2 / d1)) / math.log(r)
        grenz = e1[schl] + d1 / (r ** pord - 1.0)
        aus.append((schl, pord, grenz, (lo, hi),
                    "im Intervall" if lo <= grenz <= hi else "ausserhalb"))
    return aus


def main():
    art = sys.argv[1] if len(sys.argv) > 1 else "schnell"
    ganz = art == "ganz"
    print("=" * 100)
    print(f" PRUEFPROTOKOLL STROEMUNGSLOESER -- Fassung {FASSUNG}   Umfang: {art}")
    print(" Jede Stufe beurteilt sich selbst; hier wird nur eingesammelt.")
    print("=" * 100)
    ergebnis = []
    for kuerzel, was, befehl, grenze, nur_ganz in STUFEN:
        if nur_ganz and not ganz:
            print(f"  [ueber ] {kuerzel}  {was}\n              nur bei Umfang „ganz\"")
            continue
        t0 = time.time()
        try:
            p = subprocess.run(befehl, cwd=HIER, capture_output=True, text=True,
                               timeout=grenze)
            ok = p.returncode == 0
            letzte = [z for z in p.stdout.strip().splitlines() if z.strip()]
            beleg = letzte[-1].strip() if letzte else "keine Ausgabe"
            if not ok and p.stderr.strip():
                beleg = p.stderr.strip().splitlines()[-1][:120]
        except subprocess.TimeoutExpired:
            ok, beleg = False, f"Zeitgrenze {grenze} s ueberschritten"
        dauer = time.time() - t0
        ergebnis.append((kuerzel, ok))
        print(f"  [{'RICHTIG' if ok else 'FALSCH ':7s}] {kuerzel}  {was}   ({dauer:.0f} s)")
        print(f"              {beleg}")
    # ---- S15: der Trend der Netzreihe ------------------------------------------------
    ok, text, zeilen = netzreihe_urteil()
    if ok is None:
        print(f"  [ueber ] S15  Netzreihe gegen Schaefer-Turek\n              {text}")
    else:
        ergebnis.append(("S15", ok))
        print(f"  [{'RICHTIG' if ok else 'FALSCH ':7s}] S15  Netzreihe laeuft auf die "
              f"Referenz zu\n              {text}")
        for schl, zpd, werte, abw, besser in zeilen:
            w = "  ".join(f"{z:.0f}Z/D {v:.4f} ({a:+.1f}%)"
                          for z, v, a in zip(zpd, werte, abw))
            print(f"                {schl:3s}  {w}   {'besser' if besser else 'SCHLECHTER'}")
    # ---- S16: Richardson-Extrapolation aus drei Netzstufen ---------------------------
    import wirbelstrasse as W          # auch von S13 unten gebraucht
    stufen = []
    for d in sorted(glob.glob(os.path.join(HIER, "ergebnis_*x*.json"))):
        if "zweiwege" in d:
            continue
        e = json.load(open(d))
        if "cD" not in e or any(not math.isfinite(float(e[k])) for k in ("cD", "cL", "St")):
            continue
        stufen.append((e["nx"] / W.LX * W.DZ, e))
    stufen.sort(key=lambda s_: s_[0])
    aus = richardson(stufen)
    if aus is None:
        print("  [ueber ] S16  Richardson-Extrapolation\n"
              "              weniger als drei Netzstufen mit gleichem Verfeinerungsschritt")
    else:
        gut = all(u == "im Intervall" for _, _, _, _, u in aus)
        ergebnis.append(("S16", gut))
        print(f"  [{'RICHTIG' if gut else 'FALSCH ':7s}] S16  der extrapolierte Grenzwert "
              f"liegt im Referenzintervall")
        for schl, pord, grenz, (lo, hi), urteil in aus:
            if pord is None:
                print(f"                {schl:3s}  {urteil}")
            else:
                print(f"                {schl:3s}  beobachtete Ordnung p = {pord:.2f}   "
                      f"Grenzwert {grenz:.4f}   Referenz {lo} … {hi}   {urteil}")
        # DIAGNOSE, nicht Kriterium: Richardson setzt EINEN beherrschenden Fehleranteil
        # voraus. Dann muessten die beobachteten Ordnungen der drei Groessen aehnlich
        # sein, denn sie stammen aus derselben Diskretisierung. Streuen sie stark, ist die
        # Voraussetzung verletzt und der extrapolierte Wert nicht belastbar -- unabhaengig
        # davon, ob er zufaellig im Intervall landet. Das Urteil oben bleibt davon
        # unberuehrt; diese Zeile erklaert es nur.
        ordnungen = [pp for _, pp, _, _, _ in aus if pp is not None]
        if len(ordnungen) >= 2:
            spanne = max(ordnungen) - min(ordnungen)
            print(f"                Voraussetzung: die drei Ordnungen muessten aehnlich "
                  f"sein. Gemessen {min(ordnungen):.2f} … {max(ordnungen):.2f}, "
                  f"Spanne {spanne:.2f}")
            print(f"                -> {'vertraeglich' if spanne < 0.3 else 'UNVERTRAEGLICH: ein einziger beherrschender Fehleranteil liegt nicht vor'}")
    # ---- Erklaerblock zu S13: die Streuung durch die Lage der Treppe -----------------
    # Dies ist KEIN Kriterium und zaehlt nicht mit. Es ist die Messung, die erklaert,
    # warum S13 scheitert: Bei unveraenderter Geometrie und gleicher Aufloesung aendert
    # allein die Lage des Zylinders im Zellraster die Beiwerte. Ein Kriterium daraus zu
    # machen hiesse, sich die Latte nach dem Sprung zurechtzulegen.
    import statistics as _st
    streu = []
    for d in sorted(glob.glob(os.path.join(HIER, "streuung_*x*.json"))):
        streu.append(json.load(open(d)))
    grund = os.path.join(HIER, "ergebnis_440x82.json")
    if streu and os.path.isfile(grund):
        streu.insert(0, json.load(open(grund)))
        print("  [Messung] Erklaerung zu S16: Streuung allein durch die Lage der Treppe")
        print(f"                {len(streu)} Netze, Geometrie unveraendert, "
              f"Zellabstand um weniger als 2,5 % verschieden")
        for schl, name in (("cD", "c_D"), ("cL", "c_L"), ("St", "St ")):
            w = [e[schl] for e in streu]
            sp = max(w) - min(w)
            print(f"                {name}  {min(w):.4f} … {max(w):.4f}   Spanne {sp:.4f} "
                  f"= {100*sp/_st.fmean(w):.1f} % des Mittelwerts")

    # ---- S17: die vierte Referenzgroesse, die Druckdifferenz -------------------------
    dps = []
    for d in sorted(glob.glob(os.path.join(HIER, "ergebnis_dp_*.json"))):
        e = json.load(open(d))
        if math.isfinite(float(e.get("dp", float("nan")))):
            dps.append((e["nx"] / W.LX * W.DZ, e))
    dps.sort(key=lambda s_: s_[0])
    if len(dps) < 2:
        print("  [ueber ] S17  Druckdifferenz dp gegen 2,46 ... 2,50\n"
              "              weniger als zwei Netzstufen gerechnet")
    else:
        lo, hi = 2.46, 2.50
        mitte = 0.5 * (lo + hi)
        abw = [abs(e["dp"] / mitte - 1.0) for _, e in dps]
        besser = all(abw[i + 1] < abw[i] for i in range(len(abw) - 1))
        ergebnis.append(("S17", besser))
        print(f"  [{'RICHTIG' if besser else 'FALSCH ':7s}] S17  Druckdifferenz dp kommt "
              f"der Referenz 2,46 … 2,50 naeher")
        w = "  ".join(f"{z:.0f}Z/D {e['dp']:.4f} ({100*(e['dp']/mitte-1):+.1f}%)"
                      for (z, e), a in zip(dps, abw))
        print(f"                {w}")

    # ---- S18: der UNABHAENGIGE Weg gegen die Referenz --------------------------------
    fem = []
    for d in sorted(glob.glob(os.path.join(HIER, "ergebnis_fem_*.json"))):
        e = json.load(open(d))
        if all(math.isfinite(float(e[k])) for k in ("cD", "cL", "St")):
            fem.append(e)
    if not fem:
        print("  [ueber ] S18  Galerkin-Weg gegen Schaefer-Turek\n"
              "              noch kein Lauf vorhanden")
    else:
        fein = min(fem, key=lambda e: e["dt"])
        drin = all(lo <= fein[k] <= hi for k, (lo, hi) in
                   (("cD", W.REFERENZ["cD"]), ("cL", W.REFERENZ["cL"]),
                    ("St", W.REFERENZ["St"])))
        ergebnis.append(("S18", drin))
        print(f"  [{'RICHTIG' if drin else 'FALSCH ':7s}] S18  der unabhaengige Galerkin-Weg "
              f"trifft ALLE DREI Referenzintervalle")
        print(f"                Netz {fein['knoten']} Ecken, {fein['nzyl']} am Zylinder, "
              f"dt = {fein['dt']:.2e} s, {fein['unbekannte']} Unbekannte")
        for k, nm, (lo, hi) in (("cD", "c_D", W.REFERENZ["cD"]),
                                ("cL", "c_L", W.REFERENZ["cL"]),
                                ("St", "St ", W.REFERENZ["St"])):
            print(f"                {nm}  {fein[k]:.4f}   Referenz {lo} … {hi}   "
                  f"{'im Intervall' if lo <= fein[k] <= hi else 'DANEBEN'}")
        # Zeitordnung aus den vorhandenen Zeitstufen
        gleich = [e for e in fem if e["knoten"] == fein["knoten"]]
        gleich.sort(key=lambda e: -e["dt"])
        if len(gleich) >= 3:
            for k, nm in (("cD", "c_D"), ("cL", "c_L")):
                a, b, c = (e[k] for e in gleich[:3])
                d1, d2 = b - a, c - b
                if abs(d2) > 1e-12 and d1 * d2 > 0:
                    pp = math.log(abs(d1/d2)) / math.log(2)
                    print(f"                beobachtete Zeitordnung {nm}: p = {pp:.2f} "
                          f"(Euler erster Ordnung verlangt 1)")

    # ---- S19: die beiden Browserkerne gegen die Python-Module ------------------------
    # Ab 01.11.2026 rechnet die Seite BEIDE Wege selbst. Damit gehoeren ihre Kerne in
    # dieses Protokoll: Ein Weg, der nur im Browser laeuft und hier nicht vorkommt,
    # waere ein zweiter, ungepruefter Stand desselben Verfahrens.
    labor = os.path.join(os.path.dirname(HIER), "labor")
    for kuerzel, probe, was in (("S19a", "pruefe_kern.js", "Netzwerkkern"),
                                ("S19b", "pruefe_kern_fem.js", "Galerkin-Kern")):
        datei = os.path.join(labor, probe)
        if not os.path.isfile(datei):
            print(f"  [ueber ] {kuerzel}  {was} im Browser\n"
                  f"              {probe} nicht vorhanden")
            continue
        r = subprocess.run(["node", probe], cwd=labor, capture_output=True,
                           text=True, timeout=1800)
        aus = (r.stdout or "") + (r.stderr or "")
        m = re.search(r"== (\d+)/(\d+)", aus)
        if not m:
            ergebnis.append((kuerzel, False))
            print(f"  [FALSCH ] {kuerzel}  {was}: die Probe lieferte keine Zahl")
            continue
        gut, alle = int(m.group(1)), int(m.group(2))
        ergebnis.append((kuerzel, gut == alle))
        print(f"  [{'RICHTIG' if gut == alle else 'FALSCH ':7s}] {kuerzel}  der {was} im "
              f"Browser rechnet dasselbe wie die Python-Module")
        print(f"                {gut} von {alle} Posten, knotenweise verglichen "
              f"({probe})")

    n = sum(1 for _, o in ergebnis if o)
    print("\n" + "=" * 100)
    print(f" {n} von {len(ergebnis)} Stufen RICHTIG")
    print("=" * 100)
    return 0 if n == len(ergebnis) else 1


if __name__ == "__main__":
    raise SystemExit(main())
