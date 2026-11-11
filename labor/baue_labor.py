#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Baut CFD_Labor_<Fassung>.html -- eine einzelne Datei, die ohne Netz laeuft.

Die Fassungsnummer steht NUR hier (stehende Regel). Eingebaut werden:

    kern_cfd.js      der Rechenkern des NETZWERKWEGS, gegen die Python-Module geprueft
    kern_fem.js      der Rechenkern des GALERKIN-WEGS, ebenso geprueft
    netz_browser.js  der Vernetzer (Randpunkte, Streuung, Delaunay)
    triangle/        Triangle von Shewchuk als WebAssembly, als base64 eingebettet
    labor.js         die Bedienung
    vorlage.html     Geruest und Gestaltung
    die Manuskripte  ueber pandoc nach HTML mit MathML (kein Netz, kein CDN)
    die Bilder       als data:-URI
    die Zeitdaten    aus surrogat_daten.npz, int16 mit eigener Skala je Bild

Kein Verweis nach aussen: Die Seite laedt nichts nach. Das ist Bedingung, nicht Zier --
sie soll auf einem Rechner ohne Netz genauso laufen.
"""
FASSUNG = "1.6"
STAND = "11.11.2026"

import base64
import io
import glob
import json
import os
import re
import subprocess
import sys

import numpy as np

HIER = os.path.dirname(os.path.abspath(__file__))
WURZEL = os.path.dirname(HIER)
STROEM = os.path.join(WURZEL, "stroemung")


def entkasten(t):
    """\\boxed{...} und \\Bigl/\\Bigr kennt der MathML-Schreiber nicht."""
    while r"\boxed{" in t:
        i = t.index(r"\boxed{"); j = i + len(r"\boxed{"); tiefe = 1
        while j < len(t) and tiefe:
            if t[j] == "{": tiefe += 1
            elif t[j] == "}": tiefe -= 1
            j += 1
        t = t[:i] + t[i + len(r"\boxed{"):j - 1] + t[j:]
    return (t.replace(r"\Bigl[", r"\left[").replace(r"\Bigr]", r"\right]")
             .replace(r"\Bigl(", r"\left(").replace(r"\Bigr)", r"\right)")
             .replace(r"\bigl[", r"\left[").replace(r"\bigr]", r"\right]"))


def kernfassung(text):
    """Die Fassung des Rechenkerns steht in SEINEM Kopf, nicht hier.

    Befund 31.10.2026: Sie war im Bauskript fest eingetragen und blieb auf 1.0 stehen,
    als der Kern auf 1.1 ging (Befund B11). Eine Zahl an zwei Stellen ist eine Zahl zu
    viel -- genau das verbietet die stehende Regel zu Fassungsnummern.
    """
    m = re.search(r"Fassung (\d+\.\d+)", text[:400])
    if not m:
        raise SystemExit("FEHLER: keine Fassung im Kopf von kern_cfd.js gefunden")
    return m.group(1)


def md_zu_html(pfade, titel_weg=True):
    text = "\n\n".join(open(os.path.join(WURZEL, p), encoding="utf-8").read() for p in pfade)
    text = re.sub(r"^---\n.*?\n---\n", "", text, flags=re.S)          # YAML-Kopf weg
    text = entkasten(text)
    r = subprocess.run(["pandoc", "-f", "markdown", "-t", "html5", "--mathml",
                        "--no-highlight"], input=text, capture_output=True,
                       text=True, timeout=600)
    if r.returncode:
        raise SystemExit("pandoc: " + r.stderr[:400])
    return r.stdout


def bild_uri(pfad):
    with open(pfad, "rb") as f:
        return "data:image/png;base64," + base64.b64encode(f.read()).decode()


def zeitdaten():
    """Momentaufnahmen und Surrogat in die Form bringen, die der Browser liest."""
    q = os.path.join(HIER, "surrogat_daten.npz")
    if not os.path.isfile(q):
        return None
    d = np.load(q)
    XA, XB = d["XA"], d["XB"]                 # (Punkte, Zeiten)
    nxv, nyv = (int(x) for x in d["form"])
    n = XA.shape[1]
    # Gespeichert werden Weg A UND die Differenz, jeweils mit EIGENER Skala je Bild.
    # Zwei getrennt quantisierte Felder haetten den Unterschied (1e-14) zerstoert;
    # aus A und der Differenz laesst sich B dagegen in voller Genauigkeit bilden.
    Diff = XA - XB
    sA = np.abs(XA).max(axis=0); sA[sA == 0] = 1.0
    sD = np.abs(Diff).max(axis=0); sD[sD == 0] = 1e-300
    qA = np.rint((XA / sA) * 32767).astype("<i2").T.copy()
    qD = np.rint((Diff / sD) * 32767).astype("<i2").T.copy()
    # Das Zeitgesetz kommt als HARMONISCHE Anpassung, nicht als DMD-Operator: Auf dem
    # Grenzzyklus ist sie um den Faktor 39 genauer (0,187 % gegen 7,3 %), weil sie
    # exakt periodisch ist und sich deshalb keine Phase davonstehlen kann. Beide
    # Zahlen stehen im Bericht und auf der Seite.
    h = np.load(os.path.join(HIER, "surrogat_harmonisch.npz"))
    sur = {
        "r": int(d["U"].shape[1]), "K": int(h["K"][0]), "f": float(h["f"][0]),
        "t0": float(d["zeiten"][0]),
        "mittel": base64.b64encode(d["mittel"].ravel().astype("<f4").tobytes()).decode(),
        "moden": base64.b64encode(np.ascontiguousarray(d["U"].T).astype("<f4").tobytes()).decode(),
        "BA": h["BA"].ravel().tolist(), "BB": h["BB"].ravel().tolist(),
    }
    bericht = json.load(open(os.path.join(HIER, "surrogat_bericht.json")))
    bericht.update(json.load(open(os.path.join(HIER, "harmonisch_bericht.json"))))

    # ---- der UNABHAENGIGE Weg, ueber die PHASE zugeordnet ---------------------------
    # Netzwerk- und Galerkin-Rechnung sind nicht phasengleich: verschiedene Verfahren,
    # verschiedene Anlaufvorgaenge. Sie bei derselben absoluten Zeit nebeneinander zu
    # legen waere sinnlos -- die Wirbel staenden an verschiedenen Stellen.
    #
    # Die Phase kommt aus einer HARMONISCHEN ANPASSUNG an c_L, nicht aus dessen
    # Ableitung. Fassung 1.2, Befund 01.11.2026: Zuerst stand hier der Zustand
    # (c_L, dc_L/dt), beides auf das eigene Maximum bezogen. Das rohe c_L ist aber von
    # Schritt zu Schritt verrauscht; seine Ableitung wird vom Rauschen beherrscht, und
    # die Normierung macht das Rauschen zum Massstab. Gemessen: mittlerer
    # Zuordnungsabstand 0,31 statt rund 0,01, siebzig Prozent der Bilder schlecht
    # zugeordnet. Mit der angepassten Grundschwingung ist die Phase glatt und eindeutig.
    def phase(t, cl, T0):
        """Phase der Grundschwingung: c_L ~ A cos(w t + p0), angepasst ueber das Fenster."""
        bestes, bw, bp = np.inf, 2 * np.pi / T0, 0.0
        for w in np.linspace(2*np.pi/T0*0.9, 2*np.pi/T0*1.1, 401):
            M = np.c_[np.cos(w*t), np.sin(w*t), np.ones_like(t)]
            k, *_ = np.linalg.lstsq(M, cl, rcond=None)
            r = float(((M @ k - cl) ** 2).sum())
            if r < bestes:
                bestes, bw, bp = r, w, np.arctan2(-k[1], k[0])
        return (bw * t + bp) % (2 * np.pi), bw

    gal = os.path.join(HIER, "..", "stroemung", "fem_schnapp.npz")
    netz_fein = os.path.join(HIER, "..", "stroemung", "netz_schnapp_880.npz")
    galerkin = None
    if os.path.isfile(gal):
        G = np.load(gal)
        BG, clG, tG = G["bilder"], G["cl"], G["zeiten"]
        # Das Netzwerkfeld fuer DIESEN Vergleich kommt vom FEINSTEN Netz, sonst
        # vergleicht man eine grobe mit einer feinen Loesung und sieht die Aufloesung,
        # nicht das Verfahren (Befund 01.11.2026: |omega| 116 gegen 268).
        if os.path.isfile(netz_fein):
            Nf = np.load(netz_fein)
            BN, tNf, vN = Nf["bilder"], Nf["zeiten"], Nf["verlauf"]
            dtNf = float(Nf["dt"])
            clNf = vN[np.clip(np.round(tNf/dtNf).astype(int), 0, len(vN)-1), 2]
            TN = 0.1 / float(Nf["St"]) if float(Nf["St"]) > 0 else 0.34
            pN, _ = phase(tNf, clNf, TN)
            XN = BN.T.astype(float)
            netzquelle = f"880 x 164, St = {float(Nf['St']):.4f}"
        else:
            XN, pN, netzquelle = None, None, None
        if XN is not None:
            TG = 0.1 / float(G["St"])
            pG, _ = phase(tG, clG, TG)
            wahl = np.array([np.argmin(np.abs(((pG - q + np.pi) % (2*np.pi)) - np.pi))
                             for q in pN])
            passung = float(np.abs(((pG[wahl] - pN + np.pi) % (2*np.pi)) - np.pi).max())
            XG = BG[wahl].T.astype(float)
            sN = np.abs(XN).max(axis=0); sN[sN == 0] = 1.0
            sG = np.abs(XG).max(axis=0); sG[sG == 0] = 1.0
            qN = np.rint((XN / sN) * 32767).astype("<i2").T.copy()
            qG = np.rint((XG / sG) * 32767).astype("<i2").T.copy()
            # Punktweiser Vergleich beziffern -- Integralgroessen und Felder
            # konvergieren verschieden schnell, und das gehoert gesagt.
            Df = XN - XG
            xs_ = (np.arange(nxv) + 0.5) * 2.2 / nxv
            ys_ = (np.arange(nyv) + 0.5) * 0.41 / nyv
            Xg, Yg = np.meshgrid(xs_, ys_, indexing="ij")
            rr = np.hypot(Xg - 0.2, Yg - 0.2).ravel()
            eff = lambda a: float(np.sqrt((a ** 2).mean()))
            nachl = Xg.ravel() > 0.3
            nahe = rr < 0.09
            kmax = int(np.abs(Df).max(axis=1).argmax())
            jmax = int(np.abs(Df[:, :]).max(axis=0).argmax())
            masse = {"gesamt": 100 * eff(Df) / eff(XN),
                     "nachlauf": 100 * eff(Df[nachl]) / eff(XN[nachl]),
                     "nahe": 100 * eff(Df[nahe]) / eff(XN[nahe]),
                     "ort_mm": float(rr[jmax] * 1000)}
            galerkin = {"datN": base64.b64encode(qN.tobytes()).decode(), "skalaN": sN.tolist(),
                        "masse": masse,
                        "dat": base64.b64encode(qG.tobytes()).decode(), "skala": sG.tolist(),
                        "frames": int(XN.shape[1]), "t": tNf.tolist(),
                        "cD": float(G["cD"]), "cL": float(G["cL"]), "St": float(G["St"]),
                        "netzquelle": netzquelle,
                        "passung_grad": float(np.degrees(passung))}

    return {"nx": nxv, "ny": nyv, "frames": n, "t": d["zeiten"].tolist(),
            "datA": base64.b64encode(qA.tobytes()).decode(), "skalaA": sA.tolist(),
            "datD": base64.b64encode(qD.tobytes()).decode(), "skalaD": sD.tolist(),
            "surrogat": sur, "bericht": bericht, "galerkin": galerkin}


def wirbelseite():
    """Der Reiter „Wirbelstrasse": das Bild und die gemessenen Zahlen."""
    sch = os.path.join(STROEM, "Wirbelstrasse_Schaubild_1.0.png")
    feld = os.path.join(STROEM, "Wirbelstrasse_Feld_1.0.png")
    erg = {}
    # Nur die Netzreihe des NETZWERKWEGS: ergebnis_<nx>x<ny>.json und
    # ergebnis_dp_<nx>x<ny>.json. Die Ergebnisse des Galerkin-Wegs heissen
    # ergebnis_fem_*, die des Lattice-Boltzmann-Moduls ergebnis_lbm_* -- sie haben
    # andere Felder und gehoeren nicht in diese Reihe. Das Pruefprotokoll filtert
    # ueber dasselbe Muster.
    import fnmatch
    for f in sorted(os.listdir(STROEM)):
        if fnmatch.fnmatch(f, "ergebnis_*x*.json") and "zweiwege" not in f:
            erg[f] = json.load(open(os.path.join(STROEM, f)))
    reihe, dps = [], {}
    for f, e in erg.items():
        if "dp" in e:
            dps[e["nx"]] = e["dp"]
        elif "cD" in e and all(np.isfinite(e[k]) for k in ("cD", "cL", "St")):
            reihe.append(e)
    reihe.sort(key=lambda e: e["nx"])
    zpd = lambda e: e["nx"] / 2.2 * 0.1
    z = ['<div class="bildkasten"><img src="%s" alt="Kármánsche Wirbelstraße">' % bild_uri(sch),
         '<figcaption>Gerechnet mit <code>wirbelstrasse.py</code> auf 880 × 164 Zellen, '
         '36 000 Zeitschritten. Oben das Wirbelstärkefeld, unten links der Ablösetakt, '
         'unten rechts die Netzreihe.</figcaption></div>']
    z.append('<div class="kasten" style="margin-top:14px"><h2>Gemessen gegen Schäfer-Turek 2D-2</h2>')
    z.append('<table class="tafel"><tr><th>Zellen über D</th><th>Netz</th>'
             '<th>c_D max</th><th>c_L max</th><th>St</th><th>Δp</th></tr>')
    for e in reihe:
        dp = dps.get(e["nx"])
        z.append(f'<tr><td>{zpd(e):.0f}</td><td>{e["nx"]} × {e["ny"]}</td>'
                 f'<td class="zahl">{e["cD"]:.4f}</td><td class="zahl">{e["cL"]:.4f}</td>'
                 f'<td class="zahl">{e["St"]:.4f}</td>'
                 f'<td class="zahl">{dp:.4f}</td></tr>'.replace(".", ",") if dp else
                 f'<tr><td>{zpd(e):.0f}</td><td>{e["nx"]} × {e["ny"]}</td>'
                 f'<td class="zahl">{e["cD"]:.4f}</td><td class="zahl">{e["cL"]:.4f}</td>'
                 f'<td class="zahl">{e["St"]:.4f}</td><td>—</td></tr>')
    z.append('<tr><td><b>Referenz</b></td><td>—</td><td><b>3,22 … 3,24</b></td>'
             '<td><b>0,99 … 1,01</b></td><td><b>0,295 … 0,305</b></td>'
             '<td><b>2,46 … 2,50</b></td></tr></table>')
    z.append('<div class="merk">Auf dem feinsten Netz liegt die <b>Strouhalzahl im '
             'Referenzintervall</b>. Jede der vier Größen kommt mit jeder Verfeinerung '
             'näher — das ist der Nachweis, nicht ein einzelner Lauf, der zufällig trifft.</div>')
    z.append('</div>')
    # Treppenstreuung
    streu = []
    for f in sorted(os.listdir(STROEM)):
        if f.startswith("streuung_") and f.endswith(".json"):
            streu.append(json.load(open(os.path.join(STROEM, f))))
    if streu:
        g440 = erg.get("ergebnis_440x82.json")
        if g440:
            streu.insert(0, dict(g440, dy=0.41 / g440["ny"], lage_x=0.0, lage_y=0.0))
        z.append('<div class="kasten" style="margin-top:14px">'
                 '<h2>Woher der Rest kommt — gemessen, nicht vermutet</h2>'
                 '<p style="font:13.5px system-ui;margin:0 0 8px">Vier Netze, <b>unveränderte '
                 'Geometrie</b>, Zellabstand um weniger als 2,5 % verschieden. Nur die Lage '
                 'des Zylinders im Zellraster ist anders.</p>'
                 '<table class="tafel"><tr><th>Netz</th><th>Δy</th><th>Mittelpunkt im Raster</th>'
                 '<th>c_D max</th><th>c_L max</th><th>St</th></tr>')
        for e in streu:
            z.append(f'<tr><td>{e["nx"]} × {e["ny"]}</td><td class="zahl">{e["dy"]:.6f}</td>'
                     f'<td class="zahl">x {e["lage_x"]:.2f} Δx, y {e["lage_y"]:.2f} Δy</td>'
                     f'<td class="zahl">{e["cD"]:.4f}</td><td class="zahl">{e["cL"]:.4f}</td>'
                     f'<td class="zahl">{e["St"]:.4f}</td></tr>')
        z.append('</table>')
        for s_, nm in (("cD", "c_D"), ("cL", "c_L"), ("St", "St")):
            w = [e[s_] for e in streu]
            sp = max(w) - min(w)
            z.append(f'<p style="margin:6px 0 0;font:13px system-ui"><b>{nm}</b>: '
                     f'{min(w):.4f} … {max(w):.4f}, Spanne {sp:.4f} '
                     f'({100*sp/np.mean(w):.1f} % vom Mittelwert)</p>')
        z.append('<div class="merk warn">Für c_D und c_L ist diese Streuung <b>so groß wie der '
                 'Fehler selbst</b>. Auf diesem Netz bestimmt nicht die Auflösung den Fehler, '
                 'sondern wo die Treppe hinfällt. Deshalb trägt auch keine '
                 'Richardson-Extrapolation: Sie setzt einen festen Vorfaktor voraus, und den '
                 'gibt es hier nicht.</div></div>')
    # ---- der unabhaengige Weg ------------------------------------------------------
    fem = []
    for f in sorted(os.listdir(STROEM)):
        if f.startswith("ergebnis_fem_") and f.endswith(".json"):
            fem.append(json.load(open(os.path.join(STROEM, f))))
    if fem:
        fem.sort(key=lambda e: (-e["dt"], e["knoten"]))
        z.append('<div class="kasten" style="margin-top:14px">'
                 '<h2>Der unabhängige zweite Weg — Galerkin auf körperangepasstem Netz</h2>'
                 '<p style="font:13.5px system-ui;margin:0 0 8px">Andere Herleitung (schwache '
                 'Form statt Flussbilanz), andere Unbekannte (Knotenwerte eines Polynomansatzes '
                 'statt Zellmittel), <b>keine Treppe</b> — der Zylinderrand liegt auf Knoten. '
                 'Zwei Wege, die verschieden scheitern können.</p>'
                 '<table class="tafel"><tr><th>Netz / Zeitschritt</th><th>Knoten</th>'
                 '<th>am Zylinder</th><th>c_D max</th><th>c_L max</th><th>St</th></tr>')
        for e in fem:
            drin = all(lo <= e[k] <= hi for k, (lo, hi) in
                       (("cD", (3.22, 3.24)), ("cL", (0.99, 1.01)), ("St", (0.295, 0.305))))
            fett = ' style="font-weight:600"' if drin else ''
            z.append(f'<tr{fett}><td>h={e["h_zyl"]}, dt={e["dt"]:.2e}</td>'
                     f'<td class="zahl">{e["knoten"]}</td><td class="zahl">{e["nzyl"]}</td>'
                     f'<td class="zahl">{e["cD"]:.4f}</td><td class="zahl">{e["cL"]:.4f}</td>'
                     f'<td class="zahl">{e["St"]:.4f}</td></tr>'.replace(".", ","))
        z.append('<tr><td><b>Referenz</b></td><td>—</td><td>—</td><td><b>3,22 … 3,24</b></td>'
                 '<td><b>0,99 … 1,01</b></td><td><b>0,295 … 0,305</b></td></tr></table>')
        bester = min(fem, key=lambda e: e["dt"])
        z.append('<div class="merk">Auf dem feinsten Zeitschritt liegen <b>alle drei Größen '
                 f'im Referenzintervall</b> — c_D {bester["cD"]:.4f}, c_L {bester["cL"]:.4f}, '
                 f'St {bester["St"]:.4f} — <b>direkt gerechnet, ohne Extrapolation</b>. '
                 'Die beobachtete Zeitordnung ist 1,00: Das Zeitschema ist erster Ordnung, '
                 'die Messung trifft genau das. Zum Vergleich blieb der Netzwerkweg mit '
                 '144 320 Zellen bei +1,8 % stehen — der Rest dort kam von der '
                 'Treppenkontur, nicht vom Ansatz.</div></div>'.replace(".", ","))
    z.append('<div class="bildkasten" style="margin-top:14px"><img src="%s" alt="Wirbelstärkefeld">'
             '<figcaption>Dasselbe Feld ohne Beiwerk.</figcaption></div>' % bild_uri(feld))
    return "\n".join(z)


def aufgabe(roh):
    """Die Aufgabe eines Programms aus SEINEM EIGENEN Kopf lesen, nicht aus einer Liste.

    Jede Datei beginnt mit einer Zeile, die sagt, was sie tut -- der Docstring bei Python,
    der Kopfkommentar bei JavaScript. Von dort kommt der Text. Eine zweite Liste neben den
    Programmen zu fuehren hiesse, zwei Staende zu pflegen; der eine veraltet dann.
    """
    for z in roh.split("\n")[:14]:
        s = z.strip()
        for weg in ('r"""', '"""', "r, ", "/* ", "// ", "# ", "* "):
            if s.startswith(weg):
                s = s[len(weg):].strip()
        if (len(s) > 18 and not s.startswith(("#!", "# -*-", "import", "from", "FASSUNG",
                                              "const", "(function", "'use strict'", "*", "/*"))
                and "coding:" not in s):
            # Die JavaScript-Koepfe lauten „kern_fem.js -- der unabhaengige Weg...";
            # dort steht die Aufgabe RECHTS vom Trenner, bei Python links. Entschieden
            # wird nicht nach Dateityp, sondern danach, welche Haelfte ein Dateiname ist.
            for trenn in (" -- ", " \u2014 ", ": "):
                if trenn in s:
                    links, rechts = s.split(trenn, 1)
                    s = rechts if re.fullmatch(r"\S+\.(js|py)", links.strip()) else s
                    break
            for trenn in (" -- ", " \u2014 ", ". "):
                if trenn in s:
                    s = s.split(trenn)[0]
                    break
            return s.rstrip(".,;:").strip()
    return ""


def quelltexte():
    """Alle Programme als lesbaren Text in die Seite legen.

    Auf Freds Einwand vom 01.11.2026: „Alle Pruefprotokolle verweisen darauf,
    nachrechnen kann ich sie ohne die Quellen nicht." Eine Seite, die auf Programme
    verweist, die ihr nicht beiliegen, ist fuer den Leser eine Sackgasse -- er kann
    glauben oder bleiben lassen, aber nicht nachrechnen. Deshalb liegen sie bei:
    jedes Programm mit seiner Groesse und seiner Pruefsumme, damit auch nachvollziehbar
    ist, dass es dasselbe ist wie das gelieferte.
    """
    import hashlib
    aus = []
    for ordner, muster in (("stroemung", ("*.py",)), ("labor", ("*.py", "*.js"))):
        wo = os.path.join(WURZEL, ordner)
        dateien = []
        for m in muster:
            dateien += sorted(glob.glob(os.path.join(wo, m)))
        for d in dateien:
            roh = open(d, encoding="utf-8").read()
            aus.append({
                "name": f"{ordner}/{os.path.basename(d)}",
                "aufgabe": aufgabe(roh),
                "sha256": hashlib.sha256(roh.encode("utf-8")).hexdigest(),
                "bytes": len(roh.encode("utf-8")),
                "zeilen": roh.count("\n") + 1,
                "text": roh,
            })
    return aus


def main():
    vorlage = open(os.path.join(HIER, "vorlage.html"), encoding="utf-8").read()
    kern = open(os.path.join(HIER, "kern_cfd.js"), encoding="utf-8").read()
    kernfem = open(os.path.join(HIER, "kern_fem.js"), encoding="utf-8").read()
    netzjs = open(os.path.join(HIER, "netz_browser.js"), encoding="utf-8").read()
    labor = open(os.path.join(HIER, "labor.js"), encoding="utf-8").read()
    # Triangle (Shewchuk) als WebAssembly. Die Emscripten-Huelle holt das wasm sonst
    # ueber das Netz nach -- in einer Datei, die ohne Netz laufen soll, geht das nicht.
    # Deshalb kommen die Bytes als base64 mit und werden ueber Module.wasmBinary
    # uebergeben; dann wird nichts geladen.
    tri_js = open(os.path.join(HIER, "triangle", "triangle.out.js"), encoding="utf-8").read()
    tri_huelle = open(os.path.join(HIER, "triangle", "index.js"), encoding="utf-8").read()
    with open(os.path.join(HIER, "triangle", "triangle.out.wasm"), "rb") as f:
        tri_wasm = base64.b64encode(f.read()).decode("ascii")
    # Die Huelle ist ein CommonJS-Modul. Drei Dinge muessen weg bzw. dazu:
    #   1. ihr `require` -- sonst steht createModule zweimal da, und das ist ein
    #      PARSE-Fehler: der ganze Skriptblock laeuft dann gar nicht erst an
    #      (Befund 01.11.2026, gefunden von der Browserprobe)
    #   2. module.exports -> ein Namensraum am Fenster
    #   3. init() muss das wasm als Bytes bekommen, sonst holt Emscripten es ueber
    #      das Netz nach -- in einer Datei, die ohne Netz laufen soll, geht das nicht
    tri_huelle = re.sub(r"^const createModule = require\([^)]*\);\s*$", "",
                        tri_huelle, flags=re.M)
    tri_huelle = tri_huelle.replace("module.exports = {", "window.TRIANGLE = {")
    tri_huelle = re.sub(r"createModule\(\{\s*locateFile:[^}]*\}\)",
                        "createModule({ wasmBinary: window.__triWasm() })", tri_huelle)
    if "wasmBinary" not in tri_huelle:
        raise SystemExit("FEHLER: init() der Triangle-Huelle nicht gefunden -- "
                         "ohne wasmBinary wuerde die Seite nachladen wollen")
    # Beide Teile bekommen einen EIGENEN Geltungsbereich. Das Emscripten-Erzeugnis
    # deklariert `var Module`, die Huelle `let Module` -- nebeneinander im selben
    # Geltungsbereich ist das wieder ein Parse-Fehler, und wieder laeuft dann der ganze
    # Block nicht. Nach aussen sichtbar bleibt nur window.TRIANGLE.
    # Das Erzeugnis von Emscripten weist seine Fabrik nur zu, wenn BEIDE da sind:
    # `module` UND `exports`. Mit nur `module` lief es durch, ohne etwas zu setzen --
    # window.TRIANGLE war dann ein Objekt, dessen init() ins Leere griff.
    tri = ("window.TRIANGLE_WASM='" + tri_wasm + "';\n"
           "window.__triWasm=function(){var s=atob(window.TRIANGLE_WASM),"
           "a=new Uint8Array(s.length);for(var i=0;i<s.length;i++)a[i]=s.charCodeAt(i);"
           "return a;};\n"
           "(function(){var module={exports:{}},exports=module.exports;\n" + tri_js + "\n"
           "window.__triFactory=module.exports;})();\n"
           "(function(){var createModule=window.__triFactory;\n" + tri_huelle + "\n})();")
    # Wie viele Posten die Kernproben haben -- aus ihren Berichten, nicht von Hand.
    proben = 0
    for name in ("kernprobe.txt", "kernprobe_fem.txt"):
        bericht = os.path.join(HIER, name)
        if os.path.isfile(bericht):
            m = re.search(r"== (\d+)/(\d+)", open(bericht, encoding="utf-8").read())
            if m:
                proben += int(m.group(2))
    qt = quelltexte()
    zd = zeitdaten()
    daten = ("window.LABORDATEN = "
             + json.dumps({"zeit": zd, "quellen": qt}, separators=(",", ":")) + ";")
    ersetzungen = {
        "@@FASSUNG@@": FASSUNG,
        "@@STAND@@": STAND,
        "@@KERNFASSUNG@@": kernfassung(kern),
        "@@KERNPROBEN@@": str(proben),
        "@@WIRBEL@@": wirbelseite(),
        "@@KURZFASSUNG@@": md_zu_html(["KURZFASSUNG.md"]),
        "@@VERFAHREN@@": md_zu_html(["METHODE_UND_NACHWEIS.md"]),
        "@@BEFUNDE@@": md_zu_html(["BEFUNDE.md"]),
        "@@ANLEITUNG@@": md_zu_html(["LIESMICH.md", "PRUEFPLAN.md", "FUER_FRED.md"]),
        "@@LIZENZ@@": md_zu_html(["LIZENZ.md"]),
        "@@KERN@@": kern,
        "@@KERNFEM@@": kernfem,
        "@@NETZJS@@": netzjs,
        "@@TRIANGLE@@": tri,
        "@@DATEN@@": daten,
        "@@LABOR@@": labor,
    }
    # Befund 01.11.2026: Hier liefen die Ersetzungen NACHEINANDER, jede ueber das
    # Ergebnis der vorigen. Seit die Quelltexte mitgeliefert werden, steht
    # `baue_labor.py` selbst in den Daten -- und darin die Marken @@KERN@@, @@LABOR@@,
    # @@DATEN@@. Der naechste Durchgang ersetzte sie MITTEN IM DATENSATZ und zerschlug
    # ihn: „Invalid or unexpected token", der ganze LABORDATEN-Block fiel aus, und mit
    # ihm die Zeitdaten und die Quelltexte. Die Seite sah dabei fast normal aus.
    #
    # Richtig ist EIN Durchgang: Jede Marke wird genau einmal angesehen, und was
    # eingesetzt wird, wird nicht noch einmal durchsucht.
    html = re.sub(r"@@[A-Z]+@@",
                  lambda m: ersetzungen.get(m.group(0), m.group(0)),
                  vorlage)
    ziel = os.path.join(WURZEL, f"CFD_Labor_{FASSUNG}.html")
    open(ziel, "w", encoding="utf-8").write(html)
    print(f"  {os.path.basename(ziel)} gebaut: {len(html)/1024/1024:.2f} MB")
    print(f"  Zeitdaten: {'%d Bilder zu %dx%d Punkten' % (zd['frames'], zd['nx'], zd['ny']) if zd else 'KEINE'}")
    # Befund 01.11.2026: Diese Probe suchte Platzhalter in der FERTIGEN Seite. Sobald ein
    # Manuskript den Namen eines Platzhalters im Fliesstext nennt -- Befund B17 tut das --,
    # schlaegt sie Alarm fuer einen Reiter, der laengst gefuellt ist. Geprueft wird jetzt,
    # was richtig ist: ob die VORLAGE einen Platzhalter fuehrt, fuer den es keine
    # Ersetzung gibt. Das ist der Fehler, den die Probe finden soll.
    offen = sorted(set(re.findall(r"@@[A-Z]+@@", vorlage)) - set(ersetzungen))
    if offen:
        raise SystemExit("FEHLER: Vorlage hat Platzhalter ohne Ersetzung: " + ", ".join(offen))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
