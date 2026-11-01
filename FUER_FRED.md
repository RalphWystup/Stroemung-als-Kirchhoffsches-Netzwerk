# Kurz für Fred (und seine KI) — Fassung 1.5, 01.11.2026

**Was wir bauen:** Eine Strömungsrechnung (CFD) als **Kirchhoffsches Netzwerk** — nicht als
Analogie, sondern weil die Finite-Volumen-Methode schon eines ist: Zelle = Knoten,
Fläche = Zweig, und die Kontinuitätsgleichung $\oint \mathbf{u}\cdot\mathbf{n}\,\mathrm{d}A = 0$ **ist**
 wörtlich die Knotenregel.
Zähigkeit = Widerstand, Trägheit = Kondensator, **Konvektion = gesteuerte Quelle**,
Druckkorrektur = Poisson-Netz.

**Wozu:** Wir brauchen für die Strömung einen **zweiten unabhängigen Weg** — so wie es ihn
für die Feldrechnung schon gibt (Netzwerk nach Kirchhoff *und* FEM nach Ritz, im
`Feld_Labor` nebeneinander mit Sprachprobe). Demonstriert an einem plakativen Bild, aber
**bewiesen, nicht geschönt**.

**Das Beweisstück:** Kármánsche Wirbelstraße am Zylinder, Re = 100. Bei Re = 100 ist der
Nachlauf **echt 2D** (3D erst ab Re ≈ 190) — also keine Vereinfachung. Dafür gibt es den
**Schäfer-Turek-Benchmark (DFG 1996)** mit eng eingegrenzten Referenzwerten:
c_D = 3,22…3,24 · c_L = 0,99…1,01 · St = 0,295…0,305 · Δp = 2,46…2,50. Dimensionslos, von
mehreren Gruppen bestätigt — nicht anpassbar.

## Stand der Beweiskette

| Stufe | Referenz | Stand |
|:--|:--|:--|
| 1 Zweigstempel | geschlossene Lösung der 1D-Konvektion-Diffusion | **bestanden**, Fehler ≤ 3,4·10⁻¹⁵ bei Pe bis 10⁶ |
| 2 Kanalströmung | Hagen-Poiseuille (exakte NS-Lösung) | **bestanden**, $u_{\max}$ exakt auf 2,2·10⁻¹⁶, Ordnung 3× exakt 2,00 |
| 3 Wirbelstraße | Schäfer-Turek 2D-2 | Netzreihe mit drei Stufen, siehe unten |
| 4 Das Bild | — | **gerechnet**, 880 × 164, `wirbelstrasse_schaubild.png` |

Zu Stufe 1, falls deine KI es einordnen will: Der Stempel ist **Scharfetter-Gummel** aus
der Halbleitersimulation (1969). Er ist auf gleichmäßigem Gitter knotenweise exakt, bei
*jeder* Péclet-Zahl — die zentrale Differenz liegt bei Pe_h = 50 um Faktor 1000 daneben,
der Aufwind um 13 bis 16 Prozent bei Pe_h ≈ 1…5.

### Die Netzreihe — Abweichung von der Mitte des Referenzintervalls

| Zellen über D | c_D | c_L | St | Δp |
|--:|--:|--:|--:|--:|
| 10 | +16,1 % | −39,0 % | −15,8 % | −2,9 % |
| 20 | +7,0 % | −11,1 % | −6,4 % | −1,8 % |
| **40** | **+1,8 %** | **−2,6 %** | **−1,3 %** | **−1,0 %** |

Auf dem feinsten Netz liegt **St = 0,2962 im Referenzintervall** 0,295 … 0,305. Jede der
vier Größen kommt mit jeder Verfeinerung näher. Das ist der Nachweis — nicht ein einzelner
Lauf, der zufällig trifft.

### Was NICHT aufgeht, und warum

Die Richardson-Extrapolation aus den drei Stufen läuft über die Referenz hinaus
(c_D → 3,08 statt 3,22…3,24). Das Kriterium ist damit **nicht erfüllt**, und wir lockern es
nicht nachträglich. Stattdessen haben wir gemessen, woran es liegt:

Vier Netze, **gleiche Geometrie**, Zellabstand um weniger als 2,5 % verschieden — nur die
Lage des Zylinders im Zellraster ist anders. Allein dadurch streut

- c_D um **4,7 %** (Abstand zur Referenz auf diesem Netz: 7,0 %)
- c_L um **14,6 %** (Abstand: 11,1 %)
- St um nur **1,3 %** (Abstand: 6,4 %)

Für c_D und c_L ist die Streuung also so groß wie der Fehler selbst. Auf diesem Netz
bestimmt nicht die Auflösung den Fehler, sondern **wo die Treppe hinfällt** — und weil das
beim Verfeinern springt, gibt es keinen festen Vorfaktor, den eine Extrapolation
voraussetzen müsste. Der Rest kommt von der **Körperdarstellung**, nicht vom
Netzwerkansatz. Abhilfe wären geschnittene Zellen oder ein körperangepasstes Netz; das ist
der nächste Schritt und ausdrücklich noch nicht getan.

## Der unabhängige zweite Weg — das Wichtigste zuerst

Der Netzwerkansatz wird **gegengerechnet mit einem zweiten Löser aus anderer Herleitung**:
Galerkin mit Taylor-Hood-Elementen auf einem körperangepassten Dreiecksnetz. Schwache Form
statt Flussbilanz, Knotenwerte eines Polynomansatzes statt Zellmittel, Zylinderrand auf
Knoten statt Treppe, andere Kraftermittlung. Zwei Verfahren, die **verschieden scheitern
können** — das ist das Maß für Unabhängigkeit, nicht ob sie verschieden aussehen.

Er trifft **alle drei Benchmarkgrößen direkt im Intervall**, ohne Extrapolation:

| | c_D | c_L | St |
|:--|--:|--:|--:|
| Netzwerk, feinstes Netz (144 320 Zellen) | 3,2894 | 0,9744 | 0,2962 |
| **Galerkin, 2097 Knoten, dt = 1,25·10⁻⁴** | **3,2318** | **1,0092** | **0,3002** |
| Referenz Schäfer-Turek 2D-2 | 3,22 … 3,24 | 0,99 … 1,01 | 0,295 … 0,305 |

Die beobachtete Zeitordnung ist 1,00 — das Zeitschema ist erster Ordnung, die Messung
trifft genau das. Vorher, zur Probe: Hagen-Poiseuille kommt bei diesem Weg auf
4,9·10⁻¹⁴ heraus, weil die Parabel im P2-Ansatzraum **exakt enthalten** ist. Diese Schärfe
hat prompt einen Randfehler gefunden, den eine Probe mit der Schranke „klein" durchgelassen
hätte.

**Warum Galerkin und nicht Ritz:** Ritz minimiert ein Funktional; das setzt Symmetrie
voraus. Die Konvektion ist nicht symmetrisch, und die Inkompressibilität macht schon aus
Stokes einen Sattelpunkt — es gibt nichts zu minimieren. Galerkin braucht kein Funktional
und fällt beim Zähigkeitsterm wörtlich mit Ritz zusammen.

## Zwei Aufbauten desselben Verfahrens — und was das beweist

Kein gemeinsamer Rechenkern; gemeinsam sind nur Geometrie, Stoffwerte und Zeitschritt.

| | Weg A | Weg B |
|:--|:--|:--|
| Aufbau | jeder Zweig einzeln gestempelt | Flächenflüsse als ganze Felder |
| Darstellung | dünnbesetzte Knotenadmittanzmatrix | gar keine Matrix |
| Drucknetz | direkt gelöst (LU) | matrixfrei, konjugierte Gradienten |

Geprüft wird am **vollständigen Löser**: 6000 Zeitschritte mitten durch die ausgebildete
Wirbelablösung, größte Abweichung 5,2·10⁻¹⁴ im Feld. Wichtig zur Einordnung: Das belegt,
dass der **Aufbau** fehlerfrei ist — nicht, dass das Verfahren richtig ist. Beide teilen
Zeitschritt, Randbehandlung und Stempel, können also nur gemeinsam scheitern. Dafür ist der
Galerkin-Weg oben da.

Die Unsymmetrie der Impulsmatrix **ist** die gesteuerte Quelle; lässt man die
Geschwindigkeit gegen null gehen, wird die Matrix symmetrisch — das passive
Widerstandsnetz. Gemessen: die Unsymmetrie geht linear mit der Geschwindigkeit.

## Das Labor: beide Wege rechnen im Browser

`CFD_Labor_1.3.html` ist **eine Datei**, die ohne Netz läuft — Doppelklick genügt. Darin
rechnet nicht nur der Netzwerkweg, sondern seit dem 1. November auch der **Galerkin-Weg**
samt eigener Vernetzung: Triangle (Shewchuk) steckt als WebAssembly in der Datei,
derselbe Vernetzer, den FEMM benutzt.

Im Reiter *Vergleichsverfahren → Galerkin* laufen beide Wege nebeneinander, bei gleichem
Knotenabstand und gleichen Auswertepunkten. Was dabei zu sehen ist:

| | c_D auf grobem Browsernetz | umgerechnet (Python) | Schäfer-Turek |
|:--|--:|--:|--:|
| Galerkin, Taylor-Hood | 2,85 | **3,2318** | 3,22 … 3,24 |
| Netzwerk, Kirchhoff | 3,63 | 3,29 (+1,8 %) | |

Beide sind auf dem groben Netz daneben, **in entgegengesetzte Richtungen** — der
Benchmarkwert liegt dazwischen. Der Galerkin-Weg unterschätzt, der Netzwerkweg
überschätzt wegen der Treppenkontur. Das ist kein Widerspruch der Verfahren, sondern der
sichtbare Netzfehler beider, und er zeigt in verschiedene Richtungen, weil die
Herleitungen verschieden sind.

Der Galerkin-Browserkern ist **knotenweise** gegen die Python-Module geprüft, 47 Posten:
Hagen-Poiseuille auf 1,0e-12, die Beiwerte über 40 Zeitschritte auf 1e-9.

## Auf deine Anmerkungen zu 1.3

Alle drei Punkte sind in **1.5** erledigt:

* **Die Quellen fehlten.** Jetzt liegen alle Programme im Wortlaut in der Seite, Reiter
  **Programme** — 34 Dateien aus `stroemung/` und `labor/`, jede mit Zeilenzahl, Größe
  und sha256. Damit lässt sich nachrechnen, worauf das Prüfprotokoll verweist.
* **Dieses Blatt stand nur im Verzeichnis.** Jetzt steht es im Reiter **Anleitung**.
* **Lizenz für den eigenen Code:** **MIT**, wie in Ralphs übrigen Arbeiten. Eigener
  Reiter **Lizenz**, dazu `LICENSE` und `LIZENZ.md`.

**Zu deiner Werkzeugregel — die Sperre ist enger, als sie aussieht.** Triangle blockiert
nicht den Galerkin-Weg, sondern nur den Vernetzer:

| | hängt an Triangle? |
|:--|:--|
| `kern_fem.js` — der ganze Löser | **nein**, kein einziger Verweis |
| `netz_browser.js` | nur in **fünf Zeilen** von `netz()` |

Die Schnittstelle ist eine Punktliste hinein, eine Dreiecksliste heraus. Wenn ihr gmsh
ohnehin benutzt, sind das fünf Zeilen Arbeit — und ihr habt den vollständigen,
MIT-lizenzierten Löser ohne Triangle-Bindung. Nachprüfbar im Reiter Programme:
`kern_fem.js` enthält die Zeichenfolge „Triangle" kein einziges Mal.

**Zu deiner fachlichen Einschätzung: du hast recht.** Das Labor belegt laminare
2D-Strömung bei Re = 100 gegen eine fremde Referenz. Für Laufrad #44 — dreidimensional,
turbulent, drehend — rechnet es nichts nach und behauptet das nirgends. Der Weg dorthin
ist der, den du nennst.

## Was wir von dir brauchen

- **Wie heißt deine Fundstelle?** Name oder Adresse genügt, ich hole es selbst.
- Dein code-server verlangt Anmeldung — da komme ich nicht hinein und versuche es nicht.
- Zur Einordnung: Eine fremde CFD-Rechnung ist für uns **Gegenprobe**, keine Referenz
  (sie trägt deren Netz und Turbulenzmodell mit). Als Gegenprobe trotzdem willkommen.
- **Gemessene** Lüfterdaten (Kennlinie nach ISO 5801 / AMCA 210, PIV oder LDA) wären mehr
  wert als jede Rechnung.

## Kein Lüfterrad als Beweisstück

Das braucht 3D und ein Turbulenzmodell. Ein Turbulenzmodell ist eine Modellannahme mit
abgestimmten Konstanten — damit kann man abgleichen, nicht beweisen. Für die spätere
*Anwendung* liegt **FAN-01** hier (EAA-Benchmark, Niederdruck-Axiallüfter, Zenodo
CC-BY-4.0: gemessene Kennlinie 20 Punkte, 718…7483 m³/h bei 253…7 Pa, dazu LDA-Felder und
Geometrie).
