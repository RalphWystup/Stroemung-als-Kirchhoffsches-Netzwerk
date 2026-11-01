---
title: "Strömungsrechnung als Kirchhoffsches Netzwerk"
subtitle: "Kurzfassung — Fassung 1.3"
author: "Prof. Dr.-Ing. Ralph Wystup"
date: "1. November 2026"
lang: de
---

## Aufgabe

Eine Strömungsrechnung nach dem **Kirchhoffschen Netzwerkansatz** aufbauen — demselben, der
im `Feld_Labor` schon die Feldrechnung trägt — und an einem plakativen Beispiel gegen eine
Referenz prüfen, die sich nicht beeinflussen lässt. Analytisch und numerisch belegt, nichts
geschätzt.

## Zielsetzung

Für Strömungen **zwei Rechenwege** haben, die verschieden scheitern können — so wie für die
Feldrechnung Kirchhoff und Ritz nebeneinanderstehen. Beweisstück: die **Kármánsche
Wirbelstraße am Kreiszylinder, Re = 100** (Schäfer-Turek 2D-2). Dort ist der Nachlauf echt
zweidimensional — dreidimensionale Instabilitäten setzen erst ab Re ≈ 190 ein —, die
zweidimensionale Rechnung ist also keine Vereinfachung, sondern vollständig.

## Methode — wie die Navier-Stokes-Gleichung zum Netzwerk wird

Nicht als Analogie, sondern als Identität. Integriert man die Erhaltungssätze über eine
Zelle, so wird **jeder Term zu einem Netzwerkelement**:

| Term | integriert über die Zelle | im Netzwerk |
|:--|:--|:--|
| $\nabla\!\cdot\mathbf{u}=0$ | $\sum_\text{Flächen}\mathbf{u}\!\cdot\!\mathbf{n}\,A=0$ | **Knotenregel** |
| $\mu\nabla^2\mathbf{u}$ | $\sum \frac{\mu A}{d}(u_\text{Nachbar}-u_\text{Zelle})$ | **Leitwert** $G=\mu A/d$ |
| $\rho\,\partial\mathbf{u}/\partial t$ | $\rho V\,\dot u$ | **Kapazität** $C=\rho V$ |
| $\nabla p$ | Potentialdifferenz je Zweig | $p$ ist ein **Potential** |
| $(\mathbf{u}\!\cdot\!\nabla)\mathbf{u}$ | gerichteter Fluss durch die Fläche | **gesteuerte Quelle** |

Die ersten vier Zeilen ergeben ein **passives** Netz aus Widerständen und Kondensatoren —
genau das Netz des Feld-Labors. Die fünfte ist der Teil, den ein passives Netz nicht kann:
Konvektion ist gerichtet, hängt von der Lösung ab und macht die Matrix unsymmetrisch. Das
ist im Netzwerk nichts Fremdes, sondern eine **gesteuerte Quelle**, wie sie in jedem
Transistormodell steht.

Der entscheidende Schritt ist der **Zweigstempel**. Für einen einzelnen Zweig hat die
Konvektions-Diffusions-Gleichung eine geschlossene Lösung; man nimmt sie, statt zu nähern
(Scharfetter und Gummel 1969, aus der Halbleitersimulation):

$$J = \frac{D}{h}\left[B(-\mathrm{Pe}_h)c_i - B(\mathrm{Pe}_h)c_{i+1}\right],\quad
B(x)=\frac{x}{e^x-1},\quad g_\text{hin}-g_\text{rück} = u_f A .$$

Die **Differenz der beiden Leitwerte ist der Volumenstrom** — die gesteuerte Quelle steht
also wörtlich in der Unsymmetrie des Stempels. Für $u\to0$ wird die Matrix symmetrisch, und
übrig bleibt das Widerstandsnetz. Der Stempel ist knotenweise **exakt bei jeder
Péclet-Zahl**; übliche CFD *nähert* hier (zentral, Aufwind, QUICK, TVD).

Dazu das **versetzte Gitter**: Potentiale auf Knoten, Flüsse auf Zweigen. Das ist keine
Reparatur des Schachbrettproblems, sondern die topologisch richtige Zuordnung — das
Schachbrett entsteht umgekehrt gerade dadurch, dass man eine Zweiggröße auf einen Knoten
legt. Der Druck wird über ein zweites, rein passives Netz korrigiert (Poisson); die
Divergenz des Zwischenfelds ist der eingeprägte Strom.

## Beide Wege rechnen in der Seite

Seit dem 1. November 2026 rechnet das Browser-Labor **beide** Wege selbst — nicht nur
der Netzwerkweg, auch der Galerkin-Weg samt eigener Vernetzung. Dafür steckt Triangle
(Shewchuk) als WebAssembly in der Datei, derselbe Vernetzer, den auch FEMM benutzt.

Der Gleichungslöser ist dort ein anderer: In Python wird das Sattelsystem direkt
zerlegt, im Browser löst MINRES mit einem Blockvorkonditionierer. **Das ist kein anderer
Ansatz** — der Ansatz ist die schwache Form mit Taylor-Hood-Elementen, und der ist
derselbe. Belegt wird das mit 47 Prüfposten, die den Browserkern knotenweise gegen die
Python-Module halten: Hagen-Poiseuille auf $1{,}0\cdot10^{-12}$, die Beiwerte über 40
Zeitschritte auf $10^{-9}$.

Was die Seite dabei offen sagt: Die feinen Netze des Python-Laufs sind im Browser nicht
in Minuten zu rechnen. Auf dem groben Netz liegt der Galerkin-Weg **unter**, der
Netzwerkweg **über** dem Benchmarkwert — beide noch nicht umgerechnet, und zwar in
entgegengesetzte Richtungen. Der Benchmarkwert liegt dazwischen.

## Wie das geprüft wird

Stufenweise, jede Stufe gegen etwas, das sich nicht beeinflussen lässt: der Zweigstempel
gegen die geschlossene Lösung der 1D-Konvektion-Diffusion ($3{,}4\cdot10^{-15}$ bis
$\mathrm{Pe}=10^6$), die Kanalströmung gegen Hagen-Poiseuille (Ordnung $h^2$, dreimal exakt
2,00), die Knotenregel am **Ende** jedes Zeitschritts ($4{,}3\cdot10^{-15}$), die
Wirbelstraße gegen Schäfer-Turek.

Dazu zwei Proben verschiedener Art, die nicht verwechselt werden dürfen. **Die Sprachprobe**
stellt zwei *Aufbauten* derselben Diskretisierung gegeneinander — gestempelte Matrix gegen
matrixfreie Feldoperation —; über 6000 Zeitschritte stimmen sie auf $5{,}2\cdot10^{-14}$
überein. Das beweist, dass der Aufbau fehlerfrei ist, und **nicht**, dass das Verfahren
richtig ist: Beide teilen Zeitschritt, Randbehandlung und Stempel, können also nur
gemeinsam scheitern.

**Der unabhängige zweite Weg** ist deshalb ein anderer: **Galerkin** mit
Taylor-Hood-Elementen auf körperangepasstem Netz — schwache Form statt Flussbilanz,
Knotenwerte eines Polynomansatzes statt Zellmittel, Zylinderrand auf Knoten statt Treppe,
andere Kraftermittlung. Dass beide am Ende ein dünnbesetztes System lösen, ist eine Frage
des Werkzeugs, nicht des Ansatzes. *Ritz* ist hier nicht anwendbar: Er minimiert ein
Funktional, und die Konvektion ist nicht symmetrisch — es gibt keines. Galerkin ist die
Fortsetzung desselben Gedankens dorthin, wo die Symmetrie verlorengeht, und fällt beim
Zähigkeitsterm wörtlich mit Ritz zusammen.

## Wozu das Surrogat

Der vollständige Löser kostet je Fall Minuten bis Stunden; ein angelerntes Ersatzmodell —
POD für die Raumformen, harmonisch angepasstes Zeitgesetz — wertet in Mikrosekunden aus.
Es trägt den Zeitverlauf im Browser-Labor (eine volle Rechnung je Bildwechsel wäre
unmöglich); es erlaubt später, **Betriebspunkte zu durchfahren**, statt jeden einzeln
nachzurechnen — das ist der Unterschied zwischen Nachrechnen und Auslegen; und es macht
Struktur sichtbar: Dass 16 Moden 100,0000 % der Energie erfassen und die Singulärwerte
**paarweise** kommen, ist eine physikalische Aussage — der Grenzzyklus ist
niederdimensional. Brauchbar ist es nur, weil sein Fehler **gemessen** wird: 0,187 % der
Schwankungsbreite, gegen 7,33 % für einen angepassten linearen Zeitoperator.

## Fazit

**Der Netzwerkweg**: In der Netzreihe (10/20/40 Zellen über $D$) laufen alle vier
Benchmarkgrößen monoton auf die Referenz zu; auf dem feinsten Netz $c_D$ $+1{,}8\,\%$,
$c_L$ $-2{,}6\,\%$, $\Delta p$ $-1{,}0\,\%$, $\mathrm{St}$ **im Intervall**. Der Rest ist
**gemessen**: Vier Netze mit unveränderter Geometrie und gleicher Auflösung, nur anderer
Lage des Zylinders im Raster, streuen in $c_D$ um 4,7 % und in $c_L$ um 14,6 % — so groß
wie der Fehler selbst. Nicht die Auflösung bestimmt ihn, sondern die **Treppenkontur**.

**Der unabhängige Weg** schließt die Lücke. Er gibt zuerst Hagen-Poiseuille **exakt**
wieder ($4{,}9\cdot10^{-14}$ im Profil), weil die Parabel im Ansatzraum liegt. Am
Benchmark, 2097 Knoten, 81 davon auf dem Zylinderrand:

| $\Delta t$ | $c_D^{\max}$ | $c_L^{\max}$ | $\mathrm{St}$ |
|--:|--:|--:|--:|
| $5\cdot10^{-4}$ | 3,2456 | 1,0494 | 0,3008 |
| $2{,}5\cdot10^{-4}$ | 3,2364 | 1,0225 | 0,3005 |
| $\mathbf{1{,}25\cdot10^{-4}}$ | **3,2318** | **1,0092** | **0,3002** |
| Referenz | 3,22 … 3,24 | 0,99 … 1,01 | 0,295 … 0,305 |

**Alle drei Größen liegen im Referenzintervall — direkt gerechnet, ohne Extrapolation.**
Die beobachtete Konvergenzordnung in der Zeit beträgt 1,00 beziehungsweise 1,02; das
Zeitschema ist explizit-impliziter Euler, also erster Ordnung. Der Zeitfehler ist damit
nicht nur klein, sondern **verstanden**.

**Was damit belegt ist.** Zwei Wege, die verschieden scheitern können — Flussbilanz am
Kontrollvolumen mit Treppenkörper gegen schwache Form mit körperangepasstem Netz —, treffen
dieselbe gemessene Referenz und einander.

**Einschränkungen und nächste Schritte.** Das Zeitschema ist erster Ordnung; der Zeitfehler
ist damit der größte verbliebene. Ein turbulenter Fall verlangt ein Turbulenzmodell — damit
lässt sich abgleichen, nicht beweisen. Als Nächstes: Zeitschema zweiter Ordnung,
geschnittene Zellen auch für den Netzwerkweg, dann 3D und das Lüfterrad (FAN-01 liegt bereit).
