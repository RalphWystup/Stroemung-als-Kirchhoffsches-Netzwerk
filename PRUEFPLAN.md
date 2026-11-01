---
title: "Prüfplan aus den Aussagen des Auftraggebers"
subtitle: "Strömungsrechnung als Kirchhoffsches Netzwerk"
author: "Prof. Dr.-Ing. Ralph Wystup"
date: "31. Oktober 2026"
lang: de
---

Fassung 1.3 — gegenüber 1.2 neu: K15 (beide Wege rechnen in der HTML) und K16
(Prüfung über die kennzeichnende Zeit hinaus; Anlass Befund B20).
Fassung 1.2: K13 und K14 (die Oberfläche selbst; Anlass Befund B17).
Fassung 1.1: K7a und K7b (unabhängiger zweiter Weg), K8 berichtigt.

Nach Grundsatz 40 wird jedes Projekt gegen die **wörtlichen** Aussagen des Auftraggebers
geprüft, nicht gegen meine Auslegung davon. Links steht sein Wortlaut, rechts, woran sich
die Erfüllung ablesen lässt und wo der Beleg liegt.

| Nr. | Sein Wortlaut | Kriterium | Beleg |
|:--|:--|:--|:--|
| K1 | „ich will aber ein schönes Bild berechnen" | Es gibt ein Feldbild des Falls, und ich habe es **angesehen**, nicht nur erzeugt | `schaubild.py`, `wirbelstrasse_schaubild.png` |
| K2 | „eine CFD simulation dazu mit unserem Feldmodell zu vergleichen" | Derselbe Fall läuft über das Netzwerkmodell **und** über den klassischen Finite-Volumen-Aufbau | `wirbelstrasse_zweiwege.py` |
| K3 | „um eine zweite Methode zu finden" | Die beiden Wege teilen sich keinen Rechenkern — nur Geometrie, Stoffwerte und Zeitschritt | Weg A: Zweigstempel → Matrix → direkter Löser · Weg B: Feldoperationen → matrixfrei → CG |
| K4 | „anhand einen plakativen beispiels" | Der Fall zeigt das Gemeinte ohne Erklärung: Wirbel sind als Wirbel zu sehen | Wirbelstärkefeld mit abwechselnder Drehrichtung |
| K5 | „ohne Fake reine analytisch und nummerisch bewiesen!" | Keine Zahl ist geschätzt. Jede Aussage hat entweder eine Herleitung oder eine Messung hinter sich | `METHODE_UND_NACHWEIS.md`, `pruefprotokoll.py` |
| K6 | „Kármánsche Wirbelstraße … mit dem Störkörper und den Wirbeln" | Der Störkörper ist im Bild als Körper erkennbar; die Ablösung ist nicht angestoßen, sondern wächst von selbst | Grenzzyklus stellt sich aus der Ruhe ein |
| K7 | „gerechnet nach beiden verfahren" | Nicht nur die Operatoren, sondern der **vollständige Lauf** geht beide Wege | S10/S11 im Prüfprotokoll |
| K7a | „ich möchte die komplette berechnung unabhängig … komplett unabhängig“ | Ein zweiter Weg, der **verschieden scheitern kann** — nicht ein zweiter Aufbau desselben Verfahrens | Galerkin auf körperangepasstem Netz, S13/S14/S18 |
| K7b | „die herleitung der berechnungsansätze … die werkzeuge können dieselben sein“ | Die Unabhängigkeit wird an der **Herleitung** gemessen, nicht am Gleichungslöser | schwache Form gegen Flussbilanz; beide lösen dünnbesetzt |
| K8 | „Genau wie im Feld-Labor: zwei getrennt aufgebaute Löser … genau so!“ | Im Feld-Labor stehen Kirchhoff **und Ritz** nebeneinander — also zwei Herleitungen, nicht zwei Aufbauten | §5.6 in `METHODE_UND_NACHWEIS.md` |
| K9 | stehend: „rechnen, nie raten" | Offene Fragen werden entschieden, nicht vermutet — auch wenn das einen weiteren Lauf kostet | Befund B5: die ausbleibende Ablösung wurde durch das feinere Netz entschieden |
| K10 | stehend: „bei Fehlern sagen, warum sie damals nicht auffielen" | Jeder Befund nennt die Lücke in der bestehenden Prüfung und schließt sie | `BEFUNDE.md`, sechs Befunde |
| K11 | stehend: „alles auf C: sichern" | Projekt liegt vollständig auf dem Rechner des Auftraggebers | Ablage über die Brücke |
| K12 | „den jetson kannst du auch nehmen wenn du ihn brauchst" | Der Jetson wird genommen, **wenn** er gebraucht wird — nicht, um ihn benutzt zu haben | begründete Entscheidung, siehe unten |
| K13 | „ich breauche noch eine kurzfassung, 1 seite in dr html, Aufgabe, Zielsetzung, methode und fazit" · „die kurzfassung kann auch 2 Seiten lang sein" | Die Kurzfassung steht **als Reiter in der HTML**, mit diesen vier Abschnitten, dazu Netzwerk→Navier-Stokes, die Überprüfung und der Zweck des Surrogats | Reiter „Kurzfassung" in `CFD_Labor_1.2.html`; geprüft gegen die **Sollliste** in `pruefe_labor.js` |
| K14 | „geht das auch wieder als standalone html mit dokumentation ales reiter, vor allem die unterschiedlichen vergleichsverfahren" | Eine einzige Datei ohne Nachladen; die Reiterleiste führt **alle** geforderten Seiten, und die Prüfung vergleicht sie gegen eine Sollliste **außerhalb** der Seite | `baue_labor.py` bettet alles ein; `pruefe_labor.js` prüft Soll gegen Ist |
| K15 | „die html soll immer beide wege rechnen können" · „wenn alles gerechnet wird und miteinander verglichen und zwar online ist das doch ok" | **Beide** Wege rechnen in der Seite selbst — Netzwerk *und* Galerkin, nicht abgespielt. Dazu gehört der Vernetzer: das Netz entsteht im Browser | Reiter „Vergleichsverfahren → Galerkin"; `kern_fem.js`, `netz_browser.js`, Triangle als eingebettetes WebAssembly |
| K16 | stehend: „alles nachrechnen, nichts als geprüft ausgeben, was nicht gelaufen ist" | Ein Zeitschrittverfahren wird mindestens über die **kennzeichnende Zeit des Vorgangs** geprüft, nicht nur über ein paar Schritte | Befund B20; Dauerlauf bis t = 1,2 s in `pruefe_kern_fem.js` |

# Zum Jetson

Er wird nicht gebraucht — aber das ist gemessen, nicht angenommen. `zeitaufteilung.py`
zerlegt einen Zeitschritt auf dem Netz 440 × 82:

| Teil | je Zeitschritt | Anteil |
|:--|--:|--:|
| Konvektion und Zähigkeit (Feldoperationen) | 6,82 ms | 37,5 % |
| Drucknetz lösen (Einsetzen in die LU-Zerlegung) | 10,40 ms | 57,2 % |
| Rest (Divergenz, Korrektur, Kraftbilanz) | 0,97 ms | 5,3 % |
| **zusammen** | **18,19 ms** | |

Einmalig davor: Netz stempeln 462 ms, LU-Zerlegung 139 ms — bei 18 000 Zeitschritten
fällt das nicht ins Gewicht.

Nur der Löseanteil käme für eine Grafikkarte überhaupt in Frage, und selbst wenn er auf
null fiele, würde ein Zeitschritt höchstens **2,34-mal** schneller. Dafür müsste man vom
direkten Löser (einmal zerlegen, dann nur noch einsetzen) auf ein iteratives Verfahren
wechseln; das Einsetzen in eine Dreiecksmatrix ist von Natur aus sequentiell und wird auf
einer GPU nicht schneller. Bei den übrigen knapp 43 % — reinen Feldoperationen — ist der Jetson
Orin Nano diesem Server deutlich unterlegen. Ein Gewinn ist damit nicht zu erwarten, ein
Verlust wahrscheinlich.

Das wird hier festgehalten, damit die Entscheidung nachvollziehbar ist und nicht als
Versäumnis erscheint. Sollte später ein Fall kommen, bei dem das Drucknetz nicht mehr
direkt zerlegt werden kann — in drei Dimensionen wird die Zerlegung zu groß —, kehrt sich
die Rechnung um: Dann trägt ein iterativer Löser, und dann ist AmgX auf dem Jetson genau
das richtige Werkzeug. Es ist aus dem REM-Projekt bereits gebaut und geprüft.
