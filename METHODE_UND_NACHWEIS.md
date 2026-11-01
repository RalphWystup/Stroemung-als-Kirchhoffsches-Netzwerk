---
title: "Strömungsrechnung als Kirchhoffsches Netzwerk"
subtitle: "Verfahren und Nachweis — Fassung 1.3"
author: "Prof. Dr.-Ing. Ralph Wystup"
date: "31. Oktober 2026"
lang: de
---

Fassung 1.7, 01.11.2026. Gegenüber 1.6 neu: die Programme liegen im Wortlaut in der
Seite (Reiter Programme), neuer Abschnitt Quelltexte — auf Einwand von Fred.
Fassung 1.6, 01.11.2026. Gegenüber 1.5 neu: Danksagung an den Urheber des Vernetzers
und seine beiden Arbeiten im Schrifttum, dazu Cahouet/Chabard und Paige/Saunders.
Fassung 1.5, 01.11.2026. Gegenüber 1.4 neu: die Dateitafeln in §6 nachgezogen
(beide Browserkerne, Vernetzer, Triangle) und die Prüfstufen S19a/S19b.
Fassung 1.4, 01.11.2026. Gegenüber 1.3 neu: **§5.6a** — beide Wege rechnen jetzt im
Browser, einschließlich Vernetzung; der Gleichungslöser dort ist MINRES mit
Blockvorkonditionierer, geprüft knotenweise gegen die Python-Module.
Fassung 1.3, 01.11.2026. Gegenüber 1.2 neu: **der unabhängige zweite Weg** (§5.6) — Galerkin mit Taylor-Hood-Elementen auf körperangepasstem Netz, der alle drei Benchmarkintervalle direkt trifft. Gegenüber 1.1 neu: das Browser-Labor `CFD_Labor_1.0.html` samt seinem gegen die Python-Module geprüften Rechenkern (§6). Gegenüber 1.0 neu: der Impulsoperator als Knotenadmittanzmatrix
(§3.4), die Sprachprobe auf dem vollständigen Löser statt nur auf einzelnen Operatoren
(§5.5), die Druckdifferenz als vierte, unabhängige Referenzgröße (§5.6), die
Netzverfeinerungsreihe mit beobachteter Konvergenzordnung (§5.4) und das Prüfprotokoll
(§5.8). Berichtigt: §5.3 beschrieb ein Rohr statt des gerechneten ebenen Kanals und
behauptete eine Prüfung der Druckkopplung, die dort nicht stattfindet (Befund B9).

# 1 Aufgabe

Gesucht ist die Strömung eines inkompressiblen zähen Fluids:

$$\rho\left(\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u}\right)
= -\nabla p + \mu\nabla^2\mathbf{u} + \mathbf{f},
\qquad \nabla\cdot\mathbf{u} = 0 .$$

Behauptung dieses Textes: Diese Gleichungen lassen sich **vollständig** als elektrisches
Netzwerk aus Widerständen, Kondensatoren und gesteuerten Quellen schreiben — nicht als
Analogie, sondern als dieselbe Gleichung in anderer Schreibweise. Der Nachweis erfolgt in
Abschnitt 5 an geschlossenen Lösungen und an einem Benchmark mit veröffentlichten
Referenzwerten.

# 2 Vom Erhaltungssatz zum Netzwerk

## 2.1 Die Kontinuitätsgleichung ist die Knotenregel

Man integriere $\nabla\cdot\mathbf{u}=0$ über ein Kontrollvolumen $V$ und wende den
Gaußschen Satz an:

$$\int_V \nabla\cdot\mathbf{u}\,\mathrm{d}V = \oint_S \mathbf{u}\cdot\mathbf{n}\,\mathrm{d}A
= \sum_{\text{Flächen } f} u_f A_f = \sum_f Q_f = 0 .$$

Das ist Kirchhoffs Knotenregel, Wort für Wort: Die Summe der Zweigströme an einem Knoten
ist null. **Zelle = Knoten, Fläche = Zweig, Volumenstrom = Strom.** Hier wird nichts
übertragen; es ist dieselbe Aussage.

## 2.2 Die Impulsgleichung als Zweigbilanz

Für eine Geschwindigkeitskomponente, über die Zelle integriert:

$$\underbrace{\rho V \frac{\mathrm{d}u_P}{\mathrm{d}t}}_{\text{Trägheit}}
+ \underbrace{\sum_f \dot m_f\, u_f}_{\text{Konvektion}}
= \underbrace{\sum_f \frac{\mu A_f}{d_f}\,(u_N - u_P)}_{\text{Zähigkeit}}
- \underbrace{V\frac{\partial p}{\partial x}}_{\text{Druck}}
+ \underbrace{V f_x}_{\text{Volumenkraft}}$$

Term für Term:

| Term | Bauelement | Wert |
|:--|:--|:--|
| $\rho V\,\mathrm{d}u/\mathrm{d}t$ | Kondensator | $C=\rho V$, Energie $\tfrac12 C u^2 = \tfrac12\rho V u^2$ — die kinetische Energie |
| $\mu A_f/d_f\,(u_N-u_P)$ | Widerstand zum Nachbarn | $G=\mu A_f/d_f$ |
| $\dot m_f u_f$ | **gesteuerte Stromquelle** | Stärke aus $\dot m_f$, das selbst von Knotenwerten abhängt |
| $-V\,\partial p/\partial x$ | Quelle aus dem Drucknetz | Kopplung |
| $V f_x$ | eingeprägte Stromquelle | Schwerkraft, Coriolis |

Knotenpotential ist hier die Geschwindigkeitskomponente. In drei Dimensionen sind das drei
solche Netze auf demselben Knotensatz.

Der entscheidende Unterschied zur Feldrechnung: Ein passiver Widerstand hat **einen**
Leitwert, in beide Richtungen denselben; die zugehörige Matrix ist symmetrisch. Konvektion
ist gerichtet und erzeugt eine **unsymmetrische** Matrix. Genau diese Unsymmetrie ist die
gesteuerte Quelle (Abschnitt 3.3).

## 2.3 Warum der Druck ein eigenes Netz braucht

Löst man die Impulsnetze mit einem geratenen Druckfeld $p^*$, so erfüllt das Ergebnis
$\mathbf{u}^*$ die Knotenregel nicht. Die diskrete Impulsgleichung lautet

$$a_P u_P = \sum_{nb} a_{nb} u_{nb} + b - V(\nabla p)_P .$$

Mit $\mathbf{u}=\mathbf{u}^*+\mathbf{u}'$ und $p=p^*+p'$, Subtraktion der Gleichung für
$\mathbf{u}^*$ und Vernachlässigung von $\sum a_{nb}u'_{nb}$ (das ist die SIMPLE-Näherung,
sie betrifft nur den *Weg* zur Lösung, nicht die Lösung selbst):

$$a_P u'_P = -V(\nabla p')_P
\qquad\Longrightarrow\qquad
\mathbf{u}' = -\frac{V}{a_P}\nabla p' .$$

Einsetzen in $\nabla\cdot\mathbf{u}=0$ liefert

$$\nabla\cdot\left(\frac{V}{a_P}\nabla p'\right) = \nabla\cdot\mathbf{u}^* .$$

Das ist ein **reines Widerstandsnetz**: Leitwert $V/a_P$, Knotenpotential $p'$, und als
eingeprägter Strom an jedem Knoten genau der Betrag, der dort an der Knotenregel noch
fehlt. Es ist dieselbe Poisson-Struktur, die das `Feld_Labor` bereits rechnet.

Dann wird korrigiert und wiederholt. Beide Nichtlinearitäten — die gesteuerte Quelle hängt
von der Lösung ihres eigenen Netzes ab, und Impuls und Druck hängen wechselseitig — werden
durch diese äußere Iteration behandelt. Alternativ linearisiert ein Newton-Schritt das
ganze gekoppelte System auf einmal; das konvergiert schneller, verlangt aber einen
Vorkonditionierer für die Sattelpunktmatrix.

## 2.4 Die duale Zuordnung — wo die Induktivität sitzt

Abschnitt 2.2 wählt die Geschwindigkeit als Knotenpotential; dort trägt der Kondensator die
Trägheit. Die klassische hydraulische Zuordnung ist die duale: Knotenpotential ist der
**Druck**, Zweigstrom der **Volumenstrom**. Für eine Rohrstrecke der Länge $l$ und des
Querschnitts $A$ gilt dann aus dem Impulssatz der Fluidsäule

$$\Delta p\,A = \rho A l\,\frac{\mathrm{d}u}{\mathrm{d}t}
\qquad\Longrightarrow\qquad
\Delta p = \underbrace{\frac{\rho l}{A}}_{L}\,\frac{\mathrm{d}Q}{\mathrm{d}t} ,$$

also eine **Induktivität** $L=\rho l/A$ mit $\tfrac12 LQ^2 = \tfrac12 (\rho l A) u^2$, wieder
der kinetischen Energie. Die Kompressibilität liefert $C=V/(\rho c^2)$. Mit den Belägen
$L'=\rho/A$ und $C'=A/(\rho c^2)$ ergibt sich

$$\frac{1}{\sqrt{L'C'}} = c, \qquad \sqrt{\frac{L'}{C'}} = \frac{\rho c}{A} ,$$

also Schallgeschwindigkeit und Wellenwiderstand: **ein Rohr ist eine Leitung**, beschrieben
durch die Telegraphengleichung. Im inkompressiblen Grenzfall $c\to\infty$ geht $C\to 0$; die
Leitung entartet, der Druck wird von einer dynamischen Größe zu einer Zwangsbedingung.
Das ist genau die Sattelpunktstruktur aus Abschnitt 2.3.

Nach Zeitdiskretisierung mit der Trapezregel werden $L$ und $C$ zu Begleitmodellen —
$G=\Delta t/(2L)$ bzw. $G=2C/\Delta t$, jeweils parallel zu einer Stromquelle aus der
Vorgeschichte. In jedem Zeitschritt bleibt ein reines Widerstandsnetz mit gesteuerten
Quellen übrig; Trägheit, Kompressibilität und Konvektion fügen sich damit in dieselbe
Maschinerie.

# 3 Der Zweigstempel

## 3.1 Warum die naheliegenden Stempel versagen

Auf einem Zweig ist der Fluss $J = u\,c - D\,\mathrm{d}c/\mathrm{d}x$. Maßgeblich ist die
**Gitter-Pécletzahl** $\mathrm{Pe}_h = u h/D$.

Die *zentrale Differenz* ist zweiter Ordnung, erzeugt aber für $\mathrm{Pe}_h>2$ negative
Leitwerte und damit Schwingungen von Zelle zu Zelle. Der *Aufwind* ist bedingungslos
stabil, aber nur erster Ordnung; sein Abbruchfehler wirkt wie eine zusätzliche Diffusion
mit $D_{\text{num}} = u h/2$. Bei Farbstoff in Wasser ($D\approx 10^{-9}\,$m²/s,
$u\approx 0{,}1\,$m/s, $h\approx 10^{-4}\,$m) ist $D_{\text{num}} \approx 5\cdot10^{-6}$ —
das **Fünftausendfache** des physikalischen Werts. Das Bild sähe dann glatt aus und zeigte
doch nur das Gitter.

## 3.2 Scharfetter-Gummel: die Zweiglösung als Leitwert

Statt die Ableitung zu nähern, löst man die Gleichung auf dem Zweig exakt. Mit konstantem
$u$, $D$ und $J$ zwischen $x=0$ (Wert $c_i$) und $x=h$ (Wert $c_{i+1}$):

$$\frac{\mathrm{d}c}{\mathrm{d}x} = \frac{u c - J}{D}
\qquad\Longrightarrow\qquad
c(x) = \frac{J}{u} + \left(c_i-\frac{J}{u}\right)e^{ux/D} .$$

Einsetzen von $c(h)=c_{i+1}$ und Auflösen nach $J$ liefert

$$J = u\,\frac{c_i e^{\mathrm{Pe}_h} - c_{i+1}}{e^{\mathrm{Pe}_h}-1} .$$

Mit der Bernoulli-Funktion $B(x)=x/(e^x-1)$, für die $B(-x)=x e^{x}/(e^{x}-1)$ gilt, wird
daraus die Stempelform

$$\boxed{\;J_{i+1/2} = \frac{D}{h}\Bigl[\,B(-\mathrm{Pe}_h)\,c_i \;-\; B(\mathrm{Pe}_h)\,c_{i+1}\,\Bigr]\;}$$

Die Grenzfälle bestätigen die Bauart: $\mathrm{Pe}_h\to 0$ gibt $B\to 1$ und damit reine
Diffusion $J=(D/h)(c_i-c_{i+1})$; $\mathrm{Pe}_h\to+\infty$ gibt $J=u\,c_i$, also Aufwind
von links; $\mathrm{Pe}_h\to-\infty$ gibt $J=u\,c_{i+1}$, Aufwind von rechts. Der Stempel
geht also von selbst stetig zwischen den beiden Grenzverhalten über.

Numerisch ist $B$ in vier Bereichen auszuwerten — Reihenentwicklung nahe null, `expm1`
im Mittelfeld, und zwei Grenzfälle gegen Über- und Unterlauf. Das ist keine Kosmetik:
Ohne diese Fallunterscheidung ist der Stempel genau dort falsch, wo er gebraucht wird.

## 3.3 Die Unsymmetrie ist die gesteuerte Quelle

Die beiden Vorfaktoren $B(-\mathrm{Pe}_h)$ und $B(\mathrm{Pe}_h)$ sind **verschieden**. Der
Zweig leitet in Strömungsrichtung anders als dagegen. Ein passiver Widerstand kann das
nicht; die Differenz

$$\frac{D}{h}\bigl[B(-\mathrm{Pe}_h) - B(\mathrm{Pe}_h)\bigr] = u$$

ist genau die gesteuerte Quelle, deren Stärke vom Strömungszustand abhängt. Da
$\mathrm{Pe}_h$ seinerseits von der Lösung abhängt, ist der Stempel nichtlinear und wird
durch Newton nachgeführt.

## 3.4 Die Impulsgleichung als Knotenadmittanzmatrix

Bis hierher steht der Stempel für einen einzelnen Zweig. Setzt man ihn für **jeden** Zweig
und bilanziert an jedem Knoten, so entsteht eine Matrix — und zwar genau die
Knotenadmittanzmatrix der Netzwerkanalyse. Für den Zweig zwischen den Knoten $a$ und $b$
trägt der Stempel vier Einträge ein:

$$M_{aa} \mathrel{-}= g_\text{hin},\quad M_{ab} \mathrel{+}= g_\text{rück},\quad
  M_{bb} \mathrel{-}= g_\text{rück},\quad M_{ba} \mathrel{+}= g_\text{hin} .$$

Bei $g_\text{hin} = g_\text{rück}$ wäre das die symmetrische Matrix eines reinen
Widerstandsnetzes. Weil der Scharfetter-Gummel-Stempel die beiden Leitwerte um den
Volumenstrom auseinanderzieht, ist sie **unsymmetrisch** — und diese Unsymmetrie ist die
gesteuerte Quelle, jetzt nicht mehr an einem Zweig, sondern als Eigenschaft der ganzen
Matrix:

$$M - M^{\mathsf T} \;\neq\; 0 \quad\Longleftrightarrow\quad \text{es strömt.}$$

Das ist nachprüfbar und wird geprüft: Lässt man die Geschwindigkeit gegen null gehen, muss
$\|M - M^{\mathsf T}\|/\|M\|$ mit ihr verschwinden. Gemessen wird $9{,}2\cdot10^{-9}$ bei
einer Geschwindigkeit von $10^{-9}$ — die Unsymmetrie geht **linear** mit der
Geschwindigkeit, wie es der Stempel verlangt.

### Was mit den Randwerten geschieht

Ein Geisterwert außerhalb des Gebiets ist entweder eine Kopie eines inneren Unbekannten
(freier Ausstrom), seine Spiegelung mit Faktor $-1$ (Haftbedingung, damit die Wand
zwischen Zelle und Spiegel bei $u=0$ liegt) oder ein echter Vorgabewert (Einströmprofil).
Nur der dritte Fall ist eine Konstante; die ersten beiden gehören als Matrixeintrag auf
den Knoten, von dem sie abhängen. Wer sie stattdessen in den Lastvektor schreibt, bekommt
zwar dasselbe Residuum, aber nicht dieselbe Matrix — und merkt den Unterschied erst, wenn
er sie implizit lösen will (Befund B6).

# 4 Das versetzte Gitter

## 4.1 Das Schachbrettproblem

Liegen Druck und Geschwindigkeit auf denselben Knoten, so wird der Druckgradient zentral
gebildet. Für ein Schachbrettfeld $p_i = (-1)^i A$ folgt

$$\left.\frac{\partial p}{\partial x}\right|_i = \frac{p_{i+1}-p_{i-1}}{2h}
= \frac{A-A}{2h} = 0 .$$

Die Oszillation ist für die Gleichungen unsichtbar und kann beliebig anwachsen.

## 4.2 Warum das Netzwerkdenken das versetzte Gitter vorschreibt

Sitzt die Geschwindigkeit dagegen auf der Fläche zwischen zwei Druckknoten, so ist

$$\left.\frac{\partial p}{\partial x}\right|_{i+1/2} = \frac{p_{i+1}-p_i}{h} = -\frac{2A}{h} ,$$

also groß — die Oszillation wird sofort bestraft. Wesentlich ist der **Versatz um eine halbe
Zelle**, nicht die Nähe zweier Knoten.

Im Netzwerkbild ist das keine Reparatur, sondern die richtige Zuordnung: **Potentiale leben
auf Knoten, Flüsse auf Zweigen.** Der Druck ist ein Potential (Zellmitte), der Durchfluss
ein Fluss (Fläche). Das Schachbrettmuster entsteht umgekehrt gerade dadurch, dass man eine
Zweiggröße fälschlich auf einen Knoten legt. Das versetzte Gitter nach Harlow und Welch
(1965) ist damit die topologisch korrekte Aufstellung.

Zwei Alternativen bestehen und sind erwähnenswert: gemischte Elemente nach Taylor und Hood
(gleiche Orte, aber reicherer Geschwindigkeitsansatz, der die inf-sup-Bedingung erfüllt)
sowie stabilisierte gleichordnige Verfahren (PSPG), die einen zusätzlichen Druckterm
einführen. Letztere bringen einen frei zu wählenden Parameter mit und werden hier
vermieden.

# 5 Der Nachweis

## 5.1 Grundsatz: was als Referenz zählt

Referenz ist eine **geschlossene Lösung** oder eine **Messung**. Die Ausgabe eines anderen
Rechenprogramms ist eine Gegenprobe, keine Referenz: Sie trägt dessen Netz, Randbedingungen
und — bei turbulenten Rechnungen — dessen Turbulenzmodell mit. Übereinstimmung bedeutet dann
nur, dass dieselben Annahmen getroffen wurden.

Jede Stufe wird gegen eine Referenz geprüft, die von der Rechnung unabhängig ist, und keine
Stufe wird übersprungen: Ein Fehler im Zweigstempel fällt bei $\mathrm{Pe}=10^6$ auf die
zehnte Stelle auf, ginge aber in einer Turbulenzrechnung spurlos unter.

## 5.2 Stufe 1 — der Zweigstempel gegen die geschlossene Lösung

Stationäre eindimensionale Konvektion-Diffusion mit $c(0)=0$, $c(L)=1$ hat die Lösung

$$c(x) = \frac{e^{\mathrm{Pe}\,x/L}-1}{e^{\mathrm{Pe}}-1},
\qquad \mathrm{Pe} = \frac{uL}{D} .$$

Scharfetter-Gummel ist auf gleichmäßigem Gitter **knotenweise exakt**; das ist eine
beweisbare Eigenschaft des Stempels, keine Erfahrungstatsache. Gemessen auf 20 Zellen,
größter Knotenfehler:

| $\mathrm{Pe}$ | $\mathrm{Pe}_h$ | Scharfetter-Gummel | zentral | Aufwind |
|--:|--:|--:|--:|--:|
| 1 | 0,05 | $3{,}4\cdot10^{-15}$ | $2{,}5\cdot10^{-5}$ | $2{,}9\cdot10^{-3}$ |
| 5 | 0,25 | $5{,}0\cdot10^{-16}$ | $1{,}8\cdot10^{-3}$ | $3{,}9\cdot10^{-2}$ |
| 20 | 1 | $3{,}5\cdot10^{-18}$ | $3{,}5\cdot10^{-2}$ | $1{,}3\cdot10^{-1}$ |
| 100 | 5 | $6{,}8\cdot10^{-21}$ | $4{,}4\cdot10^{-1}$ | $1{,}6\cdot10^{-1}$ |
| $10^3$ | 50 | $1{,}2\cdot10^{-35}$ | $1{,}4$ | $2{,}0\cdot10^{-2}$ |
| $10^6$ | $5\cdot10^4$ | $0$ | $1{,}25\cdot10^{3}$ | $2{,}0\cdot10^{-5}$ |

**Urteil: bestanden.** Der größte Fehler über alle Fälle beträgt $3{,}4\cdot10^{-15}$, also
Maschinengenauigkeit.

Zur Einordnung des Aufwinds: Sein kleiner Fehler bei $\mathrm{Pe}=10^6$ ist kein Verdienst.
Die exakte Lösung ist dort eine Stufe am rechten Rand, und eine Stufe trifft der Aufwind
trivial. Sein tatsächlicher Fehler liegt bei $\mathrm{Pe}_h\approx 1\ldots5$ und beträgt
dort 13 bis 16 Prozent.

## 5.3 Stufe 2 — Kanalströmung gegen Hagen-Poiseuille

Berichtigt am 31.10.2026. Dieser Abschnitt beschrieb vorher ein **Rohr** mit
$u(r)=2\bar u(1-r^2/R^2)$ und $\Delta p = 8\mu L Q/(\pi R^4)$, während das Programm einen
**ebenen Kanal** rechnet. Und er behauptete, hier werde die Druckkopplung mitgeprüft — das
ist nachweislich nicht so (Befund B2). Beides steht jetzt richtig da.

Die ausgebildete laminare Strömung im ebenen Kanal der Höhe $h$ unter dem Antrieb $G$ ist
eine **exakte** Lösung der Navier-Stokes-Gleichungen:

$$u(y) = \frac{G}{2\nu}\,y\,(h-y),
\qquad u_{\max} = \frac{G h^2}{8\nu},
\qquad \bar u = \frac{2}{3}\,u_{\max} .$$

Gemessen wird gegen diese Formel, nicht gegen eine andere Rechnung:

| $n_y$ | $u_{\max}$ | Fehler | $\bar u$ | Fehler | Profil $L_2$ | Ordnung |
|--:|--:|--:|--:|--:|--:|--:|
| 8 | 0,125000000 | $2{,}2\cdot10^{-16}$ | 0,085937500 | $2{,}60\cdot10^{-3}$ | $1{,}95\cdot10^{-3}$ | |
| 16 | 0,125000000 | $2{,}4\cdot10^{-14}$ | 0,083984375 | $6{,}51\cdot10^{-4}$ | $4{,}88\cdot10^{-4}$ | 2,00 |
| 32 | 0,125000000 | $2{,}4\cdot10^{-13}$ | 0,083496094 | $1{,}63\cdot10^{-4}$ | $1{,}22\cdot10^{-4}$ | 2,00 |
| 64 | 0,125000000 | $1{,}1\cdot10^{-12}$ | 0,083374023 | $4{,}07\cdot10^{-5}$ | $3{,}05\cdot10^{-5}$ | 2,00 |

Zwei Dinge sind daran bemerkenswert.

**Der Scheitelwert ist auf Maschinengenauigkeit exakt**, und zwar schon bei acht Zellen.
Das ist kein Glück: Die Zweitableitung eines Polynoms zweiten Grades wird vom
Zweipunktstempel exakt wiedergegeben. Der verbleibende Fehler in $\bar u$ und im
$L_2$-Maß ist allein der **Quadraturfehler** der Mittelung über die Zellen — und der
fällt sauber mit $h^2$, dreimal hintereinander exakt 2,00.

**Was diese Stufe NICHT prüft.** In der ausgebildeten Kanalströmung ist die Divergenz
identisch null. Die rechte Seite des Drucknetzes ist damit null, gleich welcher Vorfaktor
und welches Vorzeichen davorstehen. Die Druckkopplung wird hier also **nicht** geprüft —
sie wird erst in Stufe 3 geprüft, wo die Divergenz im Zwischenfeld von null verschieden
ist. Genau diese Lücke hat einen Fehler mit dem Faktor 2560 und falschem Vorzeichen
durchgelassen (Befund B2). Stufe 2 prüft die **Zähigkeitswiderstände** und die
**Haftbedingung**, und nur die.

**Noch offen.** Die Taylor-Aris-Dispersion — ein Farbstoff verteilt sich längs des Rohres
nicht mit $D$, sondern mit $D_\text{eff} = D + a^2\bar u^2/(48 D)$ — wäre eine weitere
geschlossene Lösung, die empfindlich auf numerische Diffusion reagiert und damit genau die
Schwachstelle aus Abschnitt 3.1 prüft. Sie ist **noch nicht gerechnet** und steht hier als
Vorhaben, nicht als Nachweis.

## 5.4 Stufe 3 — die Wirbelstraße gegen Schäfer und Turek

Der Benchmark „Flow around a cylinder" (Schäfer und Turek, DFG-Prioritätsprogramm 1996)
ist der Standardfall der zweidimensionalen instationären Strömungsrechnung. Kanal
$2{,}2 \times 0{,}41$, Zylinder mit $D=0{,}1$ bei $(0{,}2;\,0{,}2)$ — leicht außermittig,
damit die Ablösung anläuft —, parabolisches Einströmprofil mit $\bar u = 1$, $\nu=10^{-3}$,
also $\mathrm{Re} = \bar u D/\nu = 100$.

Veröffentlichte Referenzintervalle, von mehreren Gruppen unabhängig bestätigt:

| Größe | Referenz |
|:--|:--|
| Widerstandsbeiwert $c_D^{\max}$ | 3,22 … 3,24 |
| Auftriebsbeiwert $c_L^{\max}$ | 0,99 … 1,01 |
| Strouhalzahl $\mathrm{St}=fD/\bar u$ | 0,295 … 0,305 |
| Druckdifferenz $\Delta p(t_0+\tfrac{1}{2f})$ | 2,46 … 2,50 |

Diese Zahlen sind dimensionslos und stammen aus Rechnung *und* Messung. Sie lassen sich
nicht anpassen — ein Verfahren trifft sie oder nicht.

Warum gerade dieser Fall: Bei $\mathrm{Re}=100$ ist der Nachlauf **echt zweidimensional**;
dreidimensionale Instabilitäten setzen erst ab $\mathrm{Re}\approx 190$ ein. Die
zweidimensionale Rechnung ist hier also keine Vereinfachung, sondern vollständig.

Einschränkung, die vorab zu nennen ist: Auf kartesischem Gitter ist der runde Zylinder nur
näherungsweise darstellbar. Die Strouhalzahl ist davon wenig betroffen, die Beiwerte $c_D$
und $c_L$ dagegen empfindlich. Es wird deshalb eine Netzverfeinerungsreihe mitgeführt und
der erreichte Wert mit seinem Trend berichtet — nicht ein einzelner Lauf, der zufällig
im Intervall liegt.

### Die Netzverfeinerungsreihe

<!-- TABELLE-NETZREIHE -->

| Zellen über $D$ | Netz | $c_D^{\max}$ | $c_L^{\max}$ | $\mathrm{St}$ | $\Delta p$ |
|--:|:--|--:|--:|--:|--:|
| 10 | 220 × 41 | 3,7514 | 0,6103 | 0,2525 | 2,4075 |
| 20 | 440 × 82 | 3,4546 | 0,8888 | 0,2809 | 2,4343 |
| 40 | 880 × 164 | 3,2894 | 0,9744 | 0,2962 | 2,4547 |
| **Referenz** | **—** | **3,22 … 3,24** | **0,99 … 1,01** | **0,295 … 0,305** | **2,46 … 2,50** |

Abweichung von der Mitte des jeweiligen Referenzintervalls:

| Zellen über $D$ | $c_D$ | $c_L$ | $\mathrm{St}$ | $\Delta p$ |
|--:|--:|--:|--:|--:|
| 10 | +16,1 % | −39,0 % | −15,8 % | −2,9 % |
| 20 | +7,0 % | −11,1 % | −6,4 % | −1,8 % |
| 40 | +1,8 % | −2,6 % | −1,3 % | −1,0 % |

<!-- ENDE-TABELLE-NETZREIHE -->

### Wie die Ordnung bestimmt wird — und was dabei mitgerechnet wird

Mit **zwei** Netzstufen lässt sich nur extrapolieren, wenn man die Konvergenzordnung
vorgibt. Das wäre Raten. Mit **drei** Stufen bei gleichem Verfeinerungsverhältnis $r$
folgt sie aus den Daten selbst:

$$p = \frac{\ln\bigl|(f_{4h}-f_{2h})/(f_{2h}-f_h)\bigr|}{\ln r},
\qquad f_\text{exakt} = f_h + \frac{f_h - f_{2h}}{r^{\,p}-1} .$$

Zwei Dinge sind dabei ehrlich zu benennen.

**Erstens** ist die Voraussetzung, dass die drei Werte im asymptotischen Bereich liegen.
Das zeigt sich daran, dass die Änderungen sich mit festem Faktor verkleinern — der Quotient
muss größer als eins sein. Ist er es nicht, wird **kein** Grenzwert ausgegeben, sondern
„nicht im asymptotischen Bereich". Ein erfundener Wert wäre schlimmer als keiner.

**Zweitens** ist der Zeitschritt über die CFL-Bedingung an den Zellabstand gekoppelt,
$\Delta t \propto h$. Die Verfeinerung verfeinert also Raum **und** Zeit zugleich. Die
beobachtete Ordnung $p$ ist damit eine **gemeinsame** Ordnung beider Diskretisierungen,
nicht die räumliche allein. Für die Extrapolation ist das zulässig, solange ein
Fehleranteil überwiegt; es wäre aber falsch, $p$ als Aussage über das Raumschema zu lesen.

**Drittens** — und das ist die eigentliche Einschränkung dieses Falls — ist der runde
Zylinder auf kartesischem Gitter eine Treppe. Die Treppe verfälscht vor allem die
Wandschubspannung und die Lage der Ablösepunkte. Das erklärt, warum die integrierten
Kraftbeiwerte empfindlicher reagieren als die Strouhalzahl und viel empfindlicher als die
Druckdifferenz aus Abschnitt 5.6. Ein Verfahren mit körperangepasstem Netz oder mit
geschnittenen Zellen würde hier schneller konvergieren; das ist eine Aussage über die
**Körperdarstellung**, nicht über den Netzwerkansatz.

### Was die Extrapolation ergibt — und warum sie hier nicht trägt

| Größe | beobachtete Ordnung $p$ | extrapolierter Grenzwert | Referenz | Urteil |
|:--|--:|--:|:--|:--|
| $c_D^{\max}$ | 0,85 | 3,082 | 3,22 … 3,24 | außerhalb |
| $c_L^{\max}$ | 1,70 | 1,012 | 0,99 … 1,01 | außerhalb |
| $\mathrm{St}$ | 0,90 | 0,314 | 0,295 … 0,305 | außerhalb |

**Das Kriterium ist damit nicht erfüllt**, und es wird nicht nachträglich gelockert. Die
Zahlen stehen hier so, wie sie herausgekommen sind.

Sie sagen aber mehr als „nicht erfüllt". Die Richardson-Extrapolation setzt voraus, dass
**ein** Fehleranteil den Rest beherrscht, $f_h = f_\text{exakt} + C h^{p}$ mit festem $C$
und festem $p$. Dann müssten die drei Größen, die aus **derselben** Diskretisierung
stammen, auch ähnliche Ordnungen zeigen. Gemessen wird 0,85 … 1,70, eine Spanne von 0,86.
Die Voraussetzung ist also nachweislich verletzt — und ein extrapolierter Wert, dessen
Voraussetzung verletzt ist, ist keine Aussage, gleichgültig wo er landet. Das Prüfprotokoll
gibt diese Spanne deshalb als eigene Zeile aus, damit das Urteil nicht nackt dasteht.

Woher die Verletzung kommt, ist nicht zu raten, sondern zu messen; das geschieht im
nächsten Abschnitt.

### Die Treppe, gemessen statt vermutet

<!-- TABELLE-TREPPE -->

| Netz | $\Delta y$ | Mittelpunkt im Raster | $c_D^{\max}$ | $c_L^{\max}$ | $\mathrm{St}$ |
|:--|--:|:--|--:|--:|--:|
| 440 × 82 | 0,005000 | x 0,00 Δx, y 0,00 Δy | 3,4546 | 0,8888 | 0,2809 |
| 441 × 83 | 0,004940 | x 0,09 Δx, y 0,49 Δy | 3,4808 | 0,9098 | 0,2812 |
| 442 × 84 | 0,004881 | x 0,18 Δx, y 0,98 Δy | 3,3194 | 0,7841 | 0,2846 |
| 443 × 82 | 0,005000 | x 0,27 Δx, y 0,00 Δy | 3,4089 | 0,8589 | 0,2828 |

- **$c_D$**: 3,3194 … 3,4808, Spanne 0,1615 (4,7 % vom Mittelwert)
- **$c_L$**: 0,7841 … 0,9098, Spanne 0,1257 (14,6 % vom Mittelwert)
- **$St$**: 0,2809 … 0,2846, Spanne 0,0037 (1,3 % vom Mittelwert)

<!-- ENDE-TABELLE-TREPPE -->

**So ist die Probe gebaut.** Die Geometrie bleibt unverändert — Kanal, Zylinder,
Stoffwerte, Einströmprofil, alles wie im Benchmark. Verändert wird allein die Zellzahl, um
eins, zwei und drei. Der Zellabstand ändert sich dadurch um weniger als zweieinhalb
Prozent; die Auflösung ist also dieselbe. Was sich stark ändert, ist die **Lage** des
Zylindermittelpunkts im Zellraster und damit das Treppenmuster. Die bedeckte Fläche liegt
in allen vier Fällen bei $100\pm1\,\%$ der wahren Kreisfläche — es geht also nicht um einen
zu großen oder zu kleinen Körper, sondern um die zackige Kontur.

**Was dabei herauskommt.** Allein durch die Lage der Treppe streut

- $c_D$ um 4,7 %, während der Abstand zur Referenz auf diesem Netz 7,0 % beträgt,
- $c_L$ um 14,6 %, bei einem Abstand von 11,1 %,
- $\mathrm{St}$ dagegen nur um 1,3 %, bei einem Abstand von 6,4 %.

Für $c_D$ und $c_L$ ist die Streuung also **so groß wie der Fehler selbst**. Damit ist
belegt, was vorher nur eine Vermutung war: Auf diesem Netz wird der Fehler nicht von der
Auflösung bestimmt, sondern davon, **wo die Treppe hinfällt**. Und weil diese Lage beim
Verfeinern springt und nicht stetig zuläuft, gibt es keinen festen Vorfaktor $C$ in
$f_h = f_\text{exakt} + C h^{p}$ — die Voraussetzung der Richardson-Extrapolation ist
verletzt, und zwar nachweisbar und nicht bloß behauptet.

Die Rangfolge ist dabei kein Zufall und bestätigt die Deutung: Am stärksten betroffen ist
$c_L$, das ganz von der Lage der Ablösepunkte lebt; dann $c_D$; am wenigsten $\mathrm{St}$,
das vom Takt der Nachlaufinstabilität bestimmt wird und nicht von der Wandkontur. Und
$\Delta p$, eine Druckablesung in einigem Abstand von der Wand, liegt schon auf dem
gröbsten Netz nur 2,9 % daneben.

**Was daraus folgt — und was nicht.** Es folgt *nicht*, dass das Verfahren ungenau ist.
Die Netzreihe läuft in allen vier Größen monoton auf die Referenz zu und erreicht auf dem
feinsten Netz $+1{,}8\,\%$ ($c_D$), $-2{,}6\,\%$ ($c_L$), $-1{,}3\,\%$ ($\mathrm{St}$) und
$-1{,}0\,\%$ ($\Delta p$); die Strouhalzahl liegt **im** Referenzintervall. Es folgt, dass
der verbleibende Rest von der **Körperdarstellung** herrührt und nicht vom Netzwerkansatz.
Wer ihn beseitigen will, braucht keine feineren Zellen, sondern eine andere Darstellung des
Körpers: geschnittene Zellen mit Flächenanteilen oder ein körperangepasstes Netz. Das ist
der nächste Schritt und ausdrücklich noch nicht getan.

## 5.5 Die Sprachprobe: zwei getrennt aufgebaute Löser

Beide Wege lösen dasselbe Problem, werden aber unabhängig voneinander aufgebaut. Sie haben
**keinen gemeinsamen Rechenkern** — nur Geometrie, Stoffwerte, Randbedingungen und
Zeitschritt sind dieselben, weil sie zur Aufgabe gehören und nicht zum Verfahren.

| | Weg A — Netzwerk nach Kirchhoff | Weg B — klassische Finite Volumen |
|:--|:--|:--|
| Aufbau | jeder Zweig **einzeln** gestempelt | Flächenflüsse als **ganze Felder** |
| Darstellung | dünnbesetzte Knotenadmittanzmatrix | gar keine Matrix |
| Zeitschritt | Matrix-Vektor-Produkt | Feldoperation mit verschobenen Schnitten |
| Drucknetz | **direkt** gelöst (LU-Zerlegung) | **matrixfrei**, konjugierte Gradienten |
| Randwerte | als Matrixeinträge beziehungsweise Lastvektor | als Geisterschicht um das Feld |

Geprüft wird auf drei Ebenen, von innen nach außen:

**Ebene 1 — die einzelnen Operatoren.** Für zufällig belegte Felder muss die Matrixwirkung
mit der Feldoperation übereinstimmen.

| Operator | größte relative Abweichung | Programm |
|:--|--:|:--|
| Drucknetz, nackt | $2\cdot10^{-16}$ | `ns_loeser.py` |
| Drucknetz mit Störkörper und Ausstrom | $1{,}7\cdot10^{-16}$ | `drucknetz_koerper_probe.py` |
| Impulsoperator $u$, mit und ohne vorgegebenes Einströmprofil | $1{,}4\cdot10^{-16}$ | `impuls_matrix.py` |
| Impulsoperator $v$ | $2{,}5\cdot10^{-16}$ | `impuls_matrix.py` |

Das ist Rundungsgenauigkeit — mehr kann in Gleitkommaarithmetik nicht erreicht werden.

**Ebene 2 — die Eigenschaften, die ohne Rechnung gelten müssen.** Sie prüfen nicht die
Übereinstimmung zweier Programme, sondern die Übereinstimmung mit der Physik:

- Gleichförmige Strömung bleibt unverändert (keine Scheinkraft aus der Diskretisierung).
- Verschwindende Geschwindigkeit führt auf reine Diffusion, und die Impulsmatrix wird
  **symmetrisch** — das passive Widerstandsnetz.
- Konstanter Druck treibt im geschlossenen Netz **keinen** Strom; mit Ausstromrand treibt
  er genau den einen Zweig ans Bezugspotential, mit dem Leitwert $\Delta y/(\Delta x/2)$.
- $g_\text{hin} - g_\text{rück} = u_f A$ — die Unsymmetrie **ist** der Volumenstrom,
  geprüft über acht Zehnerpotenzen der Geschwindigkeit.

**Ebene 3 — der vollständige Löser.** Beide Wege werden vom selben Anfangszustand über
Hunderte bis Tausende von Zeitschritten gerechnet und Schritt für Schritt verglichen. Das
ist die eigentliche Probe: Einzelne Operatoren können übereinstimmen und der
zusammengesetzte Zeitschritt trotzdem auseinanderlaufen, weil Randbehandlung, Reihenfolge
und Krafterfassung dazwischenliegen.

| Netz | Zeitschritte | $\max|u_A-u_B|$ rel. | $\max|v_A-v_B|$ rel. | $\max|c_{D,A}-c_{D,B}|$ | $\max|c_{L,A}-c_{L,B}|$ |
|:--|--:|--:|--:|--:|--:|
| $110\times21$ | 192 | $1{,}4\cdot10^{-14}$ | $1{,}5\cdot10^{-14}$ | $8{,}4\cdot10^{-13}$ | $1{,}4\cdot10^{-12}$ |
| $220\times41$, bis in die Ablösung | 6000 | $5{,}2\cdot10^{-14}$ | $5{,}5\cdot10^{-14}$ | $3{,}5\cdot10^{-12}$ | $4{,}0\cdot10^{-12}$ |

Der zweite Lauf ist der aussagekräftige: Er läuft 6000 Zeitschritte lang durch die
**ausgebildete Wirbelablösung**, also durch den Teil, in dem kleine Unterschiede anwachsen
würden, wenn es welche gäbe. Die Abweichung wächst von $1{,}4\cdot10^{-14}$ auf
$5{,}2\cdot10^{-14}$ — knapp das Vierfache über das Dreißigfache an Schritten. Das ist das
Anwachsen angesammelter Rundungsfehler und **nicht** das Auseinanderlaufen zweier
verschiedener Verfahren; letzteres würde sich im Grenzzyklus exponentiell zeigen.

Zusätzlich werden beide Bahnen getrennt ausgewertet — eingeschwungener Teil, Maximum
darin, Strouhalzahl aus den Nulldurchgängen —, so dass nicht nur dieselben Felder, sondern
dieselben **Beiwerte** belegt sind.

**Dritte Übereinstimmung, ungeplant.** `wirbelstrasse_zweiwege.py` ist eine eigenständige
Nachbildung des Zeitschritts, unabhängig von `wirbelstrasse.py` geschrieben. Beide müssen
auf demselben Netz dieselben Beiwerte liefern — das ist eine Gegenprobe auf die
Zeitschrittführung selbst und nicht nur auf die Operatoren. Netz $220\times41$,
$t_\text{Ende}=8$ s:

| | $c_D^{\max}$ | $c_L^{\max}$ | $\mathrm{St}$ |
|:--|--:|--:|--:|
| `wirbelstrasse.py` | 3,751423999 | 0,610308625 | 0,252525253 |
| `wirbelstrasse_zweiwege.py`, Weg A | 3,751423999 | 0,610308625 | 0,252525253 |

Alle ausgegebenen Stellen stimmen. Beim ersten Vergleich standen 3,751424 gegen 3,751430 —
sechs gleiche Stellen und ein Rest, für den es eine naheliegende Erklärung gab (verschiedene
Auswertefenster, 12 s gegen 8 s). Naheliegend ist nicht nachgewiesen: Der Produktionslöser
wurde deshalb mit demselben Fenster noch einmal gerechnet, und der Unterschied verschwand
vollständig. Die Erklärung stimmte — aber sie wurde geprüft, nicht geglaubt.

Dieses Vorgehen ist aus dem `Feld_Labor` übernommen, wo Netzwerkaufbau nach Kirchhoff und
Elementaufbau nach Ritz bereits so nebeneinanderstehen.

### Was die Sprachprobe nicht zeigt

Sie zeigt, dass zwei Aufbauwege dieselbe **Diskretisierung** meinen. Sie zeigt **nicht**,
dass diese Diskretisierung die Navier-Stokes-Gleichungen richtig wiedergibt — zwei Wege
können einträchtig dasselbe falsch rechnen. Dafür sind die Stufen 1 bis 3 da, die gegen
geschlossene Lösungen und gegen gemessene Referenzwerte prüfen. Die Sprachprobe ergänzt
sie, sie ersetzt sie nicht.

## 5.6 Der unabhängige zweite Weg — Galerkin auf körperangepasstem Netz

### Warum die Sprachprobe nicht genügt

Das Maß für Unabhängigkeit ist nicht, ob zwei Wege verschieden *aussehen*, sondern ob sie
verschieden **scheitern** können. Weg A und Weg B teilen Zeitschritt, Randbehandlung,
Körperbehandlung, Kraftbilanz und den Zweigstempel. Wäre einer davon begrifflich falsch,
lieferten beide einträchtig dasselbe falsche Ergebnis — auf $5\cdot10^{-14}$ genau. Die
Sprachprobe beweist, dass der **Aufbau** fehlerfrei ist. Mehr nicht, und das ist nicht
wenig; aber es ist nicht das, was im `Feld_Labor` Kirchhoff und Ritz leisten.

### Galerkin, nicht Ritz

Ritz **minimiert ein Funktional**. Das setzt voraus, dass es eines gibt: Der Operator muss
symmetrisch und positiv definit sein. Für die Poisson-Gleichung ist das so — die Energie
ist das Funktional, ihr Minimum ist die Lösung. Deshalb steht im Feld-Labor Ritz.

Navier-Stokes hat kein solches Funktional. Der konvektive Term $(\mathbf{u}\cdot\nabla)
\mathbf{u}$ ist **nicht symmetrisch** und aus keinem Potential ableitbar; dazu macht die
Inkompressibilität schon aus dem Stokes-Problem einen **Sattelpunkt**, kein Minimum. Es
gibt nichts zu minimieren.

Galerkin braucht kein Funktional: Man multipliziert die Gleichung mit Ansatzfunktionen und
integriert — die schwache Form. Für symmetrische Probleme liefert Galerkin dieselben
Gleichungen wie Ritz (Ritz-Galerkin-Äquivalenz); beim Zähigkeitsterm hier ist das wörtlich
der Fall. **Galerkin ist nicht ein anderes Verfahren als Ritz, sondern die Fortsetzung
desselben Gedankens dorthin, wo die Symmetrie verlorengeht.**

### Was dieser Weg anders macht

| | Netzwerkweg | Galerkin-Weg |
|:--|:--|:--|
| Herleitung | Flussbilanz am Kontrollvolumen | schwache Form, Ansatzfunktionen |
| Unbekannte | Zellmittel und Flächenflüsse | Knotenwerte eines Polynomansatzes |
| Ansatz | zellweise konstant | $P_2$ für $\mathbf{u}$, $P_1$ für $p$ (Taylor-Hood) |
| Netz | kartesisch, Zylinder als **Treppe** | Dreiecke, Zylinderrand **auf Knoten** |
| Druck | eigenes Netz, Projektion | im selben Sattelpunktsystem |
| Kraft | Impulsbilanz der abgedeckten Zweige | **konsistent** aus der schwachen Form |

Taylor-Hood erfüllt die inf-sup-Bedingung und braucht deshalb **keinen abgestimmten
Stabilisierungsparameter** — genau der Grund, aus dem Abschnitt 4.2 gleichordnige
Verfahren mit PSPG ablehnt. Die Kraft kommt aus der schwachen Form selbst: Testet man die
Impulsgleichung mit einer Funktion, die auf dem Zylinder gleich eins und auf allen anderen
Rändern null ist, so steht auf der einen Seite das Randintegral der Spannungen. Dieser Weg
ist deutlich genauer als das Aufsummieren von Wandspannungen und hat mit der Impulsbilanz
des Netzwerkwegs nichts gemein.

Der Zeitschritt ist halbimplizit: Zähigkeit und Druck implizit, Konvektion explizit. Damit
ist die Systemmatrix **konstant** und wird nur einmal zerlegt (0,2 s); je Zeitschritt
bleibt ein Einsetzen, rund 21 ms bei 18 503 Unbekannten.

### Die Probe, die nicht „klein" verlangt, sondern null

Hagen-Poiseuille ist für diesen Weg eine besonders scharfe Probe: Die Parabel liegt im
$P_2$-Raum, der lineare Druck im $P_1$-Raum — die exakte Lösung liegt also **im
Ansatzraum** und muss bis auf Rundung herauskommen.

| Größe | gemessen |
|:--|--:|
| Geschwindigkeitsprofil | $4{,}9\cdot10^{-14}$ relativ |
| Querkomponente | $4{,}5\cdot10^{-14}$ |
| Druck, Abweichung von der Geraden | $1{,}0\cdot10^{-14}$ |
| Druckgefälle gegen $-8\nu U_m/H^2$ | auf zehn Stellen |

Genau diese Schärfe hat einen Fehler gefunden, den eine Probe mit der Schranke „klein"
durchgelassen hätte: Zwei Eckknoten am Ausstrom waren von der Haftbedingung ausgenommen
(Befund B15). Der Fehler betrug 3,3 %.

### Das Ergebnis

Netz mit 2097 Ecken, davon 81 auf dem Zylinderrand, 18 503 Unbekannte:

| $\Delta t$ | $c_D^{\max}$ | $c_L^{\max}$ | $\mathrm{St}$ |
|--:|--:|--:|--:|
| $5\cdot10^{-4}$ | 3,2456 | 1,0494 | 0,3008 |
| $2{,}5\cdot10^{-4}$ | 3,2364 | 1,0225 | 0,3005 |
| $1{,}25\cdot10^{-4}$ | **3,2318** | **1,0092** | **0,3002** |
| Referenz | 3,22 … 3,24 | 0,99 … 1,01 | 0,295 … 0,305 |

**Alle drei Größen liegen im Referenzintervall — direkt gerechnet, ohne Extrapolation.**

Die beobachtete Konvergenzordnung in der Zeit beträgt **1,00** ($c_D$) beziehungsweise
**1,02** ($c_L$). Das Zeitschema ist explizit-impliziter Euler, also erster Ordnung: Die
Messung trifft genau das, was das Verfahren verlangt. Der Zeitfehler ist damit nicht nur
klein, sondern **verstanden** — anders als die Treppenstreuung des Netzwerkwegs, die keiner
festen Ordnung folgt. Bei $\mathrm{St}$ liegt kein Ordnungswert vor: Die Nulldurchgänge
liegen auf dem Zeitraster, und die Unterschiede sind bereits an dessen Auflösung angelangt.

Die Netzreihe daneben (40, 81, 121 Knoten auf dem Zylinderrand, Zeitschritt jeweils
mitskaliert) gibt $c_D$ = 3,2054 / 3,2456 / 3,2446 und $\mathrm{St}$ = 0,2976 / 0,3008 /
0,3027 — **alle drei Stufen mit der Strouhalzahl im Intervall**, und zwischen der zweiten
und dritten Stufe ändert sich $c_D$ nur noch um 0,03 %.

### Was damit belegt ist

Zwei Wege, die **verschieden scheitern können** — Flussbilanz am Kontrollvolumen mit
Treppenkörper gegen schwache Form mit körperangepasstem Netz —, treffen dieselbe gemessene
Referenz und einander. Das ist die Absicherung, die der Netzwerkansatz für die Feldrechnung
längst hat und die der Strömung bis hierher fehlte.

Nebenbei bestätigt dieser Weg die Deutung aus Abschnitt 5.4: Ohne Treppe erreicht man die
Referenzintervalle mit 2097 Knoten, während der Netzwerkweg mit 144 320 Zellen bei
$+1{,}8\,\%$ stehenbleibt. Der verbliebene Fehler dort kam von der Körperdarstellung, nicht
vom Ansatz.

### Was ein Benchmark bescheinigt — und was nicht

Legt man die **Felder** beider Wege übereinander (phasenrichtig zugeordnet, gleiche
Auflösung, gleiche Auswertung), so stimmen sie bei weitem nicht so gut überein wie die
Beiwerte:

| | Effektivwert der Differenz |
|:--|--:|
| über das ganze Gebiet | 27 % des Feldes |
| nur im Nachlauf ($x>0{,}3$ m) | 20 % |
| nahe am Zylinder ($r<90$ mm) | 40 % |

Die größte Abweichung sitzt 55 mm vom Zylindermittelpunkt, also **5 mm vor der Wand** —
mitten in der Grenzschicht, die hier rund $D/\sqrt{\mathrm{Re}}=10$ mm dick ist. Dort kann
die Treppenkontur den Wandgradienten nicht darstellen.

Das ist kein Widerspruch zu den $1\ldots2\,\%$ in den Beiwerten, sondern eine Lehre:
**Integralgrößen konvergieren erheblich schneller als die Felder selbst.** Ein Benchmark
wie Schäfer-Turek bescheinigt $c_D$, $c_L$, $\mathrm{St}$ und $\Delta p$ — er bescheinigt
**nicht** das Feld an jeder Stelle. Wer ein Verfahren für eine örtliche Größe einsetzen
will — Wandschubspannung, Erosion, örtlicher Wärmeübergang —, muss genau diese Größe
prüfen und darf sich nicht auf die Kraftbeiwerte verlassen.

## 5.6a Beide Wege im Browser — gerechnet, nicht abgespielt

Bis hierher lief der Galerkin-Weg in Python, und im Browser-Labor standen nur seine
Ergebnisse. Auf ausdrückliche Forderung des Auftraggebers — *„die html soll immer beide
wege rechnen können"* — rechnet die Seite jetzt **beide** Wege selbst.

Dafür waren drei Stücke nötig:

| Stück | was es tut |
|:--|:--|
| `netz_browser.js` | Randpunkte auf Kanal und Zylinderkreis, Innenpunkte nach dem Poisson-Scheiben-Verfahren mit ortsabhängigem Abstand |
| Triangle (Shewchuk) | die Delaunay-Triangulierung, als WebAssembly eingebettet — derselbe Vernetzer, den auch FEMM benutzt |
| `kern_fem.js` | P2/P1-Aufbau, Konvektion als Vektor, Dirichlet-Elimination, Sattelsystem, konsistente Kraft |

### Der Gleichungslöser — und warum er ein anderer sein darf

In Python wird das Sattelsystem **direkt zerlegt**: Die Matrix ist konstant, also fällt
die Zerlegung einmal an, und je Zeitschritt bleibt ein Einsetzen. Eine dünnbesetzte
LU-Zerlegung mit Pivotsuche im Browser nachzubauen wäre ein eigenes Vorhaben. Stattdessen
wird dort **MINRES** benutzt — das Verfahren von Paige und Saunders für symmetrische,
indefinite Systeme — mit einem Blockvorkonditionierer:

$$\hat P \;=\; \mathrm{diag}\!\left(\hat A,\; \hat A,\; \hat S\right),
\qquad
\hat S^{-1} \;\approx\; \nu\,M_p^{-1} \;+\; \frac{1}{\Delta t}\,K_p^{-1}$$

$\hat A$ ist die unvollständige Cholesky-Zerlegung von $A = M/\Delta t + \nu K$;
$M_p$ und $K_p$ sind Massenmatrix und Laplace-Operator im **Druckraum** P1. Die Form von
$\hat S$ stammt von Cahouet und Chabard und ist nicht geraten: Bei großem Zeitschritt
verhält sich der Schur-Anteil wie der zähe Term $\nu/M$, bei kleinem wie $\Delta t\,K_p$;
die Summe deckt beide Enden ab.

Damit das System überhaupt symmetrisch ist, wird die dritte Blockzeile mit $-1$
durchmultipliziert:

$$\begin{bmatrix} A & 0 & -B_x^{\mathsf T}\\ 0 & A & -B_y^{\mathsf T}\\
-B_x & -B_y & 0\end{bmatrix}$$

Dieselbe Gleichung, dieselbe Lösung — nur so geschrieben, dass MINRES anwendbar ist.
**Ein anderer Gleichungslöser ist kein anderer Ansatz.** Der Ansatz ist die schwache Form
mit Taylor-Hood-Elementen, und der ist in beiden Fassungen derselbe; das belegen die
45 Prüfposten von `pruefe_kern_fem.js`, die den Browserkern **knotenweise** gegen
`fem_galerkin.py` halten.

### Was dabei gemessen wurde

| Probe | Browser gegen Python |
|:--|--:|
| Formfunktionen, Zerlegung der Eins | $0{,}00$ |
| Matrizen $M$, $K$, $B_x$, $B_y$ | $1{,}8\cdot10^{-14}$ |
| Massenmatrix = Gebietsfläche | $0{,}894500000000$ exakt |
| Konvektionsvektor, knotenweise | $1{,}1\cdot10^{-15}$ |
| $c_D$ über 40 Zeitschritte | $9{,}5\cdot10^{-11}$ |
| $c_L$ über 40 Zeitschritte | $1{,}6\cdot10^{-9}$ |
| Felder $u_x$, $u_y$, $p$ knotenweise | $3{,}5\cdot10^{-12}\ldots7{,}6\cdot10^{-10}$ |
| **Hagen-Poiseuille** | $1{,}04\cdot10^{-12}$ |

Die letzte Zeile ist die schärfste: Im nackten Kanal ist die Lösung geschlossen bekannt,
und ein P2-Ansatz gibt eine Parabel **exakt** wieder — nicht genähert. Der verbleibende
Rest ist allein die Abbruchschranke des iterativen Lösers; Python mit direkter Zerlegung
kommt auf $1{,}6\cdot10^{-14}$.

### Die Grenze der Browserrechnung — offen benannt

Der Zeitschritt ist durch die **explizite Konvektion** begrenzt. Gemessen (Befund B20):
Bei $\mathrm{CFL} = 0{,}31$ platzt die Rechnung bei $t = 0{,}372\,$s, bei
$\mathrm{CFL} = 0{,}15$ läuft sie durch. Der Kern bestimmt seinen Zeitschritt deshalb
selbst aus dem Netz, mit $C = 0{,}15$ — ein Sicherheitsabstand von zwei gegen das
Gemessene.

Daraus folgt, was die Seite **nicht** leistet: Die feinen Netze des Python-Laufs
(18 503 Unbekannte, $\Delta t = 5\cdot10^{-4}\,$s, 24 000 Zeitschritte) sind im Browser
nicht in Minuten zu rechnen. Die Beiwerte dort liegen deshalb **nicht** auf dem
Benchmarkwert. Was die Seite zeigt, ist der **Vergleich der beiden Wege bei gleichem
Knotenabstand** — und das ist die Aussage, um die es geht. Die Benchmarkzahlen stehen in
§5.4 und §5.6 und kommen aus den vollen Läufen.

## 5.7 Die vierte Größe: die Druckdifferenz

Die Beiwerte $c_D$ und $c_L$ entstehen hier aus dem **Impulsstrom** — gezählt wird, wieviel
Impuls das Verfahren je Zeitschritt aus den abgedeckten Zweigen nimmt. Das ist ein
einziger Auswerteweg, und ein Verfahren, das ihn falsch auswertet, könnte zwei Zahlen
liefern, die zueinander passen und beide falsch sind.

Schäfer und Turek geben deshalb eine vierte Größe an, die auf einem **anderen** Weg
entsteht:

$$\Delta p = p(0{,}15;\,0{,}20) - p(0{,}25;\,0{,}20)
\quad\text{bei}\quad t_0 + \tfrac{1}{2f},\qquad
t_0 = \text{Zeitpunkt des größten } c_L ,$$

also eine reine **Druckablesung** an zwei Punkten, ohne jede Integration und ohne
Impulsbilanz. Referenz: 2,46 … 2,50.

Die beiden Punkte liegen genau auf der Zylinderoberfläche; dort liegt auf dem kartesischen
Gitter keine Zellmitte, und die nächstgelegene Zelle kann bereits im Körper liegen, wo der
Druck bedeutungslos ist. Deshalb wird aus den beiden nächsten **Fluidzellen** auf der
Achse linear auf die Oberfläche extrapoliert — ein Fehler erster Ordnung im Zellabstand,
der mit dem Netz kleiner wird und in der Reihe sichtbar ist. Geschätzt wird nichts.

| Zellen über $D$ | $\Delta p$ | Abweichung von 2,48 |
|--:|--:|--:|
| 10 | 2,4075 | −2,9 % |
| 20 | 2,4343 | −1,8 % |

Bemerkenswert ist, **wie** gut diese Größe schon auf dem groben Netz liegt: Während $c_L$
dort noch 39 % danebenliegt, sind es bei $\Delta p$ 2,9 %. Das ist kein Zufall, sondern
sagt etwas über die Fehlerquelle aus. Die Treppenkontur verfälscht vor allem die
**Wandschubspannung** und die Lage der Ablösepunkte — beides geht voll in die integrierten
Kraftbeiwerte ein. Der Druck in einigem Abstand von der Wand ist davon weit weniger
betroffen. Die Netzreihe prüft also nicht einfach „das Verfahren", sondern zeigt, welcher
Teil der Diskretisierung den Fehler trägt.

Die Periodendauer ist nebenbei eine Gegenprobe zur Strouhalzahl: gemessen 356,0 ms auf dem
Netz mit 20 Zellen über dem Durchmesser, und $\mathrm{St} = D/(\bar u\,T) = 0{,}1/0{,}356
= 0{,}2809$ — genau der Wert, den die Auswertung der Zeitreihe liefert.

## 5.8 Was nicht bewiesen werden kann

Ein turbulenter Fall — etwa die Strömung in einem Lüfterrad — verlangt ein
Turbulenzmodell. Ein Turbulenzmodell ist eine **Modellannahme**, keine Herleitung; es
enthält abgestimmte Konstanten. Rechnungen damit lassen sich gegen Messungen abgleichen,
aber nicht beweisen. Deshalb wird der Nachweis des Verfahrens am laminaren,
zweidimensionalen Fall geführt und die turbulente Anwendung davon getrennt gehalten.

Für die spätere Anwendung am Lüfterrad liegt der Datensatz **FAN-01** vor (EAA-Benchmark,
Niederdruck-Axiallüfter, Zenodo, CC-BY-4.0): gemessene Kennlinie mit 20 Betriebspunkten,
Volumenstrom 718 bis 7483 m³/h, Druckerhöhung 253 bis 7 Pa, dazu LDA-Geschwindigkeitsfelder
und die Geometrie. Das ist eine Messung und damit eine zulässige Referenz — für die
Anwendung, nicht für das Verfahren.

## 5.9 Das Prüfprotokoll

Jede Stufe ist ein eigenständiges Programm, das sich **selbst** beurteilt und mit
Rückgabewert 0 oder 1 endet. `pruefprotokoll.py` ruft sie der Reihe nach auf und sammelt
nur ein — es rechnet nichts nach, sonst wäre es eine zweite Quelle der Wahrheit.

| Stufe | Was geprüft wird | Woran gemessen |
|:--|:--|:--|
| S01 | Zweigstempel nach Scharfetter-Gummel | geschlossene Lösung der 1D-Konvektion-Diffusion |
| S02 | Knotenregel auf dem versetzten Gitter | divergenzfreies Feld $u=x,\;v=-y$ |
| S03 | Drucknetz: Matrix gegen Feldoperation | die beiden Aufbauwege gegeneinander |
| S04 | dasselbe mit Störkörper und Ausstromrand | dazu drei Eigenschaften ohne Rechnung |
| S05 | ebener Kanal | Hagen-Poiseuille, exakte Lösung, Ordnung $h^2$ |
| S06 | Konvektionsstempel | sieben Eigenschaften, die ohne Rechnung gelten |
| S07 | Impulsoperator: Matrix gegen Feldoperation | dazu Symmetrie im Grenzfall $u\to0$ |
| S08 | Knotenregel am **Ende** jedes Zeitschritts | Rundungsschwelle |
| S09 | Wirbelstärke im Bild | fünf Felder mit bekanntem $\omega$, zwei davon ortsabhängig |
| S10 | vollständiger Löser, beide Wege, kurz | die beiden Wege gegeneinander |
| S11 | vollständiger Löser, beide Wege, bis in die Ablösung | dieselben, über 6000 Zeitschritte |
| S12 | Netzreihe | läuft jede Größe auf die Referenz zu? |
| S13 | Richardson-Extrapolation | liegt der Grenzwert im Referenzintervall? |
| S14 | Druckdifferenz $\Delta p$ | kommt sie der Referenz näher? |

<!-- PROTOKOLL-ERGEBNIS -->

**Ergebnis vom 31.10.2026, Umfang „ganz": 13 von 14 Stufen richtig.**

| Stufe | Urteil | Beleg |
|:--|:--|:--|
| S01 | richtig | größter Fehler über alle Fälle $3{,}4\cdot10^{-15}$ |
| S02 | richtig | Knotenregel auf dem divergenzfreien Feld |
| S03 | richtig | beide Wege bauen denselben Operator |
| S04 | richtig | 6 von 6 |
| S05 | richtig | feinstes Netz $3{,}05\cdot10^{-5}$, Ordnung 3× exakt 2,00 |
| S06 | richtig | 7 von 7 |
| S07 | richtig | 12 von 12 |
| S08 | richtig | 2 von 2 |
| S09 | richtig | 5 von 5 |
| S10 | richtig | 4 von 4 |
| S11 | richtig | 4 von 4, 6000 Zeitschritte (1354 s) |
| S12 | richtig | jede Größe kommt mit jeder Verfeinerung näher |
| **S13** | **falsch** | Grenzwerte außerhalb; Ordnungen 0,85 … 1,70 unverträglich |
| S14 | richtig | $\Delta p$: $-2{,}9\,\% \to -1{,}8\,\% \to -1{,}0\,\%$ |

Die eine falsche Stufe ist kein Programmfehler, sondern ein **Messergebnis**: Die
Voraussetzung der Richardson-Extrapolation ist auf diesem Fall nicht erfüllt, und das
Protokoll weist das nach, statt die Hürde zu senken. Die Erklärung dazu steht als eigener
Messblock unmittelbar darunter — die Streuung durch die Lage der Treppe (§5.4).

Ein Hinweis zur Laufzeit: S11 braucht 1354 s, weil dort in jedem einzelnen der 6000
Zeitschritte **beide** Impulsmatrizen neu gestempelt werden — Zweig für Zweig, in Python.
Das ist langsam und soll es sein: Weg A ist als Nachweis gebaut, nicht als Produktionsweg.
Gerechnet wird mit `wirbelstrasse.py`, das denselben Zeitschritt in Feldoperationen macht
und auf demselben Netz dieselben neun Stellen liefert.

# 6 Beteiligte Dateien

Alle Programme liegen in `stroemung/`. Jedes beurteilt sich **selbst** und endet mit
Rückgabewert 0 oder 1; `pruefprotokoll.py` ruft sie der Reihe nach auf und sammelt nur ein.

> **Die Programme liegen im Wortlaut bei.** Im Browser-Labor stehen sie im Reiter
> **Programme** — alle 34 Dateien aus `stroemung/` und `labor/`, jede mit Zeilenzahl,
> Größe und sha256-Prüfsumme. Damit lässt sich nachrechnen, worauf hier verwiesen wird,
> und nachweisen, dass es dasselbe ist wie das Gelieferte. Auf Einwand von Fred
> (01.11.2026): „Alle Prüfprotokolle verweisen darauf, nachrechnen kann ich sie ohne die
> Quellen nicht." Eine Seite, die auf Programme verweist, die ihr nicht beiliegen, ist
> für den Leser eine Sackgasse — er kann glauben oder bleiben lassen, aber nicht
> nachrechnen.

| Datei | Fassung | Aufgabe | Stufe |
|:--|:--|:--|:--|
| `sg_stempel.py` | 1.0 | Scharfetter-Gummel-Zweigstempel, $B(x)$ in vier Zahlenbereichen | S01 |
| `gitter.py` | 1.0 | versetztes Gitter, Knotenregel je Zelle | S02 |
| `ns_loeser.py` | 1.2 | Drucknetz als Matrix **und** als Feldoperation, Zähigkeitsterm | S03 |
| `drucknetz_koerper_probe.py` | 1.0 | dasselbe mit Störkörper und Ausstromrand | S04 |
| `poiseuille.py` | 1.0 | ebener Kanal gegen die geschlossene Lösung, Ordnung $h^2$ | S05 |
| `konvektion.py` | 1.1 | Konvektion und Zähigkeit als Feldoperation (Weg B) | S06 |
| `impuls_matrix.py` | 1.1 | derselbe Operator als Knotenadmittanzmatrix (Weg A) | S07 |
| `knotenregel_probe.py` | 1.0 | Divergenz am **Ende** jedes vollständigen Zeitschritts | S08 |
| `bild.py` | 1.2 | Wirbelstärke auf der Gitterecke, gegen bekannte Felder geprüft | S09 |
| `druckdifferenz.py` | 1.0 | die vierte Referenzgröße $\Delta p$ an zwei Punkten | S17 |
| `netz_fem.py` | 1.0 | körperangepasstes Dreiecksnetz, Zylinderrand auf Knoten | S13 |
| `fem_galerkin.py` | 1.1 | Taylor-Hood-Ansatzraum, Matrizen, Hagen-Poiseuille exakt | S14 |
| `fem_wirbelstrasse.py` | 1.0 | **der unabhängige zweite Weg** am Benchmark | S18 |
| `lbm.py` | 1.0 | Lattice-Boltzmann, zwei analytische Proben (nicht weiterverfolgt) | S12 |
| `wirbelstrasse.py` | 1.4 | der Benchmarkfall, Produktionslöser | — |
| `wirbelstrasse_zweiwege.py` | 1.2 | derselbe Fall, beide Wege parallel gerechnet | S10/S11 |
| `pruefprotokoll.py` | 1.1 | alle Stufen in einem Aufruf | — |
| `tafel_netzreihe.py` | 1.0 | erzeugt die Tafeln in §5.4 aus den Ergebnisdateien | — |
| `treppenstreuung.py` | 1.0 | Streuung der Beiwerte durch die Lage der Treppe | Erklärung zu S13 |
| `zeitaufteilung.py` | 1.0 | wo die Zeit eines Zeitschritts bleibt (Jetson-Entscheidung) | — |
| `schaubild.py` | 1.2 | Feld, Ablösetakt und Netzreihe auf einem Blatt | — |

Das Browser-Labor liegt in `labor/`. Seine **beiden** Rechenkerne sind zweite Umsetzungen
derselben Verfahren und deshalb gegen die Python-Module geprüft — eine zweite Umsetzung
ohne Prüfung wäre keine zweite Meinung, sondern eine zweite Fehlerquelle. Seit Fassung 1.1
des Prüfprotokolls sind sie als Stufen **S19a** und **S19b** Teil desselben Laufs: Ein
Weg, der nur im Browser läuft und im Protokoll nicht vorkommt, wäre ein zweiter,
ungeprüfter Stand desselben Verfahrens.

| Datei | Fassung | Aufgabe |
|:--|:--|:--|
| `labor/kern_cfd.js` | 1.1 | der Rechenkern des NETZWERKWEGS im Browser |
| `labor/kern_fem.js` | 1.0 | der Rechenkern des GALERKIN-WEGS im Browser: P2/P1, MINRES, konsistente Kraft |
| `labor/netz_browser.js` | 1.0 | Vernetzer: Randpunkte, Poisson-Streuung, Delaunay über Triangle |
| `labor/triangle/` | — | Triangle (Shewchuk) als WebAssembly, mit Herkunfts- und Lizenzblatt |
| `labor/vergleichsdaten.py` | 1.0 | Vergleichsdaten für den Netzwerkkern aus den Python-Modulen |
| `labor/vergleichsdaten_fem.py` | 1.0 | Vergleichsdaten für den Galerkin-Kern, mit Netz |
| `labor/pruefe_kern.js` | 1.0 | hält den Netzwerkkern dagegen — **45 Posten** (Stufe S19a) |
| `labor/pruefe_kern_fem.js` | 1.0 | hält den Galerkin-Kern dagegen — **47 Posten**, knotenweise (S19b) |
| `labor/labor.js` | 1.0 | die Bedienung: fünf Vergleiche, Zeitverlauf, Zeichnen |
| `labor/vorlage.html` | 1.0 | Gerüst und Gestaltung, sieben Reiter |
| `labor/surrogat.py` | 1.1 | lernt das Ersatzmodell an (POD) und misst seinen Fehler |
| `labor/harmonisch.py` | 1.0 | das harmonische Zeitgesetz, gegen DMD gemessen |
| `labor/baue_labor.py` | 1.3 | baut die eine Datei; die Fassungsnummer steht nur hier |
| `labor/pruefe_labor.js` | 1.0 | Knopfprobe im echten Browser — **54 Posten** |

Daneben:

| Datei | Aufgabe |
|:--|:--|
| `referenz/characteristic/`, `referenz/lda/` | FAN-01, gemessene Lüfterkennlinie und Geschwindigkeitsfelder |
| `Feld_Labor_2.6.html` | Poisson-Netz, beide Aufbauwege — Grundlage des Drucknetzes |
| `BEFUNDE.md` | jeder gefundene Fehler mit der Lücke, die ihn durchgelassen hat |
| `PRUEFPLAN.md` | die Kriterien, abgeleitet aus dem Wortlaut des Auftraggebers |

# Danksagung

Die körperangepassten Dreiecksnetze dieser Arbeit wurden mit **Triangle** erzeugt, dem
Netzgenerator und Delaunay-Triangulierer von **Jonathan Richard Shewchuk** (Carnegie
Mellon University). Mit Dank an den Autor. Triangle ist frei verfügbar, aber
urheberrechtlich geschützt und darf ohne Lizenz des Autors nicht verkauft oder in
kommerzielle Erzeugnisse aufgenommen werden. Die WebAssembly-Hülle `triangle-wasm` 1.0.0
von Bruno Imbrizi steht unter MIT. Einzelheiten in `labor/triangle/HERKUNFT.md`.

# Quelltexte

Die Programme dieser Arbeit sind Teil des Nachweises und liegen deshalb im Wortlaut vor,
nicht als Verweis:

| wo | was |
|:--|:--|
| Reiter **Programme** im Browser-Labor | alle 34 Dateien, mit Zeilenzahl, Größe und sha256 |
| `stroemung/` | 22 Python-Programme: Löser, Proben, Prüfprotokoll |
| `labor/` | 12 Dateien: beide Browserkerne, Vernetzer, Vergleichsdaten, Proben, Bauskript |

Die Prüfsummen in der Seite erlauben den Abgleich mit den ausgelieferten Dateien: Was
dort steht, ist zeichengleich mit dem, was gerechnet hat.

# Schrifttum

- D. L. Scharfetter, H. K. Gummel: *Large-signal analysis of a silicon Read diode
  oscillator.* IEEE Trans. Electron Devices 16 (1969) 64–77.
- F. H. Harlow, J. E. Welch: *Numerical calculation of time-dependent viscous
  incompressible flow of fluid with free surface.* Phys. Fluids 8 (1965) 2182–2189.
- J. R. Shewchuk: *Triangle: Engineering a 2D Quality Mesh Generator and Delaunay
  Triangulator.* In: Applied Computational Geometry: Towards Geometric Engineering
  (M. C. Lin, D. Manocha, Hrsg.), Lecture Notes in Computer Science 1148, S. 203–222,
  Springer-Verlag, Berlin, Mai 1996.
- J. R. Shewchuk: *Delaunay Refinement Algorithms for Triangular Mesh Generation.*
  Computational Geometry: Theory and Applications 22(1–3):21–74, Mai 2002.
- J. Cahouet, J.-P. Chabard: *Some fast 3D finite element solvers for the generalized
  Stokes problem.* International Journal for Numerical Methods in Fluids 8 (1988)
  869–895.
- C. C. Paige, M. A. Saunders: *Solution of sparse indefinite systems of linear
  equations.* SIAM Journal on Numerical Analysis 12 (1975) 617–629.
- M. Schäfer, S. Turek: *Benchmark computations of laminar flow around a cylinder.*
  In: Flow Simulation with High-Performance Computers II, Vieweg 1996, 547–566.
- G. I. Taylor: *Dispersion of soluble matter in solvent flowing slowly through a tube.*
  Proc. R. Soc. London A 219 (1953) 186–203.
- R. Aris: *On the dispersion of a solute in a fluid flowing through a tube.*
  Proc. R. Soc. London A 235 (1956) 67–77.
- S. Schoder, F. Czwielong: *Dataset FAN-01: Revisiting the EAA Benchmark for a
  low-pressure axial fan.* arXiv:2211.12014.
