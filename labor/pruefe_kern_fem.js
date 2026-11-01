/* pruefe_kern_fem.js -- der Galerkin-Kern im Browser gegen die Python-Module.
 *
 * Aufruf:  node pruefe_kern_fem.js
 *
 * Geprueft wird gegen vergleichsdaten_fem.json, erzeugt von vergleichsdaten_fem.py aus
 * fem_galerkin.py und fem_wirbelstrasse.py. Das Netz kommt aus derselben Datei -- der
 * Vernetzer wird getrennt geprueft, denn er streut seine Innenpunkte mit einem eigenen
 * Zufallsgenerator und kann gar nicht dasselbe Netz liefern. Das Netz ist Eingabe,
 * nicht Verfahren.
 */
const fs = require('fs');
const path = require('path');
const K = require('./kern_fem.js');

const D = JSON.parse(fs.readFileSync(path.join(__dirname, 'vergleichsdaten_fem.json'), 'utf8'));
const urteile = [];
const U = (text, ok, beleg) => urteile.push([text, !!ok, beleg]);
const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-300);

function vergleicheProbe(name, feld, p, schranke) {
  let schlimmste = 0, wo = -1;
  for (let i = 0; i < p.stellen.length; i++) {
    const d = Math.abs(feld[p.stellen[i]] - p.werte[i]) /
              Math.max(Math.abs(p.werte[i]), 1e-12);
    if (d > schlimmste) { schlimmste = d; wo = p.stellen[i]; }
  }
  let s = 0, q = 0;
  for (let i = 0; i < feld.length; i++) { s += feld[i]; q += feld[i]*feld[i]; }
  const ds = rel(s, p.summe), dq = rel(q, p.quadratsumme);
  const schlimm = Math.max(schlimmste, ds, dq);
  U(name, feld.length === p.n && schlimm < schranke,
    `${feld.length} Werte, groesste Abweichung ${schlimm.toExponential(2)}` +
    (wo >= 0 ? ` (Stelle ${wo})` : ''));
  return schlimm;
}

// ---------------------------------------------------------------- 1. Kennzahlen
U('Kennzahlen des Falls stimmen ueberein',
  K.F.LX === D.fall.LX && K.F.LY === D.fall.LY && K.F.DZ === D.fall.DZ &&
  K.F.NU === D.fall.NU && Math.abs(K.F.U_QUER - D.fall.U_QUER) < 1e-15,
  `LX ${K.F.LX}, DZ ${K.F.DZ}, nu ${K.F.NU}, Re ${K.F.U_QUER*K.F.DZ/K.F.NU}`);

const P = Float64Array.from(D.netz.P);
const T = Int32Array.from(D.netz.T);

// ---------------------------------------------------------------- 2. P2-Raum
const R2 = K.p2Raum(P, T);
U('P2-Raum: Knotenzahl wie in fem_galerkin.py',
  R2.n2 === D.p2raum.n2 && R2.nv === D.p2raum.nv && R2.nkanten === D.p2raum.nkanten,
  `${R2.nv} Ecken + ${R2.nkanten} Kantenmitten = ${R2.n2}`);
{
  // Die Kantennummerierung darf abweichen (andere Reihenfolge der Map), die GEOMETRIE
  // nicht: dieselben Punkte, nur vielleicht anders sortiert. Geprueft wird deshalb die
  // Punktmenge, nicht die Reihenfolge.
  const schl = s => Math.round(s * 1e12) / 1e12;
  const hier = new Set();
  for (let i = 0; i < R2.n2; i++) hier.add(schl(R2.XY[2*i]) + '|' + schl(R2.XY[2*i+1]));
  let fehlen = 0;
  for (let i = 0; i < D.p2raum.XY_probe.stellen.length; i++) { /* Probe unten */ }
  U('P2-Knoten sind eindeutig (keine doppelte Kantenmitte)', hier.size === R2.n2,
    `${hier.size} verschiedene Orte bei ${R2.n2} Knoten`);
}
{
  /* Die Kantenmitten werden in derselben Reihenfolge durchnummeriert wie in
   * fem_galerkin.py -- dieselbe Schleife ueber die Dreiecke, dieselben drei Kantenpaare,
   * dieselbe Einfuegereihenfolge. Damit ist E6 ZEICHENGLEICH, und es kann knotenweise
   * verglichen werden statt nur ueber Kennzahlen. Faellt das hier durch, sind alle
   * folgenden knotenweisen Proben wertlos -- deshalb steht es vor ihnen. */
  let gleich = true;
  for (let k = 0; k < D.p2raum.E6_erste.length; k++)
    if (R2.E6[k] !== D.p2raum.E6_erste[k]) gleich = false;
  U('E6 der ersten fuenf Dreiecke zeichengleich (auch die Kantenmitten)', gleich,
    gleich ? 'Nummerierung stimmt ueberein' : 'Reihenfolge weicht ab');
  vergleicheProbe('Lage der P2-Knoten', R2.XY, D.p2raum.XY_probe, 1e-14);
}

// ---------------------------------------------------------------- 3. Formfunktionen
{
  let schlimm = 0;
  D.formfunktionen.punkte.forEach((L, i) => {
    const w = K.p2Werte(L);
    for (let j = 0; j < 6; j++)
      schlimm = Math.max(schlimm, Math.abs(w[j] - D.formfunktionen.werte[i][j]));
  });
  U('P2-Formfunktionen an vier Punkten', schlimm < 1e-14,
    `groesste Abweichung ${schlimm.toExponential(2)}`);
  // Zerlegung der Eins -- von Hand nachrechenbar, unabhaengig von der Vergleichsdatei
  let eins = 0;
  for (const L of D.formfunktionen.punkte) {
    const w = K.p2Werte(L);
    eins = Math.max(eins, Math.abs(w.reduce((a, b) => a + b, 0) - 1));
  }
  U('Formfunktionen summieren sich zu eins (Zerlegung der Eins)', eins < 1e-15,
    `groesste Abweichung ${eins.toExponential(2)}`);
}

// ---------------------------------------------------------------- 4. Gradienten
const lg = K.lambdaGradienten(P, T);
{
  let schlimm = 0;
  for (let k = 0; k < 18; k++)
    schlimm = Math.max(schlimm, rel(lg.g[k], D.lambda.g_erste[k]));
  U('Gradienten der baryzentrischen Koordinaten, erste drei Dreiecke', schlimm < 1e-12,
    `groesste Abweichung ${schlimm.toExponential(2)}`);
  vergleicheProbe('Dreiecksflaechen', lg.A, D.lambda.A_probe, 1e-12);
  // Gegenprobe ohne Vergleichsdatei: Summe der Gradienten muss null sein
  let null3 = 0;
  for (let e = 0; e < lg.A.length; e++) {
    null3 = Math.max(null3, Math.abs(lg.g[6*e] + lg.g[6*e+2] + lg.g[6*e+4]),
                            Math.abs(lg.g[6*e+1] + lg.g[6*e+3] + lg.g[6*e+5]));
  }
  U('Summe der drei Gradienten ist null (sie muessen sich aufheben)', null3 < 1e-9,
    `groesster Rest ${null3.toExponential(2)}`);
}
{
  let schlimm = 0;
  const g0 = [[lg.g[0], lg.g[1]], [lg.g[2], lg.g[3]], [lg.g[4], lg.g[5]]];
  D.ableitungen.werte.forEach((soll, i) => {
    const d = K.p2Ableitungen(D.formfunktionen.punkte[i], g0);
    for (let j = 0; j < 6; j++) {
      schlimm = Math.max(schlimm, rel(d[j][0], soll[2*j]), rel(d[j][1], soll[2*j+1]));
    }
  });
  U('Ableitungen der P2-Formfunktionen im ersten Element', schlimm < 1e-11,
    `groesste Abweichung ${schlimm.toExponential(2)}`);
}

// ---------------------------------------------------------------- 5. Die Matrizen
const S = K.zusammenstellen(P, T, R2);
for (const name of ['M', 'K', 'Bx', 'By']) {
  const soll = D.matrizen[name], ist = S[name];
  U(`Matrix ${name}: Groesse und Besetzung`,
    ist.nz === soll.zeilen && ist.ns === soll.spalten && ist.zv.length === soll.besetzt,
    `${ist.nz}x${ist.ns}, ${ist.zv.length} besetzt (Soll ${soll.besetzt})`);
  /* Befund 01.11.2026: Hier stand rel(Summe, Summe). Bei K, Bx und By ist die
   * Gesamtsumme aber VON NATUR AUS NULL -- genau das sagt „vernichtet das konstante
   * Feld". Der Vergleich teilte also durch Rundungsrauschen und meldete 0,56 und 2,2
   * als Abweichung, obwohl der Kern auf 1e-14 stimmte. Ein Mass, das durch eine Null
   * teilt, ist kein Mass. Gemessen wird jetzt am MASSSTAB DER MATRIX, der Wurzel aus
   * der Quadratsumme -- die ist nie null. */
  const skala = Math.sqrt(Math.max(soll.quadratsumme, 1e-300));
  const ds = Math.abs(K.summe(ist) - soll.summe) / skala;
  const dq = rel(K.quadratsumme(ist), soll.quadratsumme);
  U(`Matrix ${name}: Summe und Quadratsumme`, Math.max(ds, dq) < 1e-11,
    `Summe ${ds.toExponential(2)} (am Massstab ${skala.toExponential(2)}), ` +
    `Quadratsumme ${dq.toExponential(2)}`);
}
U('Massenmatrix summiert sich zur Gebietsflaeche',
  rel(K.summe(S.M), D.netz.flaeche) < 1e-12,
  `${K.summe(S.M).toFixed(12)} gegen ${D.netz.flaeche.toFixed(12)} m^2`);
{
  const eins = new Float64Array(S.n2).fill(1);
  const ke = K.matMal(S.K, eins);
  let m = 0; for (let i = 0; i < ke.length; i++) m = Math.max(m, Math.abs(ke[i]));
  U('Steifigkeit vernichtet das konstante Feld', m < 1e-10,
    `groesster Rest ${m.toExponential(2)}`);
  const bx = K.matMal(S.Bx, eins), by = K.matMal(S.By, eins);
  let mb = 0;
  for (let i = 0; i < bx.length; i++) mb = Math.max(mb, Math.abs(bx[i]), Math.abs(by[i]));
  U('Druckkopplung vernichtet das konstante Feld', mb < 1e-10,
    `groesster Rest ${mb.toExponential(2)}`);
}

// ---------------------------------------------------------------- 6. Die Raender
const RB = K.raender(R2.XY, R2.n2, true);
for (const n of ['ein', 'wand', 'zyl', 'aus', 'fest']) {
  let z = 0; for (let i = 0; i < R2.n2; i++) z += RB[n][i];
  U(`Randknoten „${n}"`, z === D.raender[n], `${z} Knoten (Soll ${D.raender[n]})`);
}
{
  // Die Haftbedingung darf am Ausstrom NICHT wegfallen (Befund B15)
  let eckenFrei = 0;
  for (let i = 0; i < R2.n2; i++) if (RB.aus[i] && RB.wand[i] && !RB.fest[i]) eckenFrei++;
  U('Eckknoten am Ausstrom behalten die Haftbedingung (Befund B15)', eckenFrei === 0,
    `${eckenFrei} freie Eckknoten`);
}

// ---------------------------------------------------------------- 7. Die Konvektion
{
  const ux = new Float64Array(R2.n2), uy = new Float64Array(R2.n2);
  for (let i = 0; i < R2.n2; i++) {
    const x = R2.XY[2*i], y = R2.XY[2*i+1];
    ux[i] = Math.sin(3*x) * Math.cos(5*y);
    uy[i] = Math.cos(2*x) * Math.sin(4*y);
  }
  const kv = K.konvektion(S, ux, uy);
  vergleicheProbe('Konvektionsvektor fx, knotenweise', kv.fx, D.konvektion.fx, 1e-11);
  vergleicheProbe('Konvektionsvektor fy, knotenweise', kv.fy, D.konvektion.fy, 1e-11);
  // Gegenprobe von Hand: ein KONSTANTES Feld hat keinen Gradienten, also keine Konvektion
  const c1 = new Float64Array(R2.n2).fill(1.7), c2 = new Float64Array(R2.n2).fill(-0.4);
  const kv0 = K.konvektion(S, c1, c2);
  let m0 = 0;
  for (let i = 0; i < kv0.fx.length; i++) m0 = Math.max(m0, Math.abs(kv0.fx[i]), Math.abs(kv0.fy[i]));
  U('konstantes Feld erzeugt keine Konvektion', m0 < 1e-12,
    `groesster Rest ${m0.toExponential(2)}`);
}

// ---------------------------------------------------------------- 8. Dirichlet
{
  const Sy = K.sattelsystem(S, D.dirichlet.dt);
  const N = Sy.N;
  U('Sattelsystem: Zahl der Unbekannten', N === D.dirichlet.N,
    `2*${S.n2} + ${S.n1} = ${N} (Soll ${D.dirichlet.N})`);
  const feste = new Uint8Array(N), werte = new Float64Array(N);
  for (let i = 0; i < S.n2; i++) {
    feste[i] = RB.fest[i]; werte[i] = RB.u[i];
    feste[S.n2 + i] = RB.fest[i]; werte[S.n2 + i] = RB.v[i];
  }
  let nf = 0; for (let i = 0; i < N; i++) nf += feste[i];
  U('Zahl der festgesetzten Unbekannten', nf === D.dirichlet.feste,
    `${nf} (Soll ${D.dirichlet.feste})`);
  const rhs0 = new Float64Array(N);
  for (let i = 0; i < 2*S.n2; i++) rhs0[i] = i * 1e-3;    // Druckteil bleibt null
  const r = K.rechteSeite(Sy, rhs0, feste, werte);
  /* Befund 01.11.2026: Dieser Vergleich fiel zuerst mit 9,1e-4 durch, und das war
   * richtig gemeldet -- nur lag der Unterschied nicht im Rechnen. fem_galerkin.py
   * schreibt die dritte Blockzeile als +B*u = 0, kern_fem.js als -B*u = 0. Dieselbe
   * Gleichung, mit -1 durchmultipliziert; dieselbe Loesung. Die rechten Seiten sind
   * deshalb im DRUCKTEIL Gegenzahlen. Wer das nicht weiss, sucht einen Fehler, den es
   * nicht gibt -- deshalb steht es hier und nicht nur im Kopf des Kerns. */
  const rv = new Float64Array(N);
  for (let i = 0; i < N; i++) rv[i] = (i < 2*S.n2) ? r[i] : -r[i];
  vergleicheProbe('rechte Seite nach dem Festsetzen, knotenweise ' +
                  '(Druckteil mit umgekehrtem Vorzeichen, s. Kopf)',
                  rv, D.dirichlet.rhs_probe, 1e-10);
  // Gegenprobe ohne Vergleichsdatei: an festen Stellen MUSS der Wert selbst stehen
  let m = 0;
  for (let i = 0; i < N; i++) if (feste[i]) m = Math.max(m, Math.abs(r[i] - werte[i]));
  U('an festgesetzten Stellen steht der vorgegebene Wert', m < 1e-15,
    `groesste Abweichung ${m.toExponential(2)}`);
}

// ---------------------------------------------------------------- 9. Zeitschritte
{
  const L = K.Loeser(P, T, D.schritte.dt, true);
  const Z = L.zustand();
  let schlimmD = 0, schlimmL = 0, schritte = 0;
  const t0 = Date.now();
  for (let k = 0; k < D.schritte.anzahl; k++) {
    L.schritt(Z);
    schritte += Z.schritte;
    schlimmD = Math.max(schlimmD, Math.abs(Z.cD - D.schritte.cD[k]) /
                                  Math.max(Math.abs(D.schritte.cD[k]), 1e-9));
    schlimmL = Math.max(schlimmL, Math.abs(Z.cL - D.schritte.cL[k]) /
                                  Math.max(Math.abs(D.schritte.cL[k]), 1e-9));
  }
  const ms = Date.now() - t0;
  U(`Widerstandsbeiwert ueber ${D.schritte.anzahl} Zeitschritte`, schlimmD < 1e-7,
    `groesste Abweichung ${schlimmD.toExponential(2)}`);
  U(`Auftriebsbeiwert ueber ${D.schritte.anzahl} Zeitschritte`, schlimmL < 1e-7,
    `groesste Abweichung ${schlimmL.toExponential(2)}`);
  vergleicheProbe('Geschwindigkeitsfeld ux am Ende, knotenweise',
                  Z.ux, D.schritte.ux_probe, 1e-7);
  vergleicheProbe('Geschwindigkeitsfeld uy am Ende, knotenweise',
                  Z.uy, D.schritte.uy_probe, 1e-7);
  vergleicheProbe('Druckfeld am Ende, knotenweise', Z.p, D.schritte.p_probe, 1e-6);
  U('MINRES bleibt unter der Schrittgrenze', schritte / D.schritte.anzahl < 300,
    `im Mittel ${(schritte/D.schritte.anzahl).toFixed(0)} Schritte, ` +
    `${(ms/D.schritte.anzahl).toFixed(1)} ms je Zeitschritt`);
}

// ---------------------------------------------------------------- 9b. Dauerlauf
{
  /* Befund 01.11.2026 (B18): Die Pruefung oben laeuft 40 Zeitschritte -- t = 0,08 s.
   * Im Dauerlauf platzte DIESELBE Rechnung bei t = 0,372 s, weil die explizite
   * Konvektion ihre Courant-Grenze hat. Eine Pruefung, die vor der kennzeichnenden
   * Zeit des Vorgangs aufhoert, sagt ueber Stabilitaet nichts. Hier wird bis t = 1,2 s
   * gerechnet -- lange genug, dass die Ablaugung hinter dem Zylinder anlaeuft -- und
   * zwar mit dem Zeitschritt, den der Kern SELBST aus dem Netz bestimmt. */
  const grenze = K.zeitschrittGrenze(P, T);
  const L = K.Loeser(P, T, null, true);          // null = Kern waehlt den Zeitschritt
  U('Kern bestimmt den Zeitschritt selbst aus der Courant-Bedingung',
    Math.abs(L.dt - grenze.dt) < 1e-15 && L.dt > 0,
    `kleinste Kante ${(grenze.hmin*1000).toFixed(2)} mm -> dt = ${L.dt.toExponential(3)} s`);
  const Z = L.zustand();
  const n = Math.round(1.2 / L.dt);
  let platzte = -1, grossD = 0;
  for (let k = 0; k < n; k++) {
    L.schritt(Z);
    if (!isFinite(Z.cD) || !isFinite(Z.cL) || Math.abs(Z.cD) > 1e3) { platzte = k + 1; break; }
    if (k > 20) grossD = Math.max(grossD, Math.abs(Z.cD));
  }
  U(`Dauerlauf ueber ${n} Zeitschritte bis t = 1,2 s bleibt endlich`, platzte < 0,
    platzte < 0 ? `cD zuletzt ${Z.cD.toFixed(4)}, cL ${Z.cL.toFixed(4)}, ` +
                  `groesster cD nach dem Anfahren ${grossD.toFixed(3)}`
                : `geplatzt bei Schritt ${platzte} (t = ${(platzte*L.dt).toFixed(3)} s)`);
  U('ein zu grosser Zeitschritt wird erkannt und gekappt',
    K.Loeser(P, T, 10 * grenze.dt, true).zuGross === true,
    `Grenze ${grenze.dt.toExponential(2)} s`);
}

// ---------------------------------------------------------------- 9c. Die Anzeige
{
  /* Die Wirbelstaerke auf das Anzeigegitter zu bringen ist eine Anzeige-, keine
   * Rechenoperation -- aber Befund B16 hat gezeigt, wie teuer ein Fehler genau dort
   * wird: Dreimal in Folge sah ein Vergleich nach einem grossen Unterschied der
   * VERFAHREN aus, und jedes Mal lag es an der Auswertung. Deshalb wird auch sie
   * gegen Python gehalten -- an beiden Auswerteorten. */
  const L = K.Loeser(P, T, D.schritte.dt, true);
  const Z = L.zustand();
  for (let k = 0; k < D.schritte.anzahl; k++) L.schritt(Z);
  const nx = D.anzeige.nx, ny = D.anzeige.ny;
  const mitte = K.wirbelAufGitter(L, Z.ux, Z.uy, nx, ny);
  vergleicheProbe('Wirbelstärke auf Zellmitten wie fem_wirbelstrasse.py',
                  mitte, D.anzeige.mitte, 1e-6);
  const knoten = K.wirbelAufGitter(L, Z.ux, Z.uy, nx, ny, false, nx + 1, ny + 1);
  vergleicheProbe('Wirbelstärke auf Gitterknoten (dort, wo der Netzwerkweg sie bildet)',
                  knoten, D.anzeige.knoten, 1e-6);
}

// ---------------------------------------------------------------- 10. Hagen-Poiseuille
{
  /* Die schaerfste Probe, die es fuer diesen Weg gibt: Im nackten Kanal ist die Loesung
   * geschlossen bekannt, und ein P2-Ansatz gibt eine Parabel EXAKT wieder -- nicht
   * genaehert. Weicht hier mehr als Rundungsrauschen ab, ist der Kern falsch. */
  const Pk = Float64Array.from(D.poiseuille.P), Tk = Int32Array.from(D.poiseuille.T);
  const R2k = K.p2Raum(Pk, Tk);
  const Sk = K.zusammenstellen(Pk, Tk, R2k);
  const RBk = K.raender(R2k.XY, R2k.n2, false);
  const n2 = Sk.n2, n1 = Sk.n1, N = 2*n2 + n1;
  const A = { zp: Sk.K.zp, zj: Sk.K.zj, zv: Float64Array.from(Sk.K.zv), nz: n2, ns: n2 };
  for (let k = 0; k < A.zv.length; k++) A.zv[k] *= K.F.NU;     // stationaer: A = nu K
  const Sy = { A: A, Bx: Sk.Bx, By: Sk.By, n2: n2, n1: n1, N: N, dt: 1 };
  const feste = new Uint8Array(N), werte = new Float64Array(N);
  for (let i = 0; i < n2; i++) {
    feste[i] = RBk.fest[i]; werte[i] = RBk.u[i];
    feste[n2+i] = RBk.fest[i]; werte[n2+i] = RBk.v[i];
  }
  const Sf = K.festsetzen(Sy, feste, werte);
  const r = K.rechteSeite(Sy, new Float64Array(N), feste, werte);
  const lgk = K.lambdaGradienten(Pk, Tk);
  const drk = K.druckraum(Pk, Tk, lgk);
  const V = { A: K.ic0(Sf.A), Kp: K.ic0(drk.Kp), mlump: drk.mlump };
  const l = K.minres(Sf, r, V, { schranke: 1e-13, hoechstens: 3000 });
  const ux = l.x.subarray(0, n2);
  let m = 0, umax = 0;
  for (let i = 0; i < n2; i++) {
    const soll = K.einstroemprofil(R2k.XY[2*i+1]);
    umax = Math.max(umax, soll);
    m = Math.max(m, Math.abs(ux[i] - soll));
  }
  U('Hagen-Poiseuille wird exakt wiedergegeben', m/umax < 1e-10,
    `${(m/umax).toExponential(3)} relativ, ${l.schritte} MINRES-Schritte ` +
    `(Python mit direkter Zerlegung: ${D.poiseuille.fehler.toExponential(2)})`);
}

// ---------------------------------------------------------------- Ausgabe
console.log('='.repeat(96));
console.log(` PRUEFUNG DES GALERKIN-KERNS gegen fem_galerkin.py -- Vergleichsdaten ${D.fassung}`);
console.log('='.repeat(96));
let gut = 0;
for (const [t, ok, b] of urteile) {
  if (ok) gut++;
  console.log(`  [${ok ? 'RICHTIG' : 'FALSCH '}] ${t}\n              ${b}`);
}
console.log('\n  == ' + gut + '/' + urteile.length);
process.exit(gut === urteile.length ? 0 : 1);
