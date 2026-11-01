/* netz_browser.js -- das koerperangepasste Dreiecksnetz im Browser, Fassung 1.0
 *
 * Dieselben drei Schritte wie netz_fem.py:
 *   1. Randpunkte setzen: Kanalrand und Zylinderkreis, mit vorgegebenem Abstand
 *   2. Innenpunkte streuen, fein am Zylinder und im Nachlauf, grob im Rest
 *   3. Delaunay-Triangulierung, dann die Dreiecke im Zylinder entfernen
 *
 * Die Triangulierung macht TRIANGLE von Jonathan Shewchuk, als WebAssembly -- derselbe
 * Vernetzer, den auch FEMM benutzt. Einen Vernetzer baut man nicht selbst; der Kern
 * sind exakte geometrische Praedikate mit adaptiver Arithmetik, und das ist eine eigene
 * Doktorarbeit. Shewchuk hat sie geschrieben.
 *
 * Das Netz ist NICHT dasselbe wie in Python: Die Innenpunkte werden mit einem eigenen
 * Zufallsgenerator gestreut, und der kann numpys PCG64 nicht nachbilden. Das ist kein
 * Mangel -- das Netz ist Eingabe, nicht Verfahren. Geprueft wird es deshalb an seinen
 * EIGENSCHAFTEN: positive Flaechen, kein Dreieck im Koerper, Gesamtflaeche gleich Kanal
 * minus Vieleck, Randknoten auf dem Kreis. Genau diese Proben fuehrt auch netz_fem.py.
 */
(function (global) {
'use strict';
const F = { LX: 2.2, LY: 0.41, XC: 0.2, YC: 0.2, DZ: 0.1 };

/* Ein kleiner, nachrechenbarer Zufallsgenerator (xorshift128+ mit fester Saat).
 * Fest heisst: Dasselbe Netz kommt bei jedem Aufruf wieder -- eine Seite, die bei
 * jedem Druck auf den Knopf andere Zahlen zeigt, ist nicht pruefbar. */
function wuerfel(saat) {
  let a = saat >>> 0 || 1, b = 362436069, c = 521288629, d = 88675123;
  return function () {
    const t = a ^ (a << 11);
    a = b; b = c; c = d;
    d = (d ^ (d >>> 19)) ^ (t ^ (t >>> 8));
    return (d >>> 0) / 4294967296;
  };
}

/* Punkte auf dem Kanalrand und auf dem Zylinderkreis. */
function randpunkte(hWand, hZyl) {
  const R = F.DZ / 2;
  const nx = Math.max(2, Math.round(F.LX / hWand));
  const ny = Math.max(2, Math.round(F.LY / hWand));
  const rand = [];
  for (let i = 0; i <= nx; i++) { rand.push(i * F.LX / nx, 0); }
  for (let i = 0; i <= nx; i++) { rand.push(i * F.LX / nx, F.LY); }
  for (let j = 1; j < ny; j++) { rand.push(0, j * F.LY / ny); }
  for (let j = 1; j < ny; j++) { rand.push(F.LX, j * F.LY / ny); }
  const nz = Math.max(12, Math.round(2 * Math.PI * R / hZyl));
  const zyl = [];
  for (let k = 0; k < nz; k++) {
    const w = 2 * Math.PI * k / nz;
    zyl.push(F.XC + R * Math.cos(w), F.YC + R * Math.sin(w));
  }
  return { rand: rand, zyl: zyl, nz: nz };
}

/* Gewuenschter Punktabstand am Ort (x,y) -- Zahl fuer Zahl wie _dichte() in netz_fem.py. */
function dichte(x, y, hFern, hZyl, o) {
  const R = F.DZ / 2;
  const d = Math.hypot(x - F.XC, y - F.YC) - R;
  let h = hZyl + (hFern - hZyl) * Math.min(Math.max(d / (o.weit * R), 0), 1);
  const nach = (x > F.XC) && (x < F.XC + o.nachlauf) && (Math.abs(y - F.YC) < o.band * R);
  if (nach) h = Math.min(h, hZyl * o.hNach);
  return h;
}

/* Randpunkte plus gestreute Innenpunkte. Ablehnungsverfahren auf einem feinen
 * Hilfsgitter: Ein Punkt wird genommen, wenn er weit genug von allen schon genommenen
 * entfernt ist (Poisson-Scheibe). Die Nachbarsuche laeuft ueber ein Kaestchenverzeichnis
 * -- sonst waere es quadratisch in der Punktzahl. */
function punkte(hFern, hZyl, mitZylinder, o) {
  o = Object.assign({ weit: 3.0, nachlauf: 1.0, hNach: 4.0, band: 2.0, saat: 7 }, o || {});
  const rp = randpunkte(hFern, hZyl);
  const P = mitZylinder ? rp.rand.concat(rp.zyl) : rp.rand.slice();
  const nRand = P.length / 2;
  const R = F.DZ / 2;
  const g = hZyl / 1.6;
  // Kaestchenverzeichnis mit Kantenlaenge hFern -- der groesste vorkommende Abstand
  const cw = hFern;
  const mx = Math.ceil(F.LX / cw), my = Math.ceil(F.LY / cw);
  const kaesten = new Array(mx * my);
  const eintragen = (x, y, i) => {
    const a = Math.min(mx-1, Math.floor(x / cw)), b = Math.min(my-1, Math.floor(y / cw));
    const k = a * my + b; (kaesten[k] || (kaesten[k] = [])).push(i);
  };
  for (let i = 0; i < nRand; i++) eintragen(P[2*i], P[2*i+1], i);
  const nah = (x, y, r) => {
    const a0 = Math.max(0, Math.floor((x - r) / cw)), a1 = Math.min(mx-1, Math.floor((x + r) / cw));
    const b0 = Math.max(0, Math.floor((y - r) / cw)), b1 = Math.min(my-1, Math.floor((y + r) / cw));
    const r2 = r * r;
    for (let a = a0; a <= a1; a++) for (let b = b0; b <= b1; b++) {
      const liste = kaesten[a * my + b];
      if (!liste) continue;
      for (const i of liste) {
        const dx = P[2*i] - x, dy = P[2*i+1] - y;
        if (dx*dx + dy*dy < r2) return true;
      }
    }
    return false;
  };
  // Hilfsgitter, gemischt
  const xs = [], ys = [];
  for (let x = g; x < F.LX; x += g) xs.push(x);
  for (let y = g; y < F.LY; y += g) ys.push(y);
  const kandidaten = new Float64Array(xs.length * ys.length * 2);
  let n = 0;
  for (const x of xs) for (const y of ys) { kandidaten[2*n] = x; kandidaten[2*n+1] = y; n++; }
  const w = wuerfel(o.saat);
  for (let i = n - 1; i > 0; i--) {                     // Mischen nach Fisher und Yates
    const j = Math.floor(w() * (i + 1));
    let t = kandidaten[2*i]; kandidaten[2*i] = kandidaten[2*j]; kandidaten[2*j] = t;
    t = kandidaten[2*i+1]; kandidaten[2*i+1] = kandidaten[2*j+1]; kandidaten[2*j+1] = t;
  }
  for (let i = 0; i < n; i++) {
    const x = kandidaten[2*i], y = kandidaten[2*i+1];
    if (mitZylinder && Math.hypot(x - F.XC, y - F.YC) <= R + 0.35 * hZyl) continue;
    const h = mitZylinder ? dichte(x, y, hFern, hZyl, o) : hFern;
    if (!nah(x, y, 0.78 * h)) { eintragen(x, y, P.length / 2); P.push(x, y); }
  }
  return { P: Float64Array.from(P), nRand: nRand, nzyl: mitZylinder ? rp.nz : 0 };
}

/* Delaunay durch Triangle, dann die Dreiecke im Zylinder entfernen und alle
 * gegen den Uhrzeigersinn drehen. */
function netz(Tri, hFern, hZyl, mitZylinder, o) {
  mitZylinder = (mitZylinder !== false);
  const pk = punkte(hFern, hZyl, mitZylinder, o);
  const ein = Tri.makeIO({ pointlist: Array.from(pk.P) });
  const aus = Tri.makeIO();
  Tri.triangulate({ pslg: false, quality: false, area: false, quiet: true }, ein, aus);
  const roh = Array.from(aus.trianglelist);
  Tri.freeIO(ein, true); Tri.freeIO(aus);
  const P = pk.P, R = F.DZ / 2;
  const T = [], A = [];
  for (let e = 0; e < roh.length / 3; e++) {
    let a = roh[3*e], b = roh[3*e+1], c = roh[3*e+2];
    const sx = (P[2*a] + P[2*b] + P[2*c]) / 3, sy = (P[2*a+1] + P[2*b+1] + P[2*c+1]) / 3;
    if (mitZylinder && Math.hypot(sx - F.XC, sy - F.YC) < R) continue;
    let A2 = (P[2*b]-P[2*a])*(P[2*c+1]-P[2*a+1]) - (P[2*c]-P[2*a])*(P[2*b+1]-P[2*a+1]);
    if (A2 < 0) { const t = b; b = c; c = t; A2 = -A2; }
    T.push(a, b, c); A.push(A2 / 2);
  }
  return { P: P, T: Int32Array.from(T), A: Float64Array.from(A),
           nRand: pk.nRand, nzyl: pk.nzyl };
}

global.NETZ_BROWSER = { F: F, wuerfel: wuerfel, randpunkte: randpunkte, dichte: dichte,
                        punkte: punkte, netz: netz };
})(typeof window !== 'undefined' ? window : globalThis);
if (typeof module !== 'undefined') module.exports = (typeof window !== 'undefined' ? window : globalThis).NETZ_BROWSER;
