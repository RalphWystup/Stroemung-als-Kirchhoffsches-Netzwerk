/* labor.js -- die Bedienung des CFD-Labors. Fassung 1.0
 * Rechnen tut kern_cfd.js (gegen Python geprueft); hier wird nur bedient und gezeichnet. */
'use strict';
const K = window.KernCFD;          // Netzwerkweg nach Kirchhoff
const KF = window.KERN_FEM;        // Galerkin-Weg, schwache Form
const NB = window.NETZ_BROWSER;    // Vernetzer (Triangle als WebAssembly)
const D = window.LABORDATEN || {};

// ---------------------------------------------------------------- Reiter
function zeige(w) {
  for (const a of window.PRUEFSTAND.reiter) {
    document.getElementById('ansicht-' + a).hidden = (w !== a);
    document.getElementById('r-' + a).classList.toggle('an', w === a);
  }
  if (w !== 'labor') window.scrollTo(0, 0);
  if (w === 'zeitverlauf') zeigeZeit(null);
  try { localStorage.setItem('cfd-reiter', w); } catch (e) {}
}

// ---------------------------------------------------------------- Zeichenhilfe
function schirm(c) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const b = c.clientWidth || c.width;
  if (c.width !== Math.round(b * dpr)) { c.width = Math.round(b * dpr); c.height = Math.round(c.height / (c.width / (b * dpr) || 1)); }
  return c.getContext('2d');
}
function farbe(n) { return getComputedStyle(document.documentElement).getPropertyValue(n).trim(); }

/* Ein kleiner Liniendiagrammzeichner -- eine Skala, die alle Marken traegt. */
function diagramm(c, reihen, opt) {
  opt = opt || {};
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const B = c.clientWidth || 900, H = opt.hoehe || 400;
  c.width = Math.round(B * dpr); c.height = Math.round(H * dpr);
  c.style.height = H + 'px';
  const x = c.getContext('2d'); x.setTransform(dpr, 0, 0, dpr, 0, 0);
  x.clearRect(0, 0, B, H);
  const L = 58, R = 14, O = 12, U = 34, bw = B - L - R, bh = H - O - U;
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
  for (const r of reihen) for (let k = 0; k < r.x.length; k++) {
    x0 = Math.min(x0, r.x[k]); x1 = Math.max(x1, r.x[k]);
    if (isFinite(r.y[k])) { y0 = Math.min(y0, r.y[k]); y1 = Math.max(y1, r.y[k]); }
  }
  if (opt.y0 !== undefined) y0 = opt.y0;
  if (opt.y1 !== undefined) y1 = opt.y1;
  if (!(y1 > y0)) { y1 = y0 + 1; }
  const sp = (y1 - y0) * 0.06; y0 -= sp; y1 += sp;
  const px = v => L + (v - x0) / (x1 - x0) * bw, py = v => O + (1 - (v - y0) / (y1 - y0)) * bh;
  const tinte = farbe('--tinte'), ged = farbe('--gedeckt'), linie = farbe('--linie');
  x.strokeStyle = linie; x.lineWidth = 1; x.font = '11px ui-monospace,Consolas,monospace';
  x.fillStyle = ged;
  for (let t = 0; t <= 5; t++) {
    const v = y0 + (y1 - y0) * t / 5, yy = Math.round(py(v)) + .5;
    x.beginPath(); x.moveTo(L, yy); x.lineTo(L + bw, yy); x.stroke();
    x.textAlign = 'right'; x.textBaseline = 'middle';
    x.fillText(opt.yfmt ? opt.yfmt(v) : v.toFixed(2), L - 7, yy);
  }
  for (let t = 0; t <= 5; t++) {
    const v = x0 + (x1 - x0) * t / 5, xx = Math.round(px(v)) + .5;
    x.textAlign = 'center'; x.textBaseline = 'top';
    x.fillText(opt.xfmt ? opt.xfmt(v) : v.toFixed(2), xx, O + bh + 7);
  }
  x.strokeStyle = ged; x.beginPath();
  x.moveTo(L, O); x.lineTo(L, O + bh); x.lineTo(L + bw, O + bh); x.stroke();
  for (const r of reihen) {
    x.strokeStyle = r.farbe; x.lineWidth = r.dick || 1.8;
    if (r.strich) x.setLineDash(r.strich); else x.setLineDash([]);
    x.beginPath();
    let an = false;
    for (let k = 0; k < r.x.length; k++) {
      const yy = r.y[k];
      if (!isFinite(yy)) { an = false; continue; }
      const X = px(r.x[k]), Y = py(Math.max(y0, Math.min(y1, yy)));
      if (!an) { x.moveTo(X, Y); an = true; } else x.lineTo(X, Y);
    }
    x.stroke();
    if (r.punkte) {
      x.fillStyle = r.farbe;
      for (let k = 0; k < r.x.length; k++)
        if (isFinite(r.y[k]) && r.y[k] >= y0 && r.y[k] <= y1) {
          x.beginPath(); x.arc(px(r.x[k]), py(r.y[k]), 2.6, 0, 7); x.fill();
        }
    }
  }
  x.setLineDash([]);
  if (opt.titel) { x.fillStyle = tinte; x.textAlign = 'left'; x.textBaseline = 'top';
    x.font = '600 12px system-ui'; x.fillText(opt.titel, L, 0); }
}

/* Feldbild: Wirbelstaerke als Farbfeld, blau-weiss-rot um null. */
function feldbild(c, w, nx, ny, skala, kreis) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const B = c.clientWidth || 900, H = Math.round(B * K.FALL.LY / K.FALL.LX);
  c.width = Math.round(B * dpr); c.height = Math.round(H * dpr);
  c.style.height = H + 'px';
  const x = c.getContext('2d'); x.setTransform(1, 0, 0, 1, 0, 0);
  const bd = x.createImageData(nx, ny);
  for (let j = 0; j < ny; j++)
    for (let i = 0; i < nx; i++) {
      const t = Math.max(-1, Math.min(1, w[i * ny + j] / (skala || 1)));
      let r, g, b;
      if (t >= 0) { r = 255; g = Math.round(255 - 200 * t); b = Math.round(255 - 215 * t); }
      else { const a = -t; r = Math.round(255 - 215 * a); g = Math.round(255 - 150 * a); b = 255; }
      const o = ((ny - 1 - j) * nx + i) * 4;
      bd.data[o] = r; bd.data[o + 1] = g; bd.data[o + 2] = b; bd.data[o + 3] = 255;
    }
  const h = document.createElement('canvas'); h.width = nx; h.height = ny;
  h.getContext('2d').putImageData(bd, 0, 0);
  x.imageSmoothingEnabled = true; x.imageSmoothingQuality = 'high';
  x.setTransform(dpr, 0, 0, dpr, 0, 0);
  x.clearRect(0, 0, B, H);
  x.drawImage(h, 0, 0, B, H);
  if (kreis !== false) {
    x.fillStyle = '#333a45'; x.strokeStyle = '#11151b'; x.lineWidth = 1;
    x.beginPath();
    x.arc(K.FALL.XC / K.FALL.LX * B, H - K.FALL.YC / K.FALL.LY * H,
          K.FALL.DZ / 2 / K.FALL.LX * B, 0, 7);
    x.fill(); x.stroke();
  }
}

// ===================================================================================
//  Vergleich 1 -- der Zweigstempel gegen die geschlossene Loesung
// ===================================================================================
const ERKLAERUNG = {
  stempel:
`Stationäre Konvektion-Diffusion in einer Dimension. Sie hat eine geschlossene
Lösung, also lässt sich jeder Stempel gegen die WAHRHEIT halten, nicht nur
gegen einen anderen Stempel.

Drei Stempel bilden denselben Zweigfluss verschieden:
  Scharfetter-Gummel   g = (D/h)·B(∓Pe_h), B(x)=x/(e^x−1)
  zentral              g = D/h ± u/2
  Aufwind              g = D/h + max(±u,0)

Scharfetter-Gummel ist knotenweise EXAKT, bei jeder Péclet-Zahl.`,
  druck:
`Sprachprobe: Derselbe Operator wird zweimal aufgebaut und die Wirkung verglichen.

  Weg A   Zweig für Zweig in eine dünnbesetzte Knotenadmittanzmatrix gestempelt
  Weg B   Flächenflüsse als ganze Felder, ohne jede Matrix

Mit Störkörper werden Zweige, deren Fläche im Körper liegt, AUFGETRENNT; die
letzte Zellspalte hängt über dy/(dx/2) am Bezugspotential.`,
  impuls:
`Dieselbe Sprachprobe für den Impulsoperator — den Teil, den ein passives Netz
nicht kann.

Die Unsymmetrie der Matrix IST die gesteuerte Quelle: g_hin − g_rück = u_f·A.
Geht die Geschwindigkeit gegen null, muss die Matrix symmetrisch werden. Dann
bleibt das reine Widerstandsnetz übrig — das Feld-Labor.`,
  loeser:
`Der vollständige Löser, beide Wege, im Browser gerechnet.

  Weg A   Impulsmatrix + Drucknetzmatrix, beide Zweig für Zweig gestempelt
  Weg B   alles als Feldoperation, Drucknetz matrixfrei

Beide lösen mit konjugierten Gradienten (einen direkten Löser gibt es im
Browser nicht) — verschieden ist der AUFBAU des Operators, nicht der Löser.
In der Python-Fassung unterscheidet sich zusätzlich der Löser.`
};
ERKLAERUNG.galerkin =
`Der ZWEITE Weg, hier im Browser gerechnet — nicht abgespielt.

  Netzwerk    Flussbilanz am Kontrollvolumen. Zelle = Knoten, Fläche = Zweig,
              Kontinuität = Knotenregel. Unbekannte sind Zellmittel, der
              Zylinder wird zur Treppe.
  Galerkin    schwache Form. Die Gleichungen werden mit Testfunktionen
              multipliziert und über das Gebiet integriert. Unbekannte sind
              Knotenwerte eines Polynomansatzes; der Zylinderrand LIEGT auf
              Knoten, es gibt keine Treppe.

      Geschwindigkeit  P2 (Ecken UND Kantenmitten)
      Druck            P1 (nur Ecken)

Diese Paarung erfüllt die inf-sup-Bedingung und braucht keinen abgestimmten
Stabilisierungsparameter. Vernetzt wird mit Triangle (Shewchuk), demselben
Vernetzer, den auch FEMM benutzt.

Dass beide am Ende ein dünnbesetztes System lösen, schadet nichts — das
Werkzeug ist nicht der Ansatz.`;

const QUELLE = {
  stempel: 'stroemung/sg_stempel.py — Stufe S01 des Prüfprotokolls\nProbe: python3 sg_stempel.py',
  druck: 'stroemung/ns_loeser.py, drucknetz_koerper_probe.py — S03/S04\nProbe: python3 drucknetz_koerper_probe.py',
  impuls: 'stroemung/impuls_matrix.py, konvektion.py — S07\nProbe: python3 impuls_matrix.py',
  loeser: 'stroemung/wirbelstrasse_zweiwege.py — S10/S11\nProbe: python3 wirbelstrasse_zweiwege.py 220 41 8.0',
  galerkin: 'stroemung/netz_fem.py, fem_galerkin.py, fem_wirbelstrasse.py — S13/S14\n' +
            'Kern hier: kern_fem.js, netz_browser.js\n' +
            'Probe: node pruefe_kern_fem.js (45 Posten gegen die Python-Module)'
};
let vergleich = 'stempel', mitKoerper = true;

function waehleVergleich(a) {
  vergleich = a;
  const kn = document.querySelectorAll('#wahl-vergleich button');
  ['stempel','druck','impuls','loeser','galerkin'].forEach((n,i)=>kn[i] && kn[i].classList.toggle('an',n===a));
  document.getElementById('k-stempel').hidden = (a !== 'stempel');
  document.getElementById('k-netz').hidden = !(a === 'druck' || a === 'impuls');
  document.getElementById('k-loeser').hidden = (a !== 'loeser');
  document.getElementById('k-galerkin').hidden = (a !== 'galerkin');
  document.getElementById('f-stempel').hidden = (a !== 'stempel');
  document.getElementById('f-netz').hidden = !(a === 'druck' || a === 'impuls');
  document.getElementById('f-loeser').hidden = (a !== 'loeser');
  document.getElementById('f-galerkin').hidden = (a !== 'galerkin');
  // Der Impulsoperator kennt KEINE Koerpermasken -- in Python wie hier. Der Koerper
  // wirkt erst im Zeitschritt, wo die abgedeckten Zweige auf null gesetzt werden.
  // Ein Knopf, der dort nichts tut, darf nicht dastehen, als taete er etwas
  // (Hausregel „gesperrte Knoepfe erklaeren"; gefunden von der Knopfprobe).
  const ohneWahl = (a === 'impuls' || a === 'galerkin');
  document.getElementById('wahl-koerper').hidden = ohneWahl;
  const hk = document.getElementById('h-kein-koerper');
  hk.hidden = !ohneWahl;
  hk.textContent = !ohneWahl ? ''
    : (a === 'impuls'
      ? 'Für den Impulsoperator gibt es hier nichts zu wählen: Er kennt keine\n' +
        'Körpermasken — weder hier noch in impuls_matrix.py. Der Störkörper wirkt\n' +
        'erst im Zeitschritt, wo die abgedeckten Zweige auf null gesetzt werden.'
      : 'Hier gibt es nichts zu wählen: Der Galerkin-Weg bekommt ein\n' +
        'körperangepasstes Netz, der Zylinder ist immer da — als Rand, nicht als\n' +
        'Maske. Ein nacktes Netz wäre ein anderer Fall (Hagen-Poiseuille), und der\n' +
        'steht in der Prüfung, nicht hier.');
  document.getElementById('erklaerung').textContent = ERKLAERUNG[a];
  document.getElementById('h-quelle').textContent = QUELLE[a];
  if (a === 'stempel') rechneStempel();
  if (a === 'druck' || a === 'impuls') rechneNetz();
}

function rechneStempel() {
  const sp = +document.getElementById('s-pe').value;
  const Pe = Math.pow(10, sp / 10), n = +document.getElementById('s-n').value;
  document.getElementById('w-pe').textContent = Pe >= 1000 ? Pe.toExponential(1) : Pe.toFixed(Pe < 10 ? 2 : 0);
  document.getElementById('w-n').textContent = n;
  const L = 1, Dif = 1, u = Pe * Dif / L;
  const xs = Array.from({length: n + 1}, (_, i) => i * L / n);
  const ex = xs.map(x => K.genau1d(x, L, Pe));
  const fein = Array.from({length: 401}, (_, i) => i * L / 400);
  const reihen = [{x: fein, y: fein.map(x => K.genau1d(x, L, Pe)),
                   farbe: farbe('--gedeckt'), dick: 3, name: 'geschlossene Lösung'}];
  const art = [['sg','Scharfetter-Gummel','--wegA'],['zentral','zentral','--wegB'],
               ['aufwind','Aufwind','--diff']];
  const zeilen = [];
  for (const [st, name, f] of art) {
    const c = Array.from(K.loese1d(n, L, u, Dif, st));
    let d = 0; for (let k = 0; k <= n; k++) d = Math.max(d, Math.abs(c[k] - ex[k]));
    reihen.push({x: xs, y: c, farbe: farbe(f), dick: 1.7, punkte: n <= 40, name: name});
    zeilen.push([name, d]);
  }
  diagramm(document.getElementById('c-stempel'), reihen,
    {hoehe: 400, y0: -0.25, y1: 1.25, titel: `Pe = ${Pe.toExponential(2)}   Pe_h = ${(Pe/n).toFixed(3)}   ${n} Zellen`,
     xfmt: v => v.toFixed(2), yfmt: v => v.toFixed(2)});
  document.getElementById('l-stempel').innerHTML =
    reihen.map(r => `<span><i style="background:${r.farbe}"></i>${r.name}</span>`).join('');
  document.getElementById('tafel-stempel').innerHTML =
    '<table class="tafel"><tr><th>Stempel</th><th>größter Knotenfehler</th><th>Urteil</th></tr>' +
    zeilen.map(([nm, d]) => `<tr><td>${nm}</td><td class="zahl">${d.toExponential(3)}</td>` +
      `<td class="${d < 1e-12 ? 'gut' : 'schlecht'}">${d < 1e-12 ? 'knotenweise exakt' : (100*d).toFixed(1) + ' % daneben'}</td></tr>`).join('') +
    '</table>';
  document.getElementById('h-stempel').textContent =
    `Pe_h = ${(Pe/n).toFixed(4)} ist die Péclet-Zahl JE ZELLE. Oberhalb von etwa 2\n` +
    `schwingt die zentrale Differenz; der Aufwind verschmiert stattdessen.\n` +
    `Scharfetter-Gummel tut weder das eine noch das andere.`;
}

// ===================================================================================
//  Vergleich 2 und 3 -- Sprachprobe fuer Drucknetz und Impulsoperator
// ===================================================================================
function zeigeNetzmass() {
  document.getElementById('w-nx').textContent = document.getElementById('s-nx').value;
  document.getElementById('w-ny').textContent = document.getElementById('s-ny').value;
}
function setzeKoerper(an) {
  mitKoerper = an;
  const kn = document.querySelectorAll('#wahl-koerper button');
  kn[0].classList.toggle('an', an); kn[1].classList.toggle('an', !an);
  rechneNetz();
}
function zufall(n, saat) {
  let z = saat || 12345;
  const a = new Float64Array(n);
  for (let k = 0; k < n; k++) { z = (z * 1103515245 + 12345) & 0x7fffffff; a[k] = z / 0x7fffffff - 0.5; }
  return a;
}
function rechneNetz() {
  const nx = +document.getElementById('s-nx').value, ny = +document.getElementById('s-ny').value;
  const g = new K.Gitter(K.FALL.LX, K.FALL.LY, nx, ny);
  const mk = K.koerperMasken(g);
  const mu = mitKoerper ? mk.mu : null, mv = mitKoerper ? mk.mv : null;
  const zeilen = []; let M = null;
  const t0 = performance.now();
  if (vergleich === 'druck') {
    const dn = K.drucknetzMatrix(g, mitKoerper, mu, mv); M = dn.A;
    const p = zufall(nx * ny, 7);
    const a = dn.A.mal(p), b = K.drucknetzWirkungVoll(g, p, mitKoerper ? dn.leer : null,
                                                      mitKoerper, mu, mv);
    // Im reinen Neumann-Fall haelt Weg A EINEN Knoten fest (sonst ist das Netz
    // stromlos aufgehaengt). Weg B kennt keine Einheitszeile. Beide unterscheiden
    // sich daher genau an diesem Knoten -- dort wird nicht verglichen, wie in der
    // Sprachprobe von ns_loeser.py.
    zeilen.push([mitKoerper ? 'Drucknetz: Weg A gegen Weg B'
                            : 'Drucknetz: Weg A gegen Weg B (ohne den Bezugsknoten)',
                 abw(a, b, mitKoerper ? -1 : 0), 1e-12]);
    const eins = new Float64Array(nx * ny).fill(1);
    if (mitKoerper) for (let k = 0; k < eins.length; k++) if (dn.leer[k]) eins[k] = 0;
    const r1 = dn.A.mal(eins);
    let innen = 0;
    for (let i = 0; i < nx - 1; i++) for (let j = 0; j < ny; j++) {
      const k = i*ny + j;
      if (k === 0 && !mitKoerper) continue;          // festgehaltener Bezugsknoten
      if (mitKoerper && dn.leer[k]) continue;        // Zellen ganz im Koerper
      innen = Math.max(innen, Math.abs(r1[k]));
    }
    zeilen.push(['konstanter Druck treibt im Innern keinen Strom', innen, 1e-9]);
    let uns = 0, mx = 0;
    for (let r = (mitKoerper ? 0 : 1); r < M.n; r++) for (let k = M.indptr[r]; k < M.indptr[r+1]; k++) {
      if (!mitKoerper && M.indices[k] === 0) continue;   // Spalte des Bezugsknotens
      mx = Math.max(mx, Math.abs(M.data[k]));
      let g2 = 0;
      for (let q = M.indptr[M.indices[k]]; q < M.indptr[M.indices[k]+1]; q++)
        if (M.indices[q] === r) g2 = M.data[q];
      uns = Math.max(uns, Math.abs(M.data[k] - g2));
    }
    zeilen.push([mitKoerper ? 'die Drucknetzmatrix ist symmetrisch'
                            : 'die Drucknetzmatrix ist symmetrisch (ohne Bezugszeile)',
                 uns / mx, 1e-14]);
  } else {
    const u = zufall((nx+1)*ny, 3), v = zufall(nx*(ny+1), 9);
    const ein = new Float64Array(ny);
    for (let j = 0; j < ny; j++) ein[j] = K.einstroemprofil(g.yu[j]);
    for (const [welch, einst, nm] of [['u', ein, 'u-Zweige, Einströmprofil vorgegeben'],
                                       ['u', null, 'u-Zweige, freier Einstrom'],
                                       ['v', null, 'v-Zweige']]) {
      const A = new K.ImpulsA(g, welch, einst);
      const x = welch === 'u' ? u : v;
      const rA = A.fuellen(u, v, K.FALL.NU).wirkung(x);
      const rB = welch === 'u' ? K.konvZaehU(g, u, v, K.FALL.NU, einst)
                               : K.konvZaehV(g, u, v, K.FALL.NU);
      zeilen.push([nm + ': Weg A gegen Weg B', abw(rA, rB), 1e-12]);
      if (!M) M = A.M;
    }
    const A0 = new K.ImpulsA(g, 'u', null);
    const u0 = Float64Array.from(u, z => z * 1e-9), v0 = Float64Array.from(v, z => z * 1e-9);
    const M0 = A0.fuellen(u0, v0, K.FALL.NU).M;
    let uns = 0, mx = 0;
    for (let r = 0; r < M0.n; r++) for (let k = M0.indptr[r]; k < M0.indptr[r+1]; k++) {
      mx = Math.max(mx, Math.abs(M0.data[k]));
      let g2 = 0;
      for (let q = M0.indptr[M0.indices[k]]; q < M0.indptr[M0.indices[k]+1]; q++)
        if (M0.indices[q] === r) g2 = M0.data[q];
      uns = Math.max(uns, Math.abs(M0.data[k] - g2));
    }
    zeilen.push(['Geschwindigkeit → 0 macht die Matrix symmetrisch', uns / mx, 1e-7]);
  }
  const ms = performance.now() - t0;
  const n_ok = zeilen.filter(([, d, s]) => d < s).length;
  document.getElementById('tafel-netz').innerHTML =
    '<table class="tafel"><tr><th>Probe</th><th>gemessen</th><th>Schranke</th><th>Urteil</th></tr>' +
    zeilen.map(([nm, d, s]) => `<tr><td>${nm}</td><td class="zahl">${d.toExponential(3)}</td>` +
      `<td class="zahl">${s.toExponential(0)}</td><td class="${d<s?'gut':'schlecht'}">` +
      `${d<s?'RICHTIG':'FALSCH'}</td></tr>`).join('') +
    `</table><p style="margin:10px 0 0;font:13px system-ui">` +
    `<b>${n_ok} von ${zeilen.length}</b> · ${nx}×${ny} Zellen · ${ms.toFixed(0)} ms</p>`;
  zeichneMuster(M);
  document.getElementById('h-netz').textContent =
    `Die Zufallsfelder werden aus einem festen Startwert erzeugt — derselbe\n` +
    `Druck, dieselbe Geschwindigkeit für beide Wege. Eigener Zufall je Weg\n` +
    `wäre kein Vergleich.`;
}
function abw(a, b, aus) {
  let d = 0, m = 0;
  for (let k = 0; k < a.length; k++) {
    if (k === aus) continue;
    d = Math.max(d, Math.abs(a[k]-b[k])); m = Math.max(m, Math.abs(b[k]));
  }
  return d / Math.max(m, 1e-300);
}
function zeichneMuster(M) {
  const c = document.getElementById('c-muster');
  if (!M) return;
  const n = M.n, S = Math.min(520, c.clientWidth || 520);
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  c.width = Math.round(S*dpr); c.height = Math.round(S*dpr); c.style.height = S + 'px';
  const x = c.getContext('2d'); x.setTransform(dpr,0,0,dpr,0,0);
  x.fillStyle = farbe('--flaeche'); x.fillRect(0,0,S,S);
  const s = S / n;
  let unsym = 0;
  for (let r = 0; r < n; r++)
    for (let k = M.indptr[r]; k < M.indptr[r+1]; k++) {
      const col = M.indices[k];
      let g2 = 0;
      for (let q = M.indptr[col]; q < M.indptr[col+1]; q++) if (M.indices[q] === r) g2 = M.data[q];
      const u = Math.abs(M.data[k] - g2) > 1e-12 * Math.abs(M.data[k]);
      if (u) unsym++;
      x.fillStyle = u ? farbe('--wegA') : farbe('--akzent');
      x.fillRect(col*s, r*s, Math.max(s,1), Math.max(s,1));
    }
  document.getElementById('b-nnz').textContent = `${n} Knoten · ${M.data.length} Einträge`;
  document.getElementById('h-muster').textContent = unsym
    ? `Rot: ${unsym} Einträge, deren Gegenstück sich unterscheidet. Diese Unsymmetrie\n` +
      `IST die gesteuerte Quelle — sie ist der Volumenstrom durch den Zweig.`
    : `Alle Einträge haben ihr symmetrisches Gegenstück: ein passives Widerstandsnetz\n` +
      `ohne gesteuerte Quellen.`;
}

// ===================================================================================
//  Vergleich 4 -- der vollstaendige Loeser, beide Wege, im Browser
// ===================================================================================
let lauf = null;
function zeigeTende() {
  document.getElementById('w-tende').textContent =
    (+document.getElementById('s-tende').value / 10).toFixed(1).replace('.', ',') + ' s';
}
function haltLoeser() { if (lauf) lauf.halt = true; }
function starteLoeser() {
  const [nx, ny] = document.getElementById('s-netz').value.split(',').map(Number);
  const t_ende = +document.getElementById('s-tende').value / 10;
  const L = new K.Loeser(nx, ny, 0.4);
  const ZA = L.zustand('A'), ZB = L.zustand('B');
  lauf = { L: L, ZA: ZA, ZB: ZB, t_ende: t_ende, halt: false, du: 0, dcD: 0,
           t0: performance.now(), bahn: [] };
  document.getElementById('b-start').disabled = true;
  document.getElementById('b-stop').disabled = false;
  haeppchen();
}
function haeppchen() {
  if (!lauf) return;
  const L = lauf.L, bis = performance.now() + 120;
  while (performance.now() < bis && lauf.ZA.t < lauf.t_ende && !lauf.halt) {
    L.schritt(lauf.ZA); L.schritt(lauf.ZB);
    let d = 0, m = 0;
    for (let k = 0; k < lauf.ZA.u.length; k++) {
      d = Math.max(d, Math.abs(lauf.ZA.u[k] - lauf.ZB.u[k]));
      m = Math.max(m, Math.abs(lauf.ZA.u[k]));
    }
    lauf.du = Math.max(lauf.du, d / Math.max(m, 1e-30));
    lauf.dcD = Math.max(lauf.dcD, Math.abs(lauf.ZA.cD - lauf.ZB.cD));
    lauf.bahn.push([lauf.ZA.t, lauf.ZA.cD, lauf.ZA.cL]);
  }
  zeichneLauf();
  if (lauf.ZA.t < lauf.t_ende && !lauf.halt) { setTimeout(haeppchen, 0); return; }
  document.getElementById('b-start').disabled = false;
  document.getElementById('b-stop').disabled = true;
  lauf.fertig = true; zeichneLauf();
}
function zeichneLauf() {
  const L = lauf.L, g = L.g, Z = lauf.ZA;
  const om = K.wirbelstaerke(g, Z.u, Z.v);
  let s = 0; for (let k = 0; k < om.length; k++) s = Math.max(s, Math.abs(om[k]));
  feldbild(document.getElementById('c-loeser'), om, g.nx - 1, g.ny - 1, s * 0.22);
  const sek = (performance.now() - lauf.t0) / 1000;
  const fertig = lauf.fertig || lauf.halt;
  document.getElementById('tafel-loeser').innerHTML =
    '<table class="tafel"><tr><th>Größe</th><th>Weg A</th><th>Weg B</th><th>|A − B|</th></tr>' +
    `<tr><td>Widerstandsbeiwert c_D</td><td class="zahl">${Z.cD.toFixed(6)}</td>` +
    `<td class="zahl">${lauf.ZB.cD.toFixed(6)}</td><td class="zahl">${Math.abs(Z.cD-lauf.ZB.cD).toExponential(2)}</td></tr>` +
    `<tr><td>Auftriebsbeiwert c_L</td><td class="zahl">${Z.cL.toFixed(6)}</td>` +
    `<td class="zahl">${lauf.ZB.cL.toFixed(6)}</td><td class="zahl">${Math.abs(Z.cL-lauf.ZB.cL).toExponential(2)}</td></tr>` +
    `<tr><td>CG-Schritte im Drucknetz</td><td class="zahl">${Z.cg}</td><td class="zahl">${lauf.ZB.cg}</td><td>—</td></tr>` +
    '</table>' +
    `<p style="margin:10px 0 0;font:13px system-ui">t = <b>${Z.t.toFixed(3)} s</b> von ${lauf.t_ende.toFixed(1)} s · ` +
    `${Z.k} Zeitschritte · ${sek.toFixed(1)} s Rechenzeit` + (fertig ? ' · <b>fertig</b>' : '') + '</p>' +
    `<p style="margin:4px 0 0;font:13px system-ui">Größte Abweichung über den ganzen Lauf: ` +
    `Geschwindigkeitsfeld <b class="${lauf.du<1e-9?'gut':'schlecht'}">${lauf.du.toExponential(3)}</b> relativ, ` +
    `c_D <b class="${lauf.dcD<1e-6?'gut':'schlecht'}">${lauf.dcD.toExponential(3)}</b> absolut.</p>`;
  document.getElementById('h-loeser').textContent =
    `Netz ${g.nx}×${g.ny}, dt = ${L.dt.toExponential(3)} s.\n` +
    `Beide Wege rechnen vom selben Anfangszustand. Was auseinanderläuft,\n` +
    `läuft an der Diskretisierung auseinander — nicht am Zufall.`;
}

// ===================================================================================
//  Zeitverlauf -- echte Rechnung und Surrogat
// ===================================================================================
let zeitArt = 'unabhaengig', spielt = 0, zIx = 0;
function rohBytes(b64) {
  const roh = atob(b64), a = new Uint8Array(roh.length);
  for (let i = 0; i < roh.length; i++) a[i] = roh.charCodeAt(i);
  return a;
}
function entpacke(b64) { return new Int16Array(rohBytes(b64).buffer); }
function entpackeF32(b64) { return new Float32Array(rohBytes(b64).buffer); }
let ZD = null;
function zeitDaten() {
  if (ZD || !D.zeit) return ZD;
  const z = D.zeit;
  ZD = { nx: z.nx, ny: z.ny, n: z.frames, t: z.t,
         A: entpacke(z.datA), Dd: entpacke(z.datD), sA: z.skalaA, sD: z.skalaD,
         N: z.nx * z.ny, sur: z.surrogat, bericht: z.bericht, gal: z.galerkin };
  if (ZD.gal) { ZD.gal.F = entpacke(ZD.gal.dat); ZD.gal.FN = entpacke(ZD.gal.datN); }
  if (ZD.sur) {
    ZD.sur.mittelF = entpackeF32(ZD.sur.mittel);
    ZD.sur.modenF = entpackeF32(ZD.sur.moden);
    ZD.sur.bahnA = []; ZD.sur.bahnB = [];
    for (let k = 0; k < ZD.n; k++) {
      const tt = ZD.t[k] - ZD.sur.t0;
      ZD.sur.bahnA.push(harmonisch(ZD.sur.BA, ZD.sur.f, ZD.sur.K, ZD.sur.r, tt));
      ZD.sur.bahnB.push(harmonisch(ZD.sur.BB, ZD.sur.f, ZD.sur.K, ZD.sur.r, tt));
    }
  }
  return ZD;
}
/* a(t) = a_0 + Summe_k [ c_k cos(2 pi k f t) + s_k sin(2 pi k f t) ].
   B liegt zeilenweise als (1+2K) x r. Exakt periodisch -- es kann sich keine Phase
   davonstehlen, wie es der angepasste lineare Zeitoperator tut. */
function harmonisch(B, f, K, r, t) {
  const a = new Float64Array(r);
  for (let m = 0; m < r; m++) {
    let v = B[m];
    for (let k = 1; k <= K; k++) {
      const w = 2 * Math.PI * k * f * t;
      v += B[(2*k-1)*r + m] * Math.cos(w) + B[(2*k)*r + m] * Math.sin(w);
    }
    a[m] = v;
  }
  return a;
}
function ausModen(sur, a) {
  const N = sur.mittelF.length, f = new Float64Array(N);
  f.set(sur.mittelF);
  for (let m = 0; m < sur.r; m++) {
    const c = a[m], o = m * N;
    if (c === 0) continue;
    for (let k = 0; k < N; k++) f[k] += c * sur.modenF[o + k];
  }
  return f;
}
function perzentil(f, q) {
  const a = Float64Array.from(f, Math.abs);
  a.sort();
  return a[Math.min(a.length - 1, Math.floor(q * a.length))] || 1;
}
const ZEITARTEN = ['unabhaengig', 'wege', 'surrogat'];
function zeitLaenge(z) {
  return (zeitArt === 'unabhaengig' && z.gal) ? z.gal.frames : z.n;
}
function zeitAchse(z) {
  return (zeitArt === 'unabhaengig' && z.gal) ? z.gal.t : z.t;
}
function waehleZeit(a) {
  zeitArt = a;
  const kn = document.querySelectorAll('#wahl-zeit button');
  ZEITARTEN.forEach((n, i) => kn[i] && kn[i].classList.toggle('an', n === a));
  const z = zeitDaten();
  if (z) {
    const L = zeitLaenge(z);
    document.getElementById('s-zeit').max = String(L - 1);
    if (zIx >= L) zIx = 0;
  }
  zeigeZeit(zIx);
}
function spielen() {
  if (spielt) { clearInterval(spielt); spielt = 0; document.getElementById('b-spiel').textContent = 'Abspielen'; return; }
  document.getElementById('b-spiel').textContent = 'Anhalten';
  spielt = setInterval(() => { const z = zeitDaten(); if (!z) return;
    zeigeZeit((zIx + 1) % zeitLaenge(z)); }, 70);
}
function zeigeZeit(ix) {
  const z = zeitDaten();
  if (!z) { document.getElementById('h-zeit').textContent = 'Keine Zeitdaten eingebaut.'; return; }
  const LAENGE = zeitLaenge(z);
  document.getElementById('s-zeit').max = String(LAENGE - 1);
  if (ix !== null) { zIx = Math.min(ix, LAENGE - 1);
    document.getElementById('s-zeit').value = String(zIx); }
  if (zIx >= LAENGE) zIx = LAENGE - 1;
  const lupe = Math.pow(10, +document.getElementById('s-lupe').value / 10);
  document.getElementById('w-lupe').textContent = '×' + (lupe >= 1000 ? lupe.toExponential(0) : lupe.toFixed(0));
  document.getElementById('w-zeit').textContent =
    (zeitAchse(z)[zIx]).toFixed(3).replace('.', ',') + ' s';
  const N = z.N, o = zIx * N;
  let f1, f2, f3, t1, t2, t3;
  if (zeitArt === 'unabhaengig') {
    const G = z.gal;
    if (!G) { document.getElementById('h-zeit').textContent =
      'Der unabhängige Weg ist in dieser Fassung nicht eingebaut.'; return; }
    f1 = new Float64Array(N); f2 = new Float64Array(N); f3 = new Float64Array(N);
    const sn = G.skalaN[zIx] / 32767, sg = G.skala[zIx] / 32767;
    for (let k = 0; k < N; k++) {
      f1[k] = G.FN[o + k] * sn;
      f2[k] = G.F[o + k] * sg;
      f3[k] = f1[k] - f2[k];
    }
    t1 = 'Netzwerk nach Kirchhoff — Flussbilanz, Zylinder als Treppe (' + G.netzquelle + ')';
    t2 = 'Galerkin — schwache Form, körperangepasstes Netz';
    t3 = 'Differenz der beiden Herleitungen, gleiche Skala';
  } else if (zeitArt === 'wege') {
    f1 = new Float64Array(N); f2 = new Float64Array(N); f3 = new Float64Array(N);
    const sa = z.sA[zIx] / 32767, sd = z.sD[zIx] / 32767;
    for (let k = 0; k < N; k++) {
      f1[k] = z.A[o + k] * sa;
      f3[k] = z.Dd[o + k] * sd;
      f2[k] = f1[k] - f3[k];
    }
    t1 = 'Weg A — Netzwerk nach Kirchhoff (gerechnet)';
    t2 = 'Weg B — klassische Finite Volumen (gerechnet)';
    t3 = 'Differenz A − B — bei gleicher Skala ist nichts zu sehen';
    t1 = t1; t2 = t2;
  } else {
    const sur = z.sur;
    if (!sur) { document.getElementById('h-zeit').textContent = 'Kein Surrogat eingebaut.'; return; }
    f1 = ausModen(sur, sur.bahnA[zIx]);
    f2 = new Float64Array(N);
    const sa = z.sA[zIx] / 32767;
    for (let k = 0; k < N; k++) f2[k] = z.A[o + k] * sa;
    f3 = new Float64Array(N);
    for (let k = 0; k < N; k++) f3[k] = f1[k] - f2[k];
    t1 = 'Surrogat — ' + sur.r + ' POD-Moden, Zeitgesetz harmonisch angepasst';
    t2 = 'echte Rechnung, Weg A';
    t3 = 'Differenz Surrogat − Rechnung';

  }
  let s1 = 0, s3 = 0;
  for (let k = 0; k < N; k++) { s1 = Math.max(s1, Math.abs(f1[k])); s3 = Math.max(s3, Math.abs(f3[k])); }
  // Die Farbskala nimmt das 99. Perzentil von |omega|, nicht das Maximum: Dicht an der
  // Wand ist die Wirbelstaerke um ein Vielfaches groesser als im Nachlauf, und eine
  // Skala nach dem Maximum liesse die Wirbel verblassen. Genau so haelt es auch das
  // matplotlib-Bild in bild.py.
  const sk1 = perzentil(f1, 0.99);   // dasselbe Perzentil wie bild.py
  // Der dritte Schirm liegt anfangs auf DERSELBEN Skala wie die beiden oberen. Dann ist
  // er leer -- und genau das ist die Aussage: Bei gleicher Skala gibt es keinen
  // Unterschied zu sehen. Erst die Verstaerkung macht sichtbar, was uebrig bleibt.
  // Eine von vornherein eigene Skala haette Rundungsrauschen wie ein Feld aussehen lassen.
  const sk3 = sk1 / Math.max(lupe, 1);
  feldbild(document.getElementById('c-z1'), f1, z.nx, z.ny, sk1);
  feldbild(document.getElementById('c-z2'), f2, z.nx, z.ny, sk1);
  feldbild(document.getElementById('c-z3'), f3, z.nx, z.ny, sk3);
  document.getElementById('t-z1').textContent = t1;
  document.getElementById('t-z2').textContent = t2;
  document.getElementById('t-z3').textContent = t3;
  const bez = 'Skala ±' + sk1.toFixed(1) + ' 1/s · größter Wert ±' + s1.toFixed(0);
  document.getElementById('b-z1').textContent = bez;
  document.getElementById('b-z2').textContent = bez;
  const noetig = sk1 / Math.max(s3, 1e-300);
  document.getElementById('b-z3').textContent =
    'Skala ±' + sk3.toExponential(1) + ' 1/s · größter Wert ±' + s3.toExponential(1) +
    (noetig >= 1.5 ? ' · nötige Verstärkung ×' + noetig.toExponential(0)
                   : ' · übersteigt die Skala stellenweise');
  const be = z.bericht || {};
  const B = z.bericht || {}, G2 = z.gal;
  document.getElementById('h-zeit').textContent = zeitArt === 'unabhaengig'
    ? `Zwei verschiedene HERLEITUNGEN desselben Falls: oben die Flussbilanz am\n` +
      `Kontrollvolumen (Zylinder als Treppe), in der Mitte die schwache Form mit\n` +
      `Ansatzfunktionen auf körperangepasstem Netz. Unten ihre Differenz, auf\n` +
      `derselben Skala — sie ist sichtbar, und das ist richtig so: Zwei Verfahren,\n` +
      `die verschieden scheitern können, stimmen nicht auf Rundung überein.\n` +
      `Galerkin trifft alle drei Benchmarkgrößen im Intervall: c_D ` +
      `${(G2 ? G2.cD : 0).toFixed(4)}, c_L ${(G2 ? G2.cL : 0).toFixed(4)}, ` +
      `St ${(G2 ? G2.St : 0).toFixed(4)}.\n` +
      `Zugeordnet wird über die PHASE des Grenzzyklus, nicht über die absolute Zeit —\n` +
      `die beiden Rechnungen laufen nicht im Gleichtakt (Restfehler ` +
      `${(G2 ? G2.passung_grad : 0).toFixed(2)}°).\n\n` +
      `PUNKTWEISE stimmen die Felder weit weniger überein als die Beiwerte: Effektivwert\n` +
      `der Differenz ${(G2 && G2.masse ? G2.masse.gesamt : 0).toFixed(0)} % des Feldes — im Nachlauf ` +
      `${(G2 && G2.masse ? G2.masse.nachlauf : 0).toFixed(0)} %, nahe am Zylinder ` +
      `${(G2 && G2.masse ? G2.masse.nahe : 0).toFixed(0)} %,\n` +
      `die größte Abweichung ${(G2 && G2.masse ? G2.masse.ort_mm : 0).toFixed(0)} mm vom Mittelpunkt, also dicht an der Wand (Radius 50 mm).\n` +
      `Das ist kein Widerspruch, sondern eine Lehre: Integralgrößen konvergieren viel\n` +
      `schneller als die Felder selbst. Ein Benchmark bescheinigt c_D, c_L und St —\n` +
      `nicht das Feld an jeder Stelle.`
    : zeitArt === 'wege'
    ? `Oben und in der Mitte dieselbe Strömung, auf zwei getrennt aufgebauten Wegen\n` +
      `gerechnet. Der dritte Schirm liegt auf DERSELBEN Skala wie die beiden oberen und ist\n` +
      `deshalb leer. Drehen Sie die Verstärkung auf, bis etwas erscheint: Es braucht\n` +
      `rund das ${(1/(be.unterschied_echt||1e-14)).toExponential(0)}-fache — dann sieht man Rundungsrauschen, kein Feld.`
    : `Oben das angelernte Surrogat (${(z.sur||{}).r} POD-Moden, Zeitgesetz harmonisch\n` +
      `angepasst), in der Mitte die\n` +
      `echte Rechnung, unten die Differenz. Das Surrogat ist aus diesem Fall gelernt\n` +
      `und gilt nur für ihn — nicht für andere Reynoldszahlen oder Geometrien.`;
  const b = z.bericht || {};
  document.getElementById('h-zeitfuss').textContent =
    `Surrogat: ${(z.sur||{}).r} POD-Moden, ${(100*(b.energie||0)).toFixed(4)} % der Energie erfasst;\n` +
    `Zeitgesetz als harmonische Anpassung mit ${(z.sur||{}).K} Oberwellen der Grundfrequenz ` +
    `${(b.f||0).toFixed(4)} Hz.\n` +
    `Gemessener Fehler gegen den echten Löser über den ganzen Verlauf: ` +
    `${(100*(b.fehler_harmonisch||0)).toFixed(3)} % der Schwankungsbreite,\n` +
    `davon ${(100*(b.fehler_trunkierung||0)).toFixed(3)} % allein durch das Abschneiden auf ${(z.sur||{}).r} Moden.\n` +
    `Zum Vergleich: ein angepasster linearer Zeitoperator (DMD) kommt auf demselben\n` +
    `Datensatz auf ${(100*(b.fehler_dmd||0)).toFixed(2)} % — er ist nicht exakt periodisch, und über 6 Perioden\n` +
    `läuft daraus ein Phasenfehler auf. Beide Zahlen sind gemessen, keine ist geschätzt.\n` +
    `Die beiden Surrogate (Weg A und Weg B) unterscheiden sich um ${(b.unterschied_surrogat||0).toExponential(2)} —\n` +
    `sie erben also die Übereinstimmung der echten Rechnung (${(b.unterschied_echt||0).toExponential(2)}).\n\n`  +
    `Gespeichert sind Weg A und die DIFFERENZ, jeweils mit eigener Skala; Weg B wird\n` +
    `daraus als A − Differenz gebildet. So bleibt der Unterschied in voller relativer\n` +
    `Genauigkeit erhalten — zwei getrennt quantisierte Felder hätten ihn zerstört.`;
  const unbenutzt =
    `Gespeichert sind Weg A und die DIFFERENZ, jeweils mit eigener Skala; Weg B wird\n` +
    `daraus als A − Differenz gebildet. So bleibt der Unterschied in voller relativer\n` +
    `Genauigkeit erhalten — zwei getrennt quantisierte Felder hätten ihn zerstört.`;
}

// ---------------------------------------------------------------- Start
// ================================================================ Galerkin: zweiter Weg
/* Beide Wege laufen hier NEBENEINANDER im Browser:
 *   Weg 1  kern_cfd.js   -- Netzwerk nach Kirchhoff auf kartesischen Zellen
 *   Weg 2  kern_fem.js   -- Galerkin, Taylor-Hood, koerperangepasstes Dreiecksnetz
 * Sie teilen nur Geometrie, Stoffwerte und die Zeit -- keine Zeile Rechenweg.
 */
let gLauf = null, gBereit = null;

function zeigeGTende() {
  document.getElementById('w-gtende').textContent =
    (+document.getElementById('s-gtende').value / 10).toFixed(1).replace('.', ',') + ' s';
}

function netzBild(c, M) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const B = c.clientWidth || 900, H = Math.round(B * K.FALL.LY / K.FALL.LX);
  c.width = Math.round(B * dpr); c.height = Math.round(H * dpr);
  c.style.height = H + 'px';
  const x = c.getContext('2d');
  x.setTransform(dpr, 0, 0, dpr, 0, 0);
  x.fillStyle = '#ffffff'; x.fillRect(0, 0, B, H);
  const sx = B / K.FALL.LX, sy = H / K.FALL.LY;
  x.strokeStyle = '#9aa3ad'; x.lineWidth = 0.4;
  x.beginPath();
  for (let e = 0; e < M.T.length / 3; e++) {
    const a = M.T[3*e], b = M.T[3*e+1], c3 = M.T[3*e+2];
    x.moveTo(M.P[2*a]*sx, H - M.P[2*a+1]*sy);
    x.lineTo(M.P[2*b]*sx, H - M.P[2*b+1]*sy);
    x.lineTo(M.P[2*c3]*sx, H - M.P[2*c3+1]*sy);
    x.closePath();
  }
  x.stroke();
  // der Zylinderrand, damit man sieht, dass er auf Knoten liegt
  x.strokeStyle = '#b23a3a'; x.lineWidth = 1.4;
  x.beginPath();
  x.arc(K.FALL.XC*sx, H - K.FALL.YC*sy, K.FALL.DZ/2*sx, 0, 2*Math.PI);
  x.stroke();
}

function gBahnBild(c, bahnG, bahnN, dtG, dtN) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const B = c.clientWidth || 900, H = 300;
  c.width = Math.round(B*dpr); c.height = Math.round(H*dpr);
  const x = c.getContext('2d'); x.setTransform(dpr, 0, 0, dpr, 0, 0);
  x.fillStyle = '#fff'; x.fillRect(0, 0, B, H);
  const L = 54, Rr = 12, O = 14, Uu = 26;
  let tmax = 0, ymin = 1e9, ymax = -1e9;
  const reihen = [[bahnG, '#2e6f9e'], [bahnN, '#b23a3a']];
  /* Der Anfahrstoss darf den Massstab nicht bestimmen: Bei t = 0 springt c_D auf ueber
   * 200, weil die Stroemung aus der Ruhe schlagartig anliegt. Das ist richtig gerechnet
   * und trotzdem nichts, was man sehen will -- danach waere alles andere ein Strich auf
   * der Null. Der Massstab kommt deshalb aus der Zeit NACH einem Zehntel des Laufs. */
  for (const [bahn] of reihen) for (const [tt] of (bahn || [])) tmax = Math.max(tmax, tt);
  if (!(tmax > 0)) { x.fillStyle = '#5b636e'; x.font = '12px system-ui';
    x.fillText('noch nichts gerechnet', L, H/2); return; }
  const ab = 0.1 * tmax;
  for (const [bahn] of reihen) for (const [tt, cd, cl] of (bahn || [])) {
    if (tt < ab) continue;
    if (isFinite(cd)) { ymin = Math.min(ymin, cd); ymax = Math.max(ymax, cd); }
    if (isFinite(cl)) { ymin = Math.min(ymin, cl); ymax = Math.max(ymax, cl); }
  }
  if (!(ymax > ymin)) { ymin = -1; ymax = 4; }
  const luft = 0.12 * (ymax - ymin);
  ymin -= luft; ymax += luft;
  const px = tt => L + (B - L - Rr) * tt / tmax;
  const py = v => O + (H - O - Uu) * (ymax - v) / (ymax - ymin);
  x.strokeStyle = '#d9dad4'; x.lineWidth = 1;
  x.beginPath();
  for (let k = 0; k <= 4; k++) { const v = ymin + k*(ymax-ymin)/4; x.moveTo(L, py(v)); x.lineTo(B-Rr, py(v)); }
  x.stroke();
  x.fillStyle = '#5b636e'; x.font = '11px ui-monospace,monospace'; x.textAlign = 'right';
  for (let k = 0; k <= 4; k++) { const v = ymin + k*(ymax-ymin)/4; x.fillText(v.toFixed(1), L-6, py(v)+3); }
  x.textAlign = 'center';
  for (let k = 0; k <= 4; k++) { const tt = tmax*k/4; x.fillText(tt.toFixed(1)+' s', px(tt), H-8); }
  for (const [bahn, farbe] of reihen) {
    if (!bahn || !bahn.length) continue;
    for (const [feld, strich] of [[1, []], [2, [4, 3]]]) {
      x.strokeStyle = farbe; x.lineWidth = 1.5; x.setLineDash(strich);
      x.beginPath();
      let angesetzt = false;
      bahn.forEach(b => {
        const v = b[feld];
        if (!isFinite(v) || v < ymin || v > ymax) { angesetzt = false; return; }
        const X = px(b[0]), Y = py(v);
        if (!angesetzt) { x.moveTo(X, Y); angesetzt = true; } else x.lineTo(X, Y);
      });
      x.stroke(); x.setLineDash([]);
    }
  }
  x.textAlign = 'left'; x.font = '11px system-ui';
  x.fillStyle = '#2e6f9e'; x.fillText('Galerkin  — c_D   -- c_L', L + 6, O + 12);
  x.fillStyle = '#b23a3a'; x.fillText('Netzwerk  — c_D   -- c_L', L + 6, O + 26);
}

/* Triangle einmal hochfahren. init() ist asynchron -- das WebAssembly muss erst
 * uebersetzt werden. Befund 01.11.2026: Ohne dieses Warten griff der erste Aufruf auf
 * Module._malloc, das es noch nicht gab; die Seite meldete „Netz wird erzeugt …" und
 * blieb stehen. Ein Knopf, der nichts tut und nicht sagt warum, ist ein Fehler
 * (Hausregel „gesperrte Knoepfe erklaeren"). */
function triBereit() {
  if (!gBereit) gBereit = window.TRIANGLE.init();
  return gBereit;
}

function starteGalerkin() {
  const hz = +document.getElementById('s-hzyl').value;
  const tEnde = +document.getElementById('s-gtende').value / 10;
  const h = document.getElementById('h-galerkin');
  h.textContent = 'Vernetzer wird hochgefahren …';
  document.getElementById('b-gstart').disabled = true;
  document.getElementById('b-gstop').disabled = false;
  triBereit().then(() => {
    h.textContent = 'Netz wird erzeugt …';
    const t0 = performance.now();
    const M = NB.netz(window.TRIANGLE, hz * 3.5, hz, true);
    const LG = KF.Loeser(M.P, M.T, null, true);      // Zeitschritt aus der Courant-Bedingung
    const tNetz = performance.now() - t0;
    netzBild(document.getElementById('c-gnetz'), M);
    document.getElementById('b-gnetz').textContent =
      `${M.P.length/2} Ecken · ${M.T.length/3} Dreiecke · ${LG.N} Unbekannte · ${tNetz.toFixed(0)} ms`;
    /* Der Netzwerkweg auf GLEICHEM KNOTENABSTAND -- nicht auf gleicher Kantenlänge.
     * Befund B16: Zwei Dinge nebeneinanderzulegen ist erst dann ein Vergleich, wenn
     * alles andere gleich ist. Taylor-Hood P2 trägt Knoten auf Ecken UND Kantenmitten,
     * der Knotenabstand ist also die halbe Kantenlänge. Ein Netzwerkgitter mit der
     * vollen Kantenlänge hätte die halbe Auflösung, und der Unterschied im Bild wäre
     * die Auflösung, nicht das Verfahren. */
    const dxSoll = hz / 2;
    const nx = Math.max(60, Math.round(K.FALL.LX / dxSoll));
    const ny = Math.max(12, Math.round(K.FALL.LY / dxSoll));
    const LN = new K.Loeser(nx, ny, 0.4);
    const ZN = LN.zustand('A');
    /* Der Netzwerkweg fuehrt seine Bahn am LAUF, nicht am Zustand (so ist es im
     * Reiter „Ganzer Loeser" gebaut). Hier wird sie deshalb selbst mitgeschrieben --
     * sonst zeichnet das Diagramm nur eine Kurve, waehrend die Legende zwei verspricht.
     * Befund 01.11.2026: genau das war der Fall, und die Knopfprobe hat es nicht
     * gemeldet, weil sie nur Farbtoene zaehlte. */
    gLauf = { M: M, LG: LG, ZG: LG.zustand(), LN: LN, ZN: ZN, tEnde: tEnde,
              halt: false, t0: performance.now(), nx: nx, ny: ny, bahnN: [] };
    h.textContent =
      `Galerkin: ${LG.N} Unbekannte, Knoten alle ${(LG.grenze.hmin*500).toFixed(1)} mm ` +
      `(P2: halbe Kantenlänge), dt = ${(LG.dt*1000).toFixed(2)} ms aus der Courant-Bedingung\n` +
      `Netzwerk: ${nx} × ${ny} Zellen, dx = ${(K.FALL.LX/nx*1000).toFixed(1)} mm, ` +
      `dt = ${(LN.dt*1000).toFixed(2)} ms — gleicher Knotenabstand, damit der Vergleich ` +
      `das Verfahren zeigt und nicht die Auflösung`;
    gHaeppchen();
  }).catch(e => {
    h.textContent = 'Der Vernetzer ließ sich nicht starten: ' + (e && e.message);
    document.getElementById('b-gstart').disabled = false;
    document.getElementById('b-gstop').disabled = true;
  });
}

function haltGalerkin() { if (gLauf) gLauf.halt = true; }

function gHaeppchen() {
  if (!gLauf) return;
  const bis = performance.now() + 150;
  while (performance.now() < bis && gLauf.ZG.t < gLauf.tEnde && !gLauf.halt) {
    gLauf.LG.schritt(gLauf.ZG);
    while (gLauf.ZN.t < gLauf.ZG.t && gLauf.ZN.t < gLauf.tEnde) {
      gLauf.LN.schritt(gLauf.ZN);
      gLauf.bahnN.push([gLauf.ZN.t, gLauf.ZN.cD, gLauf.ZN.cL]);
    }
  }
  const Z = gLauf.ZG, ZN = gLauf.ZN;
  /* Beide Felder an DENSELBEN Punkten und auf DERSELBEN Skala (Befund B16). Der
   * Netzwerkweg bildet die Wirbelstaerke auf den Gitterknoten (i*dx, j*dy); genau
   * dort wird auch das Galerkin-Feld ausgewertet, und zwar linear interpoliert --
   * nicht blockgemittelt, denn die Mittelung glaettet die Wandspitze weg, die der
   * Vergleich zeigen soll. */
  const g = gLauf.LN.g, nxA = g.nx - 1, nyA = g.ny - 1;
  const wN = K.wirbelstaerke(g, ZN.u, ZN.v);
  const wG = KF.wirbelAufGitter(gLauf.LG, Z.ux, Z.uy, nxA, nyA, false, g.nx, g.ny);
  const sk = Math.max(perzentil(wG, 0.99), perzentil(wN, 0.99));
  feldbild(document.getElementById('c-gfeld'), wG, nxA, nyA, sk);
  feldbild(document.getElementById('c-nfeld'), wN, nxA, nyA, sk);
  document.getElementById('b-gfeld').textContent = `t = ${Z.t.toFixed(3)} s · Skala ±${sk.toFixed(0)} 1/s`;
  document.getElementById('b-nfeld').textContent = `t = ${ZN.t.toFixed(3)} s · gleiche Skala`;
  gBahnBild(document.getElementById('c-gbahn'), Z.bahn, gLauf.bahnN);
  document.getElementById('b-gbahn').textContent =
    `${Z.k} Zeitschritte Galerkin · ${ZN.k} Netzwerk · ` +
    `${((performance.now()-gLauf.t0)/1000).toFixed(1)} s Rechenzeit`;
  const zn = v => (v === undefined || !isFinite(v)) ? '—' : v.toFixed(4);
  document.getElementById('tafel-galerkin').innerHTML =
    '<table class="tafel"><tr><th>Größe</th><th>Galerkin</th><th>Netzwerk</th><th>Unterschied</th></tr>' +
    `<tr><td>Widerstandsbeiwert c_D</td><td class="zahl">${zn(Z.cD)}</td>` +
    `<td class="zahl">${zn(ZN.cD)}</td><td class="zahl">` +
    `${(isFinite(Z.cD)&&isFinite(ZN.cD)) ? (100*Math.abs(Z.cD-ZN.cD)/Math.max(Math.abs(Z.cD),1e-9)).toFixed(1)+' %' : '—'}</td></tr>` +
    `<tr><td>Auftriebsbeiwert c_L</td><td class="zahl">${zn(Z.cL)}</td>` +
    `<td class="zahl">${zn(ZN.cL)}</td><td class="zahl">` +
    `${(isFinite(Z.cL)&&isFinite(ZN.cL)) ? Math.abs(Z.cL-ZN.cL).toFixed(4) : '—'}</td></tr>` +
    `<tr><td>Unbekannte</td><td class="zahl">${gLauf.LG.N}</td>` +
    `<td class="zahl">${(g.nx)*(g.ny)*3}</td><td>—</td></tr>` +
    `<tr><td>Zeitschritt</td><td class="zahl">${(gLauf.LG.dt*1000).toFixed(2)} ms</td>` +
    `<td class="zahl">${(gLauf.LN.dt*1000).toFixed(2)} ms</td><td>—</td></tr>` +
    `<tr><td>Löserschritte je Zeitschritt</td><td class="zahl">${Z.schritte} MINRES</td>` +
    `<td class="zahl">${ZN.cg === undefined ? '—' : ZN.cg + ' CG'}</td><td>—</td></tr>` +
    '</table>';
  const fertig = gLauf.ZG.t >= gLauf.tEnde || gLauf.halt;
  if (fertig) {
    const g1 = Z.cD, n1 = ZN.cD, mitte = (3.22 + 3.24) / 2;
    const beide = isFinite(g1) && isFinite(n1);
    const dazwischen = beide && Math.min(g1, n1) < mitte && mitte < Math.max(g1, n1);
    document.getElementById('h-gfazit').textContent =
      'Zwei verschieden hergeleitete Rechnungen, beide hier im Browser gerechnet —\n' +
      'gleicher Knotenabstand, gleiche Auswertepunkte, gleiche Farbskala.\n\n' +
      (beide
        ? `Auf diesem groben Netz liegt der Galerkin-Weg bei c_D = ${g1.toFixed(3)}, der\n` +
          `Netzwerkweg bei ${n1.toFixed(3)}. Der Benchmarkwert ist 3,22 … 3,24` +
          (dazwischen
            ? ' — er liegt DAZWISCHEN.\nBeide Wege sind noch nicht umgerechnet, und zwar in entgegengesetzte\n' +
              'Richtungen: Der Galerkin-Weg unterschätzt auf grobem Netz, der Netzwerkweg\n' +
              'überschätzt wegen der Treppenkontur. Das ist kein Widerspruch der Verfahren,\n' +
              'sondern der sichtbare Netzfehler beider — und er zeigt in verschiedene\n' +
              'Richtungen, weil die Herleitungen verschieden sind.\n'
            : '.\n') +
          '\nDie umgerechneten Werte stehen im Reiter „Verfahren & Nachweis": Galerkin auf\n' +
          'feinem Netz c_D 3,2318 (im Intervall), Netzwerk mit 144 320 Zellen +1,8 %.\n' +
          'Die feinen Netze sind im Browser nicht in Minuten zu rechnen — was hier zu\n' +
          'sehen ist, ist der VERGLEICH der beiden Wege bei gleicher Auflösung.'
        : 'Eine der beiden Rechnungen hat keinen endlichen Wert geliefert.');
  } else {
    document.getElementById('h-gfazit').textContent =
      `läuft … t = ${Z.t.toFixed(2)} s von ${gLauf.tEnde.toFixed(1)} s. Die Ablösung ` +
      `braucht mehrere Sekunden,\nbis sie als Wirbelstraße zu sehen ist.`;
  }
  if (!fertig) { requestAnimationFrame(gHaeppchen); return; }
  gLauf.halt = true;
  document.getElementById('b-gstart').disabled = false;
  document.getElementById('b-gstop').disabled = true;
}

// ================================================================ Die Programme
/* Auf Freds Einwand vom 01.11.2026: Eine Seite, die auf Programme verweist, die ihr
 * nicht beiliegen, ist fuer den Leser eine Sackgasse -- er kann glauben oder bleiben
 * lassen, aber nicht nachrechnen. Hier liegen sie im Wortlaut. */
function fuelleQuellen() {
  const q = (D.quellen || []);
  const s = document.getElementById('s-quelle');
  if (!s || !q.length) return;
  s.innerHTML = q.map((f, i) =>
    `<option value="${i}">${f.name} — ${f.zeilen} Zeilen</option>`).join('');
  /* Befund 01.11.2026: Die Auswahlliste zeigte immer nur EINE Datei. Wer sie nicht
   * aufklappte, musste denken, es liege nur ein Programm bei. „Alle Programme sind da"
   * muss man SEHEN, nicht aufklappen müssen -- deshalb stehen sie hier alle, nach
   * Ordner gegliedert, und ein Druck darauf zeigt den Quelltext. */
  const liste = document.getElementById('quellliste');
  const gruppen = {};
  q.forEach((f, i) => {
    const o = f.name.split('/')[0];
    (gruppen[o] || (gruppen[o] = [])).push([f, i]);
  });
  const titel = { stroemung: 'stroemung/ — die Rechnung in Python',
                  labor: 'labor/ — die Seite, ihre Rechenkerne und deren Proben' };
  liste.innerHTML = Object.keys(gruppen).map(o =>
    `<h4>${titel[o] || o + '/'}</h4>` +
    '<table class="quelltafel"><tr><th>Datei</th><th>Aufgabe</th>' +
    '<th class="zahl">Zeilen</th><th class="zahl">Byte</th><th>sha256</th></tr>' +
    gruppen[o].map(([f, i]) =>
      `<tr data-i="${i}" onclick="waehleQuelle(${i})" title="${f.sha256}">` +
      `<td><code>${f.name.split('/')[1]}</code></td>` +
      `<td>${f.aufgabe || ''}</td>` +
      `<td class="zahl">${f.zeilen}</td>` +
      `<td class="zahl">${f.bytes.toLocaleString('de-DE')}</td>` +
      `<td class="zahl pruefsumme">${f.sha256.slice(0, 12)}…</td></tr>`).join('') +
    '</table>').join('');
  let zeilen = 0, bytes = 0;
  q.forEach(f => { zeilen += f.zeilen; bytes += f.bytes; });
  document.getElementById('b-quellzahl').textContent =
    `${q.length} Dateien · ${zeilen.toLocaleString('de-DE')} Zeilen · ` +
    `${(bytes / 1024).toFixed(0)} kB`;
  zeigeQuelle();
}
function waehleQuelle(i) {
  document.getElementById('s-quelle').value = String(i);
  zeigeQuelle();
  document.getElementById('quelltext').scrollIntoView({ behavior: 'smooth', block: 'start' });
}
function zeigeQuelle() {
  const q = (D.quellen || []);
  const s = document.getElementById('s-quelle');
  if (!s || !q.length) return;
  const f = q[+s.value || 0];
  document.getElementById('quelltext').textContent = f.text;
  document.getElementById('h-quelle-kopf').textContent =
    `${f.name} · ${f.zeilen} Zeilen · ${f.bytes.toLocaleString('de-DE')} Byte\n` +
    `sha256 ${f.sha256}`;
  for (const r of document.querySelectorAll('#quellliste tr[data-i]'))
    r.classList.toggle('an', r.dataset.i === s.value);
}

window.addEventListener('DOMContentLoaded', () => {
  let r = 'labor';
  try { r = localStorage.getItem('cfd-reiter') || 'labor'; } catch (e) {}
  if (!window.PRUEFSTAND.reiter.includes(r)) r = 'labor';
  waehleVergleich('stempel');
  zeigeNetzmass(); zeigeTende(); zeigeGTende(); fuelleQuellen();
  zeige(r);
});
window.addEventListener('resize', () => {
  if (vergleich === 'stempel') rechneStempel();
  if (!document.getElementById('ansicht-zeitverlauf').hidden) zeigeZeit(null);
});
