/* pruefe_kern.js -- der Browserkern gegen die geprueften Python-Module.
 *
 * Aufruf:  node pruefe_kern.js
 *
 * Verglichen wird gegen vergleichsdaten.json, das vergleichsdaten.py mit den Modulen
 * aus ../stroemung/ erzeugt hat. Die Zufallsfelder stehen in der Datei, beide Seiten
 * rechnen also wirklich dasselbe. Eigener Zufall auf jeder Seite waere kein Vergleich.
 */
const fs = require('fs');
const path = require('path');
require(path.join(__dirname, 'kern_cfd.js'));
const K = globalThis.KernCFD;
const D = JSON.parse(fs.readFileSync(path.join(__dirname, 'vergleichsdaten.json'), 'utf8'));

const urteile = [];
function pruefe(text, a, b, schranke) {
  let d = 0, m = 0;
  for (let k = 0; k < a.length; k++) {
    d = Math.max(d, Math.abs(a[k] - b[k]));
    m = Math.max(m, Math.abs(b[k]));
  }
  const rel = d / Math.max(m, 1e-300);
  urteile.push([text, rel < schranke, rel, d]);
}

for (const f of D.faelle) {
  if (f.art === 'bernoulli') {
    pruefe('B(x) ueber 16 Stuetzstellen von 0 bis 1e6',
           f.x.map(K.bernoulli), f.y, 1e-14);
  } else if (f.art === 'loese1d') {
    const c = K.loese1d(f.n, f.L, f.u, f.D, f.stempel);
    pruefe(`1D-Loesung, Stempel ${f.stempel.padEnd(7)} Pe=${(f.u * f.L / f.D).toExponential(0)}`,
           c, f.c, 1e-11);
  } else if (f.art === 'genau1d') {
    pruefe(`geschlossene Loesung           Pe=${f.Pe.toExponential(0)}`,
           f.x.map(x => K.genau1d(x, f.L, f.Pe)), f.c, 1e-13);
  } else if (f.art === 'drucknetz_nackt') {
    const g = new K.Gitter(f.Lx, f.Ly, f.nx, f.ny);
    const p = Float64Array.from(f.p);
    const dn = K.drucknetzMatrix(g, false, null, null);
    pruefe(`${f.nx}x${f.ny} Drucknetz OHNE Ausstrom, Weg A (Bezugsknoten festgehalten)`,
           dn.A.mal(p), f.A, 1e-12);
    pruefe(`${f.nx}x${f.ny} Drucknetz OHNE Ausstrom, Weg B`,
           K.drucknetzWirkungVoll(g, p, dn.leer, false, null, null), f.B, 1e-12);
  } else if (f.art === 'felder') {
    const g = new K.Gitter(f.Lx, f.Ly, f.nx, f.ny);
    const u = Float64Array.from(f.u), v = Float64Array.from(f.v), p = Float64Array.from(f.p);
    const ein = Float64Array.from(f.ein);
    const mu = Uint8Array.from(f.mu), mv = Uint8Array.from(f.mv), leer = Uint8Array.from(f.leer);
    const kennz = `${f.nx}x${f.ny}`;
    pruefe(`${kennz} Divergenz (Knotenregel)`, g.divergenz(u, v), f.divergenz, 1e-12);
    pruefe(`${kennz} Weg B: Konvektion+Zaehigkeit u`, K.konvZaehU(g, u, v, f.nu, ein), f.konv_u, 1e-12);
    pruefe(`${kennz} Weg B: Konvektion+Zaehigkeit v`, K.konvZaehV(g, u, v, f.nu), f.konv_v, 1e-12);
    const Au = new K.ImpulsA(g, 'u', ein), Av = new K.ImpulsA(g, 'v', null);
    pruefe(`${kennz} Weg A: Impulsmatrix u`, Au.fuellen(u, v, f.nu).wirkung(u), f.impuls_u, 1e-12);
    pruefe(`${kennz} Weg A: Impulsmatrix v`, Av.fuellen(u, v, f.nu).wirkung(v), f.impuls_v, 1e-12);
    const dn = K.drucknetzMatrix(g, true, mu, mv);
    pruefe(`${kennz} Weg A: Drucknetzmatrix`, dn.A.mal(p), f.druck_A, 1e-12);
    pruefe(`${kennz} Weg B: Drucknetz als Feldoperation`,
           K.drucknetzWirkungVoll(g, p, leer, true, mu, mv), f.druck_B, 1e-12);
    let gl = 0;
    for (let k = 0; k < leer.length; k++) if (leer[k] !== f.leer[k]) gl++;
    urteile.push([`${kennz} dieselben Zellen gelten als „ganz im Koerper"`, gl === 0, gl, gl]);
  }
}

console.log('='.repeat(96));
console.log(` BROWSERKERN GEGEN PYTHON -- kern_cfd.js gegen ../stroemung/, Vergleichsdaten ${D.fassung}`);
console.log('='.repeat(96));
let n = 0;
for (const [text, ok, rel] of urteile) {
  if (ok) n++;
  console.log(`  [${ok ? 'RICHTIG' : 'FALSCH '}] ${text}   ${typeof rel === 'number' ? rel.toExponential(3) : rel}`);
}
console.log(`\n  == ${n}/${urteile.length}`);
process.exit(n === urteile.length ? 0 : 1);
