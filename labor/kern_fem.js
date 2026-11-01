/* kern_fem.js -- der UNABHAENGIGE Weg im Browser: Galerkin mit Taylor-Hood, Fassung 1.0
 *
 * Der zweite Rechenweg, im Browser gerechnet -- nicht abgespielt. Die Unabhaengigkeit
 * sitzt in der HERLEITUNG:
 *
 *   kern_cfd.js   Flussbilanz am Kontrollvolumen. Zelle = Knoten, Flaeche = Zweig,
 *                 Kontinuitaet = Knotenregel. Unbekannte sind Zellmittel.
 *   kern_fem.js   schwache Form. Die Gleichungen werden mit Testfunktionen
 *                 multipliziert und ueber das Gebiet integriert; Unbekannte sind
 *                 Knotenwerte eines Polynomansatzes.
 *
 *     Geschwindigkeit  P2 (quadratisch, Knoten auf Ecken UND Kantenmitten)
 *     Druck            P1 (linear, Knoten nur auf den Ecken)
 *
 * Diese Paarung erfuellt die inf-sup-Bedingung und braucht keinen abgestimmten
 * Stabilisierungsparameter. Dazu ein koerperangepasstes Netz: Der Zylinderrand LIEGT
 * auf Knoten, es gibt keine Treppe.
 *
 * Dass beide Wege am Ende ein duennbesetztes System loesen, ist eine Frage des
 * Werkzeugs, nicht des Ansatzes.
 *
 * Geprueft gegen fem_galerkin.py und fem_wirbelstrasse.py ueber
 * vergleichsdaten_fem.json (pruefe_kern_fem.js).
 */
(function (global) {
'use strict';

// ---------------------------------------------------------------- Kennzahlen des Falls
const F = { LX: 2.2, LY: 0.41, XC: 0.2, YC: 0.2, DZ: 0.1, UM: 1.5,
            NU: 1e-3, RHO: 1.0, U_QUER: 1.0 };
function einstroemprofil(y) { return 4.0 * F.UM * y * (F.LY - y) / (F.LY * F.LY); }

// Gausspunkte auf dem Dreieck, Grad 5, sieben Punkte -- dieselben wie in fem_galerkin.py
const _a = 0.470142064105115, _b = 0.101286507323456;
const QUAD = [
  [1/3, 1/3, 1/3, 0.225],
  [1 - 2*_a, _a, _a, 0.132394152788506],
  [_a, 1 - 2*_a, _a, 0.132394152788506],
  [_a, _a, 1 - 2*_a, 0.132394152788506],
  [1 - 2*_b, _b, _b, 0.125939180544827],
  [_b, 1 - 2*_b, _b, 0.125939180544827],
  [_b, _b, 1 - 2*_b, 0.125939180544827],
];

// ---------------------------------------------------------------- Formfunktionen
/* Die sechs P2-Formfunktionen an baryzentrischen Koordinaten L = (l0,l1,l2).
 * Nummerierung im Element: 0,1,2 = Ecken; 3 = Mitte(1,2), 4 = Mitte(2,0), 5 = Mitte(0,1). */
function p2Werte(L) {
  const l0 = L[0], l1 = L[1], l2 = L[2];
  return [l0*(2*l0-1), l1*(2*l1-1), l2*(2*l2-1), 4*l1*l2, 4*l2*l0, 4*l0*l1];
}
/* Gradienten der sechs Formfunktionen; g[i] = grad(l_i), je Element konstant. */
function p2Ableitungen(L, g) {
  const l0 = L[0], l1 = L[1], l2 = L[2];
  const d = [];
  d.push([(4*l0-1)*g[0][0], (4*l0-1)*g[0][1]]);
  d.push([(4*l1-1)*g[1][0], (4*l1-1)*g[1][1]]);
  d.push([(4*l2-1)*g[2][0], (4*l2-1)*g[2][1]]);
  d.push([4*(l2*g[1][0] + l1*g[2][0]), 4*(l2*g[1][1] + l1*g[2][1])]);
  d.push([4*(l0*g[2][0] + l2*g[0][0]), 4*(l0*g[2][1] + l2*g[0][1])]);
  d.push([4*(l1*g[0][0] + l0*g[1][0]), 4*(l1*g[0][1] + l0*g[1][1])]);
  return d;
}

/* grad(l_i) je Dreieck und die Flaechen -- beides je Element konstant (gerade Kanten). */
function lambdaGradienten(P, T) {
  const ne = T.length / 3;
  const g = new Float64Array(ne * 6);          // [e][i][d]
  const A = new Float64Array(ne);
  for (let e = 0; e < ne; e++) {
    const a = T[3*e], b = T[3*e+1], c = T[3*e+2];
    const x0 = P[2*a], y0 = P[2*a+1], x1 = P[2*b], y1 = P[2*b+1], x2 = P[2*c], y2 = P[2*c+1];
    const A2 = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0);
    A[e] = 0.5 * A2;
    g[6*e+0] = (y1-y2)/A2; g[6*e+1] = (x2-x1)/A2;
    g[6*e+2] = (y2-y0)/A2; g[6*e+3] = (x0-x2)/A2;
    g[6*e+4] = (y0-y1)/A2; g[6*e+5] = (x1-x0)/A2;
  }
  return { g: g, A: A };
}

/* Knoten des P2-Raums: erst die Ecken, dann die Kantenmitten -- genau wie p2_raum(). */
function p2Raum(P, T) {
  const nv = P.length / 2, ne = T.length / 3;
  const kanten = new Map();
  const E6 = new Int32Array(ne * 6);
  const paare = [[1, 2], [2, 0], [0, 1]];
  for (let e = 0; e < ne; e++) {
    for (let j = 0; j < 3; j++) {
      const a = T[3*e + paare[j][0]], b = T[3*e + paare[j][1]];
      const k = (a < b) ? (a * 4294967296 + b) : (b * 4294967296 + a);
      if (!kanten.has(k)) kanten.set(k, kanten.size);
    }
  }
  const nk = kanten.size;
  const XY = new Float64Array((nv + nk) * 2);
  XY.set(P.subarray ? P.subarray(0, 2*nv) : P.slice(0, 2*nv));
  for (const [k, i] of kanten) {
    const a = Math.floor(k / 4294967296), b = k % 4294967296;
    XY[2*(nv+i)]     = 0.5 * (P[2*a]   + P[2*b]);
    XY[2*(nv+i) + 1] = 0.5 * (P[2*a+1] + P[2*b+1]);
  }
  for (let e = 0; e < ne; e++) {
    E6[6*e] = T[3*e]; E6[6*e+1] = T[3*e+1]; E6[6*e+2] = T[3*e+2];
    for (let j = 0; j < 3; j++) {
      const a = T[3*e + paare[j][0]], b = T[3*e + paare[j][1]];
      const k = (a < b) ? (a * 4294967296 + b) : (b * 4294967296 + a);
      E6[6*e + 3 + j] = nv + kanten.get(k);
    }
  }
  return { XY: XY, E6: E6, nv: nv, nkanten: nk, n2: nv + nk };
}

// ---------------------------------------------------------------- duennbesetzte Matrix
/* Aufbau ueber Tripel, dann zu CSR mit zusammengefassten Doppelten -- wie coo_matrix. */
function csrAusTripeln(zeilen, spalten, werte, nz, ns) {
  const zahl = new Int32Array(nz + 1);
  for (let k = 0; k < zeilen.length; k++) zahl[zeilen[k] + 1]++;
  for (let i = 0; i < nz; i++) zahl[i + 1] += zahl[i];
  const ji = new Int32Array(zeilen.length), vi = new Float64Array(zeilen.length);
  const platz = zahl.slice(0, nz);
  for (let k = 0; k < zeilen.length; k++) {
    const p = platz[zeilen[k]]++;
    ji[p] = spalten[k]; vi[p] = werte[k];
  }
  // je Zeile sortieren und Doppelte addieren
  const zp = new Int32Array(nz + 1), zj = [], zv = [];
  for (let i = 0; i < nz; i++) {
    const a = zahl[i], b = zahl[i + 1];
    const ord = [];
    for (let k = a; k < b; k++) ord.push(k);
    ord.sort((p, q) => ji[p] - ji[q]);
    let letzte = -1;
    for (const k of ord) {
      if (ji[k] === letzte) zv[zv.length - 1] += vi[k];
      else { zj.push(ji[k]); zv.push(vi[k]); letzte = ji[k]; }
    }
    zp[i + 1] = zj.length;
  }
  return { zp: zp, zj: new Int32Array(zj), zv: new Float64Array(zv), nz: nz, ns: ns };
}
function matMal(S, x, y) {
  y = y || new Float64Array(S.nz);
  for (let i = 0; i < S.nz; i++) {
    let s = 0;
    for (let k = S.zp[i]; k < S.zp[i + 1]; k++) s += S.zv[k] * x[S.zj[k]];
    y[i] = s;
  }
  return y;
}
/* y += alpha * S^T x -- gebraucht fuer die Druckkopplung und die Kraft. */
function matTMal(S, x, y, alpha) {
  alpha = (alpha === undefined) ? 1 : alpha;
  y = y || new Float64Array(S.ns);
  for (let i = 0; i < S.nz; i++) {
    const xi = alpha * x[i];
    if (xi === 0) continue;
    for (let k = S.zp[i]; k < S.zp[i + 1]; k++) y[S.zj[k]] += S.zv[k] * xi;
  }
  return y;
}
function zeilensummen(S) {
  const r = new Float64Array(S.nz);
  for (let i = 0; i < S.nz; i++) { let s = 0; for (let k = S.zp[i]; k < S.zp[i+1]; k++) s += S.zv[k]; r[i] = s; }
  return r;
}
function summe(S) { let s = 0; for (let k = 0; k < S.zv.length; k++) s += S.zv[k]; return s; }
function quadratsumme(S) { let s = 0; for (let k = 0; k < S.zv.length; k++) s += S.zv[k]*S.zv[k]; return s; }

// ---------------------------------------------------------------- Zusammenstellen
/* Massenmatrix, Steifigkeit und die Druck-Geschwindigkeits-Kopplung.
 *
 *   (du/dt, v) + c(u,u,v) + nu (grad u, grad v) - (p, div v) = 0
 *   (q, div u) = 0
 *
 * Die Konvektion c wird NICHT als Matrix gebaut, sondern je Zeitschritt als Vektor
 * ausgewertet -- damit bleibt die Systemmatrix konstant. */
function zusammenstellen(P, T, R2) {
  const ne = T.length / 3, n2 = R2.n2, n1 = R2.nv, E6 = R2.E6;
  const lg = lambdaGradienten(P, T), g = lg.g, A = lg.A;
  const Me = new Float64Array(ne * 36), Ke = new Float64Array(ne * 36);
  const Bxe = new Float64Array(ne * 18), Bye = new Float64Array(ne * 18);
  const dN = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]];
  for (const q of QUAD) {
    const L = [q[0], q[1], q[2]], w = q[3];
    const N = p2Werte(L);
    for (let e = 0; e < ne; e++) {
      const go = 6*e;
      const g0x = g[go], g0y = g[go+1], g1x = g[go+2], g1y = g[go+3], g2x = g[go+4], g2y = g[go+5];
      dN[0][0] = (4*L[0]-1)*g0x; dN[0][1] = (4*L[0]-1)*g0y;
      dN[1][0] = (4*L[1]-1)*g1x; dN[1][1] = (4*L[1]-1)*g1y;
      dN[2][0] = (4*L[2]-1)*g2x; dN[2][1] = (4*L[2]-1)*g2y;
      dN[3][0] = 4*(L[2]*g1x + L[1]*g2x); dN[3][1] = 4*(L[2]*g1y + L[1]*g2y);
      dN[4][0] = 4*(L[0]*g2x + L[2]*g0x); dN[4][1] = 4*(L[0]*g2y + L[2]*g0y);
      dN[5][0] = 4*(L[1]*g0x + L[0]*g1x); dN[5][1] = 4*(L[1]*g0y + L[0]*g1y);
      const wA = w * A[e];
      const mo = 36*e, bo = 18*e;
      for (let i = 0; i < 6; i++) {
        for (let j = 0; j < 6; j++) {
          Me[mo + 6*i + j] += wA * N[i] * N[j];
          Ke[mo + 6*i + j] += wA * (dN[i][0]*dN[j][0] + dN[i][1]*dN[j][1]);
        }
      }
      for (let k = 0; k < 3; k++) {
        for (let i = 0; i < 6; i++) {
          Bxe[bo + 6*k + i] += wA * L[k] * dN[i][0];
          Bye[bo + 6*k + i] += wA * L[k] * dN[i][1];
        }
      }
    }
  }
  // streuen
  const zi = new Int32Array(ne*36), sj = new Int32Array(ne*36);
  for (let e = 0; e < ne; e++) for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) {
    zi[36*e + 6*i + j] = E6[6*e + i]; sj[36*e + 6*i + j] = E6[6*e + j];
  }
  const bi = new Int32Array(ne*18), bj = new Int32Array(ne*18);
  for (let e = 0; e < ne; e++) for (let k = 0; k < 3; k++) for (let i = 0; i < 6; i++) {
    bi[18*e + 6*k + i] = T[3*e + k]; bj[18*e + 6*k + i] = E6[6*e + i];
  }
  return {
    M: csrAusTripeln(zi, sj, Me, n2, n2),
    K: csrAusTripeln(zi, sj, Ke, n2, n2),
    Bx: csrAusTripeln(bi, bj, Bxe, n1, n2),
    By: csrAusTripeln(bi, bj, Bye, n1, n2),
    g: g, A: A, E6: E6, n2: n2, n1: n1, ne: ne,
  };
}

/* c(u,u,v) als Vektor: integral (u . grad)u . v. Je Zeitschritt neu, ohne Matrix. */
function konvektion(S, ux, uy) {
  const ne = S.ne, E6 = S.E6, g = S.g, A = S.A;
  const fx = new Float64Array(S.n2), fy = new Float64Array(S.n2);
  const dN = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]];
  for (const q of QUAD) {
    const L = [q[0], q[1], q[2]], w = q[3];
    const N = p2Werte(L);
    for (let e = 0; e < ne; e++) {
      const go = 6*e, eo = 6*e;
      const g0x = g[go], g0y = g[go+1], g1x = g[go+2], g1y = g[go+3], g2x = g[go+4], g2y = g[go+5];
      dN[0][0] = (4*L[0]-1)*g0x; dN[0][1] = (4*L[0]-1)*g0y;
      dN[1][0] = (4*L[1]-1)*g1x; dN[1][1] = (4*L[1]-1)*g1y;
      dN[2][0] = (4*L[2]-1)*g2x; dN[2][1] = (4*L[2]-1)*g2y;
      dN[3][0] = 4*(L[2]*g1x + L[1]*g2x); dN[3][1] = 4*(L[2]*g1y + L[1]*g2y);
      dN[4][0] = 4*(L[0]*g2x + L[2]*g0x); dN[4][1] = 4*(L[0]*g2y + L[2]*g0y);
      dN[5][0] = 4*(L[1]*g0x + L[0]*g1x); dN[5][1] = 4*(L[1]*g0y + L[0]*g1y);
      let uq = 0, vq = 0, dxdx = 0, dxdy = 0, dydx = 0, dydy = 0;
      for (let i = 0; i < 6; i++) {
        const n = E6[eo + i], a = ux[n], b = uy[n];
        uq += a * N[i]; vq += b * N[i];
        dxdx += a * dN[i][0]; dxdy += a * dN[i][1];
        dydx += b * dN[i][0]; dydy += b * dN[i][1];
      }
      const wA = w * A[e];
      const cx = wA * (uq*dxdx + vq*dxdy), cy = wA * (uq*dydx + vq*dydy);
      for (let i = 0; i < 6; i++) { const n = E6[eo + i]; fx[n] += cx * N[i]; fy[n] += cy * N[i]; }
    }
  }
  return { fx: fx, fy: fy };
}

/* Welche P2-Knoten liegen wo. Der Ausstrom bleibt FREI (do-nothing).
 * Wie in fem_galerkin.py Fassung 1.1: Wand bleibt Wand bis zum letzten Knoten --
 * die Eckknoten am Ausstrom liegen auf der Wand und behalten die Haftbedingung. */
function raender(XY, n2, mitZylinder) {
  const e = 1e-9;
  const ein = new Uint8Array(n2), wand = new Uint8Array(n2), zyl = new Uint8Array(n2);
  const aus = new Uint8Array(n2), fest = new Uint8Array(n2);
  const u = new Float64Array(n2), v = new Float64Array(n2);
  for (let i = 0; i < n2; i++) {
    const x = XY[2*i], y = XY[2*i+1];
    ein[i] = x < e ? 1 : 0;
    wand[i] = (y < e || y > F.LY - e) ? 1 : 0;
    zyl[i] = (mitZylinder && Math.hypot(x - F.XC, y - F.YC) < F.DZ/2 + 1e-7) ? 1 : 0;
    aus[i] = x > F.LX - e ? 1 : 0;
    fest[i] = (ein[i] || wand[i] || zyl[i]) ? 1 : 0;
    if (ein[i]) u[i] = einstroemprofil(y);
    if (wand[i] || zyl[i]) u[i] = 0.0;
  }
  return { ein: ein, wand: wand, zyl: zyl, aus: aus, fest: fest, u: u, v: v };
}

// ---------------------------------------------------------------- Das Sattelsystem
/* Die Blockmatrix des Sattelpunktproblems, SYMMETRISIERT:
 *
 *     [  A    0   -Bx^T ]
 *     [  0    A   -By^T ]
 *     [ -Bx  -By    0   ]
 *
 * A = M/dt + nu K (SPD). In fem_galerkin.py steht die letzte Blockzeile mit +B; beide
 * Systeme haben DIESELBE Loesung -- die dritte Blockzeile ist dort q*div(u) = 0, hier
 * mit -1 durchmultipliziert. Der Grund: So ist die Matrix symmetrisch und indefinit,
 * und darauf laesst sich MINRES anwenden. Ein unsymmetrisches System braeuchte GMRES
 * mit wachsendem Speicher. Das ist eine Frage des Werkzeugs, nicht des Ansatzes.
 */
function sattelsystem(S, dt) {
  const n2 = S.n2, n1 = S.n1, N = 2*n2 + n1;
  // A = M/dt + nu K, beide haben dieselbe Besetzung (gleiche Tripel) -- Werte addieren
  const A = { zp: S.M.zp, zj: S.M.zj, zv: new Float64Array(S.M.zv.length), nz: n2, ns: n2 };
  for (let k = 0; k < A.zv.length; k++) A.zv[k] = S.M.zv[k]/dt + F.NU * S.K.zv[k];
  return { A: A, Bx: S.Bx, By: S.By, n2: n2, n1: n1, N: N, dt: dt };
}

/* y = Sys * x, ohne die Blockmatrix je aufzustellen. */
function sattelMal(Sy, x, y) {
  const n2 = Sy.n2, n1 = Sy.n1;
  y = y || new Float64Array(Sy.N);
  y.fill(0);
  const xu = x.subarray(0, n2), xv = x.subarray(n2, 2*n2), xp = x.subarray(2*n2);
  const yu = y.subarray(0, n2), yv = y.subarray(n2, 2*n2), yp = y.subarray(2*n2);
  matMal(Sy.A, xu, yu);
  matMal(Sy.A, xv, yv);
  matTMal(Sy.Bx, xp, yu, -1);          // -Bx^T p
  matTMal(Sy.By, xp, yv, -1);          // -By^T p
  const t1 = matMal(Sy.Bx, xu), t2 = matMal(Sy.By, xv);
  for (let i = 0; i < n1; i++) yp[i] = -(t1[i] + t2[i]);
  return y;
}

/* Dirichlet-Knoten symmetrisch eliminieren -- dieselbe Rechnung wie festsetzen() in
 * fem_galerkin.py, nur ohne Matrizenprodukte: Zeile UND Spalte leeren, Diagonale auf
 * eins. Die bekannten Werte wandern vorher als Spaltenbeitrag in die rechte Seite. */
function festsetzen(Sy, feste, werte) {
  const n2 = Sy.n2, N = Sy.N;
  const A = Sy.A;
  // Kopie, damit das urspruengliche System fuer die Kraft unangetastet bleibt
  const Af = { zp: A.zp, zj: A.zj, zv: Float64Array.from(A.zv), nz: A.nz, ns: A.ns };
  for (let i = 0; i < n2; i++) {
    if (feste[i]) {
      for (let k = Af.zp[i]; k < Af.zp[i+1]; k++) Af.zv[k] = (Af.zj[k] === i) ? 1 : 0;
    } else {
      for (let k = Af.zp[i]; k < Af.zp[i+1]; k++) if (feste[Af.zj[k]]) Af.zv[k] = 0;
    }
  }
  const streich = M => {
    const Mf = { zp: M.zp, zj: M.zj, zv: Float64Array.from(M.zv), nz: M.nz, ns: M.ns };
    for (let i = 0; i < Mf.nz; i++)
      for (let k = Mf.zp[i]; k < Mf.zp[i+1]; k++) if (feste[Mf.zj[k]]) Mf.zv[k] = 0;
    return Mf;
  };
  return { A: Af, Bx: streich(Sy.Bx), By: streich(Sy.By),
           n2: n2, n1: Sy.n1, N: N, dt: Sy.dt, feste: feste, werte: werte };
}

/* Die rechte Seite auf dieselbe Weise behandeln: Spaltenbeitrag der bekannten Werte
 * abziehen, dann an den festen Stellen den Wert selbst eintragen. */
function rechteSeite(Sy, rhs, feste, werte) {
  const bekannt = new Float64Array(Sy.N);
  for (let i = 0; i < Sy.N; i++) if (feste[i]) bekannt[i] = werte[i];
  const bei = sattelMal(Sy, bekannt);
  const r = new Float64Array(Sy.N);
  for (let i = 0; i < Sy.N; i++) r[i] = feste[i] ? werte[i] : rhs[i] - bei[i];
  return r;
}

// ---------------------------------------------------------------- Vorkonditionierer
/* Unvollstaendige Cholesky-Zerlegung ohne Zusatzbesetzung, IC(0) -- dasselbe Mittel wie
 * im Netzwerkkern. A = M/dt + nu K ist bei kleinem dt massendominiert und damit stark
 * diagonaldominant; IC(0) ist dort sehr gut. */
function ic0(A) {
  const n = A.nz;
  const zp = A.zp, zj = A.zj;
  const L = Float64Array.from(A.zv);
  const dia = new Int32Array(n);
  for (let i = 0; i < n; i++) {
    dia[i] = -1;
    for (let k = zp[i]; k < zp[i+1]; k++) if (zj[k] === i) { dia[i] = k; break; }
  }
  for (let i = 0; i < n; i++) {
    for (let k = zp[i]; k < zp[i+1]; k++) {
      const j = zj[k];
      if (j > i) break;
      let s = L[k];
      // Summe ueber gemeinsame Spalten kleiner j
      let p = zp[i], q = zp[j];
      while (p < zp[i+1] && q < zp[j+1]) {
        const a = zj[p], b = zj[q];
        if (a >= j || b >= j) break;
        if (a === b) { s -= L[p] * L[q]; p++; q++; }
        else if (a < b) p++;
        else q++;
      }
      if (j === i) {
        L[k] = Math.sqrt(Math.max(s, 1e-300));
      } else {
        L[k] = s / L[dia[j]];
      }
    }
    // oberes Dreieck wird nicht gebraucht
    for (let k = zp[i]; k < zp[i+1]; k++) if (zj[k] > i) L[k] = 0;
  }
  return { zp: zp, zj: zj, zv: L, n: n, dia: dia };
}
/* z = (L L^T)^-1 r -- Vorwaerts-, dann Rueckwaertseinsetzen.
 *
 * L steht zeilenweise (CSR). Fuer L^T braeuchte man die SPALTE i, die in Zeilenform
 * nicht direkt dasteht. Der uebliche Griff: von unten nach oben laufen, z[i] durch die
 * Diagonale teilen und den Beitrag SOFORT auf alle j < i verteilen -- das ist genau
 * dieselbe Rechnung, nur andersherum sortiert. */
function ic0Loese(P_, r, z) {
  const n = P_.n, zp = P_.zp, zj = P_.zj, L = P_.zv, dia = P_.dia;
  z = z || new Float64Array(n);
  for (let i = 0; i < n; i++) {                     // L y = r
    let s = r[i];
    for (let k = zp[i]; k < dia[i]; k++) s -= L[k] * z[zj[k]];
    z[i] = s / L[dia[i]];
  }
  for (let i = n - 1; i >= 0; i--) {                // L^T z = y
    z[i] /= L[dia[i]];
    const zi = z[i];
    for (let k = zp[i]; k < dia[i]; k++) z[zj[k]] -= L[k] * zi;
  }
  return z;
}

/* Druckraum P1: Massenmatrix (zusammengezogen auf die Diagonale) und Laplace.
 * Beide werden NUR fuer den Vorkonditionierer gebraucht -- sie stehen in keiner
 * Gleichung des Verfahrens und aendern die Loesung nicht, nur die Zahl der Schritte
 * bis dorthin. */
function druckraum(P, T, lg) {
  const ne = T.length / 3, n1 = P.length / 2;
  const mlump = new Float64Array(n1);
  const zi = new Int32Array(ne*9), sj = new Int32Array(ne*9), v = new Float64Array(ne*9);
  for (let e = 0; e < ne; e++) {
    const A = lg.A[e], go = 6*e;
    const g = [[lg.g[go], lg.g[go+1]], [lg.g[go+2], lg.g[go+3]], [lg.g[go+4], lg.g[go+5]]];
    for (let k = 0; k < 3; k++) {
      mlump[T[3*e+k]] += A / 3;
      for (let l = 0; l < 3; l++) {
        zi[9*e + 3*k + l] = T[3*e+k]; sj[9*e + 3*k + l] = T[3*e+l];
        v[9*e + 3*k + l] = A * (g[k][0]*g[l][0] + g[k][1]*g[l][1]);
      }
    }
  }
  const Kp = csrAusTripeln(zi, sj, v, n1, n1);
  // kleine Verschiebung, damit IC(0) nicht am konstanten Nullraum scheitert
  let spur = 0;
  for (let i = 0; i < n1; i++)
    for (let k = Kp.zp[i]; k < Kp.zp[i+1]; k++) if (Kp.zj[k] === i) spur += Kp.zv[k];
  const eps = 1e-8 * spur / n1;
  for (let i = 0; i < n1; i++)
    for (let k = Kp.zp[i]; k < Kp.zp[i+1]; k++) if (Kp.zj[k] === i) Kp.zv[k] += eps;
  return { mlump: mlump, Kp: Kp };
}

/* MINRES fuer das symmetrische, indefinite Sattelsystem, mit Blockvorkonditionierer
 *
 *     Pv = diag( A^,  A^,  S^ )
 *
 * A^  = IC(0) von A  -- A ist bei kleinem dt massendominiert, das greift sehr gut
 * S^  = nu * Mp^-1 + dt * Kp^-1   (nach Cahouet und Chabard)
 *
 * Die Form von S^ ist nicht geraten: Der Schur-Anteil des Sattelsystems verhaelt sich
 * bei grossem Zeitschritt wie der zaehe Anteil (nu durch Masse), bei kleinem wie dt mal
 * Laplace. Diese Summe deckt beide Enden ab. Ein Vorkonditionierer aendert die Loesung
 * NICHT -- nur die Zahl der Schritte bis dorthin. Falsch liegen kann man damit nicht,
 * nur langsam sein; und wie langsam, steht in der Schrittzahl, die der Loeser meldet.
 *
 * Der Ablauf ist das Verfahren von Paige und Saunders: Lanczos mit Vorkonditionierer,
 * daraus eine Folge von Givens-Drehungen, die das kleinste Residuum ueber dem bisher
 * aufgespannten Raum trifft. Von Hand nachrechenbar -- jede Zeile ist ein Skalarprodukt
 * oder eine Linearkombination.
 */
function minres(Sy, rhs, V, opt) {
  opt = opt || {};
  const N = Sy.N, n2 = Sy.n2, n1 = Sy.n1;
  const schranke = (opt.schranke === undefined) ? 1e-10 : opt.schranke;
  const hoechstens = opt.hoechstens || 500;
  const punkt = (a, b) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i]*b[i]; return s; };

  const anwenden = (r, z) => {                        // z = Pv^-1 r
    z = z || new Float64Array(N);
    ic0Loese(V.A, r.subarray(0, n2), z.subarray(0, n2));
    ic0Loese(V.A, r.subarray(n2, 2*n2), z.subarray(n2, 2*n2));
    const rp = r.subarray(2*n2), zp_ = z.subarray(2*n2);
    const h = ic0Loese(V.Kp, rp);
    /* Befund 01.11.2026: Hier stand `Sy.dt * h[i]`. Richtig ist `h[i] / Sy.dt` --
     * das Verhalten bei kleinem Zeitschritt ist S ~ dt * Kp, also S^-1 ~ Kp^-1 / dt.
     * Mit dt = 2e-3 lag der Druckanteil um dt^2 = 4e-6 daneben; MINRES lief in jedem
     * Zeitschritt in die Schrittgrenze und brach UNFERTIG ab. Die Widerstandsbeiwerte
     * sahen trotzdem fast richtig aus (6e-6), der Auftrieb nicht (63 %) -- weil cL die
     * kleine Differenz zweier grosser Zahlen ist und als erstes zerfaellt. */
    for (let i = 0; i < n1; i++) zp_[i] = F.NU * rp[i] / V.mlump[i] + h[i] / Sy.dt;
    return z;
  };

  const x = opt.start ? Float64Array.from(opt.start) : new Float64Array(N);
  let vAlt = new Float64Array(N);
  let v = new Float64Array(N);
  const Ax = sattelMal(Sy, x);
  for (let i = 0; i < N; i++) v[i] = rhs[i] - Ax[i];
  let z = anwenden(v);
  let gamma = Math.sqrt(Math.max(punkt(v, z), 0));
  const gamma0 = gamma;
  if (gamma0 < 1e-300) return { x: x, schritte: 0, rest: 0 };
  let gammaAlt = 1;
  let eta = gamma, s0 = 0, s1 = 0, c0 = 1, c1 = 1;
  let w = new Float64Array(N), w1 = new Float64Array(N), w2 = new Float64Array(N);
  let k = 0, rest = 1;
  for (k = 1; k <= hoechstens; k++) {
    const zs = new Float64Array(N);
    for (let i = 0; i < N; i++) zs[i] = z[i] / gamma;
    const Az = sattelMal(Sy, zs);
    const delta = punkt(Az, zs);
    const vNeu = new Float64Array(N);
    const f1 = delta / gamma, f2 = gamma / gammaAlt;
    for (let i = 0; i < N; i++) vNeu[i] = Az[i] - f1 * v[i] - f2 * vAlt[i];
    const zNeu = anwenden(vNeu);
    const gammaNeu = Math.sqrt(Math.max(punkt(vNeu, zNeu), 0));
    const a0 = c1 * delta - c0 * s1 * gamma;
    const a1 = Math.sqrt(a0 * a0 + gammaNeu * gammaNeu);
    const a2 = s1 * delta + c0 * c1 * gamma;
    const a3 = s0 * gamma;
    c0 = c1; c1 = (a1 === 0) ? 1 : a0 / a1;
    s0 = s1; s1 = (a1 === 0) ? 0 : gammaNeu / a1;
    for (let i = 0; i < N; i++) w[i] = (zs[i] - a3 * w2[i] - a2 * w1[i]) / (a1 || 1);
    const ce = c1 * eta;
    for (let i = 0; i < N; i++) x[i] += ce * w[i];
    const tmp = w2; w2 = w1; w1 = w; w = tmp;
    eta = -s1 * eta;
    rest = Math.abs(eta) / gamma0;
    if (rest < schranke) { k++; break; }
    vAlt = v; v = vNeu; z = zNeu; gammaAlt = gamma; gamma = gammaNeu;
    if (gamma < 1e-300) { k++; break; }
  }
  return { x: x, schritte: k - 1, rest: rest };
}

/* Die groesste Zeitschrittweite, bei der das Verfahren noch stabil ist.
 *
 * Zaehigkeit und Druck sind implizit und setzen keine Grenze. Die Konvektion ist
 * EXPLIZIT, und dafuer gilt die Courant-Bedingung: In einem Zeitschritt darf die
 * Stroemung nicht weiter laufen als von einem Knoten zum naechsten. Bei P2 liegen die
 * Knoten auf halber Kantenlaenge, der Abstand ist also h/2.
 *
 *     dt <= C * (h_min/2) / u_max
 *
 * Befund 01.11.2026: Diese Grenze fehlte. Die Pruefung lief nur 40 Zeitschritte
 * (t = 0,08 s) und war gruen; im Dauerlauf platzte dieselbe Rechnung bei t = 0,372 s.
 * Gemessen: CFL 0,31 platzt, CFL 0,15 laeuft durch. Mit C = 0,15 liegt ein Sicherheits-
 * abstand von zwei dazwischen. Eine Pruefung, die vor der kennzeichnenden Zeit des
 * Vorgangs aufhoert, sagt ueber Stabilitaet nichts. */
function zeitschrittGrenze(P, T, C) {
  C = (C === undefined) ? 0.15 : C;
  let hmin = Infinity;
  const paare = [[0,1],[1,2],[2,0]];
  for (let e = 0; e < T.length/3; e++) {
    for (const pr of paare) {
      const i = T[3*e+pr[0]], j = T[3*e+pr[1]];
      const d = Math.hypot(P[2*i]-P[2*j], P[2*i+1]-P[2*j+1]);
      if (d < hmin) hmin = d;
    }
  }
  const umax = 4.0 * F.UM / 4.0 * 1.0;        // Hoechstgeschwindigkeit im Einstrom = UM
  return { dt: C * (hmin/2) / F.UM, hmin: hmin, umax: F.UM };
}

// ---------------------------------------------------------------- Der Zeitschritt
/* Halbimplizit: Zaehigkeit und Druck implizit, Konvektion explizit -- genau wie in
 * fem_wirbelstrasse.py. Damit ist die Systemmatrix konstant und wird nur EINMAL
 * aufgebaut und vorkonditioniert; je Zeitschritt bleibt ein Loesen. */
function Loeser(P, T, dt, mitZylinder) {
  const grenze = zeitschrittGrenze(P, T);
  const gekappt = (dt === undefined || dt === null || dt > grenze.dt);
  if (dt === undefined || dt === null) dt = grenze.dt;
  const R2 = p2Raum(P, T);
  const S = zusammenstellen(P, T, R2);
  const lg = { g: S.g, A: S.A };
  const RB = raender(R2.XY, R2.n2, mitZylinder !== false);
  const Sy = sattelsystem(S, dt);
  const n2 = S.n2, n1 = S.n1, N = Sy.N;
  const feste = new Uint8Array(N), werte = new Float64Array(N);
  for (let i = 0; i < n2; i++) {
    feste[i] = RB.fest[i]; werte[i] = RB.u[i];
    feste[n2 + i] = RB.fest[i]; werte[n2 + i] = RB.v[i];
  }
  const Sf = festsetzen(Sy, feste, werte);
  const dr = druckraum(P, T, lg);
  const V = { A: ic0(Sf.A), Kp: ic0(dr.Kp), mlump: dr.mlump };
  // Testfunktion fuer die Kraft: eins auf dem Zylinder, null sonst
  const vz = new Float64Array(n2);
  for (let i = 0; i < n2; i++) vz[i] = RB.zyl[i] ? 1 : 0;
  const nenner = 0.5 * F.RHO * F.U_QUER * F.U_QUER * F.DZ;
  return {
    P: P, T: T, R2: R2, S: S, RB: RB, Sy: Sy, Sf: Sf, V: V, vz: vz, nenner: nenner,
    n2: n2, n1: n1, N: N, dt: dt, feste: feste, werte: werte,
    grenze: grenze, zuGross: gekappt,
    zustand: function () {
      const ux = new Float64Array(n2), uy = new Float64Array(n2);
      for (let i = 0; i < n2; i++) ux[i] = RB.fest[i] ? RB.u[i] : 0;
      return { ux: ux, uy: uy, p: new Float64Array(n1), t: 0, k: 0,
               cD: 0, cL: 0, schritte: 0, rest: 0, bahn: [] };
    },
    schritt: function (Z) {
      const kv = konvektion(S, Z.ux, Z.uy);
      const rhs = new Float64Array(N);
      const mu = matMal(S.M, Z.ux), mv = matMal(S.M, Z.uy);
      for (let i = 0; i < n2; i++) {
        rhs[i] = mu[i] / dt - kv.fx[i];
        rhs[n2 + i] = mv[i] / dt - kv.fy[i];
      }
      const r = rechteSeite(Sy, rhs, feste, werte);
      const start = new Float64Array(N);
      start.set(Z.ux, 0); start.set(Z.uy, n2); start.set(Z.p, 2*n2);
      const l = minres(Sf, r, V, { start: start, schranke: 1e-11, hoechstens: 400 });
      const uxn = l.x.subarray(0, n2), uyn = l.x.subarray(n2, 2*n2), pr = l.x.subarray(2*n2);
      // Konsistente Kraft: Residuum der Impulsgleichung, getestet mit vz.
      // Das ist nicht das Aufsummieren von Wandspannungen, sondern das Randintegral
      // selbst -- deutlich genauer, und ohne jede Gemeinsamkeit mit der Impulsbilanz
      // der abgedeckten Zweige im Netzwerkweg.
      const dux = new Float64Array(n2), duy = new Float64Array(n2);
      for (let i = 0; i < n2; i++) { dux[i] = (uxn[i] - Z.ux[i]) / dt; duy[i] = (uyn[i] - Z.uy[i]) / dt; }
      const Mdx = matMal(S.M, dux), Mdy = matMal(S.M, duy);
      const Kux = matMal(S.K, uxn), Kuy = matMal(S.K, uyn);
      const Bxp = matTMal(S.Bx, pr, new Float64Array(n2), 1);
      const Byp = matTMal(S.By, pr, new Float64Array(n2), 1);
      let Fx = 0, Fy = 0;
      for (let i = 0; i < n2; i++) {
        if (vz[i] === 0) continue;
        Fx -= (Mdx[i] + F.NU * Kux[i] + kv.fx[i] - Bxp[i]) * vz[i];
        Fy -= (Mdy[i] + F.NU * Kuy[i] + kv.fy[i] - Byp[i]) * vz[i];
      }
      Z.ux = Float64Array.from(uxn); Z.uy = Float64Array.from(uyn); Z.p = Float64Array.from(pr);
      Z.k += 1; Z.t = Z.k * dt;
      Z.cD = Fx / nenner; Z.cL = Fy / nenner;
      Z.schritte = l.schritte; Z.rest = l.rest;
      Z.bahn.push([Z.t, Z.cD, Z.cL]);
      return Z;
    },
  };
}

/* Wirbelstaerke des Galerkin-Feldes auf ein kartesisches Anzeigegitter.
 * omega = dv/dx - du/dy ist bei P2 je Element linear und ueber die Kanten unstetig;
 * fuer die Anzeige wird es flaechengewichtet auf die Ecken gemittelt und dann im
 * Dreieck baryzentrisch interpoliert. Das ist eine Anzeigeoperation, keine
 * Rechenoperation -- genau wie in fem_wirbelstrasse.py. */
function wirbelAufGitter(L_, ux, uy, nx, ny, mitte, nxG, nyG) {
  nxG = nxG || nx; nyG = nyG || ny;
  const P = L_.P, T = L_.T, S = L_.S, E6 = S.E6, g = S.g, A = S.A, ne = S.ne;
  const Lc = [1/3, 1/3, 1/3];
  const omE = new Float64Array(ne);
  for (let e = 0; e < ne; e++) {
    const go = 6*e;
    const g0x = g[go], g0y = g[go+1], g1x = g[go+2], g1y = g[go+3], g2x = g[go+4], g2y = g[go+5];
    const dN = [
      [(4*Lc[0]-1)*g0x, (4*Lc[0]-1)*g0y], [(4*Lc[1]-1)*g1x, (4*Lc[1]-1)*g1y],
      [(4*Lc[2]-1)*g2x, (4*Lc[2]-1)*g2y],
      [4*(Lc[2]*g1x + Lc[1]*g2x), 4*(Lc[2]*g1y + Lc[1]*g2y)],
      [4*(Lc[0]*g2x + Lc[2]*g0x), 4*(Lc[0]*g2y + Lc[2]*g0y)],
      [4*(Lc[1]*g0x + Lc[0]*g1x), 4*(Lc[1]*g0y + Lc[0]*g1y)]];
    let dvdx = 0, dudy = 0;
    for (let i = 0; i < 6; i++) { const n = E6[6*e+i]; dvdx += uy[n]*dN[i][0]; dudy += ux[n]*dN[i][1]; }
    omE[e] = dvdx - dudy;
  }
  const nv = P.length / 2;
  const zahl = new Float64Array(nv), sum = new Float64Array(nv);
  for (let e = 0; e < ne; e++) for (let k = 0; k < 3; k++) {
    sum[T[3*e+k]] += omE[e] * A[e]; zahl[T[3*e+k]] += A[e];
  }
  const omV = new Float64Array(nv);
  for (let i = 0; i < nv; i++) omV[i] = sum[i] / Math.max(zahl[i], 1e-30);
  // Punktsuche ueber ein Kaestchenverzeichnis -- sonst waere es nx*ny*ne
  const mx = 120, my = 40;
  const faecher = new Array(mx*my);
  for (let e = 0; e < ne; e++) {
    let x0 = 1e9, x1 = -1e9, y0 = 1e9, y1 = -1e9;
    for (let k = 0; k < 3; k++) {
      const n = T[3*e+k], x = P[2*n], y = P[2*n+1];
      if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y;
    }
    const i0 = Math.max(0, Math.floor(x0/F.LX*mx)), i1 = Math.min(mx-1, Math.floor(x1/F.LX*mx));
    const j0 = Math.max(0, Math.floor(y0/F.LY*my)), j1 = Math.min(my-1, Math.floor(y1/F.LY*my));
    for (let i = i0; i <= i1; i++) for (let j = j0; j <= j1; j++) {
      const f = i*my + j; (faecher[f] || (faecher[f] = [])).push(e);
    }
  }
  /* Wo ausgewertet wird, entscheidet der Aufrufer -- und das ist kein Nebenpunkt.
   * Befund B16: Dreimal in Folge sah ein Vergleich nach einem grossen Unterschied der
   * VERFAHREN aus, und jedes Mal lag es an der Auswertung. Der Netzwerkweg bildet die
   * Wirbelstaerke auf den GITTERKNOTEN (i*dx, j*dy), nicht in den Zellmitten. Wer hier
   * Zellmitten nimmt und dort Knoten, vergleicht zwei verschiedene Orte. */
  const versatz = (mitte === false) ? 1.0 : 0.5;
  const Z = new Float64Array(nx*ny);
  for (let i = 0; i < nx; i++) {
    const x = (i + versatz) * F.LX / nxG;
    const fi = Math.min(mx-1, Math.max(0, Math.floor(x/F.LX*mx)));
    for (let j = 0; j < ny; j++) {
      const y = (j + versatz) * F.LY / nyG;
      const fj = Math.min(my-1, Math.max(0, Math.floor(y/F.LY*my)));
      const liste = faecher[fi*my + fj];
      let wert = 0;
      if (liste) for (const e of liste) {
        const a = T[3*e], b = T[3*e+1], c = T[3*e+2];
        const xa = P[2*a], ya = P[2*a+1], xb = P[2*b], yb = P[2*b+1], xc = P[2*c], yc = P[2*c+1];
        const d = (yb-yc)*(xa-xc) + (xc-xb)*(ya-yc);
        const l0 = ((yb-yc)*(x-xc) + (xc-xb)*(y-yc)) / d;
        const l1 = ((yc-ya)*(x-xc) + (xa-xc)*(y-yc)) / d;
        const l2 = 1 - l0 - l1;
        if (l0 >= -1e-12 && l1 >= -1e-12 && l2 >= -1e-12) {
          wert = l0*omV[a] + l1*omV[b] + l2*omV[c];
          break;
        }
      }
      Z[i*ny + j] = wert;
    }
  }
  return Z;
}

global.KERN_FEM = {
  sattelsystem: sattelsystem, sattelMal: sattelMal, festsetzen: festsetzen,
  rechteSeite: rechteSeite, ic0: ic0, ic0Loese: ic0Loese, druckraum: druckraum,
  minres: minres, Loeser: Loeser, wirbelAufGitter: wirbelAufGitter,
  zeitschrittGrenze: zeitschrittGrenze,
  F: F, QUAD: QUAD, einstroemprofil: einstroemprofil,
  p2Werte: p2Werte, p2Ableitungen: p2Ableitungen, lambdaGradienten: lambdaGradienten,
  p2Raum: p2Raum, zusammenstellen: zusammenstellen, konvektion: konvektion,
  raender: raender,
  csrAusTripeln: csrAusTripeln, matMal: matMal, matTMal: matTMal,
  zeilensummen: zeilensummen, summe: summe, quadratsumme: quadratsumme,
};
})(typeof window !== 'undefined' ? window : globalThis);
if (typeof module !== 'undefined') module.exports = (typeof window !== 'undefined' ? window : globalThis).KERN_FEM;
