/* kern_cfd.js -- Rechenkern des CFD-Labors, Fassung 1.1 (31.10.2026)
 *
 * Dieselben Verfahren wie die Python-Programme in ../stroemung/, Zeile fuer Zeile
 * nachgebaut. Der Kern liegt ABSICHTLICH in einer eigenen Datei: Er wird gegen die
 * Python-Fassung geprueft (pruefe_kern.js gegen vergleichsdaten.json), bevor er in die
 * Seite eingebaut wird. Ein zweiter, ungeprueffter Loeser im Browser waere kein zweiter
 * Weg, sondern eine zweite Fehlerquelle.
 *
 * Felder liegen zeilenweise (row-major) in Float64Array:
 *     p  (nx, ny)      p[i*ny + j]
 *     u  (nx+1, ny)    u[i*ny + j]
 *     v  (nx, ny+1)    v[i*(ny+1) + j]
 */
(function (global) {
'use strict';

// ===================================================================================
//  1  Der Zweigstempel
// ===================================================================================

/* B(x) = x/(e^x - 1), numerisch sauber ueber den ganzen Bereich.
   Vier Bereiche wie in sg_stempel.py -- ohne sie ist der Stempel genau dort falsch,
   wo er gebraucht wird. */
function bernoulli(x) {
  if (Math.abs(x) < 1e-10) return 1.0 - x / 2.0;   // Reihe 1 - x/2 + x^2/12 - ...
  if (x > 500.0) return 0.0;                        // x/e^x -> 0
  if (x < -500.0) return -x;                        // x/(0-1) = -x
  return x / Math.expm1(x);                         // expm1: keine Ausloeschung
}

/* Zweigleitwerte aus der FLAECHENGESCHWINDIGKEIT (nicht aus dem Volumenstrom!).
   Pe = u_f*d/nu ; g_hin - g_rueck = u_f*A ist die gesteuerte Quelle. */
function sgLeitwerte(uf, nu, d, A) {
  const Pe = uf * d / nu, s = nu * A / d;
  return [s * bernoulli(-Pe), s * bernoulli(Pe)];
}

// ===================================================================================
//  2  Stationaere 1D-Konvektion-Diffusion -- drei Stempel gegen die exakte Loesung
// ===================================================================================

function thomas(a, b, c, d) {               // Tridiagonalloeser, n Unbekannte
  const n = d.length, cp = new Float64Array(n), dp = new Float64Array(n);
  cp[0] = c[0] / b[0]; dp[0] = d[0] / b[0];
  for (let i = 1; i < n; i++) {
    const m = b[i] - a[i] * cp[i - 1];
    cp[i] = c[i] / m; dp[i] = (d[i] - a[i] * dp[i - 1]) / m;
  }
  const x = new Float64Array(n); x[n - 1] = dp[n - 1];
  for (let i = n - 2; i >= 0; i--) x[i] = dp[i] - cp[i] * x[i + 1];
  return x;
}

/* n Zellen auf [0,L], c(0)=0, c(L)=1. Bilanziert wird um die inneren Knoten --
   Kirchhoffs Knotenregel. Der Stempel bestimmt nur, WIE der Zweigfluss entsteht. */
function loese1d(n, L, u, D, stempel) {
  const h = L / n, N = n + 1, Pe_h = u * h / D;
  let gh, gr;
  if (stempel === 'sg')          { gh = D / h * bernoulli(-Pe_h); gr = D / h * bernoulli(Pe_h); }
  else if (stempel === 'zentral'){ gh = D / h + u / 2.0;          gr = D / h - u / 2.0; }
  else if (stempel === 'aufwind'){ gh = D / h + Math.max(u, 0);   gr = D / h + Math.max(-u, 0); }
  else throw new Error('unbekannter Stempel: ' + stempel);
  const a = new Float64Array(N), b = new Float64Array(N),
        c = new Float64Array(N), d = new Float64Array(N);
  b[0] = 1; d[0] = 0;
  for (let i = 1; i < n; i++) { a[i] = -gh; b[i] = gh + gr; c[i] = -gr; d[i] = 0; }
  b[n] = 1; d[n] = 1;
  return thomas(a, b, c, d);
}

/* c(x) = (e^{Pe x/L} - 1)/(e^{Pe} - 1), ueberlaufsicher. */
function genau1d(x, L, Pe) {
  if (Math.abs(Pe) < 1e-12) return x / L;
  if (Pe > 500) return Math.exp(Pe * (x / L - 1.0));
  if (Pe < -500) return Math.expm1(Pe * x / L);
  return Math.expm1(Pe * x / L) / Math.expm1(Pe);
}

// ===================================================================================
//  3  Das versetzte Gitter
// ===================================================================================

function Gitter(Lx, Ly, nx, ny) {
  this.Lx = Lx; this.Ly = Ly; this.nx = nx; this.ny = ny;
  this.dx = Lx / nx; this.dy = Ly / ny;
  this.xp = new Float64Array(nx); this.yp = new Float64Array(ny);
  for (let i = 0; i < nx; i++) this.xp[i] = (i + 0.5) * this.dx;
  for (let j = 0; j < ny; j++) this.yp[j] = (j + 0.5) * this.dy;
  this.xu = new Float64Array(nx + 1); this.yv = new Float64Array(ny + 1);
  for (let i = 0; i <= nx; i++) this.xu[i] = i * this.dx;
  for (let j = 0; j <= ny; j++) this.yv[j] = j * this.dy;
  this.yu = this.yp; this.xv = this.xp;
}
Gitter.prototype.leerU = function () { return new Float64Array((this.nx + 1) * this.ny); };
Gitter.prototype.leerV = function () { return new Float64Array(this.nx * (this.ny + 1)); };
Gitter.prototype.leerP = function () { return new Float64Array(this.nx * this.ny); };

/* Knotenregel je Zelle: Summe der Zweigstroeme durch die vier Flaechen. */
Gitter.prototype.divergenz = function (u, v) {
  const nx = this.nx, ny = this.ny, dx = this.dx, dy = this.dy, d = new Float64Array(nx * ny);
  for (let i = 0; i < nx; i++)
    for (let j = 0; j < ny; j++)
      d[i * ny + j] = (u[(i + 1) * ny + j] - u[i * ny + j]) / dx
                    + (v[i * (ny + 1) + j + 1] - v[i * (ny + 1) + j]) / dy;
  return d;
};

// ===================================================================================
//  4  Geisterschichten -- genau wie _rand_u / _rand_v in konvektion.py
// ===================================================================================

/* U hat die Form (nx+3, ny+2); U[(1+i)*(ny+2) + (1+j)] = u[i*ny+j]. */
function randU(g, u, einstroem) {
  const nx = g.nx, ny = g.ny, W = ny + 2, U = new Float64Array((nx + 3) * W);
  for (let i = 0; i <= nx; i++)
    for (let j = 0; j < ny; j++) U[(1 + i) * W + 1 + j] = u[i * ny + j];
  for (let i = 0; i <= nx; i++) {
    U[(1 + i) * W + 0] = -u[i * ny + 0];                 // Haftbedingung unten
    U[(1 + i) * W + ny + 1] = -u[i * ny + ny - 1];       // oben
  }
  for (let b = 0; b < W; b++) {
    U[0 * W + b] = einstroem
      ? (b === 0 || b === W - 1 ? 0 : einstroem[b - 1])  // vorgegebenes Profil
      : U[1 * W + b];                                     // frei
    U[(nx + 2) * W + b] = U[(nx + 1) * W + b];            // freier Ausstrom
  }
  return U;
}

/* V hat die Form (nx+2, ny+3); V[(1+i)*(ny+3) + (1+j)] = v[i*(ny+1)+j]. */
function randV(g, v) {
  const nx = g.nx, ny = g.ny, W = ny + 3, m = ny + 1, V = new Float64Array((nx + 2) * W);
  for (let i = 0; i < nx; i++)
    for (let j = 0; j < m; j++) V[(1 + i) * W + 1 + j] = v[i * m + j];
  for (let i = 0; i < nx; i++) {
    V[(1 + i) * W + 0] = -v[i * m + 0];
    V[(1 + i) * W + ny + 2] = -v[i * m + m - 1];
  }
  for (let b = 0; b < W; b++) {
    V[0 * W + b] = V[1 * W + b];
    V[(nx + 1) * W + b] = V[nx * W + b];
  }
  return V;
}

// ===================================================================================
//  5  Weg B -- Konvektion und Zaehigkeit als Feldoperation, ohne jede Matrix
// ===================================================================================

function konvZaehU(g, u, v, nu, einstroem) {
  const nx = g.nx, ny = g.ny, dx = g.dx, dy = g.dy;
  const U = randU(g, u, einstroem), V = randV(g, v), WU = ny + 2, WV = ny + 3;
  const Jx = new Float64Array((nx + 2) * ny), Jy = new Float64Array((nx + 1) * (ny + 1));
  for (let a = 0; a <= nx + 1; a++)
    for (let j = 0; j < ny; j++) {
      const uf = 0.5 * (U[a * WU + 1 + j] + U[(a + 1) * WU + 1 + j]);
      const gg = sgLeitwerte(uf, nu, dx, dy);
      Jx[a * ny + j] = gg[0] * U[a * WU + 1 + j] - gg[1] * U[(a + 1) * WU + 1 + j];
    }
  for (let a = 0; a <= nx; a++)
    for (let j = 0; j <= ny; j++) {
      const vf = 0.5 * (V[a * WV + 1 + j] + V[(a + 1) * WV + 1 + j]);
      const gg = sgLeitwerte(vf, nu, dy, dx);
      Jy[a * (ny + 1) + j] = gg[0] * U[(1 + a) * WU + j] - gg[1] * U[(1 + a) * WU + j + 1];
    }
  const r = new Float64Array((nx + 1) * ny), Vz = dx * dy;
  for (let i = 0; i <= nx; i++)
    for (let j = 0; j < ny; j++)
      r[i * ny + j] = (Jx[i * ny + j] - Jx[(i + 1) * ny + j]
                     + Jy[i * (ny + 1) + j] - Jy[i * (ny + 1) + j + 1]) / Vz;
  return r;
}

function konvZaehV(g, u, v, nu) {
  const nx = g.nx, ny = g.ny, dx = g.dx, dy = g.dy, m = ny + 1;
  const U = randU(g, u, null), V = randV(g, v), WU = ny + 2, WV = ny + 3;
  const Jx = new Float64Array((nx + 1) * m), Jy = new Float64Array(nx * (ny + 2));
  for (let a = 0; a <= nx; a++)
    for (let j = 0; j < m; j++) {
      const uf = 0.5 * (U[(1 + a) * WU + j] + U[(1 + a) * WU + j + 1]);
      const gg = sgLeitwerte(uf, nu, dx, dy);
      Jx[a * m + j] = gg[0] * V[a * WV + 1 + j] - gg[1] * V[(a + 1) * WV + 1 + j];
    }
  for (let i = 0; i < nx; i++)
    for (let b = 0; b <= ny + 1; b++) {
      const vf = 0.5 * (V[(1 + i) * WV + b] + V[(1 + i) * WV + b + 1]);
      const gg = sgLeitwerte(vf, nu, dy, dx);
      Jy[i * (ny + 2) + b] = gg[0] * V[(1 + i) * WV + b] - gg[1] * V[(1 + i) * WV + b + 1];
    }
  const r = new Float64Array(nx * m), Vz = dx * dy;
  for (let i = 0; i < nx; i++)
    for (let j = 0; j < m; j++)
      r[i * m + j] = (Jx[i * m + j] - Jx[(i + 1) * m + j]
                    + Jy[i * (ny + 2) + j] - Jy[i * (ny + 2) + j + 1]) / Vz;
  return r;
}

// ===================================================================================
//  6  Duennbesetzte Matrizen -- Weg A stempelt Zweig fuer Zweig
// ===================================================================================

/* Sammelt (Zeile, Spalte, Wert) und baut daraus CSR. Das Muster (welche Zweige es
   gibt) ist die TOPOLOGIE des Netzes und aendert sich nie; nur die Leitwerte haengen
   vom Stroemungszustand ab. Deshalb wird das Muster einmal gebaut und danach nur noch
   neu gefuellt -- in genau derselben Stempelreihenfolge. */
function CSR(n) { this.n = n; this.zi = []; this.zj = []; }
CSR.prototype.ein = function (i, j) { this.zi.push(i); this.zj.push(j); };
CSR.prototype.fertig = function () {
  const n = this.n, zahl = new Int32Array(n);
  const spalten = new Array(n);
  for (let r = 0; r < n; r++) spalten[r] = new Map();
  for (let k = 0; k < this.zi.length; k++) spalten[this.zi[k]].set(this.zj[k], 0);
  this.indptr = new Int32Array(n + 1);
  for (let r = 0; r < n; r++) {
    const s = Array.from(spalten[r].keys()).sort((a, b) => a - b);
    s.forEach((c, ix) => spalten[r].set(c, ix));
    zahl[r] = s.length; this.indptr[r + 1] = this.indptr[r] + s.length;
    this.sortiert = this.sortiert || [];
    this.sortiert[r] = s;
  }
  const nnz = this.indptr[n];
  this.indices = new Int32Array(nnz); this.data = new Float64Array(nnz);
  for (let r = 0; r < n; r++)
    for (let ix = 0; ix < this.sortiert[r].length; ix++)
      this.indices[this.indptr[r] + ix] = this.sortiert[r][ix];
  this.slot = new Int32Array(this.zi.length);
  for (let k = 0; k < this.zi.length; k++)
    this.slot[k] = this.indptr[this.zi[k]] + spalten[this.zi[k]].get(this.zj[k]);
  this.zi = null; this.zj = null; this.sortiert = null;
  return this;
};
CSR.prototype.leeren = function () { this.data.fill(0); this.k = 0; };
CSR.prototype.schreibe = function (w) { this.data[this.slot[this.k++]] += w; };
CSR.prototype.mal = function (x, y) {
  const n = this.n, ip = this.indptr, ic = this.indices, d = this.data;
  y = y || new Float64Array(n);
  for (let r = 0; r < n; r++) {
    let s = 0;
    for (let k = ip[r]; k < ip[r + 1]; k++) s += d[k] * x[ic[k]];
    y[r] = s;
  }
  return y;
};
CSR.prototype.diagonale = function () {
  const n = this.n, dg = new Float64Array(n);
  for (let r = 0; r < n; r++)
    for (let k = this.indptr[r]; k < this.indptr[r + 1]; k++)
      if (this.indices[k] === r) dg[r] = this.data[k];
  return dg;
};

/* Weg A fuer das Drucknetz: jeder innere Zweig traegt den Leitwert A_f/d_f.
   Zweige, deren Flaeche im Koerper liegt, werden AUFGETRENNT. Die letzte Zellspalte
   haengt ueber dy/(dx/2) am Bezugspotential -- sonst ist das Netz stromlos aufgehaengt
   und die Projektion kann den Ausstrom nicht korrigieren. */
function drucknetzMatrix(g, ausstromRechts, uMaske, vMaske) {
  const nx = g.nx, ny = g.ny, N = nx * ny, nr = (i, j) => i * ny + j;
  const A = new CSR(N), grad = new Int32Array(N), paare = [];
  const zweig = (i1, j1, i2, j2, G) => {
    const a = nr(i1, j1), b = nr(i2, j2);
    A.ein(a, a); A.ein(a, b); A.ein(b, b); A.ein(b, a);
    paare.push(G, -G, G, -G);
    grad[a]++; grad[b]++;
  };
  for (let i = 0; i < nx; i++)
    for (let j = 0; j < ny; j++) {
      if (i + 1 < nx && (!uMaske || !uMaske[(i + 1) * ny + j])) zweig(i, j, i + 1, j, g.dy / g.dx);
      if (j + 1 < ny && (!vMaske || !vMaske[i * (ny + 1) + j + 1])) zweig(i, j, i, j + 1, g.dx / g.dy);
    }
  const leer = new Uint8Array(N);
  for (let r = 0; r < N; r++) if (grad[r] === 0) leer[r] = 1;
  if (ausstromRechts) for (let j = 0; j < ny; j++) { A.ein(nr(nx - 1, j), nr(nx - 1, j)); paare.push(g.dy / (g.dx / 2)); }
  for (let r = 0; r < N; r++) if (leer[r]) { A.ein(r, r); paare.push(1.0); }
  if (!ausstromRechts) A.ein(0, 0), paare.push(0);   // Platz fuer den Bezugsknoten
  A.fertig(); A.leeren();
  for (const w of paare) A.schreibe(w);
  if (!ausstromRechts) {
    // Reines Neumann-Problem: EIN Knoten wird festgehalten, sonst ist das Netz
    // stromlos aufgehaengt und die Matrix singulaer. Genau wie drucknetz_matrix in
    // ns_loeser.py. Fassung 1.1, Befund B11 -- in Fassung 1.0 fehlte das, und die
    // Kernprobe konnte es nicht sehen, weil sie nur den Fall MIT Ausstrom pruefte.
    for (let k = A.indptr[0]; k < A.indptr[1]; k++) A.data[k] = (A.indices[k] === 0) ? 1 : 0;
  }
  return { A: A, leer: leer, bezug: ausstromRechts ? -1 : 0 };
}

/* Weg B fuer dasselbe Netz: Flaechenfluesse als ganze Felder, ohne jede Matrix.
   r[i] = Summe ueber die Zweige am Knoten i von G*(p_i - p_Nachbar) = (+A*p). */
function drucknetzWirkungVoll(g, p, leer, ausstromRechts, uMaske, vMaske, out) {
  const nx = g.nx, ny = g.ny, dx = g.dx, dy = g.dy, Gx = dy / dx, Gy = dx / dy;
  const r = out || new Float64Array(nx * ny);
  r.fill(0);
  for (let i = 1; i < nx; i++)
    for (let j = 0; j < ny; j++) {
      if (uMaske && uMaske[i * ny + j]) continue;      // Zweig im Koerper: aufgetrennt
      const F = Gx * (p[i * ny + j] - p[(i - 1) * ny + j]);
      r[i * ny + j] += F; r[(i - 1) * ny + j] -= F;
    }
  for (let i = 0; i < nx; i++)
    for (let j = 1; j < ny; j++) {
      if (vMaske && vMaske[i * (ny + 1) + j]) continue;
      const F = Gy * (p[i * ny + j] - p[i * ny + j - 1]);
      r[i * ny + j] += F; r[i * ny + j - 1] -= F;
    }
  if (ausstromRechts)
    for (let j = 0; j < ny; j++) r[(nx - 1) * ny + j] += (dy / (dx / 2)) * p[(nx - 1) * ny + j];
  if (leer) for (let k = 0; k < r.length; k++) if (leer[k]) r[k] = p[k];
  // Der Bezugsknoten des reinen Neumann-Falls wird hier ABSICHTLICH nicht gesetzt:
  // Weg B kennt keine Einheitszeile, er bilanziert nur. Beide Wege unterscheiden sich
  // daher genau an diesem einen Knoten -- und genau dort wird nicht verglichen,
  // so wie in der Sprachprobe von ns_loeser.py.
  return r;
}

/* Weg A fuer den Impuls: Muster einmal, Werte je Zeitschritt.
   Geisterwerte, die Kopien (freier Ausstrom) oder Spiegelungen (Haftbedingung, -1)
   innerer Unbekannter sind, stehen als MATRIXEINTRAG auf dem Knoten, von dem sie
   abhaengen -- nicht als Konstante in b (Befund B6). */
function stempelImpulsU(g, u, v, nu, ein, nimm, nimmB) {
  const nx = g.nx, ny = g.ny, dx = g.dx, dy = g.dy;
  const U = randU(g, u, ein), V = randV(g, v), WU = ny + 2, WV = ny + 3, iVz = 1 / (dx * dy);
  const nr = (i, j) => i * ny + j;
  const ghx = (a, j) => sgLeitwerte(0.5 * (U[a * WU + 1 + j] + U[(a + 1) * WU + 1 + j]), nu, dx, dy);
  const ghy = (a, j) => sgLeitwerte(0.5 * (V[a * WV + 1 + j] + V[(a + 1) * WV + 1 + j]), nu, dy, dx);
  for (let i = 0; i <= nx; i++)
    for (let j = 0; j < ny; j++) {
      const a = nr(i, j);
      let G = ghx(i, j), gh = G[0] * iVz, gr = G[1] * iVz;
      nimm(a, a, -gr);
      if (i >= 1) nimm(a, nr(i - 1, j), gh);
      else if (!ein) nimm(a, nr(0, j), gh);
      else nimmB(a, gh * ein[j]);
      G = ghx(i + 1, j); gh = G[0] * iVz; gr = G[1] * iVz;
      nimm(a, a, -gh);
      nimm(a, nr(i <= nx - 1 ? i + 1 : nx, j), gr);
      G = ghy(i, j); gh = G[0] * iVz; gr = G[1] * iVz;
      nimm(a, a, -gr);
      if (j >= 1) nimm(a, nr(i, j - 1), gh); else nimm(a, nr(i, 0), -gh);
      G = ghy(i, j + 1); gh = G[0] * iVz; gr = G[1] * iVz;
      nimm(a, a, -gh);
      if (j <= ny - 2) nimm(a, nr(i, j + 1), gr); else nimm(a, nr(i, ny - 1), -gr);
    }
}

function stempelImpulsV(g, u, v, nu, nimm) {
  const nx = g.nx, ny = g.ny, dx = g.dx, dy = g.dy, m = ny + 1;
  const U = randU(g, u, null), V = randV(g, v), WU = ny + 2, WV = ny + 3, iVz = 1 / (dx * dy);
  const nr = (i, j) => i * m + j;
  const ghx = (a, j) => sgLeitwerte(0.5 * (U[(1 + a) * WU + j] + U[(1 + a) * WU + j + 1]), nu, dx, dy);
  const ghy = (i, b) => sgLeitwerte(0.5 * (V[(1 + i) * WV + b] + V[(1 + i) * WV + b + 1]), nu, dy, dx);
  for (let i = 0; i < nx; i++)
    for (let j = 0; j < m; j++) {
      const a = nr(i, j);
      let G = ghx(i, j), gh = G[0] * iVz, gr = G[1] * iVz;
      nimm(a, a, -gr);
      nimm(a, nr(i >= 1 ? i - 1 : 0, j), gh);
      G = ghx(i + 1, j); gh = G[0] * iVz; gr = G[1] * iVz;
      nimm(a, a, -gh);
      nimm(a, nr(i <= nx - 2 ? i + 1 : nx - 1, j), gr);
      G = ghy(i, j); gh = G[0] * iVz; gr = G[1] * iVz;
      nimm(a, a, -gr);
      if (j >= 1) nimm(a, nr(i, j - 1), gh); else nimm(a, nr(i, 0), -gh);
      G = ghy(i, j + 1); gh = G[0] * iVz; gr = G[1] * iVz;
      nimm(a, a, -gh);
      if (j <= m - 2) nimm(a, nr(i, j + 1), gr); else nimm(a, nr(i, m - 1), -gr);
    }
}

/* Baut Muster UND Werte. Beim ersten Aufruf entsteht das Muster, danach wird nur
   noch gefuellt -- in derselben Stempelreihenfolge, also in dieselben Plaetze. */
function ImpulsA(g, welch, ein) {
  this.g = g; this.welch = welch; this.ein = ein;
  const n = welch === 'u' ? (g.nx + 1) * g.ny : g.nx * (g.ny + 1);
  this.M = new CSR(n); this.b = new Float64Array(n); this.n = n;
  const u0 = g.leerU(), v0 = g.leerV();
  const st = welch === 'u'
    ? (nimm, nimmB) => stempelImpulsU(g, u0, v0, 1, ein, nimm, nimmB)
    : (nimm) => stempelImpulsV(g, u0, v0, 1, nimm);
  st((i, j) => this.M.ein(i, j), () => {});
  this.M.fertig();
}
ImpulsA.prototype.fuellen = function (u, v, nu) {
  this.M.leeren(); this.b.fill(0);
  const schreib = (i, j, w) => this.M.schreibe(w);
  if (this.welch === 'u') stempelImpulsU(this.g, u, v, nu, this.ein, schreib, (a, w) => { this.b[a] += w; });
  else stempelImpulsV(this.g, u, v, nu, schreib);
  return this;
};
ImpulsA.prototype.wirkung = function (x, out) {
  const y = this.M.mal(x, out);
  for (let k = 0; k < y.length; k++) y[k] += this.b[k];
  return y;
};

// ===================================================================================
//  7  Konjugierte Gradienten mit Diagonalvorkonditionierung
// ===================================================================================

/* Unvollstaendige Cholesky-Zerlegung IC(0): L mit demselben Besetzungsmuster wie das
   untere Dreieck von A. Das Drucknetz aendert sich nie -- die Zerlegung wird also EINMAL
   gebaut und dient danach in jedem Zeitschritt als Vorkonditionierer. Sie ersetzt den
   Loeser nicht, sie beschleunigt ihn: Geloest wird weiter mit konjugierten Gradienten,
   das Ergebnis ist dasselbe, nur mit weniger Schritten. */
function ic0(A) {
  const n = A.n, ip = A.indptr, ic = A.indices, d = A.data;
  const lp = new Int32Array(n + 1), zahl = new Int32Array(n);
  for (let i = 0; i < n; i++) {
    let c = 0;
    for (let k = ip[i]; k < ip[i + 1]; k++) if (ic[k] <= i) c++;
    zahl[i] = c; lp[i + 1] = lp[i] + c;
  }
  const li = new Int32Array(lp[n]), lv = new Float64Array(lp[n]);
  for (let i = 0, t = 0; i < n; i++)
    for (let k = ip[i]; k < ip[i + 1]; k++)
      if (ic[k] <= i) { li[t] = ic[k]; lv[t] = d[k]; t++; }
  for (let i = 0; i < n; i++) {
    for (let k = lp[i]; k < lp[i + 1]; k++) {
      const j = li[k];
      let s = lv[k];
      if (j < i) {
        let a = lp[i], bq = lp[j];
        while (a < k && bq < lp[j + 1] - 1) {              // nur Spalten < j
          if (li[a] === li[bq]) { s -= lv[a] * lv[bq]; a++; bq++; }
          else if (li[a] < li[bq]) a++; else bq++;
        }
        const dj = lv[lp[j + 1] - 1];
        lv[k] = dj !== 0 ? s / dj : 0;
      } else {                                             // Diagonale
        let a = lp[i];
        while (a < k) { s -= lv[a] * lv[a]; a++; }
        lv[k] = s > 1e-300 ? Math.sqrt(s) : Math.sqrt(Math.abs(d[ip[i]]) + 1e-300);
      }
    }
  }
  const y = new Float64Array(n);
  return function (r, z) {                                 // loest L L^T z = r
    for (let i = 0; i < n; i++) {
      let s = r[i];
      for (let k = lp[i]; k < lp[i + 1] - 1; k++) s -= lv[k] * y[li[k]];
      y[i] = s / lv[lp[i + 1] - 1];
    }
    for (let i = 0; i < n; i++) z[i] = y[i];
    for (let i = n - 1; i >= 0; i--) {
      z[i] /= lv[lp[i + 1] - 1];
      for (let k = lp[i]; k < lp[i + 1] - 1; k++) z[li[k]] -= lv[k] * z[i];
    }
    return z;
  };
}

function cg(matvec, b, x, diag, tol, maxiter) {
  const n = b.length, r = new Float64Array(n), z = new Float64Array(n),
        p = new Float64Array(n), Ap = new Float64Array(n);
  matvec(x, Ap);
  let nb = 0;
  for (let k = 0; k < n; k++) { r[k] = b[k] - Ap[k]; nb += b[k] * b[k]; }
  nb = Math.sqrt(nb) || 1;
  const vor = (typeof diag === 'function')
    ? diag
    : (rr, zz) => { for (let k = 0; k < n; k++) zz[k] = rr[k] / (diag[k] || 1); return zz; };
  let rz = 0;
  vor(r, z);
  for (let k = 0; k < n; k++) { p[k] = z[k]; rz += r[k] * z[k]; }
  let it = 0, rest = 0;
  for (; it < maxiter; it++) {
    rest = 0; for (let k = 0; k < n; k++) rest += r[k] * r[k];
    if (Math.sqrt(rest) / nb < tol) break;
    matvec(p, Ap);
    let pAp = 0; for (let k = 0; k < n; k++) pAp += p[k] * Ap[k];
    if (!isFinite(pAp) || pAp === 0) break;
    const al = rz / pAp;
    for (let k = 0; k < n; k++) { x[k] += al * p[k]; r[k] -= al * Ap[k]; }
    let rz2 = 0;
    vor(r, z);
    for (let k = 0; k < n; k++) rz2 += r[k] * z[k];
    const be = rz2 / rz; rz = rz2;
    for (let k = 0; k < n; k++) p[k] = z[k] + be * p[k];
  }
  return { x: x, schritte: it, rest: Math.sqrt(rest) / nb };
}

// ===================================================================================
//  8  Der Benchmarkfall und der Zeitschritt -- beide Wege nebeneinander
// ===================================================================================

const FALL = { LX: 2.2, LY: 0.41, XC: 0.2, YC: 0.2, DZ: 0.1, UM: 1.5, NU: 1e-3, RHO: 1.0 };
FALL.U_QUER = 2.0 / 3.0 * FALL.UM;                 // = 1,0
FALL.RE = FALL.U_QUER * FALL.DZ / FALL.NU;         // = 100

function einstroemprofil(y) { return 4.0 * FALL.UM * y * (FALL.LY - y) / (FALL.LY * FALL.LY); }

function koerperMasken(g) {
  const nx = g.nx, ny = g.ny, r2 = (FALL.DZ / 2) * (FALL.DZ / 2);
  const mu = new Uint8Array((nx + 1) * ny), mv = new Uint8Array(nx * (ny + 1)),
        mp = new Uint8Array(nx * ny);
  const innen = (x, y) => (x - FALL.XC) * (x - FALL.XC) + (y - FALL.YC) * (y - FALL.YC) < r2;
  for (let i = 0; i <= nx; i++) for (let j = 0; j < ny; j++)
    if (innen(g.xu[i], g.yu[j])) mu[i * ny + j] = 1;
  for (let i = 0; i < nx; i++) for (let j = 0; j <= ny; j++)
    if (innen(g.xv[i], g.yv[j])) mv[i * (ny + 1) + j] = 1;
  for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++)
    if (innen(g.xp[i], g.yp[j])) mp[i * ny + j] = 1;
  return { mu: mu, mv: mv, mp: mp };
}

/* Randbedingungen -- sie gehoeren zur AUFGABE, nicht zum Weg, und sind deshalb fuer
   beide Wege dieselben. v sitzt bei y = 0 und y = Ly genau auf der Wand (dort null);
   u sitzt bei y = dy/2, also INNEN -- die Haftbedingung steckt im Spiegelwert. */
function raender(g, us, vs, ein) {
  const nx = g.nx, ny = g.ny, m = ny + 1;
  for (let i = 0; i < nx; i++) { vs[i * m] = 0; vs[i * m + ny] = 0; }
  for (let j = 0; j < m; j++) vs[j] = 0;                           // Einstroemung ohne Quer
  for (let j = 0; j < ny; j++) { us[j] = ein[j]; us[nx * ny + j] = us[(nx - 1) * ny + j]; }
  for (let j = 0; j < m; j++) vs[(nx - 1) * m + j] = vs[(nx - 2) * m + j];
}

function Loeser(nx, ny, cfl) {
  const g = new Gitter(FALL.LX, FALL.LY, nx, ny);
  this.g = g;
  const mk = koerperMasken(g);
  this.mu = mk.mu; this.mv = mk.mv; this.mp = mk.mp;
  this.ein = new Float64Array(ny);
  for (let j = 0; j < ny; j++) this.ein[j] = einstroemprofil(g.yu[j]);
  const dn = drucknetzMatrix(g, true, this.mu, this.mv);
  this.Ap = dn.A; this.leer = dn.leer; this.diag = ic0(dn.A);
  this.dt = Math.min((cfl || 0.4) * Math.min(g.dx, g.dy) / (2 * FALL.UM),
                     0.2 * Math.min(g.dx, g.dy) ** 2 / FALL.NU);
  this.ImpU = new ImpulsA(g, 'u', this.ein);
  this.ImpV = new ImpulsA(g, 'v', null);
  this.nenner = 0.5 * FALL.RHO * FALL.U_QUER * FALL.U_QUER * FALL.DZ;
  this.zustand = (weg) => {
    const u = g.leerU(), v = g.leerV();
    for (let i = 0; i <= nx; i++) for (let j = 0; j < ny; j++) u[i * ny + j] = this.ein[j];
    for (let k = 0; k < u.length; k++) if (this.mu[k]) u[k] = 0;
    for (let k = 0; k < v.length; k++) if (this.mv[k]) v[k] = 0;
    return { u: u, v: v, p: g.leerP(), t: 0, k: 0, weg: weg, cg: 0 };
  };
}

/* Ein vollstaendiger Zeitschritt. weg === 'A': Impuls als MATRIX, Drucknetz als Matrix.
   weg === 'B': Impuls als Feldoperation, Drucknetz matrixfrei. Geloest wird in beiden
   Faellen mit demselben Verfahren (CG) -- im Browser gibt es keinen direkten Loeser.
   Was sich unterscheidet, ist also der AUFBAU des Operators, nicht der Loeser. */
Loeser.prototype.schritt = function (Z) {
  const g = this.g, nx = g.nx, ny = g.ny, m = ny + 1, dt = this.dt, A = g.dx * g.dy;
  let rU, rV;
  if (Z.weg === 'A') {
    rU = this.ImpU.fuellen(Z.u, Z.v, FALL.NU).wirkung(Z.u);
    rV = this.ImpV.fuellen(Z.u, Z.v, FALL.NU).wirkung(Z.v);
  } else {
    rU = konvZaehU(g, Z.u, Z.v, FALL.NU, this.ein);
    rV = konvZaehV(g, Z.u, Z.v, FALL.NU);
  }
  const us = new Float64Array(Z.u.length), vs = new Float64Array(Z.v.length);
  for (let k = 0; k < us.length; k++) us[k] = Z.u[k] + dt * rU[k];
  for (let k = 0; k < vs.length; k++) vs[k] = Z.v[k] + dt * rV[k];
  raender(g, us, vs, this.ein);
  let fx1 = 0, fy1 = 0;
  for (let k = 0; k < us.length; k++) if (this.mu[k]) { fx1 += us[k]; us[k] = 0; }
  for (let k = 0; k < vs.length; k++) if (this.mv[k]) { fy1 += vs[k]; vs[k] = 0; }
  const d = g.divergenz(us, vs), b = new Float64Array(d.length);
  for (let k = 0; k < b.length; k++) b[k] = this.leer[k] ? 0 : -d[k] * A / dt;
  const mv = (Z.weg === 'A')
    ? (x, y) => this.Ap.mal(x, y)
    : (x, y) => drucknetzWirkungVoll(g, x, this.leer, true, this.mu, this.mv, y);
  const erg = cg(mv, b, Z.p, this.diag, 1e-12, 600);
  Z.cg = erg.schritte;
  const p = Z.p;
  for (let i = 1; i < nx; i++) for (let j = 0; j < ny; j++)
    us[i * ny + j] -= dt * (p[i * ny + j] - p[(i - 1) * ny + j]) / g.dx;
  for (let i = 0; i < nx; i++) for (let j = 1; j < ny; j++)
    vs[i * m + j] -= dt * (p[i * ny + j] - p[i * ny + j - 1]) / g.dy;
  for (let j = 0; j < ny; j++) us[nx * ny + j] -= dt * (0 - p[(nx - 1) * ny + j]) / (g.dx / 2);
  let fx2 = 0, fy2 = 0;
  for (let k = 0; k < us.length; k++) if (this.mu[k]) { fx2 += us[k]; us[k] = 0; }
  for (let k = 0; k < vs.length; k++) if (this.mv[k]) { fy2 += vs[k]; vs[k] = 0; }
  for (let j = 0; j < ny; j++) us[j] = this.ein[j];
  Z.u = us; Z.v = vs; Z.t += dt; Z.k += 1;
  const f = A * FALL.RHO / dt / this.nenner;
  Z.cD = (fx1 + fx2) * f; Z.cL = (fy1 + fy2) * f;
  return Z;
};

/* omega = dv/dx - du/dy auf den INNEREN Gitterecken -- dort haben beide Anteile
   denselben Ort (Befund B10). Form (nx-1, ny-1). */
function wirbelstaerke(g, u, v) {
  const nx = g.nx, ny = g.ny, m = ny + 1, o = new Float64Array((nx - 1) * (ny - 1));
  for (let i = 1; i < nx; i++)
    for (let j = 1; j < ny; j++)
      o[(i - 1) * (ny - 1) + (j - 1)] =
        (v[i * m + j] - v[(i - 1) * m + j]) / g.dx
      - (u[i * ny + j] - u[i * ny + j - 1]) / g.dy;
  return o;
}

global.KernCFD = {
  FALL: FALL, einstroemprofil: einstroemprofil, koerperMasken: koerperMasken,
  Loeser: Loeser, wirbelstaerke: wirbelstaerke, raender: raender, ic0: ic0,
  CSR: CSR, drucknetzMatrix: drucknetzMatrix, drucknetzWirkungVoll: drucknetzWirkungVoll,
  stempelImpulsU: stempelImpulsU, stempelImpulsV: stempelImpulsV, ImpulsA: ImpulsA, cg: cg,
  bernoulli: bernoulli, sgLeitwerte: sgLeitwerte,
  loese1d: loese1d, genau1d: genau1d, thomas: thomas,
  Gitter: Gitter, randU: randU, randV: randV,
  konvZaehU: konvZaehU, konvZaehV: konvZaehV
};
})(typeof window !== 'undefined' ? window : globalThis);
