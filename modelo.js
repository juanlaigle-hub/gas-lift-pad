/* Núcleo físico de la demo de gas lift en un pad.
 * Modelo estacionario de pozo: IPR de Vogel + columna de producción con modelo de deriva (Zuber–Findlay,
 * C0 = 1,2, velocidad de deriva de Harmathy) para la fracción de líquido, fricción de Darcy (Haaland),
 * gas ideal con z constante y perfil lineal de temperatura.
 * Unidades internas SI (Pa, m, s, kg); interfaz en m³/d de líquido y miles de Sm³/d (kSm³/d) de gas.
 * Simplificaciones documentadas en el README: un solo régimen de flujo (deriva), todo el gas de formación libre,
 * z constante, pozo vertical equivalente hasta el talón.
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.GL = factory();
})(this, function () {
  'use strict';
  const g = 9.81, R = 8.314, PSTD = 101325, TSTD = 288.15;
  const C = {
    M: 0.0217,            // kg/mol, gas de densidad relativa 0,75
    z: 0.9,               // factor de compresibilidad medio (constante)
    rhoO: 820,            // kg/m³, petróleo liviano (~41 °API)
    rhoW: 1050,           // kg/m³, agua de formación
    muL: 1.0e-3,          // Pa·s, viscosidad del líquido a temperatura de fondo
    muG: 1.5e-5,          // Pa·s
    sigma: 0.025,         // N/m, tensión superficial gas-líquido
    C0: 1.2,              // coeficiente de distribución (Zuber–Findlay)
    D: 0.062,             // m, tubing de 2 7/8" (diámetro interno)
    rug: 4.5e-5,          // m, rugosidad del acero
    Tsup: 293.15, Tfondo: 368.15, // K, 20 °C en superficie y 95 °C en el talón
    Ltalon: 2800,         // m, profundidad vertical del talón
    pwh: 1.0e6,           // Pa, presión de cabeza (línea de conducción)
    dpValv: 0.5e6,        // Pa, caída de presión en la válvula operadora
    paso: 50              // m, tramo de integración
  };
  const rhoGstd = PSTD * C.M / (R * TSTD);     // kg/Sm³
  const temp = z => C.Tsup + (C.Tfondo - C.Tsup) * z / C.Ltalon;

  // Pozos del pad (datos simulados). pr: presión media de drenaje; qmax: caudal de líquido a pwf = 0 (Vogel);
  // wct: corte de agua; glr: gas de formación por m³ de líquido; Linj: profundidad de la válvula operadora.
  const POZOS = [
    { id: 'P-1', pr: 24e6, qmax: 300, wct: 0.08, glr: 180, Linj: 2400, etapa: 'recién salido de surgencia' },
    { id: 'P-2', pr: 21e6, qmax: 230, wct: 0.15, glr: 140, Linj: 2400, etapa: 'primer año de gas lift' },
    { id: 'P-3', pr: 19e6, qmax: 190, wct: 0.25, glr: 100, Linj: 2400, etapa: 'primer año de gas lift' },
    { id: 'P-4', pr: 17e6, qmax: 150, wct: 0.40, glr: 70,  Linj: 2400, etapa: 'maduro, más agua' },
    { id: 'P-5', pr: 15.5e6, qmax: 120, wct: 0.50, glr: 55, Linj: 2400, etapa: 'maduro, más agua' },
    { id: 'P-6', pr: 20e6, qmax: 110, wct: 0.20, glr: 120, Linj: 2400, etapa: 'baja productividad' }
  ];
  // Declinación hiperbólica de la productividad (Arps): qmax(t) = qmax0 / (1 + b·Di·t)^(1/b)
  const DECL = { Di: 0.004, b: 1.0 };   // 1/d
  const factorDecl = (t, d = DECL) => Math.pow(1 + d.b * d.Di * t, -1 / d.b);

  // IPR de Vogel: caudal de líquido [m³/d] para una pwf [Pa]
  function ipr(w, pwf) {
    const x = Math.min(1, Math.max(0, pwf / w.pr));
    return w.qmax * (1 - 0.2 * x - 0.8 * x * x);
  }

  // Gradiente de un tramo con flujo homogéneo. ql: m³/d de líquido; qgStd: Sm³/d de gas; D: diámetro.
  function gradiente(p, z, ql, wct, qgStd, D) {
    const T = temp(z);
    const qlS = ql / 86400, qgS = qgStd / 86400;
    const rhoL = wct * C.rhoW + (1 - wct) * C.rhoO;
    const mL = qlS * rhoL, mG = qgS * rhoGstd;
    const qgIn = qgS * (PSTD / p) * (T / TSTD) * C.z;        // m³/s de gas en condiciones del tramo
    const A = Math.PI * D * D / 4, Qt = qlS + qgIn;
    if (Qt <= 0) return rhoL * g;
    const rhoG = p * C.M / (C.z * R * T);
    const vsl = qlS / A, vsg = qgIn / A, vm = vsl + vsg;
    // Modelo de deriva: el gas sube más rápido que la mezcla (C0·vm + vd); vd de Harmathy para burbujas
    const vd = 1.53 * Math.pow(g * C.sigma * (rhoL - rhoG) / (rhoL * rhoL), 0.25);
    const alfa = Math.min(0.98, vsg / (C.C0 * vm + vd));
    const rhoS = alfa * rhoG + (1 - alfa) * rhoL;            // densidad con deslizamiento (gravedad)
    const rhoN = (mL + mG) / Qt;                            // densidad sin deslizamiento (fricción)
    const lam = qlS / Qt;
    const mu = lam * C.muL + (1 - lam) * C.muG;
    const Re = Math.max(1, rhoN * vm * D / mu);
    const f = Re < 2300 ? 64 / Re : Math.pow(-1.8 * Math.log10(Math.pow(C.rug / D / 3.7, 1.11) + 6.9 / Re), -2);
    return rhoS * g + f * rhoN * vm * vm / (2 * D);          // Pa/m (la presión aumenta hacia abajo)
  }

  // Marcha de presión desde la cabeza hasta el talón. Devuelve presión en la válvula y en el talón.
  // anom: { parafina: m de tubing con diámetro reducido desde la superficie, Dpar, LinjReal }
  function columna(w, ql, qgi, anom = {}) {
    const Linj = anom.LinjReal || w.Linj;
    const wct = anom.wct != null ? anom.wct : w.wct;
    let p = C.pwh, z = 0, pInj = null;
    const qgForm = w.glr * ql;
    while (z < C.Ltalon) {
      const dz = Math.min(C.paso, C.Ltalon - z);
      const qg = qgForm + (z < Linj ? qgi * 1000 : 0);
      const D = anom.parafina && z < anom.parafina ? anom.Dpar : C.D;
      // Runge–Kutta de 2º orden
      const k1 = gradiente(p, z, ql, wct, qg, D);
      const k2 = gradiente(p + k1 * dz, z + dz, ql, wct, qg, D);
      const pn = p + 0.5 * (k1 + k2) * dz;
      if (pInj === null && z + dz >= Linj) pInj = p + (pn - p) * (Linj - z) / dz;
      p = pn; z += dz;
    }
    return { pTalon: p, pInj };
  }

  // Punto de operación: cruce de IPR y columna. qgi en kSm³/d. Devuelve caudales y presiones.
  function operar(w, qgi, t = 0, anom = {}) {
    const wt = Object.assign({}, w, { qmax: w.qmax * factorDecl(t, anom.decl) * (anom.prodFac || 1), pr: w.pr * (anom.prFac || 1) });
    const wct = anom.wct != null ? anom.wct : w.wct;
    const res = pwf => columna(wt, ipr(wt, pwf), qgi, anom).pTalon - pwf;
    let lo = C.pwh, hi = wt.pr * 0.999;
    if (res(hi) > 0) return sinFlujo(wt, qgi, wct, anom); // la columna pesa más que el reservorio: no fluye
    if (res(lo) < 0) lo = C.pwh;                          // (no ocurre en el rango de la demo)
    for (let i = 0; i < 50; i++) { const m = 0.5 * (lo + hi); if (res(m) > 0) lo = m; else hi = m; }
    const pwf = 0.5 * (lo + hi), ql = ipr(wt, pwf);
    const col = columna(wt, ql, qgi, anom);
    return { ql, qo: ql * (1 - wct), qw: ql * wct, wct, pwf, pInj: col.pInj, pc: presionCasing(col.pInj, anom.LinjReal || w.Linj), qgi };
  }
  // Pozo que no fluye: por el tubing sube solo el gas inyectado; la presión de casing es la de esa columna de gas.
  function sinFlujo(w, qgi, wct, anom) {
    const col = columna(w, 0, Math.max(qgi, 1e-3), anom);
    return { ql: 0, qo: 0, qw: 0, wct, pwf: w.pr, pInj: col.pInj, pc: presionCasing(col.pInj, anom.LinjReal || w.Linj), qgi };
  }

  // Presión de inyección en superficie: la del tubing en la válvula más la caída en la válvula, menos el peso del gas en el anular.
  function presionCasing(pInj, Linj) {
    if (!isFinite(pInj)) return NaN;
    const Tm = 0.5 * (C.Tsup + temp(Linj));
    return (pInj + C.dpValv) * Math.exp(-C.M * g * Linj / (C.z * R * Tm));
  }

  // Curva de rendimiento: producción de petróleo vs gas inyectado
  function curva(w, t = 0, qgMax = 80, paso = 1, anom = {}) {
    const out = [];
    for (let q = 0; q <= qgMax + 1e-9; q += paso) out.push({ qgi: q, qo: operar(w, q, t, anom).qo });
    return out;
  }

  // Reparto óptimo del gas por programación dinámica (exacta para la grilla, admite curvas no cóncavas).
  // curvas[i][k].qo con qgi = k·paso. Devuelve índices k por pozo.
  function optimizar(curvas, capacidad, paso = 1) {
    const n = curvas.length, K = Math.floor(capacidad / paso + 1e-9);
    let best = new Float64Array(K + 1); const elec = [];
    for (let i = 0; i < n; i++) {
      const nb = new Float64Array(K + 1).fill(-1e18), ch = new Int32Array(K + 1);
      for (let c = 0; c <= K; c++) for (let k = 0; k < curvas[i].length && k <= c; k++) {
        const v = best[c - k] + curvas[i][k].qo; if (v > nb[c]) { nb[c] = v; ch[c] = k; }
      }
      elec.push(ch); best = nb;
    }
    let c = 0; for (let j = 0; j <= K; j++) if (best[j] > best[c] + 1e-9) c = j;   // no hace falta usar todo el gas
    const ks = new Array(n);
    for (let i = n - 1; i >= 0; i--) { ks[i] = elec[i][c]; c -= ks[i]; }
    return ks.map(k => k * paso);
  }
  function uniforme(n, capacidad) { return new Array(n).fill(capacidad / n); }
  function producir(pozos, gases, t = 0) { return pozos.map((w, i) => operar(w, gases[i], t)); }

  // Generador de números aleatorios reproducible (para los datos simulados del ensayo)
  function rng(seed) { let s = seed >>> 0; return () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; }; }
  function normal(r) { let u = 0, v = 0; while (!u) u = r(); v = r(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); }

  // Ensayo de escalones con declinación. Historia previa a inyección constante + escalones crecientes.
  function ensayo(w, { tInicio = 30, qgBase = 25, escalones = [10, 25, 40, 55, 70], dias = 4, ruido = 0.015, seed = 7 } = {}) {
    const r = rng(seed), hist = [], test = [];
    for (let d = 0; d < tInicio; d++) hist.push({ t: d, qgi: qgBase, qo: operar(w, qgBase, d).qo * (1 + ruido * normal(r)) });
    let t = tInicio;
    for (const q of escalones) { t += dias; test.push({ t, qgi: q, qo: operar(w, q, t).qo * (1 + ruido * normal(r)) }); }
    return { hist, test, tFin: t };
  }
  // Método habitual: parábola qo = a + b·qgi − c·qgi² ajustada a los puntos del ensayo, sin considerar el tiempo.
  function ajusteParabola(test) {
    // mínimos cuadrados con ecuaciones normales 3×3
    const S = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], y = [0, 0, 0];
    for (const p of test) { const v = [1, p.qgi, p.qgi * p.qgi]; for (let i = 0; i < 3; i++) { y[i] += v[i] * p.qo; for (let j = 0; j < 3; j++) S[i][j] += v[i] * v[j]; } }
    const det = m => m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    const D0 = det(S), c = [0, 1, 2].map(k => det(S.map((row, i) => row.map((v, j) => (j === k ? y[i] : v)))) / D0);
    const qopt = c[2] < 0 ? -c[1] / (2 * c[2]) : NaN;
    return { coef: c, qopt, f: q => c[0] + c[1] * q + c[2] * q * q };
  }
  // Método propuesto: modelo físico con dos parámetros (productividad inicial f y declinación Di).
  // Se tabula qo(k) para cada caudal de inyección usado, con k = f·declinación(t), y se ajusta por grilla.
  function tablaK(w, qgis) {
    const K = []; for (let k = 0.3; k <= 1.5001; k += 0.02) K.push(k);
    const T = {}; for (const q of qgis) T[q] = K.map(k => operar(w, q, 0, { prodFac: k }).qo);
    const interp = (q, k) => { const t = T[q], x = (k - 0.3) / 0.02, i = Math.max(0, Math.min(K.length - 2, Math.floor(x))), u = x - i; return t[i] * (1 - u) + t[i + 1] * u; };
    return interp;
  }
  function ajusteModelo(w, datos) {
    const qgis = [...new Set(datos.map(p => p.qgi))], qo = tablaK(w, qgis); let mejor = null;
    for (let Di = 0; Di <= 0.0120001; Di += 0.0002) for (let f = 0.6; f <= 1.4001; f += 0.005) {
      let e = 0; for (const p of datos) { const m = qo(p.qgi, f * factorDecl(p.t, { Di, b: DECL.b })); e += (m - p.qo) ** 2; }
      if (!mejor || e < mejor.e) mejor = { e, f, Di };
    }
    return mejor;
  }
  // Óptimo económico de un pozo: inyectar mientras el petróleo extra pague el gas (relación marginal mínima).
  function curvaK(w, k, qgMax = 80, paso = 1) { const out = []; for (let q = 0; q <= qgMax + 1e-9; q += paso) out.push({ qgi: q, qo: operar(w, q, 0, { prodFac: k }).qo }); return out; }

  // Diagnóstico por residuos: compara lo medido con lo que predice el modelo calibrado.
  const ANOMALIAS = {
    normal:   { nombre: 'Operación normal', anom: {} },
    parafina: { nombre: 'Parafina en el tubing', anom: { parafina: 900, Dpar: 0.040 } },
    valvula:  { nombre: 'Inyección por una válvula superior', anom: { LinjReal: 1500 } },
    frachit:  { nombre: 'Frac hit con entrada de agua', anom: { wctMas: 0.35, prFac: 1.03 } }
  };
  // sev ∈ (0, 1]: severidad de la anomalía (parafina: estrechamiento; frac hit: agua que entra). La válvula es todo o nada.
  function anomalia(w, clave, sev = 1) {
    const a = Object.assign({}, ANOMALIAS[clave].anom);
    if (a.Dpar) a.Dpar = C.D - sev * (C.D - a.Dpar);
    if (a.wctMas) { a.wct = Math.min(0.95, w.wct + sev * a.wctMas); a.prFac = 1 + sev * (a.prFac - 1); delete a.wctMas; }
    return a;
  }
  function medir(w, qgi, t, clave, r, sev = 1) {
    const a = anomalia(w, clave, sev);
    const o = operar(w, qgi, t, a), n = () => normal(r);
    return { ql: o.ql * (1 + 0.02 * n()), pc: o.pc + 0.05e6 * n(), wct: Math.min(0.99, Math.max(0, o.wct + 0.01 * n())), qo: o.qo };
  }
  // Huellas: cada anomalía produce un patrón de residuos (Δql relativo, Δpc, Δwct) que se calcula con el mismo modelo.
  function huellas(w, qgi, t) {
    const pred = operar(w, qgi, t), H = {};
    for (const k of Object.keys(ANOMALIAS)) {
      const o = operar(w, qgi, t, anomalia(w, k, 1));
      H[k] = residuos(pred, { ql: o.ql, pc: o.pc, wct: o.wct });
    }
    return { pred, H };
  }
  function residuos(pred, m) { return [(m.ql - pred.ql) / Math.max(1, pred.ql), (m.pc - pred.pc) / 1e6, m.wct - pred.wct]; }
  const ESCALA = [0.02, 0.05, 0.01];   // desvío del ruido de medición en cada residuo
  function diagnosticar(w, qgi, t, m) {
    const { pred, H } = huellas(w, qgi, t), r = residuos(pred, m), ll = {}, sev = {};
    for (const k of Object.keys(H)) {
      if (k === 'normal') { ll[k] = -0.5 * r.reduce((s, v, i) => s + (v / ESCALA[i]) ** 2, 0); sev[k] = 0; continue; }
      // severidad que mejor explica los residuos (aproximación lineal en la severidad), acotada a [0, 1]
      let num = 0, den = 0; H[k].forEach((h, i) => { num += r[i] * h / ESCALA[i] ** 2; den += h * h / ESCALA[i] ** 2; });
      const s = Math.min(1, Math.max(0, den > 0 ? num / den : 0)); sev[k] = s;
      ll[k] = -0.5 * r.reduce((acc, v, i) => acc + ((v - s * H[k][i]) / ESCALA[i]) ** 2, 0) - 1.5; // penalización por el parámetro extra
    }
    const mx = Math.max(...Object.values(ll)); let Z = 0; const P = {};
    for (const k of Object.keys(ll)) { P[k] = Math.exp(ll[k] - mx); Z += P[k]; }
    for (const k of Object.keys(P)) P[k] /= Z;
    return { pred, r, P, H, sev };
  }

  return { C, POZOS, DECL, factorDecl, ipr, columna, operar, curva, optimizar, uniforme, producir, ensayo, ajusteParabola, ajusteModelo, curvaK, rng, normal, ANOMALIAS, anomalia, medir, diagnosticar, presionCasing };
});
