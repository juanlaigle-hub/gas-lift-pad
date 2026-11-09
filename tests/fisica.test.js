// Pruebas de la física y de los algoritmos de la demo. Se corren con:  node --test tests/fisica.test.js
const test = require('node:test');
const assert = require('node:assert/strict');
const GL = require('../modelo.js');
const W = GL.POZOS;

test('IPR de Vogel: caudal máximo con pwf = 0 y nulo con pwf = pr', () => {
  const w = W[2];
  assert.ok(Math.abs(GL.ipr(w, 0) - w.qmax) < 1e-9);
  assert.ok(Math.abs(GL.ipr(w, w.pr)) < 1e-9);
});

test('la presión aumenta con la profundidad y una columna con gas pesa menos', () => {
  const w = W[3], sinGas = GL.columna(w, 50, 0).pTalon, conGas = GL.columna(w, 50, 30).pTalon;
  assert.ok(sinGas > GL.C.pwh && conGas > GL.C.pwh);
  assert.ok(conGas < sinGas, `con gas ${conGas} < sin gas ${sinGas}`);
});

test('curva de gas lift: los pozos depletados no surgen, el gas los hace producir y el exceso de gas resta', () => {
  for (const w of W.slice(2)) {
    const c = GL.curva(w, 0, 80, 1), max = c.reduce((a, b) => (b.qo > a.qo ? b : a));
    assert.equal(c[0].qo, 0, `${w.id} no surge sin gas`);
    assert.ok(max.qgi > 5 && max.qgi < 80, `${w.id}: máximo interior (${max.qgi})`);
    assert.ok(c[80].qo < max.qo, `${w.id}: con exceso de gas la producción baja`);
  }
});

test('el punto de operación cumple a la vez la IPR y la columna', () => {
  const w = W[1], o = GL.operar(w, 30, 0);
  assert.ok(Math.abs(GL.ipr(w, o.pwf) - o.ql) < 1e-6);
  assert.ok(Math.abs(GL.columna(w, o.ql, 30).pTalon - o.pwf) / o.pwf < 1e-4);
});

test('presión de casing en el rango de inyección de gas lift (5 a 12 MPa)', () => {
  for (const w of W.slice(1)) { const pc = GL.operar(w, 35, 0).pc / 1e6; assert.ok(pc > 5 && pc < 12, `${w.id}: ${pc}`); }
});

test('la declinación reduce la producción con el tiempo', () => {
  const w = W[2]; assert.ok(GL.operar(w, 35, 365).qo < GL.operar(w, 35, 0).qo);
});

test('reparto óptimo: respeta la capacidad, nunca es peor que el uniforme y ahorra gas', () => {
  const cs = W.map(w => GL.curva(w, 0, 80, 1)), tot = g => g.reduce((a, q, i) => a + cs[i][q].qo, 0);
  for (const cap of [60, 120, 180]) {
    const go = GL.optimizar(cs, cap, 1), gu = GL.uniforme(6, cap).map(Math.round);
    assert.ok(go.reduce((a, b) => a + b, 0) <= cap);
    assert.ok(tot(go) >= tot(gu) - 1e-9);
    let g = 0; while (tot(GL.optimizar(cs, g, 1)) < tot(gu)) g++;
    assert.ok(g < cap, `cap ${cap}: el óptimo iguala al uniforme con ${g}`);
  }
});

test('reparto óptimo: con dos pozos coincide con la búsqueda exhaustiva', () => {
  const cs = [GL.curva(W[2], 0, 80, 1), GL.curva(W[4], 0, 80, 1)], cap = 50;
  let best = -1;
  for (let a = 0; a <= cap; a++) for (let b = 0; a + b <= cap; b++) best = Math.max(best, cs[0][a].qo + cs[1][b].qo);
  const go = GL.optimizar(cs, cap, 1);
  assert.ok(Math.abs(cs[0][go[0]].qo + cs[1][go[1]].qo - best) < 1e-9);
});

test('ensayo con declinación: el modelo recupera el óptimo y la declinación', () => {
  for (const i of [2, 3, 4]) {
    const w = W[i], en = GL.ensayo(w), fit = GL.ajusteModelo(w, en.hist.concat(en.test));
    const am = c => c.reduce((a, b) => (b.qo > a.qo ? b : a)).qgi;
    const real = am(GL.curvaK(w, GL.factorDecl(en.tFin)));
    const mod = am(GL.curvaK(w, fit.f * GL.factorDecl(en.tFin, { Di: fit.Di, b: GL.DECL.b })));
    assert.ok(Math.abs(mod - real) <= 2, `${w.id}: modelo ${mod} vs real ${real}`);
    assert.ok(Math.abs(fit.Di - GL.DECL.Di) < 0.0015, `${w.id}: Di ${fit.Di}`);
  }
});

test('diagnóstico: reconoce las anomalías de severidad media en al menos el 95 % de los casos', () => {
  const r = GL.rng(5); let ok = 0, n = 0;
  for (const w of W) for (const k of ['parafina', 'valvula', 'frachit']) for (let j = 0; j < 8; j++) {
    const d = GL.diagnosticar(w, 35, 0, GL.medir(w, 35, 0, k, r, 0.6));
    const top = Object.keys(d.P).reduce((a, b) => (d.P[b] > d.P[a] ? b : a)); n++; if (top === k) ok++;
  }
  assert.ok(ok / n >= 0.95, `${ok}/${n}`);
});

test('diagnóstico: huellas físicas con el signo correcto', () => {
  const w = W[2], d = GL.diagnosticar(w, 35, 0, GL.medir(w, 35, 0, 'normal', GL.rng(1)));
  assert.ok(d.H.parafina[0] < 0 && d.H.parafina[1] > 0, 'parafina: menos líquido, más presión de casing');
  assert.ok(d.H.valvula[0] < 0 && d.H.valvula[1] < 0, 'válvula superior: menos líquido, menos presión de casing');
  assert.ok(d.H.frachit[2] > 0.2, 'frac hit: mucha más agua');
});
