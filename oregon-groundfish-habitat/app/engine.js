// Scenario engine: habitat (fitted delta-GAM) x per-recruit fishing layer. Pure functions, no DOM.
(function (root) {
  const FM = 1.8288; // metres per fathom

  // ---------- life history ----------
  function lengthAt(g, a) {
    let L;
    if (g.a2 === null || g.a2 === undefined) L = g.Linf + (g.L1 - g.Linf) * Math.exp(-g.k * (a - g.a1));
    else L = g.L1 + (g.L2 - g.L1) * (1 - Math.exp(-g.k * (a - g.a1))) / (1 - Math.exp(-g.k * (g.a2 - g.a1)));
    return Math.max(L, 1);
  }
  const logistic = (x) => 1 / (1 + Math.exp(-x));

  function schedule(sp) {
    const lh = sp.lh, n = lh.amax + 1;
    const L = new Float64Array(n), Lmid = new Float64Array(n), W = new Float64Array(n), mat = new Float64Array(n);
    for (let a = 0; a < n; a++) {
      L[a] = lengthAt(lh.growth, a);
      Lmid[a] = lengthAt(lh.growth, a + 0.5);
      W[a] = lh.wl.a * Math.pow(L[a], lh.wl.b);
      const x = lh.mat.type === "age" ? a : L[a];
      mat[a] = 1 / (1 + Math.exp(lh.mat.slope * (x - lh.mat.x50)));
    }
    return { n, L, Lmid, W, mat, M: lh.M };
  }

  // spawning output per recruit given age-specific fishing mortality
  function sbpr(sch, Fa) {
    let N = 1, s = 0;
    for (let a = 0; a < sch.n; a++) {
      const Z = sch.M + (Fa ? Fa[a] : 0);
      if (a === sch.n - 1) N = N / (1 - Math.exp(-Z)); // plus group
      s += N * sch.mat[a] * sch.W[a];
      N *= Math.exp(-Z);
    }
    return s;
  }

  // ---------- bag-limit arithmetic ----------
  // legal fish encountered per angler-trip ~ NegBin(mean mu, size theta); R(b) = E[min(X, b)]
  function expectedRetained(mu, theta, b) {
    if (b <= 0) return 0;
    let p = Math.pow(theta / (theta + mu), theta), cdf = 0, r = 0;
    const q = mu / (theta + mu);
    for (let x = 0; x < b; x++) { cdf += p; r += 1 - cdf; p *= ((x + theta) / (x + 1)) * q; }
    return r;
  }
  // mortality on legal fish per unit of encounter mortality: retained + released-over-bag x dmort
  function bagFactor(rec, b) {
    const R = expectedRetained(rec.mu, rec.theta, b);
    return (R + rec.dmort * (rec.mu - R)) / rec.mu;
  }

  // age-specific F for one cell class
  function fAtAge(sch, sp, Frec, Fcom, reg) {
    const n = sch.n, out = new Float64Array(n), r = sp.rec, c = sp.com;
    const g = bagFactor(r, reg.bag);
    for (let a = 0; a < n; a++) {
      const L = sch.Lmid[a];
      const sr = logistic((L - r.sel50) / r.selw), sc = logistic((L - c.sel50) / c.selw);
      const ret = reg.sizeCm > 0 ? logistic((L - reg.sizeCm) / 0.75) : 1;
      out[a] = Frec * sr * (ret * g + (1 - ret) * r.dmort) + Fcom * sc;
    }
    return out;
  }

  // Calibrate total apical F so open-water SPR under 2026 rules matches the assessment.
  function calibrate(sp, sch) {
    const sq = { bag: sp.regs.bag, sizeCm: sp.regs.size_cm };
    const g0 = bagFactor(sp.rec, sq.bag), s = sp.rec.share, s0 = sbpr(sch, null);
    const sprAt = (th) => sbpr(sch, fAtAge(sch, sp, th * s / g0, th * (1 - s), sq)) / s0;
    let lo = 0, hi = 5;
    for (let i = 0; i < 80; i++) { const mid = (lo + hi) / 2; if (sprAt(mid) > sp.status.spr) lo = mid; else hi = mid; }
    const th = (lo + hi) / 2;
    return { theta: th, FrecEnc0: th * s / g0, Fcom0: th * (1 - s), sbpr0: s0, sprCheck: sprAt(th) };
  }

  // ---------- habitat ----------
  function tIndex(T, t0, dt, nT) {
    const i = Math.round((T - t0) / dt);
    return i < 0 ? 0 : i >= nT ? nT - 1 : i;
  }

  // density per cell for a temperature anomaly, using either point-estimate or one posterior draw of the T smooth
  function density(m, grid, dT, fp, fg, tg, out) {
    const n = grid.n;
    for (let i = 0; i < n; i++) {
      const k = tIndex(grid.tclim[i] + dT, tg.t0, tg.dt, tg.n);
      out[i] = logistic(m.eta_p[i] + fp[k]) * Math.exp(m.eta_g[i] + fg[k]);
    }
    return out;
  }

  root.HabitatEngine = { FM, lengthAt, schedule, sbpr, expectedRetained, bagFactor, fAtAge, calibrate, density, tIndex };
  if (typeof module !== "undefined") module.exports = root.HabitatEngine;
})(typeof window !== "undefined" ? window : globalThis);
