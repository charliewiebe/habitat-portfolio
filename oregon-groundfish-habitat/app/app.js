(function () {
  "use strict";
  const E = window.HabitatEngine, D = DATA, G = D.grid, N = G.n, FM = E.FM;
  const SPECIES = Object.keys(D.params);
  const $ = (s, el = document) => el.querySelector(s);
  const MINUS = "−";
  const fmtPct = (v, d = 0) => (v > 0 ? "+" : v < 0 ? MINUS : "") + Math.abs(v * 100).toFixed(d) + "%";
  const fmtPctAbs = (v, d = 0) => (v * 100).toFixed(d) + "%";
  const fmtDens = (v) => v >= 100 ? v.toFixed(0) : v >= 10 ? v.toFixed(1) : v >= 1 ? v.toFixed(2) : v >= 0.01 ? v.toFixed(3) : v.toExponential(1);
  const fmtLat = (v) => Math.abs(v).toFixed(3) + "°" + (v >= 0 ? "N" : "S");
  const fmtLon = (v) => Math.abs(v).toFixed(3) + "°" + (v < 0 ? "W" : "E");
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

  // ---------------- static per-cell arrays ----------------
  const depth = Float64Array.from(G.depth), area = Float64Array.from(G.area), tclim = Float64Array.from(G.tclim);
  const availRec = depth.map((d) => 1 / (1 + Math.exp((d - 150) / 25))); // recreational reach fades past ~150 m
  const CLOSURES = D.geo.closures;
  const DEPTH_BINS = [[0, 100], [100, 150], [150, 200], [200, 300], [300, 500], [500, 800], [800, 1300], [1300, 2000]];
  const binOf = depth.map((d) => DEPTH_BINS.findIndex(([a, b]) => d >= a && d < b));

  // ---------------- model caches ----------------
  const model = {};
  function buildModel(k, adv) {
    const base = D.params[k];
    const sp = JSON.parse(JSON.stringify(base));
    sp.rec.share = adv.share; sp.rec.mu = adv.mu; sp.rec.dmort = adv.dmort;
    const sch = E.schedule(sp), cal = E.calibrate(sp, sch);
    const m = D.species[k];
    const fpD = m.fp_draws, fgD = m.fg_draws, nd = fpD.length;
    if (!model[k] || !model[k].dens0D) {
      const dens0D = [];
      for (let d = 0; d < nd; d++) dens0D.push(E.density(m, G, 0, fpD[d], fgD[d], D.tgrid, new Float32Array(N)));
      model[k] = { dens0D };
    }
    Object.assign(model[k], { sp, sch, cal, m });
    // fixed colour domain from status-quo density
    const d0 = E.density(m, G, 0, m.fp_t, m.fg_t, D.tgrid, new Float64Array(N));
    const logs = Array.from(d0).filter((v) => v > 0).map(Math.log10).sort((a, b) => a - b);
    const q = (p) => logs[Math.min(logs.length - 1, Math.floor(p * logs.length))];
    model[k].hi = q(0.995); model[k].lo = Math.max(q(0.06), model[k].hi - 2.5);
    model[k].sq = null;
  }

  // ---------------- state ----------------
  const defaultClosureMode = (c) => c.kind === "reserve" ? "all" : c.id === "yrca_stonewall_bank" ? "rec" : "off";
  const advDefault = (k) => ({ share: D.params[k].rec.share, mu: D.params[k].rec.mu, dmort: D.params[k].rec.dmort });
  const state = {
    sp: "lingcod", layer: "habitat",
    overlays: { closures: true, contours: true, hauls: false },
    dT: 0, hypDepth: 0, hypSev: 0.5,
    regs: {}, closures: {}, band: { mode: "off", fm0: 30, fm1: 100 }, displace: true,
    win: [30, 1830], adv: {}, pinned: null, hover: null,
  };
  function regsDefault(k) { const r = D.params[k].regs; return { bag: r.bag, sizeIn: r.size_in, recEffort: 1, comMult: 1 }; }
  SPECIES.forEach((k) => { state.regs[k] = regsDefault(k); state.adv[k] = advDefault(k); });
  CLOSURES.forEach((c) => (state.closures[c.id] = defaultClosureMode(c)));

  function sqScenario(k) {
    const closures = {}; CLOSURES.forEach((c) => (closures[c.id] = defaultClosureMode(c)));
    return { k, dT: 0, hypDepth: 0, hypSev: 0, regs: regsDefault(k), closures, band: { mode: "off", fm0: 30, fm1: 100 }, displace: true };
  }
  function curScenario() {
    return { k: state.sp, dT: state.dT, hypDepth: state.hypDepth, hypSev: state.hypSev, regs: state.regs[state.sp],
      closures: state.closures, band: state.band, displace: state.displace };
  }

  // ---------------- scenario computation ----------------
  function closureMasks(sc) {
    const recC = new Uint8Array(N), comC = new Uint8Array(N);
    for (const c of CLOSURES) {
      const mode = sc.closures[c.id]; if (mode === "off") continue;
      for (const i of c.cells) { recC[i] = 1; if (mode === "all") comC[i] = 1; }
    }
    if (sc.band.mode !== "off") {
      for (let i = 0; i < N; i++) {
        const fm = depth[i] / FM;
        if (fm >= sc.band.fm0 && fm <= sc.band.fm1) { recC[i] = 1; if (sc.band.mode === "all") comC[i] = 1; }
      }
    }
    return { recC, comC };
  }

  function evaluate(sc) {
    const M = model[sc.k], m = M.m;
    const dens = E.density(m, G, sc.dT, m.fp_t, m.fg_t, D.tgrid, new Float64Array(N));
    const hyp = new Float64Array(N).fill(1);
    if (sc.hypDepth > 0 && sc.hypSev > 0) for (let i = 0; i < N; i++) if (depth[i] <= sc.hypDepth) hyp[i] = 1 - sc.hypSev;
    for (let i = 0; i < N; i++) dens[i] *= hyp[i];
    const { recC, comC } = closureMasks(sc);
    let wr = 0, wrC = 0, wc = 0, wcC = 0;
    for (let i = 0; i < N; i++) {
      const a = area[i] * dens[i], r = a * availRec[i];
      wr += r; wc += a; if (recC[i]) wrC += r; if (comC[i]) wcC += a;
    }
    const Cr = Math.min(wrC / wr, 0.95), Cc = Math.min(wcC / wc, 0.95);
    const sq = M.sq;
    const Drec = sc.displace && sq ? Math.min(2, (1 - sq.Cr) / (1 - Cr)) : 1; // at most doubles effort in open water
    const Dcom = sc.displace && sq ? Math.min(2, (1 - sq.Cc) / (1 - Cc)) : 1;
    const reg = { bag: sc.regs.bag, sizeCm: sc.regs.sizeIn * 2.54 };
    const g = E.bagFactor(M.sp.rec, reg.bag);
    const FrO = M.cal.FrecEnc0 * sc.regs.recEffort * Drec, FcO = M.cal.Fcom0 * sc.regs.comMult * Dcom;
    const sprC = [], Fapi = [], FrecC = [], FcomC = [];
    for (let k = 0; k < 4; k++) {
      const ro = k >> 1, co = k & 1, Fr = ro ? FrO : 0, Fc = co ? FcO : 0;
      sprC.push(E.sbpr(M.sch, E.fAtAge(M.sch, M.sp, Fr, Fc, reg)) / M.cal.sbpr0);
      Fapi.push(Fr * g + Fc); FrecC.push(Fr * g); FcomC.push(Fc);
    }
    const cls = new Uint8Array(N), spr = new Float64Array(N);
    for (let i = 0; i < N; i++) { cls[i] = (recC[i] ? 0 : 2) + (comC[i] ? 0 : 1); spr[i] = sprC[cls[i]]; }
    return { sc, dens, hyp, recC, comC, cls, spr, sprC, Fapi, FrecC, FcomC, Cr, Cc, Drec, Dcom, g, reg };
  }

  function summarize(R) {
    const M = model[R.sc.k], sq = M.sq, [w0, w1] = state.win;
    const t = M.m.t_range;
    let H = 0, S = 0, Hsq = 0, Ssq = 0, Hcl = 0, Hex = 0;
    const bins = new Float64Array(DEPTH_BINS.length), binsSq = new Float64Array(DEPTH_BINS.length);
    for (let i = 0; i < N; i++) {
      const h = R.dens[i] * area[i], hs = sq.dens[i] * area[i];
      bins[binOf[i]] += h; binsSq[binOf[i]] += hs;
      if (depth[i] < w0 || depth[i] > w1) continue;
      H += h; S += h * R.spr[i]; Hsq += hs; Ssq += hs * sq.spr[i];
      if (R.recC[i] || R.comC[i]) Hcl += h;
      const T = tclim[i] + R.sc.dT; if (T < t[0] || T > t[1]) Hex += h;
    }
    // temperature-effect uncertainty from posterior draws of the temperature smooth
    let ciH = null, ciS = null;
    if (R.sc.dT !== 0) {
      const m = M.m, nd = m.fp_draws.length, rh = [], rs = [], tmp = new Float32Array(N);
      for (let d = 0; d < nd; d++) {
        E.density(m, G, R.sc.dT, m.fp_draws[d], m.fg_draws[d], D.tgrid, tmp);
        const d0 = M.dens0D[d];
        let h = 0, s = 0, h0 = 0, s0 = 0;
        for (let i = 0; i < N; i++) {
          if (depth[i] < w0 || depth[i] > w1) continue;
          const a = tmp[i] * R.hyp[i] * area[i], a0 = d0[i] * area[i];
          h += a; s += a * R.spr[i]; h0 += a0; s0 += a0 * sq.spr[i];
        }
        rh.push(h / h0 - 1); rs.push(s / s0 - 1);
      }
      const qq = (arr, p) => { const a = arr.slice().sort((x, y) => x - y); return a[Math.min(a.length - 1, Math.max(0, Math.round(p * (a.length - 1))))]; };
      ciH = [qq(rh, 0.025), qq(rh, 0.975)]; ciS = [qq(rs, 0.025), qq(rs, 0.975)];
    }
    return { H, S, Hsq, Ssq, dH: H / Hsq - 1, dS: S / Ssq - 1, sprW: S / H, sprWsq: Ssq / Hsq, closedShare: Hcl / H,
      extrap: Hex / H, ciH, ciS, bins, binsSq };
  }

  function refreshSq(k) {
    const M = model[k];
    M.sq = null;
    const sq = evaluate(sqScenario(k));
    M.sq = sq;
  }

  // ---------------- colour ----------------
  let theme = {};
  function hexToRgb(h) { h = h.trim().replace("#", ""); if (h.length === 3) h = h.split("").map((c) => c + c).join(""); const n = parseInt(h, 16); return [(n >> 16) & 255, (n >> 8) & 255, n & 255]; }
  function lut(stops) {
    const rgb = stops.map(hexToRgb), out = new Array(256);
    for (let i = 0; i < 256; i++) {
      const t = (i / 255) * (rgb.length - 1), j = Math.min(rgb.length - 2, Math.floor(t)), f = t - j;
      const c = rgb[j].map((v, q) => Math.round(v + (rgb[j + 1][q] - v) * f));
      out[i] = `rgb(${c[0]},${c[1]},${c[2]})`;
    }
    return out;
  }
  function readTheme() {
    const cs = getComputedStyle(document.documentElement), v = (n) => cs.getPropertyValue(n).trim();
    theme = {
      water: v("--water"), land: v("--land"), coast: v("--coast"), reg: v("--reg"), regFill: v("--reg-fill"),
      ink: v("--ink"), ink2: v("--ink-2"), ink3: v("--ink-3"), rule: v("--rule"), panel: v("--panel"), accent: v("--accent"),
      hab: lut([0, 1, 2, 3, 4, 5, 6].map((i) => v("--hab-" + i))),
      div: lut([v("--div-neg"), v("--div-mid"), v("--div-pos")]),
      f: lut([0, 1, 2, 3, 4].map((i) => v("--f-" + i))),
      t: lut([0, 1, 2, 3, 4].map((i) => v("--t-" + i))),
      habStops: [0, 1, 2, 3, 4, 5, 6].map((i) => v("--hab-" + i)),
      divStops: [v("--div-neg"), v("--div-mid"), v("--div-pos")],
      fStops: [0, 1, 2, 3, 4].map((i) => v("--f-" + i)),
      tStops: [0, 1, 2, 3, 4].map((i) => v("--t-" + i)),
    };
  }

  // ---------------- layers ----------------
  const LAYERS = [
    { id: "habitat", label: "Habitat" },
    { id: "spawn", label: "Spawning potential" },
    { id: "change", label: "Change vs 2026" },
    { id: "fishing", label: "Fishing mortality" },
    { id: "temp", label: "Bottom temp" },
  ];
  function cellColor(i, R, M) {
    const L = state.layer;
    if (L === "habitat" || L === "spawn") {
      const v = L === "habitat" ? R.dens[i] : R.dens[i] * R.spr[i];
      const t = v > 0 ? (Math.log10(v) - M.lo) / (M.hi - M.lo) : 0;
      return theme.hab[Math.max(0, Math.min(255, Math.round(t * 255)))];
    }
    if (L === "change") {
      const base = M.sq.dens[i] * M.sq.spr[i];
      if (base <= Math.pow(10, M.lo)) return theme.div[128];
      const r = (R.dens[i] * R.spr[i]) / base - 1;
      const t = Math.max(-1, Math.min(1, r / 0.5));
      return theme.div[Math.round((t + 1) * 127.5)];
    }
    if (L === "fishing") {
      const fmax = Math.max(M.sq.Fapi[3] * 2.5, 1e-6);
      const t = Math.min(1, R.Fapi[R.cls[i]] / fmax);
      return theme.f[Math.round(t * 255)];
    }
    const T = tclim[i] + R.sc.dT, t = (T - 3) / (10 - 3);
    return theme.t[Math.max(0, Math.min(255, Math.round(t * 255)))];
  }

  // ---------------- map ----------------
  const canvas = $("#map"), stage = $("#stage"), tip = $("#tip");
  const ctx = canvas.getContext("2d");
  let view = { frames: [], W: 0, H: 0, dpr: 1 };
  const bucket = new Map();
  for (let i = 0; i < N; i++) {
    const key = Math.floor(G.x[i] / 2000) + ":" + Math.floor(G.y[i] / 2000);
    if (!bucket.has(key)) bucket.set(key, []); bucket.get(key).push(i);
  }
  function nearestCell(x, y) {
    const bx = Math.floor(x / 2000), by = Math.floor(y / 2000);
    let best = -1, bd = 1450 * 1450;
    for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) {
      const l = bucket.get(bx + dx + ":" + (by + dy)); if (!l) continue;
      for (const i of l) { const d = (G.x[i] - x) ** 2 + (G.y[i] - y) ** 2; if (d < bd) { bd = d; best = i; } }
    }
    return best;
  }

  const mapwrap = stage.parentElement;
  function layout() {
    const cs = getComputedStyle(mapwrap);
    const W = mapwrap.clientWidth - parseFloat(cs.paddingLeft) - parseFloat(cs.paddingRight), F = D.frames, gap = 14;
    const dx = F[0].x1 - F[0].x0, dy = F.map((f) => f.y1 - f.y0);
    let s, side, H;
    if (window.innerWidth >= 1240) {
      const Hs = mapwrap.clientHeight - parseFloat(cs.paddingTop) - parseFloat(cs.paddingBottom)
        - $(".toolbar").offsetHeight - $("#legend").offsetHeight - 2 * 10;
      const sSide = Math.min((W - gap) / 2 / dx, Hs / Math.max(...dy));
      const sStack = Math.min(W / dx, (Hs - gap) / (dy[0] + dy[1]));
      side = sSide >= sStack; s = side ? sSide : sStack;
      H = side ? Math.max(...dy) * s : (dy[0] + dy[1]) * s + gap;
      stage.style.height = Math.floor(H) + "px";
    } else {
      side = W >= 520;
      s = side ? (W - gap) / 2 / dx : W / dx;
      H = side ? Math.max(...dy) * s : (dy[0] + dy[1]) * s + gap;
      stage.style.height = Math.ceil(H) + "px";
    }
    const fw = dx * s, used = side ? fw * 2 + gap : fw, ox = Math.max(0, (W - used) / 2);
    view.frames = F.map((f, j) => {
      const fh = (f.y1 - f.y0) * s;
      const left = side ? ox + j * (fw + gap) : ox;
      const top = side ? 0 : (j === 0 ? 0 : (F[0].y1 - F[0].y0) * s + gap);
      return { ...f, left, top, w: fw, h: fh, s };
    });
    view.W = W; view.H = H; view.dpr = window.devicePixelRatio || 1;
    canvas.width = Math.round(W * view.dpr); canvas.height = Math.round(H * view.dpr);
    canvas.style.width = W + "px"; canvas.style.height = H + "px";
  }
  const px = (f, x) => f.left + (x - f.x0) * f.s;
  const py = (f, y) => f.top + (f.y1 - y) * f.s;
  function path(f, pts, close) {
    ctx.beginPath();
    pts.forEach((p, j) => (j ? ctx.lineTo(px(f, p[0]), py(f, p[1])) : ctx.moveTo(px(f, p[0]), py(f, p[1]))));
    if (close) ctx.closePath();
  }

  let hatch = null;
  function hatchPattern() {
    const c = document.createElement("canvas"); c.width = c.height = 8;
    const g = c.getContext("2d"); g.strokeStyle = theme.reg; g.globalAlpha = 0.55; g.lineWidth = 1.2;
    g.beginPath(); g.moveTo(0, 8); g.lineTo(8, 0); g.moveTo(-2, 2); g.lineTo(2, -2); g.moveTo(6, 10); g.lineTo(10, 6); g.stroke();
    return ctx.createPattern(c, "repeat");
  }

  function drawMap(R) {
    const M = model[state.sp];
    ctx.setTransform(view.dpr, 0, 0, view.dpr, 0, 0);
    ctx.clearRect(0, 0, view.W, view.H);
    hatch = hatchPattern();
    const [w0, w1] = state.win;
    for (const f of view.frames) {
      ctx.save();
      ctx.beginPath(); ctx.rect(f.left, f.top, f.w, f.h); ctx.clip();
      ctx.fillStyle = theme.water; ctx.fillRect(f.left, f.top, f.w, f.h);
      // graticule
      ctx.strokeStyle = theme.rule; ctx.lineWidth = 0.7; ctx.globalAlpha = 0.7;
      for (const gl of D.geo.graticule) { path(f, gl.pts); ctx.stroke(); }
      ctx.globalAlpha = 1;
      // cells
      const cs = 2000 * f.s + 0.7;
      for (let i = 0; i < N; i++) {
        const y = G.y[i]; if (y < f.y0 - 2000 || y > f.y1 + 2000) continue;
        const out = depth[i] < w0 || depth[i] > w1;
        ctx.globalAlpha = out ? 0.16 : 1;
        ctx.fillStyle = cellColor(i, R, M);
        ctx.fillRect(px(f, G.x[i]) - cs / 2, py(f, y) - cs / 2, cs, cs);
      }
      ctx.globalAlpha = 1;
      // depth-band closure hatch
      if (state.overlays.closures && state.band.mode !== "off") {
        ctx.fillStyle = hatch;
        for (let i = 0; i < N; i++) {
          const fm = depth[i] / FM; if (fm < state.band.fm0 || fm > state.band.fm1) continue;
          const y = G.y[i]; if (y < f.y0 - 2000 || y > f.y1 + 2000) continue;
          ctx.fillRect(px(f, G.x[i]) - cs / 2, py(f, y) - cs / 2, cs, cs);
        }
      }
      // contours
      if (state.overlays.contours) {
        ctx.strokeStyle = theme.ink3; ctx.lineWidth = 0.8;
        const labelled = new Set();
        for (const c of D.geo.contours) {
          ctx.globalAlpha = c.depth === 200 ? 0.9 : 0.55;
          ctx.setLineDash(c.depth >= 500 ? [3, 3] : []);
          path(f, c.pts); ctx.stroke();
        }
        ctx.setLineDash([]); ctx.globalAlpha = 1;
        ctx.font = "500 10px " + getComputedStyle(document.body).getPropertyValue("--font-data");
        ctx.fillStyle = theme.ink2; ctx.textAlign = "center"; ctx.textBaseline = "middle";
        for (const c of D.geo.contours) {
          if (labelled.has(c.depth) || c.pts.length < 30) continue;
          const yMid = f.y0 + (f.y1 - f.y0) * (c.depth === 200 ? 0.62 : c.depth === 100 ? 0.42 : c.depth === 500 ? 0.3 : 0.18);
          let bp = null, bd = 1e12;
          for (const p of c.pts) { const d = Math.abs(p[1] - yMid); if (d < bd && p[1] > f.y0 && p[1] < f.y1) { bd = d; bp = p; } }
          if (!bp || bd > 15000) continue;
          labelled.add(c.depth);
          const X = px(f, bp[0]), Y = py(f, bp[1]), txt = c.depth + " m";
          ctx.save(); ctx.fillStyle = theme.water; ctx.globalAlpha = 0.85;
          const tw = ctx.measureText(txt).width; ctx.fillRect(X - tw / 2 - 2, Y - 6, tw + 4, 12); ctx.restore();
          ctx.fillText(txt, X, Y);
        }
      }
      // closures (chart convention: magenta boundaries for regulated areas)
      if (state.overlays.closures) {
        for (const c of CLOSURES) {
          const mode = state.closures[c.id];
          path(f, c.ring, true);
          if (mode !== "off") { ctx.fillStyle = theme.regFill; ctx.fill(); }
          ctx.strokeStyle = theme.reg; ctx.lineWidth = mode === "off" ? 1 : 1.8;
          ctx.globalAlpha = mode === "off" ? 0.5 : 1;
          ctx.setLineDash(mode === "off" ? [2, 3] : [7, 3]);
          ctx.stroke(); ctx.setLineDash([]); ctx.globalAlpha = 1;
        }
      }
      // land
      ctx.fillStyle = theme.land; ctx.strokeStyle = theme.coast; ctx.lineWidth = 0.9;
      for (const r of D.geo.land) { path(f, r, true); ctx.fill(); ctx.stroke(); }
      // survey tows
      if (state.overlays.hauls) {
        const H = D.hauls, v = H[state.sp];
        for (let j = 0; j < v.length; j++) {
          const y = H.y[j]; if (y < f.y0 || y > f.y1) continue;
          const X = px(f, H.x[j]), Y = py(f, y);
          ctx.beginPath();
          if (v[j] > 0) { ctx.arc(X, Y, 1.4 + Math.min(3, Math.log10(1 + v[j]) * 1.3), 0, 6.283); ctx.fillStyle = theme.ink; ctx.globalAlpha = 0.7; ctx.fill(); }
          else { ctx.arc(X, Y, 1.3, 0, 6.283); ctx.strokeStyle = theme.ink3; ctx.globalAlpha = 0.55; ctx.lineWidth = 0.7; ctx.stroke(); }
        }
        ctx.globalAlpha = 1;
      }
      // labels
      const bodyFont = getComputedStyle(document.body).getPropertyValue("--font-body");
      const waterFont = getComputedStyle(document.body).getPropertyValue("--font-water");
      const dataFont = getComputedStyle(document.body).getPropertyValue("--font-data");
      const dispFont = getComputedStyle(document.body).getPropertyValue("--font-display");
      ctx.textBaseline = "middle"; ctx.textAlign = "left";
      ctx.font = "500 11px " + bodyFont;
      for (const p of D.geo.ports) {
        if (p.y < f.y0 + 3000 || p.y > f.y1 - 3000) continue;
        const X = px(f, p.x), Y = py(f, p.y);
        ctx.fillStyle = theme.ink; ctx.beginPath(); ctx.arc(X, Y, 2.4, 0, 6.283); ctx.fill();
        ctx.lineWidth = 3; ctx.strokeStyle = theme.land; ctx.strokeText(p.name, X + 6, Y);
        ctx.fillStyle = theme.ink2; ctx.fillText(p.name, X + 6, Y);
      }
      ctx.font = "italic 400 12px " + waterFont; ctx.fillStyle = theme.ink2; ctx.textAlign = "center";
      for (const ft of D.features) {
        if (ft.y < f.y0 || ft.y > f.y1) continue;
        const big = ft.name.startsWith("P A");
        ctx.font = (big ? "italic 400 14px " : "italic 400 12px ") + waterFont;
        ctx.globalAlpha = big ? 0.75 : 1;
        if (big) { ctx.save(); ctx.translate(px(f, ft.x), py(f, ft.y)); ctx.rotate(-Math.PI / 2); ctx.fillText(ft.name, 0, 0); ctx.restore(); }
        else { ctx.lineWidth = 3; ctx.strokeStyle = theme.water; ctx.strokeText(ft.name, px(f, ft.x), py(f, ft.y)); ctx.fillText(ft.name, px(f, ft.x), py(f, ft.y)); }
        ctx.globalAlpha = 1;
      }
      // latitude labels
      ctx.font = "400 10px " + dataFont; ctx.fillStyle = theme.ink3; ctx.textAlign = "left";
      for (const gl of D.geo.graticule) {
        if (gl.kind !== "lat") continue;
        const p = gl.pts[0], q = gl.pts[gl.pts.length - 1];
        const yAt = p[1] + (q[1] - p[1]) * ((f.x0 - p[0]) / (q[0] - p[0]));
        const Yl = py(f, yAt); if (Yl < f.top + 52 || Yl > f.top + f.h - 34) continue;
        ctx.fillText(gl.label, f.left + 6, Yl - 7);
      }
      // frame title
      ctx.textAlign = "left"; ctx.textBaseline = "alphabetic";
      ctx.font = "400 10.5px " + bodyFont;
      const subW = ctx.measureText(f.sub).width;
      ctx.font = "700 13px " + dispFont;
      const boxW = Math.max(ctx.measureText(f.label.toUpperCase()).width, subW) + 14;
      ctx.save(); ctx.globalAlpha = 0.88; ctx.fillStyle = theme.panel;
      ctx.fillRect(f.left + 3, f.top + 4, boxW, 36); ctx.restore();
      ctx.fillStyle = theme.ink;
      ctx.fillText(f.label.toUpperCase(), f.left + 10, f.top + 20);
      ctx.font = "400 10.5px " + bodyFont; ctx.fillStyle = theme.ink2;
      ctx.fillText(f.sub, f.left + 10, f.top + 34);
      // scale bar: 10 nautical miles
      const nm = 18520 * f.s, sx = f.left + 12, sy = f.top + f.h - 16;
      ctx.fillStyle = theme.ink; ctx.fillRect(sx, sy, nm / 2, 3); ctx.strokeStyle = theme.ink; ctx.lineWidth = 1; ctx.strokeRect(sx, sy, nm, 3);
      ctx.font = "400 9.5px " + dataFont; ctx.fillStyle = theme.ink2; ctx.textBaseline = "bottom";
      ctx.fillText("0", sx - 2, sy - 2); ctx.fillText("10 nmi", sx + nm - 12, sy - 2);
      // north arrow
      const nx = f.left + f.w - 18, ny = f.top + f.h - 34;
      ctx.beginPath(); ctx.moveTo(nx, ny - 9); ctx.lineTo(nx + 5, ny + 6); ctx.lineTo(nx, ny + 2); ctx.lineTo(nx - 5, ny + 6); ctx.closePath();
      ctx.fillStyle = theme.ink; ctx.fill();
      ctx.font = "600 9px " + dataFont; ctx.textAlign = "center"; ctx.textBaseline = "top"; ctx.fillText("N", nx, ny + 8);
      // pinned / hovered cell
      for (const [idx, col, lw] of [[state.hover, theme.ink, 1.2], [state.pinned, theme.accent, 2.2]]) {
        if (idx === null || idx < 0) continue;
        if (G.y[idx] < f.y0 || G.y[idx] > f.y1) continue;
        ctx.strokeStyle = col; ctx.lineWidth = lw;
        ctx.strokeRect(px(f, G.x[idx]) - cs / 2 - 1, py(f, G.y[idx]) - cs / 2 - 1, cs + 2, cs + 2);
      }
      ctx.restore();
      ctx.strokeStyle = theme.rule; ctx.lineWidth = 1; ctx.strokeRect(f.left + 0.5, f.top + 0.5, f.w - 1, f.h - 1);
    }
  }

  function frameAt(mx, my) { return view.frames.find((f) => mx >= f.left && mx <= f.left + f.w && my >= f.top && my <= f.top + f.h); }
  function cellAt(ev) {
    const r = canvas.getBoundingClientRect(), mx = ev.clientX - r.left, my = ev.clientY - r.top;
    const f = frameAt(mx, my); if (!f) return { idx: -1, mx, my };
    const x = f.x0 + (mx - f.left) / f.s, y = f.y1 - (my - f.top) / f.s;
    return { idx: nearestCell(x, y), mx, my, x, y };
  }

  function cellRows(i, R) {
    const M = model[state.sp], sp = D.params[state.sp];
    const T = tclim[i] + R.sc.dT;
    const p = 1 / (1 + Math.exp(-(M.m.eta_p[i] + M.m.fp_t[E.tIndex(T, D.tgrid.t0, D.tgrid.dt, D.tgrid.n)])));
    const names = CLOSURES.filter((c) => state.closures[c.id] !== "off" && c.cells.includes(i)).map((c) => c.name + (state.closures[c.id] === "rec" ? " (sport)" : ""));
    if (state.band.mode !== "off") { const fm = depth[i] / FM; if (fm >= state.band.fm0 && fm <= state.band.fm1) names.push(`Depth-band closure ${state.band.fm0}–${state.band.fm1} fm`); }
    const tr = M.m.t_range, ex = T < tr[0] || T > tr[1];
    return {
      title: `${fmtLat(G.lat[i])} ${fmtLon(G.lon[i])}`,
      rows: [
        ["Depth", `${depth[i].toFixed(0)} m · ${(depth[i] / FM).toFixed(0)} fm`],
        ["Bottom temp", `${T.toFixed(1)} °C${R.sc.dT ? ` (${tclim[i].toFixed(1)} ${R.sc.dT > 0 ? "+" : MINUS} ${Math.abs(R.sc.dT).toFixed(1)})` : ""}${ex ? " *" : ""}`],
        ["Encounter prob.", fmtPctAbs(p * R.hyp[i])],
        [`${sp.name} density`, `${fmtDens(R.dens[i])} kg/ha`],
        ["Spawning ratio", fmtPctAbs(R.spr[i])],
        ["Fishing mort.", `${R.Fapi[R.cls[i]].toFixed(3)} yr⁻¹`],
      ],
      closures: names, extrap: ex,
    };
  }
  function rowsHtml(o) {
    return `<h4>${o.title}</h4><dl>${o.rows.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join("")}</dl>` +
      (o.closures.length ? `<p class="hint" style="margin-top:6px;color:var(--reg)">Closed: ${o.closures.map(esc).join("; ")}</p>` : "") +
      (o.extrap ? `<p class="hint" style="margin-top:4px">* Outside the 3.2–9.0 °C range the survey sampled; temperature effect held at its edge value.</p>` : "");
  }

  canvas.addEventListener("pointermove", (ev) => {
    const h = cellAt(ev);
    if (h.idx !== state.hover) { state.hover = h.idx; scheduleDraw(); }
    if (h.idx < 0) {
      const f = frameAt(h.mx, h.my);
      if (f && ev.pointerType === "mouse") {
        tip.innerHTML = `<p class="hint" style="margin:0">No survey grid here. The trawl survey samples 55–1,280 m, so nearshore reefs and the deep slope are not mapped.</p>`;
        placeTip(h.mx, h.my); tip.hidden = false;
      } else tip.hidden = true;
      return;
    }
    if (ev.pointerType !== "mouse") return;
    tip.innerHTML = rowsHtml(cellRows(h.idx, current)) + `<p class="hint" style="margin-top:6px">Click to pin this cell.</p>`;
    placeTip(h.mx, h.my); tip.hidden = false;
  });
  function placeTip(mx, my) {
    const tw = tip.offsetWidth || 230, th = tip.offsetHeight || 160;
    let x = mx + 14, y = my + 14;
    if (x + tw > view.W) x = mx - tw - 14;
    if (y + th > view.H) y = Math.max(4, my - th - 14);
    tip.style.left = Math.max(4, x) + "px"; tip.style.top = y + "px";
  }
  canvas.addEventListener("pointerleave", () => { tip.hidden = true; state.hover = null; scheduleDraw(); });
  canvas.addEventListener("click", (ev) => {
    const h = cellAt(ev);
    state.pinned = h.idx >= 0 && h.idx !== state.pinned ? h.idx : null;
    renderInspector(); scheduleDraw();
  });
  document.addEventListener("keydown", (ev) => { if (ev.key === "Escape" && state.pinned !== null) { state.pinned = null; renderInspector(); scheduleDraw(); } });

  // ---------------- legend ----------------
  function gradientCss(stops) { return `linear-gradient(90deg, ${stops.join(", ")})`; }
  function renderLegend() {
    const M = model[state.sp], el = $("#legend"); let ramp = "";
    const tickVals = (lo, hi) => [lo, (lo + hi) / 2, hi].map((v) => fmtDens(Math.pow(10, v)));
    if (state.layer === "habitat" || state.layer === "spawn") {
      const t = tickVals(M.lo, M.hi);
      ramp = `<div class="ramp"><span>${state.layer === "habitat" ? "Expected survey catch" : "Habitat × spawning potential ratio"} (kg/ha, log scale)</span>
        <div class="bar" style="background:${gradientCss(theme.habStops)}"></div><div class="ticks"><span>≤${t[0]}</span><span>${t[1]}</span><span>≥${t[2]}</span></div></div>`;
    } else if (state.layer === "change") {
      ramp = `<div class="ramp"><span>Spawning potential vs 2026 status quo</span><div class="bar" style="background:${gradientCss(theme.divStops)}"></div>
        <div class="ticks"><span>≤${MINUS}50%</span><span>0</span><span>≥+50%</span></div></div>`;
    } else if (state.layer === "fishing") {
      const fmax = M.sq.Fapi[3] * 2.5;
      ramp = `<div class="ramp"><span>Fishing mortality on fully selected fish (yr⁻¹)</span><div class="bar" style="background:${gradientCss(theme.fStops)}"></div>
        <div class="ticks"><span>0</span><span>${(fmax / 2).toFixed(3)}</span><span>≥${fmax.toFixed(3)}</span></div></div>`;
    } else {
      ramp = `<div class="ramp"><span>Bottom temperature (°C)</span><div class="bar" style="background:${gradientCss(theme.tStops)}"></div>
        <div class="ticks"><span>≤3</span><span>6.5</span><span>≥10</span></div></div>`;
    }
    const keys = [];
    if (state.overlays.closures) keys.push(`<span class="key"><span class="sw" style="border:2px dashed var(--reg);background:var(--reg-fill)"></span>Closure in effect</span>`,
      `<span class="key"><span class="sw" style="border:1px dotted var(--reg)"></span>Defined, not in effect</span>`);
    if (state.overlays.hauls) keys.push(`<span class="key"><span class="sw" style="background:var(--ink);border-radius:50%;width:9px;height:9px"></span>Tow with ${esc(D.params[state.sp].name.toLowerCase())}</span>`,
      `<span class="key"><span class="sw" style="border:1px solid var(--ink-3);border-radius:50%;width:8px;height:8px"></span>Tow without</span>`);
    keys.push(`<span class="key"><span class="sw" style="background:var(--water);border:1px solid var(--rule)"></span>Not surveyed</span>`);
    el.innerHTML = ramp + keys.join("");
  }

  // ---------------- controls ----------------
  const controls = $("#controls");
  function sliderHtml(id, label, min, max, step, value, valText, ref, disabled) {
    return `<div class="ctl"><div class="ctl-head"><label for="${id}">${label}</label><span class="val" id="${id}-v">${valText}</span></div>
      <input type="range" id="${id}" min="${min}" max="${max}" step="${step}" value="${value}" ${disabled ? "disabled" : ""}>
      ${ref ? `<span class="ref">${ref}</span>` : ""}</div>`;
  }
  function segHtml(id, value, opts) {
    return `<div class="seg" role="group" data-seg="${id}">${opts.map(([v, l]) => `<button type="button" data-v="${v}" aria-pressed="${v === value}">${l}</button>`).join("")}</div>`;
  }
  const bagText = (b) => (b === 0 ? "No retention" : b + " fish");
  const sizeText = (s) => (s === 0 ? "None" : `${s} in · ${(s * 2.54).toFixed(0)} cm`);
  const pctText = (v) => Math.round(v * 100) + "%";
  const dTText = (v) => (v > 0 ? "+" : v < 0 ? MINUS : "±") + Math.abs(v).toFixed(1) + " °C";
  const hypText = (d) => (d === 0 ? "Off" : `to ${d} m · ${(d / FM).toFixed(0)} fm`);
  const depthText = (d) => `${d.toLocaleString()} m · ${(d / FM).toFixed(0)} fm`;

  function renderControls() {
    const k = state.sp, P = D.params[k], r = state.regs[k], a = state.adv[k];
    const anom = D.tempAnomByYear;
    const yrs = Object.keys(anom), warm = yrs.reduce((p, c) => (anom[c] > anom[p] ? c : p)), cool = yrs.reduce((p, c) => (anom[c] < anom[p] ? c : p));
    const yrca = CLOSURES.filter((c) => c.kind === "yrca"), res = CLOSURES.filter((c) => c.kind === "reserve");
    const closureRow = (c) => `<div class="closure"><span class="nm">${esc(c.name)}</span>
      ${segHtml("cl:" + c.id, state.closures[c.id], [["off", "Off"], ["rec", "Sport"], ["all", "All"]])}
      <span class="meta">${c.cite.replace('50 CFR ', '')} · ${Math.round(c.area_km2)} km² · ${c.cells.length} cells</span></div>`;
    controls.innerHTML = `
      <div class="group">
        <h2>Species</h2>
        <div class="species" role="radiogroup" aria-label="Species">
          ${SPECIES.map((s) => { const p = D.params[s]; return `<button type="button" class="sp-card" role="radio" aria-checked="${s === k}" data-sp="${s}">
            <span class="dot"></span><span><b>${p.name}</b> <i>${p.sci}</i></span><small>${p.niche}</small></button>`; }).join("")}
        </div>
      </div>
      <div class="group">
        <h2>Ocean conditions</h2>
        ${sliderHtml("dT", "Bottom temperature anomaly", -1.5, 2, 0.1, state.dT, dTText(state.dT), "0 = average 2003–2024 bottom temperature at each cell")}
        <div class="chips" aria-label="Temperature presets">
          <button type="button" class="chip" data-dt="0">Average</button>
          <button type="button" class="chip" data-dt="${anom[cool]}">${cool} (coolest survey year)</button>
          <button type="button" class="chip" data-dt="${anom[warm]}">${warm} (warmest)</button>
          <button type="button" class="chip" data-dt="1">+1 °C</button>
          <button type="button" class="chip" data-dt="2">+2 °C</button>
        </div>
        ${sliderHtml("hyp", "Low-oxygen water on the shelf", 0, 200, 5, state.hypDepth, hypText(state.hypDepth), "Scenario overlay: habitat shallower than this depth is reduced")}
        ${sliderHtml("hypSev", "Habitat lost inside low-oxygen zone", 0, 1, 0.05, state.hypSev, pctText(state.hypSev), "", state.hypDepth === 0)}
      </div>
      <div class="group">
        <h2>Fishing rules</h2>
        <p class="hint">${esc(P.regs.note)}</p>
        ${sliderHtml("bag", P.regs.bag_label, 0, P.regs.bag_max, 1, r.bag, bagText(r.bag), `2026 rule: ${P.regs.bag}${P.regs.bag_prior !== P.regs.bag ? ` · previous: ${P.regs.bag_prior}` : ""}`)}
        ${P.regs.size_max_in > 0
          ? sliderHtml("size", "Minimum size, sport", 0, P.regs.size_max_in, 1, r.sizeIn, sizeText(r.sizeIn), `2026 rule: ${P.regs.size_in} in`)
          : `<p class="hint">No minimum size for ${esc(P.name.toLowerCase())} in 2026 Oregon sport rules.</p>`}
        ${sliderHtml("recE", "Sport fishing effort", 0, 1.5, 0.05, r.recEffort, pctText(r.recEffort), "Season length and trips, relative to 2026")}
        ${sliderHtml("com", "Commercial mortality", 0, 2, 0.05, r.comMult, pctText(r.comMult), "Trip limits and allocation, relative to current")}
        <p class="sub-h">Yelloweye Rockfish Conservation Areas</p>
        <div>${yrca.map(closureRow).join("")}</div>
        <p class="hint">Stonewall Bank is closed to Oregon sport groundfish fishing. The others are defined in federal rules and can be activated in season.</p>
        <p class="sub-h">Marine reserves</p>
        <div>${res.map(closureRow).join("")}</div>
        <p class="hint">Reserves sit mostly inshore of the survey grid, so they cover few mapped cells.</p>
        <p class="sub-h">Depth-band closure (RCA style)</p>
        <div class="closure" style="border-top:0"><span class="nm">Close ${state.band.fm0}–${state.band.fm1} fm</span>
          ${segHtml("band", state.band.mode, [["off", "Off"], ["rec", "Sport"], ["all", "All"]])}</div>
        ${sliderHtml("fm0", "Shallow edge", 10, 300, 5, state.band.fm0, state.band.fm0 + " fm · " + Math.round(state.band.fm0 * FM) + " m", "", state.band.mode === "off")}
        ${sliderHtml("fm1", "Deep edge", 20, 400, 5, state.band.fm1, state.band.fm1 + " fm · " + Math.round(state.band.fm1 * FM) + " m", "", state.band.mode === "off")}
        <label class="check"><input type="checkbox" id="displace" ${state.displace ? "checked" : ""}><span>Displaced effort moves to open water<br><span class="hint">Up to double the current effort in open cells; the rest is lost. Off: closed-area effort is removed.</span></span></label>
        <button type="button" class="btn" id="reset">Reset to 2026 rules and average conditions</button>
      </div>
      <div class="group">
        <h2>Depth window</h2>
        ${sliderHtml("w0", "Shallowest", 30, 1000, 10, state.win[0], depthText(state.win[0]))}
        ${sliderHtml("w1", "Deepest", 60, 1830, 10, state.win[1], depthText(state.win[1]))}
        <p class="hint">Totals and the map use only cells inside this window.</p>
      </div>
      <div class="group">
        <details class="adv">
          <summary>Model assumptions for ${esc(P.name.toLowerCase())}</summary>
          <div class="adv-body">
            <p class="hint">Changing these refits current fishing mortality so 2026 conditions still reproduce the assessment's spawning potential ratio (${fmtPctAbs(P.status.spr)}, ${P.status.year}).</p>
            ${sliderHtml("share", "Sport share of current fishing mortality", 0, 1, 0.01, a.share, pctText(a.share))}
            ${sliderHtml("mu", "Legal fish encountered per angler-trip", 0.1, 4, 0.1, a.mu, a.mu.toFixed(1))}
            ${sliderHtml("dmort", "Release mortality", 0, 0.6, 0.01, a.dmort, pctText(a.dmort))}
            <button type="button" class="btn" id="resetAdv">Restore defaults</button>
          </div>
        </details>
      </div>`;
  }

  controls.addEventListener("click", (ev) => {
    const spb = ev.target.closest("[data-sp]");
    if (spb) { selectSpecies(spb.dataset.sp); return; }
    const chip = ev.target.closest("[data-dt]");
    if (chip) { state.dT = Math.round(parseFloat(chip.dataset.dt) * 100) / 100; syncSlider("dT", state.dT, dTText(state.dT)); update(); return; }
    const segBtn = ev.target.closest(".seg button");
    if (segBtn) {
      const key = segBtn.parentElement.dataset.seg, v = segBtn.dataset.v;
      segBtn.parentElement.querySelectorAll("button").forEach((b) => b.setAttribute("aria-pressed", b === segBtn));
      if (key === "band") { state.band.mode = v; $("#fm0").disabled = $("#fm1").disabled = v === "off"; }
      else state.closures[key.slice(3)] = v;
      update(); return;
    }
    if (ev.target.id === "reset") {
      state.dT = 0; state.hypDepth = 0; state.hypSev = 0.5; state.regs[state.sp] = regsDefault(state.sp);
      CLOSURES.forEach((c) => (state.closures[c.id] = defaultClosureMode(c)));
      state.band = { mode: "off", fm0: 30, fm1: 100 }; state.displace = true; state.win = [30, 1830];
      renderControls(); update(); return;
    }
    if (ev.target.id === "resetAdv") { state.adv[state.sp] = advDefault(state.sp); buildModel(state.sp, state.adv[state.sp]); refreshSq(state.sp); renderControls(); $("details.adv").open = true; update(); }
  });
  function syncSlider(id, v, text) { const el = $("#" + id); if (el) el.value = v; const t = $("#" + id + "-v"); if (t) t.textContent = text; }
  let advTimer = null;
  controls.addEventListener("input", (ev) => {
    const id = ev.target.id, v = parseFloat(ev.target.value), r = state.regs[state.sp], a = state.adv[state.sp];
    const set = (t) => { const el = $("#" + id + "-v"); if (el) el.textContent = t; };
    switch (id) {
      case "dT": state.dT = v; set(dTText(v)); break;
      case "hyp": state.hypDepth = v; set(hypText(v)); $("#hypSev").disabled = v === 0; break;
      case "hypSev": state.hypSev = v; set(pctText(v)); break;
      case "bag": r.bag = v; set(bagText(v)); break;
      case "size": r.sizeIn = v; set(sizeText(v)); break;
      case "recE": r.recEffort = v; set(pctText(v)); break;
      case "com": r.comMult = v; set(pctText(v)); break;
      case "fm0": state.band.fm0 = Math.min(v, state.band.fm1 - 5); set(state.band.fm0 + " fm · " + Math.round(state.band.fm0 * FM) + " m"); bandLabel(); break;
      case "fm1": state.band.fm1 = Math.max(v, state.band.fm0 + 5); set(state.band.fm1 + " fm · " + Math.round(state.band.fm1 * FM) + " m"); bandLabel(); break;
      case "w0": state.win[0] = Math.min(v, state.win[1] - 20); set(depthText(state.win[0])); break;
      case "w1": state.win[1] = Math.max(v, state.win[0] + 20); set(depthText(state.win[1])); break;
      case "displace": state.displace = ev.target.checked; break;
      case "share": a.share = v; set(pctText(v)); break;
      case "mu": a.mu = v; set(v.toFixed(1)); break;
      case "dmort": a.dmort = v; set(pctText(v)); break;
      default: return;
    }
    if (id === "share" || id === "mu" || id === "dmort") {
      clearTimeout(advTimer);
      advTimer = setTimeout(() => { buildModel(state.sp, a); refreshSq(state.sp); update(); }, 60);
      return;
    }
    update();
  });
  controls.addEventListener("change", (ev) => { if (ev.target.id === "displace") { state.displace = ev.target.checked; update(); } });
  function bandLabel() { const el = controls.querySelector('[data-seg="band"]'); if (el) el.previousElementSibling.textContent = `Close ${state.band.fm0}–${state.band.fm1} fm`; }

  function selectSpecies(k) {
    state.sp = k; state.pinned = null;
    renderControls(); update();
  }

  // ---------------- toolbar ----------------
  $("#layers").innerHTML = LAYERS.map((l) => `<button type="button" data-layer="${l.id}" aria-pressed="${l.id === state.layer}">${l.label}</button>`).join("");
  $("#layers").addEventListener("click", (ev) => {
    const b = ev.target.closest("[data-layer]"); if (!b) return;
    state.layer = b.dataset.layer;
    $("#layers").querySelectorAll("button").forEach((x) => x.setAttribute("aria-pressed", x === b));
    renderLegend(); scheduleDraw();
  });
  $("#overlays").innerHTML = [["closures", "Closures"], ["contours", "Depth contours"], ["hauls", "Survey tows"]]
    .map(([k, l]) => `<label><input type="checkbox" data-ov="${k}" ${state.overlays[k] ? "checked" : ""}>${l}</label>`).join("");
  $("#overlays").addEventListener("change", (ev) => { const k = ev.target.dataset.ov; if (!k) return; state.overlays[k] = ev.target.checked; renderLegend(); scheduleDraw(); });

  // ---------------- charts ----------------
  const chartTip = $("#chartTip");
  document.addEventListener("pointerover", (ev) => {
    const t = ev.target.closest && ev.target.closest("[data-tip]");
    if (!t) { chartTip.hidden = true; return; }
    chartTip.textContent = t.dataset.tip; chartTip.hidden = false;
  });
  document.addEventListener("pointermove", (ev) => {
    if (chartTip.hidden) return;
    const w = chartTip.offsetWidth; let x = ev.clientX + 12; if (x + w > window.innerWidth - 8) x = ev.clientX - w - 12;
    chartTip.style.left = x + "px"; chartTip.style.top = ev.clientY + 14 + "px";
  });

  function sprBullet(sum, sp) {
    const W = 360, H = 78, l = 10, r = 10, x = (v) => l + v * (W - l - r);
    const tgt = sp.status.target, cur = sum.sprW, base = sum.sprWsq;
    const ticks = [0, 0.25, 0.5, 0.75, 1];
    return `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Spawning potential ratio ${fmtPctAbs(cur)}, target ${fmtPctAbs(tgt)}">
      <rect style="fill:var(--panel-2);stroke:var(--rule)" x="${x(0)}" y="22" width="${x(1) - x(0)}" height="18" rx="3"/>
      <rect style="fill:var(--bad)" x="${x(0)}" y="22" width="${x(tgt) - x(0)}" height="18" rx="3" opacity="0.12"/>
      <rect style="fill:var(--ink)" x="${x(0)}" y="26" width="${Math.max(2, x(cur) - x(0))}" height="10" rx="2" data-tip="Scenario: ${fmtPctAbs(cur, 1)}"/>
      <line style="stroke:var(--accent)" x1="${x(base)}" x2="${x(base)}" y1="18" y2="44" stroke-width="2.5" data-tip="2026 status quo: ${fmtPctAbs(base, 1)}"/>
      <line style="stroke:var(--bad)" x1="${x(tgt)}" x2="${x(tgt)}" y1="14" y2="48" stroke-width="1.5" stroke-dasharray="3 2"/>
      <text x="${x(tgt)}" y="10" text-anchor="middle" class="lab">Target SPR${Math.round(tgt * 100)}%</text>
      ${ticks.map((t) => `<text x="${x(t)}" y="60" text-anchor="${t === 0 ? "start" : t === 1 ? "end" : "middle"}">${Math.round(t * 100)}%</text>`).join("")}
      <text x="${x(0)}" y="74" class="lab">▬ Scenario</text><text style="fill:var(--accent)" x="${x(0.3)}" y="74" class="lab">│ 2026 status quo</text>
      <text x="${x(1)}" y="74" text-anchor="end">100% = unfished</text>
    </svg>`;
  }

  function depthBars(sum) {
    const W = 360, H = 170, l = 34, r = 6, t = 10, b = 38;
    const tot = sum.binsSq.reduce((a, c) => a + c, 0);
    const sq = Array.from(sum.binsSq, (v) => v / tot), sc = Array.from(sum.bins, (v) => v / tot);
    const ymax = Math.max(0.05, ...sq, ...sc) * 1.1, n = DEPTH_BINS.length, bw = (W - l - r) / n;
    const y = (v) => t + (1 - v / ymax) * (H - t - b);
    const step = ymax > 0.4 ? 0.2 : ymax > 0.2 ? 0.1 : 0.05;
    let grid = "";
    for (let v = 0; v <= ymax + 1e-9; v += step) grid += `<line class="grid" x1="${l}" x2="${W - r}" y1="${y(v)}" y2="${y(v)}"/><text x="${l - 4}" y="${y(v) + 3}" text-anchor="end">${Math.round(v * 100)}%</text>`;
    const [w0, w1] = state.win;
    const bars = DEPTH_BINS.map(([a, bb], i) => {
      const x0 = l + i * bw, inWin = bb > w0 && a < w1;
      const lab = i === 0 ? "<100" : i === n - 1 ? ">1,300" : `${a}–${bb}`;
      return `<g opacity="${inWin ? 1 : 0.35}">
        <rect style="stroke:var(--ink-3)" x="${x0 + bw * 0.14}" y="${y(sq[i])}" width="${bw * 0.34}" height="${y(0) - y(sq[i])}" fill="none" stroke-width="1.2" rx="2"/>
        <rect style="fill:var(--hab-4)" x="${x0 + bw * 0.52}" y="${y(sc[i])}" width="${bw * 0.34}" height="${y(0) - y(sc[i])}" rx="2"/>
        <rect class="hit" x="${x0}" y="${t}" width="${bw}" height="${H - t - b}" data-tip="${lab} m: status quo ${fmtPctAbs(sq[i], 1)}, scenario ${fmtPctAbs(sc[i], 1)}"/>
        ${i > 0 ? `<line class="axis" x1="${x0}" x2="${x0}" y1="${y(0)}" y2="${y(0) + 4}"/><text x="${x0}" y="${H - b + 14}" text-anchor="middle">${a >= 1000 ? (a / 1000).toFixed(1) + "k" : a}</text>` : ""}</g>`;
    }).join("");
    return `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Share of Oregon habitat index by depth band">
      ${grid}<line class="axis" x1="${l}" x2="${W - r}" y1="${y(0)}" y2="${y(0)}"/>${bars}
      <text x="${(W + l) / 2}" y="${H - 10}" text-anchor="middle" class="lab">Depth (m)</text></svg>`;
  }

  function trendChart(sp) {
    const obs = D.species[sp].obs_by_year, yrs = Object.keys(obs).map(Number).sort((a, b) => a - b);
    const W = 360, H = 140, l = 34, r = 16, t = 10, b = 26;
    const vals = yrs.map((yv) => obs[yv][0]), ymax = Math.min(1, Math.max(...vals) * 1.15);
    const x = (yv) => l + ((yv - yrs[0]) / (yrs[yrs.length - 1] - yrs[0])) * (W - l - r);
    const y = (v) => t + (1 - v / ymax) * (H - t - b);
    let d = "", prev = null;
    yrs.forEach((yv) => { d += (prev !== null && yv - prev === 1 ? "L" : "M") + x(yv).toFixed(1) + " " + y(obs[yv][0]).toFixed(1); prev = yv; });
    const mean = vals.reduce((a, c) => a + c, 0) / vals.length;
    const step = ymax > 0.5 ? 0.2 : ymax > 0.2 ? 0.1 : 0.05;
    let grid = "";
    for (let v = 0; v <= ymax + 1e-9; v += step) grid += `<line class="grid" x1="${l}" x2="${W - r}" y1="${y(v)}" y2="${y(v)}"/><text x="${l - 4}" y="${y(v) + 3}" text-anchor="end">${Math.round(v * 100)}%</text>`;
    const xt = [2003, 2008, 2013, 2018, 2024].map((yv, j, a) => `<text x="${x(yv)}" y="${H - 8}" text-anchor="${j === 0 ? "start" : j === a.length - 1 ? "end" : "middle"}">${yv}</text>`).join("");
    const dots = yrs.map((yv) => `<circle style="fill:var(--hab-5);stroke:var(--panel)" cx="${x(yv)}" cy="${y(obs[yv][0])}" r="${yv === yrs[yrs.length - 1] ? 4.5 : 3}" stroke-width="1.5"/>
      <rect class="hit" x="${x(yv) - 8}" y="${t}" width="16" height="${H - t - b}" data-tip="${yv}: caught in ${fmtPctAbs(obs[yv][0])} of tows, mean ${fmtDens(obs[yv][1])} kg/ha"/>`).join("");
    return `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Share of Oregon survey tows that caught the species, by year">
      ${grid}<line class="axis" x1="${l}" x2="${W - r}" y1="${y(0)}" y2="${y(0)}"/>
      <line style="stroke:var(--ink-3)" x1="${l}" x2="${W - r}" y1="${y(mean)}" y2="${y(mean)}" stroke-dasharray="4 3"/>
      <text x="${W - r}" y="${y(mean) - 4}" text-anchor="end">mean ${fmtPctAbs(mean)}</text>
      <path style="stroke:var(--hab-5)" d="${d}" fill="none" stroke-width="2"/>${dots}${xt}
      <text x="${x(2020)}" y="${y(0) - 6}" text-anchor="middle">no survey</text></svg>`;
  }

  function curveChart(series, opts) {
    const W = 170, H = 130, l = 34, r = 6, t = 8, b = 30;
    const xs = series[0].c.x, xlog = opts.log;
    const tx = (v) => (xlog ? Math.log10(v) : v);
    const x0 = tx(xs[0]), x1 = tx(xs[xs.length - 1]);
    let ylo = Infinity, yhi = -Infinity;
    series.forEach((s) => { ylo = Math.min(ylo, ...s.c.lo); yhi = Math.max(yhi, ...s.c.hi); });
    const pad = (yhi - ylo) * 0.05; ylo -= pad; yhi += pad;
    const X = (v) => l + ((tx(v) - x0) / (x1 - x0)) * (W - l - r), Y = (v) => t + (1 - (v - ylo) / (yhi - ylo)) * (H - t - b);
    const band = (s) => "M" + s.c.x.map((v, i) => `${X(v).toFixed(1)} ${Y(s.c.hi[i]).toFixed(1)}`).join("L") + "L" +
      s.c.x.slice().reverse().map((v, i) => `${X(v).toFixed(1)} ${Y(s.c.lo[s.c.x.length - 1 - i]).toFixed(1)}`).join("L") + "Z";
    const line = (s) => "M" + s.c.x.map((v, i) => `${X(v).toFixed(1)} ${Y(s.c.y[i]).toFixed(1)}`).join("L");
    const yt = [ylo + pad, (ylo + yhi) / 2, yhi - pad].map((v) => `<text x="${l - 3}" y="${Y(v) + 3}" text-anchor="end">${(Math.round(v) || 0).toString().replace("-", "\u2212")}</text><line class="grid" x1="${l}" x2="${W - r}" y1="${Y(v)}" y2="${Y(v)}"/>`).join("");
    const xt = opts.ticks.map((v) => `<text x="${X(v)}" y="${H - b + 12}" text-anchor="middle">${v}</text>`).join("");
    return `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${opts.title}">${yt}
      ${series.map((s) => `<path d="${band(s)}" style="fill:${s.col}" opacity="0.16"/><path d="${line(s)}" style="fill:none;stroke:${s.col}" stroke-width="2"/>`).join("")}
      ${xt}<text x="${(W + l) / 2}" y="${H - 4}" text-anchor="middle" class="lab">${opts.xlab}</text></svg>`;
  }

  // ---------------- insights ----------------
  const insights = $("#insights");
  function renderInsights(R, S) {
    const k = state.sp, P = D.params[k], m = D.species[k];
    const isSq = Math.abs(S.dH) < 5e-4 && Math.abs(S.dS) < 5e-4;
    const ci = (c) => (c ? `95% range ${fmtPct(c[0])} to ${fmtPct(c[1])}` : "");
    const lede = isSq
      ? `This is the <b>2026 status quo</b> for ${P.name.toLowerCase()}: Oregon sport rules as of 2026 and average 2003–2024 bottom temperatures. Move a control to compare a scenario against it.`
      : `Against 2026 rules and average conditions, this scenario changes the Oregon ${P.name.toLowerCase()} habitat index by <b>${fmtPct(S.dH)}</b> and spawning potential by <b>${fmtPct(S.dS)}</b>${state.win[0] > 30 || state.win[1] < 1830 ? ` within ${state.win[0]}–${state.win[1].toLocaleString()} m` : ""}.`;
    const target = P.status.target, above = S.sprW >= target;
    const spr = S.sprW;
    insights.innerHTML = `
      <p class="lede">${lede}</p>
      <div class="tiles">
        <div class="tile"><span class="k">Habitat index</span><span class="v">${fmtPct(S.dH)}</span><span class="s">${S.ciH ? ci(S.ciH) : "vs 2026 baseline"}</span></div>
        <div class="tile"><span class="k">Spawning potential</span><span class="v">${fmtPct(S.dS)}</span><span class="s">${S.ciS ? ci(S.ciS) : "habitat × SPR, vs 2026"}</span></div>
        <div class="tile"><span class="k">Spawning ratio (SPR)</span><span class="v">${fmtPctAbs(spr)}</span>
          <span class="pill ${above ? "good" : spr >= target * 0.9 ? "warn" : "bad"}">${above ? "✓ Above" : "▼ Below"} ${Math.round(target * 100)}% target</span></div>
        <div class="tile"><span class="k">Habitat in closures</span><span class="v">${fmtPctAbs(S.closedShare, S.closedShare < 0.1 ? 1 : 0)}</span><span class="s">closed to at least one fleet</span></div>
      </div>
      ${S.extrap > 0.01 ? `<p class="hint" style="margin-top:8px"><span class="pill warn">⚠ Extrapolation</span> ${fmtPctAbs(S.extrap)} of the habitat index sits where bottom temperature leaves the 3.2–9.0 °C the survey sampled.</p>` : ""}
      <div class="card" id="inspectCard"></div>
      <div class="card">
        <h3>Spawning potential ratio</h3>
        <p>Equilibrium spawning output per recruit relative to unfished, weighted by where the fish are. The ${P.status.year} assessment estimate is ${fmtPctAbs(P.status.spr)} (${esc(P.status.basis)}).</p>
        <div class="chart">${sprBullet(S, P)}</div>
      </div>
      <div class="card">
        <h3>Where the habitat sits</h3>
        <p>Share of the 2026 Oregon habitat index by depth. Outline: status quo. Filled: this scenario.</p>
        <div class="chart">${depthBars(S)}</div>
      </div>
      <div class="card">
        <h3>Survey record</h3>
        <p>Share of Oregon tows that caught ${P.name.toLowerCase()} each year (${m.n_pos.toLocaleString()} of ${m.n_hauls.toLocaleString()} tows).</p>
        <div class="chart">${trendChart(k)}</div>
      </div>
      <div class="card">
        <h3>How well the habitat model predicts</h3>
        <table class="metrics"><thead><tr><th>Hold-out test</th><th>AUC (encounter)</th><th>Spearman (density)</th></tr></thead><tbody>
          <tr><td>Random 5-fold</td><td>${m.metrics.auc_random.toFixed(3)}</td><td>${m.metrics.rho_random.toFixed(3)}</td></tr>
          <tr><td>Latitude blocks</td><td>${m.metrics.auc_block.toFixed(3)}</td><td>${m.metrics.rho_block.toFixed(3)}</td></tr>
        </tbody></table>
        <p class="hint">Deviance explained: encounter ${fmtPctAbs(m.metrics.dev_presence)}, positive catch ${fmtPctAbs(m.metrics.dev_positive)}. Latitude blocks hold out each fifth of the coast in turn.</p>
        <div class="twocol">
          <div class="chart">${curveChart([{ c: m.curves.depth, col: "var(--cat-1)" }], { log: true, ticks: [100, 300, 1000], xlab: "Depth (m)", title: "Depth effect on encounter" })}</div>
          <div class="chart">${curveChart([{ c: m.curves.temp_p, col: "var(--cat-1)" }, { c: m.curves.temp_g, col: "var(--cat-2)" }], { log: false, ticks: [4, 6, 8], xlab: "Bottom temp (°C)", title: "Temperature effects" })}</div>
        </div>
        <p class="hint"><span class="sw" style="background:var(--cat-1);height:3px;vertical-align:middle"></span> encounter (logit) &nbsp; <span class="sw" style="background:var(--cat-2);height:3px;vertical-align:middle"></span> log catch when present. Bands are 95% intervals; the temperature effect is estimated after accounting for depth and location, so it is weakly identified.</p>
      </div>
      ${methodsHtml(k)}`;
    renderInspector();
  }

  function renderInspector() {
    const el = $("#inspectCard"); if (!el || !current) return;
    if (state.pinned === null || state.pinned < 0) {
      el.innerHTML = `<p class="hint">Click any cell on the map to pin its depth, temperature, density and fishing pressure here.</p>`;
      return;
    }
    const o = cellRows(state.pinned, current);
    el.innerHTML = `<h3>Pinned cell</h3><div class="inspector">${rowsHtml(o).replace("<dl>", '<dl class="kv">')}
      <p class="hint" style="margin-top:6px">Press Esc or click the cell again to unpin.</p></div>`;
  }

  function methodsHtml(k) {
    const P = D.params[k], lh = P.lh, g = lh.growth;
    const src = '<span class="tag src">assessment</span>', asm = '<span class="tag asm">assumed</span>';
    const growthTxt = g.a2 === null ? `L∞ ${g.Linf} cm, k ${g.k}, L(${g.a1}) ${g.L1} cm` : `L(${g.a1}) ${g.L1} cm, L(${g.a2}) ${g.L2} cm, k ${g.k}`;
    const matTxt = lh.mat.type === "age" ? `50% mature at age ${lh.mat.x50}` : `50% mature at ${lh.mat.x50} cm`;
    return `<div class="card method">
      <h3>Method, assumptions and sources</h3>
      <h4>Data</h4>
      <p>NOAA Northwest Fisheries Science Center West Coast Groundfish Bottom Trawl Survey, ${D.meta.yearMin}–${D.meta.yearMax} (no survey in 2020): ${D.meta.nHauls.toLocaleString()} tows between 42.0°N and 46.25°N at 55–1,280 m, with tow-level bottom temperature. Read from the public surveyjoin data mirror. Predictions use the survey's own 2-km sampling grid (${D.meta.nCells.toLocaleString()} Oregon cells), shown in UTM zone 10N.</p>
      <h4>Habitat model</h4>
      <p>A delta-lognormal generalized additive model per species: the chance a tow catches the species, and log catch per hectare when it does, each modelled on log depth, bottom temperature, a smooth surface over location, and a year effect. Smoothness was chosen by grid search. The baseline averages all survey years on the response scale. Bottom temperature at each grid cell comes from a separate GAM fitted to tow temperatures (deviance explained 91.5%). A temperature anomaly is added uniformly; the fitted temperature effect is held flat outside the 3.2–9.0 °C the survey sampled. The 95% ranges come from 60 posterior draws of the temperature smooth.</p>
      <h4>Fishing layer</h4>
      <p>Female equilibrium per-recruit model. Current fishing mortality is calibrated so open water under 2026 rules reproduces the assessment's spawning potential ratio. Bag limits act through a negative-binomial count of legal fish encountered per angler-trip; fish released over the bag or under the size limit die at the release mortality rate. Closures remove the chosen fleet's fishing mortality in covered cells, and displaced effort can move to open water in proportion to fish density, up to double the current effort there (sport effort fades past about 150 m). Each cell's ratio is computed locally, which assumes adults stay put; read the spawning potential map as where spawning output would come from, not as a population projection.</p>
      <div class="tablewrap"><table><thead><tr><th>${esc(P.name)}</th><th>Value</th><th></th></tr></thead><tbody>
        <tr><td>Natural mortality (female)</td><td class="n">${lh.M} yr⁻¹</td><td>${src}</td></tr>
        <tr><td>Growth</td><td class="n">${growthTxt}</td><td>${src}</td></tr>
        <tr><td>Weight-length (kg, cm)</td><td class="n">${lh.wl.a.toExponential(3)} · L^${lh.wl.b}</td><td>${src}</td></tr>
        <tr><td>Maturity</td><td class="n">${matTxt}, slope ${lh.mat.slope}</td><td>${k === "lingcod" ? `${src} ${asm}` : src}</td></tr>
        <tr><td>Status-quo SPR (${P.status.year})</td><td class="n">${fmtPctAbs(P.status.spr, 1)}</td><td>${src}</td></tr>
        <tr><td>Fraction unfished (${P.status.frac_year})</td><td class="n">${fmtPctAbs(P.status.frac_unfished, 1)}</td><td>${src}</td></tr>
        <tr><td>Sport share of fishing mortality</td><td class="n">${fmtPctAbs(state.adv[k].share)}</td><td>${asm}</td></tr>
        <tr><td>Legal fish per angler-trip (mean, size)</td><td class="n">${state.adv[k].mu}, ${P.rec.theta}</td><td>${asm}</td></tr>
        <tr><td>Release mortality</td><td class="n">${fmtPctAbs(state.adv[k].dmort)}</td><td>${asm}</td></tr>
        <tr><td>Sport / commercial selectivity, 50% at</td><td class="n">${P.rec.sel50} / ${P.com.sel50} cm</td><td>${asm}</td></tr>
      </tbody></table></div>
      <h4>Limits to keep in mind</h4>
      <ul>
        <li>The trawl survey does not sample shallower than 55 m or untrawlable rock, so nearshore reef habitat (prime adult lingcod and much canary habitat) is under-represented.</li>
        <li>Canary rockfish school; a few very large tows dominate raw survey means, so treat canary density as a relative index.</li>
        <li>Assessment estimates are coastwide (lingcod: north of 40°10′ N) and are applied to Oregon waters.</li>
        <li>The low-oxygen overlay is a what-if. Oxygen is not in this survey dataset, so its effect is not fitted.</li>
        <li>Cold bottom-water anomalies on the Oregon shelf usually mean freshly upwelled, low-oxygen water, so the fitted temperature effect partly carries an oxygen signal. Read warming results as an association, not a mechanism.</li>
        <li>Survey data run through ${D.meta.yearMax}; recruitment pulses (for example young sablefish on the shelf) are averaged into the baseline.</li>
      </ul>
      <h4>Sources</h4>
      <ul>
        <li><a href="https://www.pcouncil.org/documents/2021/12/status-of-lingcod-ophiodon-elongatus-along-the-northern-u-s-west-coast-in-2021-december-2021.pdf/" target="_blank" rel="noopener">Taylor et al. 2021, northern lingcod assessment (PFMC)</a></li>
        <li><a href="https://www.pcouncil.org/documents/2024/02/status-of-canary-rockfish-sebastes-pinniger-along-the-u-s-west-coast-in-2023.pdf" target="_blank" rel="noopener">Langseth et al. 2023, canary rockfish assessment (PFMC)</a></li>
        <li><a href="https://www.pcouncil.org/documents/2024/02/status-of-sablefish-anoplopoma-fimbria-along-the-u-s-west-coast-in-2023.pdf" target="_blank" rel="noopener">Johnson et al. 2023, sablefish assessment (PFMC)</a></li>
        <li><a href="https://myodfw.com/sport-bottomfish-seasons" target="_blank" rel="noopener">ODFW 2026 sport bottomfish regulations</a> and <a href="https://myodfw.com/news/commission-sets-2026-groundfish-seasons" target="_blank" rel="noopener">2026 season announcement</a></li>
        <li><a href="https://myodfw.com/sites/default/files/2025-10/Stonewall_Bank_Yelloweye_Rockfish_Conservation_Area.pdf" target="_blank" rel="noopener">ODFW, Stonewall Bank YRCA sport closure</a></li>
        <li><a href="https://www.ecfr.gov/current/title-50/chapter-VI/part-660/subpart-C/section-660.70" target="_blank" rel="noopener">50 CFR 660.70, groundfish conservation area coordinates</a></li>
        <li><a href="https://oregon.public.law/rules/oar_635-012-0040" target="_blank" rel="noopener">OAR 635-012-0040 and 141-142, Oregon marine reserve boundaries</a></li>
        <li><a href="https://github.com/DFO-NOAA-Pacific/surveyjoin" target="_blank" rel="noopener">surveyjoin (NOAA/DFO) survey data and WCGBTS grid</a></li>
        <li><a href="https://www.fisheries.noaa.gov/inport/item/18418" target="_blank" rel="noopener">NOAA InPort: West Coast Groundfish Bottom Trawl Survey</a></li>
      </ul>
      <p class="hint">Coastline from the GLOBE 30-arc-second land mask. Built ${D.meta.built}.</p>
    </div>`;
  }

  // ---------------- update loop ----------------
  let current = null, summary = null, drawQueued = false;
  function scheduleDraw() { if (drawQueued) return; drawQueued = true; requestAnimationFrame(() => { drawQueued = false; if (current) drawMap(current); }); }
  let updQueued = false;
  function update() {
    if (updQueued) return; updQueued = true;
    requestAnimationFrame(() => {
      updQueued = false;
      if (!model[state.sp]) { buildModel(state.sp, state.adv[state.sp]); refreshSq(state.sp); }
      current = evaluate(curScenario());
      summary = summarize(current);
      renderInsights(current, summary);
      renderLegend();
      drawMap(current);
    });
  }

  function onTheme() { readTheme(); renderLegend(); scheduleDraw(); }
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", onTheme);
  new MutationObserver(onTheme).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
  new ResizeObserver(() => { layout(); scheduleDraw(); }).observe(mapwrap);

  $("#stamp").textContent = `NOAA trawl survey ${D.meta.yearMin}–${D.meta.yearMax} · ${D.meta.nHauls.toLocaleString()} Oregon tows · 2026 ODFW rules`;
  readTheme();
  SPECIES.forEach((k) => { buildModel(k, state.adv[k]); refreshSq(k); });
  renderControls();
  layout();
  update();
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => scheduleDraw());
})();
