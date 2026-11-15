const fs = require("fs"); const E = require("../app/engine.js");
const h = fs.readFileSync("index.html", "utf8"); const m = h.match(/const DATA = (\{[\s\S]*?\});\n/); const D = JSON.parse(m[1]);
const G = D.grid, N = G.n;
for (const k of Object.keys(D.species)) {
  const s = D.species[k]; const out = [];
  for (const dT of [-1, 0.25, 1, 2, 3]) {
    const d0 = E.density(s, G, 0, s.fp_t, s.fg_t, D.tgrid, new Float64Array(N)), d1 = E.density(s, G, dT, s.fp_t, s.fg_t, D.tgrid, new Float64Array(N));
    let h0 = 0, h1 = 0, sh0 = 0, sh1 = 0; for (let i = 0; i < N; i++) { h0 += d0[i] * G.area[i]; h1 += d1[i] * G.area[i]; if (G.depth[i] < 200) { sh0 += d0[i] * G.area[i]; sh1 += d1[i] * G.area[i]; } }
    const rs = []; for (let d = 0; d < s.fp_draws.length; d++) { const a = E.density(s, G, 0, s.fp_draws[d], s.fg_draws[d], D.tgrid, new Float64Array(N)), b = E.density(s, G, dT, s.fp_draws[d], s.fg_draws[d], D.tgrid, new Float64Array(N)); let x = 0, y = 0; for (let i = 0; i < N; i++) { x += a[i] * G.area[i]; y += b[i] * G.area[i]; } rs.push(y / x - 1); }
    rs.sort((a, b) => a - b);
    out.push(`dT ${dT}: ${(100 * (h1 / h0 - 1)).toFixed(0)}% [${(100 * rs[1]).toFixed(0)}, ${(100 * rs[58]).toFixed(0)}]  shelf<200m ${(100 * (sh1 / sh0 - 1)).toFixed(0)}%`);
  }
  console.log(k, "\n  " + out.join("\n  "));
}
