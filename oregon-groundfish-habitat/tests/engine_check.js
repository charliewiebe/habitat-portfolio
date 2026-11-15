const E = require("../app/engine.js");
const P = require("../data/species_params.json");
for (const [k, sp] of Object.entries(P)) {
  const sch = E.schedule(sp), cal = E.calibrate(sp, sch);
  const spr = (bag, sizeCm, mr = 1, mc = 1) => E.sbpr(sch, E.fAtAge(sch, sp, cal.FrecEnc0 * mr, cal.Fcom0 * mc, { bag, sizeCm })) / cal.sbpr0;
  console.log(k, JSON.stringify({ theta: +cal.theta.toFixed(4), FrecEnc0: +cal.FrecEnc0.toFixed(4), Fcom0: +cal.Fcom0.toFixed(4), sprCheck: +cal.sprCheck.toFixed(4), target: sp.status.spr }));
  console.log("   L(a):", [0, 1, 2, 3, 5, 10, 20].filter(a => a <= sp.lh.amax).map(a => a + ":" + sch.L[a].toFixed(1)).join(" "), " mat50ish:", Array.from(sch.mat).findIndex(m => m >= 0.5));
  console.log("   R(b):", [0, 1, 2, 3, 5, 10].map(b => b + ":" + E.expectedRetained(sp.rec.mu, sp.rec.theta, b).toFixed(3)).join(" "));
  const b0 = sp.regs.bag, s0 = sp.regs.size_cm;
  console.log("   SPR sq", spr(b0, s0).toFixed(4), "| bag-1", spr(b0 - 1, s0).toFixed(4), "| bag+1", spr(b0 + 1, s0).toFixed(4), "| bag0", spr(0, s0).toFixed(4),
    "| size+2in", spr(b0, s0 + 5.08).toFixed(4), "| nosize", spr(b0, 0).toFixed(4), "| com0", spr(b0, s0, 1, 0).toFixed(4), "| unfished", spr(b0, s0, 0, 0).toFixed(4));
}
