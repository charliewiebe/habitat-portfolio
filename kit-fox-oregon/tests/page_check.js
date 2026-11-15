// Smoke test for the built web app: runs the page script against a stub DOM, steps through every
// training set x scenario x period x map layer and every outlook metric, and checks the baseline
// verdict against the pipeline's own viability numbers (outputs/timeseries.csv).
// Usage, from the project root: node tests/page_check.js index.html
const fs = require("fs");
const html = fs.readFileSync(process.argv[2] || "index.html", "utf8");
const scripts = [...html.matchAll(/<script>([\s\S]*?)<\/script>/g)].map(m => m[1]);
const js = scripts[scripts.length - 1];

const ctx2d = () => ({ createImageData: (w, h) => ({ data: new Uint8ClampedArray(w * h * 4) }), putImageData() {}, clearRect() {},
  drawImage() {}, beginPath() {}, moveTo() {}, lineTo() {}, stroke() {}, arc() {}, fill() {}, strokeRect() {}, setLineDash() {} });
const els = {};
const mk = id => els[id] || (els[id] = { id, innerHTML: "", textContent: "", value: "", checked: false, style: {}, className: "",
  clientWidth: 900, width: 0, height: 0, parentElement: { clientWidth: 900 }, addEventListener() {}, getContext: ctx2d,
  setAttribute() {}, getBoundingClientRect: () => ({ left: 0, top: 0, width: 900, height: 700 }) });
global.document = { getElementById: mk, createElement: () => mk("c" + Math.random()), querySelectorAll: () => [],
  documentElement: { getAttribute: () => null } };
global.getComputedStyle = () => ({ getPropertyValue: () => "#336699" });
global.matchMedia = () => ({ matches: false, addEventListener() {} });
global.requestAnimationFrame = f => { f(); return 1; }; global.cancelAnimationFrame = () => {};
global.MutationObserver = class { observe() {} }; global.addEventListener = () => {}; global.window = { devicePixelRatio: 2 };
global.atob = s => Buffer.from(s, "base64").toString("binary");

const ts = fs.readFileSync("outputs/timeseries.csv", "utf8").trim().split("\n").map(l => l.split(","));
const col = n => ts[0].indexOf(n);
const expected = Object.fromEntries(ts.slice(1).filter(r => r[col("ssp")] === "baseline")
  .map(r => [r[col("set")], { suitable: +r[col("suitable_km2")], K: +r[col("largest_network_K")] }]));

eval(js + `
let n = 0;
for (const set of SETS) for (const ssp of M.ssps) for (const period of ["baseline", ...M.periods]) for (const layer of Object.keys(LAYERS)) {
  Object.assign(S, { set, ssp, period, layer }); update(); n++;
}
for (const metric of Object.keys(METRICS)) { S.metric = metric; update(); }
S.rw = "effort"; update();
console.log("rendered " + n + " map states");

// The page quantises the habitat index to 8 bits and recomputes networks in the browser, so
// allow a small tolerance against the pipeline's float numbers.
let worst = 0;
Object.assign(S, DEF, { period: "baseline" });
for (const set of SETS) {
  const m = viability(getLayer(set, "baseline"), S).m, e = expected[set];
  const dS = Math.abs(m.suitable - e.suitable) / e.suitable, dK = Math.abs(m.largestK - e.K) / e.K;
  worst = Math.max(worst, dS, dK);
  console.log(set.padEnd(11), "suitable km2 page " + Math.round(m.suitable) + " vs pipeline " + Math.round(e.suitable),
              "| largest K page " + Math.round(m.largestK) + " vs pipeline " + Math.round(e.K));
}
if (worst > 0.05) { console.error("FAIL: page and pipeline differ by " + (worst * 100).toFixed(1) + "%"); process.exit(1); }
console.log("PASS (max difference " + (worst * 100).toFixed(1) + "%)");
`);
