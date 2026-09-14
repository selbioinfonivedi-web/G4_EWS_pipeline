/* Boot web/workstation/static/g4.js against a minimal DOM and assert what the
   default surface shows.

   This is not a browser: it checks structure, not paint, and the shim below
   implements only the handful of DOM behaviours g4.js actually uses. It exists
   because the interesting property -- "the default surface is small, and
   nothing was deleted to make it small" -- is a statement about the rendered
   mode bar and view list, which no amount of reading the source asserts.

   Run by tests/web/test_ui_surface.py. Exits non-zero on the first failure,
   printing every check either way. */
class N {
  constructor(tag){ this.tag=tag; this.kids=[]; this.attrs={}; this._text="";
    this.style={}; this.classList={ _s:new Set(),
      add:(c)=>this.classList._s.add(c), remove:(c)=>this.classList._s.delete(c),
      toggle:(c,on)=>{ on?this.classList._s.add(c):this.classList._s.delete(c); },
      contains:(c)=>this.classList._s.has(c) }; }
  set className(v){ this.attrs.class=v; }
  /* the real DOM reflects .title onto the attribute; mirror that */
  set title(v){ this.attrs.title=v; }
  get title(){ return this.attrs.title||""; }
  get className(){ return this.attrs.class||""; }
  set textContent(v){ this._text=v; this.kids=[]; }
  get textContent(){ return this._text + this.kids.map(k=>typeof k==="string"?k:k.textContent).join(""); }
  set innerHTML(v){ this._text=String(v).replace(/<[^>]*>/g,""); }
  setAttribute(k,v){ this.attrs[k]=v; }
  getAttribute(k){ return this.attrs[k]??null; }
  removeAttribute(k){ delete this.attrs[k]; }
  addEventListener(){} removeEventListener(){}
  replaceChildren(...k){ this.kids=k.filter(x=>x!=null&&x!==false); }
  append(...k){ for(const x of k) if(x!=null&&x!==false) this.kids.push(x); }
  remove(){} click(){} focus(){}
  getBoundingClientRect(){ return {top:0,left:0,width:100,height:20,bottom:20,right:100}; }
  querySelectorAll(sel){ const want = sel.replace(".",""); 
    const out=[]; const walk=(n)=>{ for(const k of n.kids) if(k instanceof N){ 
      if((k.attrs.class||"").split(/\s+/).includes(want)) out.push(k); walk(k); } }; walk(this); return out; }
  matches(){ return false; }
  /* walk every descendant of a given tag */
  all(tag, out=[]){ for(const k of this.kids){ if(k instanceof N){ if(k.tag===tag) out.push(k); k.all(tag,out); } } return out; }
}
const REG = new Map();
const node = (id) => { if(!REG.has(id)) REG.set(id, new N("div")); return REG.get(id); };
global.N = N;
global.document = {
  createElement:(t)=>new N(t), createElementNS:(ns,t)=>new N(t),
  documentElement:new N("html"), body:new N("body"), head:new N("head"),
  querySelector:(sel)=> sel.startsWith("#") ? node(sel.slice(1)) : new N("div"),
  /* the spine nav holds the stage nodes, as it does in app.html */
  querySelectorAll:()=>[], addEventListener(){}, getElementById:(id)=>node(id),
};
global.window = global;
global.addEventListener = () => {};
global.removeEventListener = () => {};
global.requestAnimationFrame = (f) => f();
global.getComputedStyle = () => ({ getPropertyValue: () => "#000" });
global.matchMedia = () => ({ matches:false, addEventListener(){} });
global.navigator = { clipboard:{ writeText:()=>Promise.resolve() } };
global.localStorage = { _m:new Map(), getItem(k){return this._m.has(k)?this._m.get(k):null;},
                        setItem(k,v){this._m.set(k,String(v));} };
global.fetch = () => Promise.reject(new Error("no network in the shim"));
global.URL = { createObjectURL:()=>"blob:x", revokeObjectURL(){} };
global.Blob = class {};
global.EventSource = class { constructor(){} close(){} };

const fs = require("fs");
let src = fs.readFileSync("web/workstation/static/g4.js","utf8");
src = src.replace(/\bboot\(\);\s*$/, "");           // don't hit the network
eval(src + "\nglobal.__x = { S, MODES, VIZ, VIZ_MORE, ALL_VIZ, renderModes, syncChrome, renderIdentity, workVisualize };");

const X = global.__x, fail = [];
const ok = (c, m) => { if(!c) fail.push(m); else console.log("  ok  " + m); };

/* default surface */
X.syncChrome(); X.renderModes();
const bar = node("modes");
const modeBtns = bar.kids.filter(k => (k.attrs.class||"").startsWith("mode"));
ok(modeBtns.length === 5, `mode bar shows 5 modes by default (got ${modeBtns.length}: ${modeBtns.map(b=>b.textContent).join("|")})`);
ok(document.body.classList.contains("no-spine"), "spine column collapsed by default");
ok(node("flow").hidden === true, "spine hidden by default");
ok(X.S.mode === "analyses", `lands on Analyses (got ${X.S.mode})`);

/* disclosure */
const disc = bar.kids.find(k => /More steps|Fewer steps/.test(k.textContent||""));
ok(!!disc, "a More steps disclosure exists");
X.S.showSteps = true; X.syncChrome(); X.renderModes();
const all = node("modes").kids.filter(k => (k.attrs.class||"").startsWith("mode"));
ok(all.length === 9, `disclosure reveals all 9 modes (got ${all.length})`);
ok(!document.body.classList.contains("no-spine"), "spine column returns with the steps");
ok(node("flow").hidden === false, "spine shown with the steps");

/* an advanced mode stays visible while selected, even folded away */
X.S.showSteps = false; X.S.mode = "configure"; X.renderModes();
const sel = node("modes").kids.filter(k => (k.attrs.class||"").startsWith("mode"));
ok(sel.some(b => /Configure/.test(b.textContent)), "the selected advanced mode does not vanish under it");

/* views */
ok(Object.keys(X.VIZ).length === 4, `4 primary views (got ${Object.keys(X.VIZ).length})`);
ok(Object.keys(X.ALL_VIZ).length === 10, `all 10 views still reachable (got ${Object.keys(X.ALL_VIZ).length})`);

/* identity: reference and period survive as the organism tooltip */
X.S.data = { identity:{ display_name:"FMDV", reference:"NC_004004", n_samples:7, n_raw:9, period:[1994,2024] } };
X.renderIdentity();
ok(/NC_004004/.test(node("idf-organism").attrs.title||"") && /1994–2024/.test(node("idf-organism").attrs.title||""),
   `reference and period kept on hover (${node("idf-organism").attrs.title})`);

if (fail.length) { console.log("\nFAILED:"); fail.forEach(f=>console.log("  x  " + f)); process.exit(1); }
console.log("\nall chrome checks passed");
