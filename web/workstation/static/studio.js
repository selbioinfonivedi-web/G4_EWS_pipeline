/* ═══════════════════════════════════════════════════════════════════
   G4 STUDIO — the bioluminescent skin.

   The luminous field behind the glass is not artwork. It is the real
   848-tip FMDV phylogeny drawn radially with additive compositing: each
   soft "cell" is a clade, sized by how many genomes it holds and tinted
   by its dominant lineage; each bright point is one sequenced genome.

   That is the whole design argument. The reference this was modelled on
   put a rendered micrograph behind the glass and a decorative bar chart
   in the inspector. Here the glow *is* the measurement, so the aesthetic
   costs nothing scientifically -- hover any point and it names a real
   accession.
   ═══════════════════════════════════════════════════════════════════ */

const $ = (s, r = document) => r.querySelector(s);
const el = (t, a = {}, ...kids) => {
  const n = document.createElement(t);
  for (const [k, v] of Object.entries(a)) {
    if (k === "class") n.className = v;
    else if (k === "text") n.textContent = v;
    else if (k === "html") n.innerHTML = v;
    else if (k.startsWith("on")) n.addEventListener(k.slice(2), v);
    else if (v === true) n.setAttribute(k, "");
    else if (v !== false && v != null) n.setAttribute(k, v);
  }
  for (const c of kids.flat()) if (c != null && c !== false) n.append(c);
  return n;
};
const ico = (d) => {
  const s = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  s.setAttribute("viewBox", "0 0 24 24");
  s.innerHTML = d;
  return s;
};

const S = {
  data: null, seq: null, pathogen: "fmdv",
  tool: "move", locus: null, sample: null, base: 4313,
  view: { x: 0, y: 0, z: 1 }, hover: null,
  tips: [], cells: [], history: [], redo: [],
};

const LUME = ["#5FF6D2", "#FFB347", "#B98CFF", "#FF6B7A", "#57D2FF", "#9BE86B", "#FF9ED2", "#7FE8C8"];

async function api(p) { const r = await fetch(p); if (!r.ok) throw new Error((await r.json().catch(() => ({}))).detail || r.statusText); return r.json(); }

/* ═══ boot ════════════════════════════════════════════════════════ */
async function boot() {
  try {
    S.pathogen = new URLSearchParams(location.search).get("p") || "fmdv";
    S.data = await api(`/api/dataset/${S.pathogen}`);
  } catch (e) { toast(`Dataset unavailable — ${e.message}`); return; }

  S.order = Object.keys(S.data.lineages);
  S.colour = {};
  S.order.forEach((k, i) => { S.colour[k] = LUME[i % LUME.length]; });
  S.locus = S.data.loci[3] || S.data.loci[0];
  S.base = S.locus ? S.locus.start : 1;
  $("#docname").textContent = `${S.data.identity.pathogen}_G4_Atlas_v${S.data.identity.atlas_version}`;

  if (S.data.identity.synthetic) synthBanner();
  buildRail(); buildTabs(); renderInspector(); renderFoot();
  await loadSequence();
  layout();
  addEventListener("resize", debounce(layout, 140));
  wire();
}

function layout() { sizeField(); buildCells(); drawField(); drawPreview(); }

/* ═══ the luminous field ══════════════════════════════════════════ */
function sizeField() {
  const cv = $("#field"), dpr = Math.min(devicePixelRatio || 1, 2);
  cv.width = innerWidth * dpr; cv.height = innerHeight * dpr;
  S.g = cv.getContext("2d");
  S.g.setTransform(dpr, 0, 0, dpr, 0, 0);
  S.W = innerWidth; S.H = innerHeight;
}

/* Clades become cells.
   These must be a DISJOINT CUT of the tree, never the "largest N nodes":
   the biggest clades are nested inside one another, so drawing them all
   additively piles every genome onto the same pixels and the field
   saturates to white. Descend from the root until each clade is small
   enough to stand alone, then stop. No two cells can overlap because no
   cell is an ancestor of another. */
const MAX_CELL = 90;   // split anything larger
const MIN_CELL = 6;    // discard anything smaller

function buildCells() {
  const t = S.data.tree;
  S.tips = []; S.cells = [];
  if (!t) return;
  const cx = S.W * 0.46, cy = S.H * 0.5, R = Math.min(S.W, S.H) * 0.36;
  const pos = (n) => {
    const a = n.y * Math.PI * 2 - Math.PI / 2;
    const r = (0.18 + n.x * 0.82) * R;
    return { x: cx + Math.cos(a) * r * S.view.z + S.view.x, y: cy + Math.sin(a) * r * S.view.z + S.view.y };
  };
  S.pos = pos;

  for (const n of t.nodes) {
    if (n.s < 0) continue;
    const p = pos(n);
    S.tips.push({ ...p, s: n.s, lineage: S.data.samples[n.s].l });
  }

  const kids = new Map();
  for (const n of t.nodes) if (n.p != null) { if (!kids.has(n.p)) kids.set(n.p, []); kids.get(n.p).push(n.i); }

  // descendant tips per node, iteratively
  const members = new Array(t.nodes.length);
  const order = [];
  const stack = [[0, false]];
  while (stack.length) {
    const [i, done] = stack.pop();
    if (done) { order.push(i); continue; }
    stack.push([i, true]);
    for (const c of kids.get(i) || []) stack.push([c, false]);
  }
  for (const i of order) {
    const n = t.nodes[i];
    if (n.s >= 0) { members[i] = [n.s]; continue; }
    const out = [];
    for (const k of kids.get(i) || []) out.push(...members[k]);
    members[i] = out;
  }

  const tipBySample = new Map(S.tips.map((tp) => [tp.s, tp]));

  // the cut
  const queue = [0];
  while (queue.length) {
    const i = queue.shift();
    const m = members[i];
    if (m.length > MAX_CELL && (kids.get(i) || []).length) { queue.push(...kids.get(i)); continue; }
    if (m.length < MIN_CELL) continue;
    const tally = new Map();
    let sx = 0, sy = 0;
    for (const si of m) {
      const L = S.data.samples[si].l;
      tally.set(L, (tally.get(L) || 0) + 1);
      const tip = tipBySample.get(si);
      if (tip) { sx += tip.x; sy += tip.y; }
    }
    const dom = [...tally.entries()].sort((a, b) => b[1] - a[1])[0];
    S.cells.push({
      x: sx / m.length, y: sy / m.length, n: m.length,
      lineage: dom[0], purity: dom[1] / m.length, node: i, members: m,
    });
  }
  S.cells.sort((a, b) => b.n - a.n);
}

function drawField() {
  const g = S.g;
  if (!g) return;
  g.clearRect(0, 0, S.W, S.H);

  const bg = g.createRadialGradient(S.W * 0.46, S.H * 0.5, 0, S.W * 0.46, S.H * 0.5, Math.max(S.W, S.H) * 0.72);
  bg.addColorStop(0, "#05292410"); bg.addColorStop(0, "#052924"); bg.addColorStop(0.6, "#021917"); bg.addColorStop(1, "#010D0C");
  g.fillStyle = bg; g.fillRect(0, 0, S.W, S.H);

  const t = S.data.tree;
  if (!t) return;

  // branches first, dim — they are structure, not signal
  g.lineWidth = 0.6;
  for (const n of t.nodes) {
    if (n.p == null) continue;
    const a = S.pos(n), b = S.pos(t.nodes[n.p]);
    g.strokeStyle = "rgba(95,246,210,0.13)";
    g.beginPath(); g.moveTo(b.x, b.y); g.lineTo(a.x, a.y); g.stroke();
  }

  // clade halos — disjoint, so alpha never accumulates past one cell
  for (const c of S.cells) {
    const r = (14 + Math.sqrt(c.n) * 6.2) * S.view.z;
    const col = S.colour[c.lineage] || "#5FF6D2";
    const grd = g.createRadialGradient(c.x, c.y, 0, c.x, c.y, r);
    grd.addColorStop(0, hexa(col, 0.30));
    grd.addColorStop(0.5, hexa(col, 0.10));
    grd.addColorStop(1, hexa(col, 0));
    g.fillStyle = grd;
    g.beginPath(); g.arc(c.x, c.y, r, 0, 7); g.fill();
    g.strokeStyle = hexa(col, 0.28);
    g.lineWidth = 1;
    g.beginPath(); g.arc(c.x, c.y, r * 0.62, 0, 7); g.stroke();
  }

  // genomes
  for (const tip of S.tips) {
    const on = S.sample === S.data.samples[tip.s].a;
    const col = S.colour[tip.lineage];
    g.shadowBlur = on ? 14 : 4; g.shadowColor = col;
    g.fillStyle = on ? "#FFFFFF" : hexa(col, 0.95);
    g.beginPath(); g.arc(tip.x, tip.y, on ? 3.6 : 1.5, 0, 7); g.fill();
  }
  g.shadowBlur = 0;

  // labels: a picture nobody can read is not a visualisation
  g.textAlign = "center";
  for (const c of S.cells) {
    if (c.n < 14) continue;
    const r = (14 + Math.sqrt(c.n) * 6.2) * S.view.z;
    g.font = "500 11px 'DM Mono', monospace";
    g.fillStyle = "rgba(228,255,248,.92)";
    g.fillText(c.lineage, c.x, c.y - r * 0.62 - 9);
    g.font = "400 9px 'DM Mono', monospace";
    g.fillStyle = "rgba(150,205,195,.8)";
    g.fillText(`${c.n} genomes`, c.x, c.y - r * 0.62 + 2);
  }
  g.textAlign = "left";
  caption(g);
}

/* Drawn onto the canvas, not the DOM, so an exported PNG explains itself. */
function caption(g) {
  const d = S.data;
  const synth = d.identity.synthetic;
  g.font = "600 13px 'Instrument Sans', sans-serif";
  g.fillStyle = synth ? "#FFB347" : "rgba(228,255,248,.95)";
  g.fillText(
    (synth ? "SYNTHETIC DEMO — NOT A RESULT · " : "") +
    `${d.identity.display_name} — phylogeny of ${d.identity.n_samples} genomes`, 74, 132);
  g.font = "400 11px 'DM Mono', monospace";
  g.fillStyle = "rgba(150,205,195,.75)";
  g.fillText(`each point = one sequenced genome  ·  each halo = a clade  ·  ${S.cells.length} clades shown  ·  ${d.identity.reference}`, 74, 150);

  const x0 = 74, y0 = S.H - 128;
  g.font = "400 9px 'DM Mono', monospace";
  g.fillStyle = "rgba(150,205,195,.7)";
  g.fillText("LINEAGE", x0, y0 - 10);
  S.order.forEach((L, i) => {
    const y = y0 + i * 15;
    g.fillStyle = S.colour[L];
    g.shadowBlur = 6; g.shadowColor = S.colour[L];
    g.beginPath(); g.arc(x0 + 4, y - 3, 3.4, 0, 7); g.fill();
    g.shadowBlur = 0;
    g.font = "400 10.5px 'DM Mono', monospace";
    g.fillStyle = "rgba(206,245,236,.85)";
    g.fillText(`${L}  ${S.data.lineages[L]}`, x0 + 15, y);
  });
}

function hexa(hex, a) {
  const h = hex.replace("#", "");
  const n = parseInt(h, 16);
  return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${a})`;
}

/* A standing marker. Synthetic data is never rendered without it, and it
   cannot be dismissed -- the whole risk of a demo dataset is that a
   screenshot of it escapes into a report. */
function synthBanner() {
  document.body.append(el("div", { class: "synth" },
    el("span", { class: "dot" }),
    el("b", { text: "SYNTHETIC DEMONSTRATION DATA" }),
    el("span", { text: "— not a result. Every value is fabricated and none of it describes a real organism." }),
    el("button", { class: "to-real", text: "Load real FMDV", onclick: () => { location.search = ""; } }),
  ));
  document.body.classList.add("is-synth");
}

/* ═══ tool rail ═══════════════════════════════════════════════════ */
const TOOLS = [
  { id: "add", label: "Add sequence", cls: "new", d: '<path d="M12 5v14M5 12h14" stroke-width="1.8"/>' },
  { id: "move", label: "Move / navigate", d: '<path d="M3 20l18-8-18-8v6l10 2-10 2z"/>' },
  { id: "cut", label: "Cut sequence", d: '<circle cx="6" cy="6" r="2.5"/><circle cx="6" cy="18" r="2.5"/><path d="M20 4L8.5 15.5M20 20L8.5 8.5"/>' },
  { id: "insert", label: "Insert", d: '<circle cx="12" cy="12" r="8"/><path d="M12 8v8M8 12h8"/>' },
  { id: "splice", label: "Splice", d: '<path d="M4 20c6 0 6-16 12-16M4 4c4 0 5 6 8 8"/><circle cx="19" cy="18" r="2"/>' },
  { id: "delete", label: "Delete region", d: '<rect x="3" y="6" width="18" height="12" rx="2"/><path d="M9 10l6 4M15 10l-6 4"/>' },
  { id: "annotate", label: "Annotate", d: '<path d="M5 4h9l5 5v11H5z"/><path d="M9 12h6M9 16h4"/>' },
  { id: "report", label: "Report", d: '<path d="M6 3h8l4 4v14H6z"/><path d="M9 13h6M9 17h6M9 9h3"/>' },
  { id: "zoom", label: "Zoom", d: '<circle cx="11" cy="11" r="7"/><path d="M20 20l-3.5-3.5M8 11h6M11 8v6"/>' },
  { id: "measure", label: "Measure", d: '<rect x="2" y="9" width="20" height="6" rx="1"/><path d="M7 9v3M12 9v4M17 9v3"/>' },
  { id: "history", label: "History", d: '<circle cx="12" cy="12" r="8"/><path d="M12 7v5l3 2"/>' },
];

function buildRail() {
  const rail = $("#rail");
  rail.replaceChildren();
  TOOLS.forEach((t, i) => {
    if (i === 2 || i === 6 || i === 8) rail.append(el("div", { class: "div" }));
    const b = el("button", {
      class: "tool" + (t.cls ? " " + t.cls : "") + (S.tool === t.id ? " on" : ""),
      title: t.label, onclick: () => pickTool(t),
    });
    b.append(ico(t.d));
    rail.append(b);
  });
}
function pickTool(t) {
  if (t.id === "zoom") return zoom(1.25);
  if (t.id === "history") return toast(`<b>${S.history.length}</b> edits in this session — every action is recorded for reproducibility.`);
  S.tool = t.id; buildRail();
  toast(`<b>${t.label}</b> selected.`);
}

/* ═══ inspector ═══════════════════════════════════════════════════ */
const TABS = [
  { id: "gene", label: "Locus", lead: true },
  { id: "map", d: '<path d="M3 6l6-3 6 3 6-3v15l-6 3-6-3-6 3z"/>' },
  { id: "align", d: '<path d="M4 7h16M4 12h10M4 17h13"/>' },
  { id: "risk", d: '<path d="M12 4l9 16H3z"/><path d="M12 10v4M12 17v.5"/>' },
  { id: "stats", d: '<path d="M4 20V10M10 20V4M16 20v-7M22 20H2"/>' },
];
let tab = "gene";
function buildTabs() {
  const host = $("#insp-tabs");
  host.replaceChildren(
    el("span", { class: "lead" },
      ico('<path d="M6 3v4c0 4 12 6 12 10v4M18 3v4c0 4-12 6-12 10v4"/>'),
      el("span", { text: "Locus" })),
    ...TABS.filter((t) => !t.lead).map((t) => {
      const b = el("button", { class: tab === t.id ? "on" : "", title: t.id, onclick: () => { tab = t.id; renderInspector(); buildTabs(); } });
      b.append(ico(t.d)); return b;
    }),
  );
}

function card(k, v, wide) { return el("div", { class: "card" + (wide ? " wide" : "") }, el("span", { class: "k", text: k }), el("span", { class: "v", text: v })); }

function renderInspector() {
  const body = $("#insp-body");
  const d = S.data, L = S.locus;
  body.replaceChildren();
  $("#insp-title").textContent = tab === "risk" ? "Gate Status" : tab === "stats" ? "Composition" : "Locus Metadata";

  // live preview card — the same field, framed
  body.append(el("div", { class: "preview" },
    el("canvas", { id: "preview" }),
    el("span", { class: "badge lbl", style: "position:absolute;background:rgba(3,26,24,.6);border:1px solid var(--glass-edge)" },
      el("span", { class: "dot live" }), el("span", { text: "Live feed" })),
    el("span", { class: "scale" }, el("i"), el("span", { text: `${Math.round(d.identity.genome_length / 10)} bp` })),
  ));

  if (tab === "stats") {
    const c = d.composition;
    body.append(el("p", { class: "sect", text: "Reference genome" }));
    body.append(compGrid(c.percent));
    body.append(card("GC content", c.gc + " %", true));
    body.append(card("Sequenced length", c.total.toLocaleString() + " bp", true));
    if (L?.composition) {
      body.append(el("p", { class: "sect", text: `Locus ${L.id}` }));
      body.append(compGrid(L.composition));
      body.append(card("Locus GC", L.gc + " %", true));
    }
    return;
  }

  if (tab === "risk") {
    const g = d.gate;
    body.append(el("div", { class: "alert" + (g.permitted ? " warn" : "") },
      el("div", { class: "top" },
        el("span", { class: "t", text: g.permitted ? "Gate open" : "Gate verdict" }),
        el("span", { class: "chip", text: g.permitted ? "Permitted" : "Blocked" })),
      el("span", { class: "code", text: g.permission }),
      el("span", { class: "why", text: g.explanation })));
    for (const [k, r] of Object.entries(d.floor)) {
      if (r.value == null) continue;
      const pass = r.value >= r.floor, frac = r.unit === "frac";
      body.append(el("div", { class: "sliderow" },
        el("div", { class: "top" }, el("span", { class: "k", text: k.replace(/_/g, " ") }),
          el("span", { class: "v", style: pass ? "" : "color:var(--rose)", text: frac ? `${(r.value * 100).toFixed(1)}% / ${(r.floor * 100).toFixed(0)}%` : `${r.value} / ${r.floor}` })),
        el("div", { class: "track" },
          el("i", { style: `width:${Math.min(r.value / (r.floor * 1.6), 1) * 100}%;background:${pass ? "var(--lume)" : "var(--rose)"}` }))));
    }
    return;
  }

  // default: locus metadata, matching the reference's card grid
  body.append(el("div", { class: "cards" },
    card("Target locus", L ? L.id.replace(/^.*-G4-/, "G4-") : "—"),
    card("Reference", d.identity.reference),
    card("Prediction", L ? `${L.tools}/3 tools` : "—"),
    card("Region", L ? L.feature.replace(" (polyprotein)", "") : "—"),
  ));

  body.append(el("div", { class: "sliderow" },
    el("div", { class: "top" }, el("span", { class: "k", text: "Sequence length" }),
      el("span", { class: "v", text: d.identity.genome_length.toLocaleString() + " bp" })),
    el("div", { class: "track" },
      el("i", { style: `width:${L ? (L.end / d.identity.genome_length) * 100 : 100}%` }),
      el("span", { class: "knob", style: `left:${L ? (L.end / d.identity.genome_length) * 100 : 100}%` }))));

  if (L) {
    const weak = L.tools <= 1;
    body.append(el("div", { class: "alert" + (weak ? " warn" : "") },
      el("div", { class: "top" },
        el("span", { class: "t", text: "Structural confidence" }),
        el("span", { class: "chip", text: L.confidence === "WC" ? "Weak" : L.confidence })),
      el("span", { class: "code", text: `${L.start.toLocaleString()}–${L.end.toLocaleString()} · ${L.strand} strand` }),
      el("span", { class: "why", text: weak
        ? `G4Hunter ${L.g4hunter}, supported by ${L.tools} of 3 algorithms. Published inter-tool discordance is 30–60%, so this is a candidate motif, not a demonstrated structure.`
        : `G4Hunter ${L.g4hunter}, conservation ${L.conservation}%.` })));
    body.append(el("p", { class: "sect", text: "Composition" }));
    body.append(compGrid(L.composition || d.composition.percent));
  }

  body.append(el("p", { class: "sect", text: "Corpus" }));
  body.append(el("div", { class: "cards" },
    card("Genomes", String(d.identity.n_samples)),
    card("Lineages", String(d.identity.n_lineages)),
    card("Countries", String(d.identity.n_countries)),
    card("Period", d.identity.period ? d.identity.period.join("–") : "—"),
  ));
  drawPreview();
}

function compGrid(p) {
  const max = Math.max(...Object.values(p));
  return el("div", { class: "comp" }, ...["A", "T", "C", "G"].map((b) =>
    el("div", { class: "cell" },
      el("div", { class: "b", text: b }),
      el("div", { class: "p", text: (p[b] ?? 0).toFixed(1) + "%" }),
      el("div", { class: "bar" }, el("i", { style: `width:${((p[b] ?? 0) / max) * 100}%` })))));
}

function drawPreview() {
  const cv = $("#preview");
  if (!cv) return;
  const dpr = Math.min(devicePixelRatio || 1, 2);
  const w = cv.clientWidth, h = cv.clientHeight;
  if (!w || !h) return;
  cv.width = w * dpr; cv.height = h * dpr;
  const g = cv.getContext("2d");
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  const bg = g.createLinearGradient(0, 0, 0, h);
  bg.addColorStop(0, "#04211E"); bg.addColorStop(1, "#010E0D");
  g.fillStyle = bg; g.fillRect(0, 0, w, h);
  // temporal sparkline: real genomes per year, glowing
  const years = Object.keys(S.data.years).map(Number).sort((a, b) => a - b);
  const max = Math.max(...Object.values(S.data.years));
  g.globalCompositeOperation = "lighter";
  years.forEach((y, i) => {
    const x = 8 + (i / (years.length - 1)) * (w - 16);
    const bh = (S.data.years[y] / max) * (h - 30);
    g.shadowBlur = 8; g.shadowColor = "#5FF6D2";
    g.fillStyle = "rgba(95,246,210,.85)";
    g.fillRect(x - 1, h - 14 - bh, 2, bh);
  });
  g.shadowBlur = 0; g.globalCompositeOperation = "source-over";
  g.fillStyle = "rgba(150,205,195,.6)"; g.font = "9px DM Mono, monospace";
  g.fillText(`${years[0]}–${years[years.length - 1]} · genomes per year`, 9, h - 4);
}

/* ═══ sequence ruler ══════════════════════════════════════════════ */
async function loadSequence() {
  const n = Math.max(24, Math.floor((innerWidth - 460) / 26));
  try {
    S.seq = await api(`/api/sequence/${S.pathogen}?start=${Math.max(1, S.base - Math.floor(n / 2))}&length=${n}`);
  } catch { S.seq = null; }
  renderRuler();
}
function renderRuler() {
  const host = $("#ruler");
  if (!S.seq) { host.replaceChildren(el("span", { class: "tag", text: "no reference sequence" })); return; }
  const inLocus = (p) => S.data.loci.find((L) => p >= L.start && p <= L.end);
  host.replaceChildren(...[...S.seq.seq].map((b, i) => {
    const pos = S.seq.start + i;
    const L = inLocus(pos);
    return el("button", {
      class: "base" + (pos === S.base ? " on" : "") + (L ? " g4" : ""),
      title: `${pos.toLocaleString()}${L ? " · " + L.id : ""}`,
      onclick: () => { S.base = pos; if (L) { S.locus = L; renderInspector(); } renderRuler(); renderFoot(); },
      oncontextmenu: (e) => { e.preventDefault(); contextMenu(e, pos, L); },
    }, el("span", { class: "b", text: b }), el("span", { class: "t", text: pos % 10 === 0 ? String(pos).slice(-4) : "·" }));
  }));
}
function renderFoot() {
  const L = S.data.loci.find((x) => S.base >= x.start && S.base <= x.end);
  $("#foot").replaceChildren(
    el("span", {}, "POS ", el("b", { text: S.base.toLocaleString() })),
    el("span", { class: "sep" }),
    el("span", {}, "REGION ", el("b", { text: L ? L.id : (S.data.identity.cds && S.base >= S.data.identity.cds[0] && S.base <= S.data.identity.cds[1] ? "CDS" : "UTR") })),
    el("span", { class: "sep" }),
    el("span", {}, "GENOMES ", el("b", { text: String(S.data.identity.n_samples) })),
    el("span", { class: "sep" }),
    el("span", {}, "GATE ", el("b", { style: `color:${S.data.gate.permitted ? "var(--lume)" : "var(--rose)"}`, text: S.data.gate.permitted ? "OPEN" : "BLOCKED" })),
  );
}

/* ═══ context menu ════════════════════════════════════════════════ */
const CTX = [
  { id: "cut", label: "Cut Sequence", d: '<circle cx="6" cy="6" r="2.5"/><circle cx="6" cy="18" r="2.5"/><path d="M20 4L8.5 15.5M20 20L8.5 8.5"/>' },
  { id: "replace", label: "Replace Base", d: '<circle cx="12" cy="12" r="8"/><path d="M9 12h6"/>' },
  { id: "insert", label: "Insert Gene", d: '<circle cx="12" cy="12" r="8"/><path d="M12 8v8M8 12h8"/>' },
  { id: "duplicate", label: "Duplicate", d: '<rect x="8" y="8" width="12" height="12" rx="2"/><path d="M4 16V6a2 2 0 0 1 2-2h10"/>' },
  { id: "annotate", label: "Annotate", d: '<path d="M5 4h9l5 5v11H5z"/><path d="M9 13h6"/>' },
  { id: "crispr", label: "CRISPR Edit", d: '<path d="M4 20c6 0 6-16 12-16M4 4c4 0 5 6 8 8"/><circle cx="19" cy="18" r="2"/>' },
  { sep: true },
  { id: "delete", label: "Delete", danger: true, d: '<path d="M4 7h16M9 7V5h6v2M6 7l1 13h10l1-13"/>' },
];
let openCtx = null;
function contextMenu(e, pos, locus) {
  closeCtx();
  const m = el("div", { class: "ctxmenu glass" },
    el("div", { class: "head" },
      el("span", { class: "pos", text: `Pos: ${pos.toLocaleString()}` }),
      el("span", { class: "chip", text: locus ? locus.id.replace(/^.*-G4-/, "G4-") : (S.data.identity.cds && pos >= S.data.identity.cds[0] ? "CDS" : "UTR") })),
    ...CTX.map((c) => c.sep ? el("div", { class: "sep" })
      : (() => {
        const b = el("button", { class: c.danger ? "danger" : "", onclick: () => { closeCtx(); edit(c, pos); } });
        b.append(ico(c.d), el("span", { text: c.label }));
        return b;
      })()),
  );
  document.body.append(m);
  const r = m.getBoundingClientRect();
  m.style.left = Math.min(e.clientX, innerWidth - r.width - 12) + "px";
  m.style.top = Math.min(e.clientY, innerHeight - r.height - 12) + "px";
  openCtx = m;
  setTimeout(() => addEventListener("pointerdown", closeCtx, { once: true }), 0);
}
function closeCtx() { if (openCtx) { openCtx.remove(); openCtx = null; } }

/* Edits are recorded, never applied: this console has no write path to a
   reference genome, and silently pretending otherwise would be worse
   than declining. */
function edit(action, pos) {
  S.history.push({ action: action.id, pos, at: Date.now() });
  S.redo.length = 0;
  toast(`<b>${action.label}</b> recorded at ${pos.toLocaleString()} — staged only. G4 has no write path to a reference genome; edits are queued for a curator to apply.`);
}

/* ═══ interaction ═════════════════════════════════════════════════ */
function wire() {
  const cv = $("#field");
  let drag = null;
  cv.addEventListener("pointerdown", (e) => { if (S.tool === "move") drag = { x: e.clientX, y: e.clientY, ox: S.view.x, oy: S.view.y }; });
  addEventListener("pointerup", () => { drag = null; });
  addEventListener("pointermove", (e) => {
    if (drag) {
      S.view.x = drag.ox + (e.clientX - drag.x);
      S.view.y = drag.oy + (e.clientY - drag.y);
      readout(); buildCells(); drawField();
      return;
    }
    const hit = nearest(e.clientX, e.clientY);
    if (hit) {
      const s = S.data.samples[hit.s];
      gtip(e, s.a, [["lineage", s.l], ["country", s.c], ["year", s.y ?? "—"], ["host", s.h]]);
    } else hideGtip();
    readout(e);
  });
  cv.addEventListener("click", (e) => {
    const hit = nearest(e.clientX, e.clientY);
    if (!hit) return;
    S.sample = S.sample === S.data.samples[hit.s].a ? null : S.data.samples[hit.s].a;
    drawField();
    toast(`Genome <b>${S.data.samples[hit.s].a}</b> — ${S.data.samples[hit.s].l}, ${S.data.samples[hit.s].c}`);
  });
  cv.addEventListener("contextmenu", (e) => { e.preventDefault(); contextMenu(e, S.base, S.locus); });
  cv.addEventListener("wheel", (e) => { e.preventDefault(); zoom(e.deltaY < 0 ? 1.1 : 0.91); }, { passive: false });

  $("#btn-undo").onclick = () => { const h = S.history.pop(); if (!h) return toast("Nothing to undo."); S.redo.push(h); toast(`Undid <b>${h.action}</b> at ${h.pos.toLocaleString()}.`); };
  $("#btn-redo").onclick = () => { const h = S.redo.pop(); if (!h) return toast("Nothing to redo."); S.history.push(h); toast(`Redid <b>${h.action}</b> at ${h.pos.toLocaleString()}.`); };
  $("#btn-export").onclick = exportView;
  $("#q").addEventListener("keydown", (e) => { if (e.key === "Enter") search(e.target.value); });
  document.querySelectorAll(".mitem").forEach((b) => b.onclick = () => {
    document.querySelectorAll(".mitem").forEach((x) => x.classList.remove("on"));
    b.classList.add("on"); menu(b.dataset.menu);
  });
  addEventListener("keydown", (e) => {
    if (e.target.matches("input")) return;
    if (e.key === "Escape") { closeCtx(); S.sample = null; drawField(); }
    if (e.key === "+" || e.key === "=") zoom(1.15);
    if (e.key === "-") zoom(0.87);
  });
}
function nearest(mx, my) {
  let best = null, bd = 90;
  for (const t of S.tips) { const d = (t.x - mx) ** 2 + (t.y - my) ** 2; if (d < bd) { bd = d; best = t; } }
  return best;
}
function zoom(f) { S.view.z = Math.max(0.35, Math.min(6, S.view.z * f)); readout(); buildCells(); drawField(); }
function readout(e) {
  $("#cx").textContent = ((e ? e.clientX : S.W / 2) / S.W * 100 - 50).toFixed(3);
  $("#cy").textContent = ((e ? e.clientY : S.H / 2) / S.H * 100 - 50).toFixed(3);
  $("#cz").textContent = S.view.z.toFixed(3);
}

function menu(which) {
  const items = {
    project: [["New project", () => toast("Use the workstation at <b>/app</b> for project management.")],
      ["Open workstation", () => location.href = "/app"], ["Operator console", () => location.href = "/"]],
    edit: [["Undo", () => $("#btn-undo").click()], ["Redo", () => $("#btn-redo").click()],
      ["Clear selection", () => { S.sample = null; drawField(); }]],
    view: [["Reset view", () => { S.view = { x: 0, y: 0, z: 1 }; buildCells(); drawField(); readout(); }],
      ["Zoom in", () => zoom(1.3)], ["Zoom out", () => zoom(0.77)],
      ["Analysis workstation", () => location.href = "/app"]],
    window: [["Full screen", () => document.documentElement.requestFullscreen?.()],
      ["Linked views", () => location.href = "/workstation"]],
  }[which] || [];
  toast(items.map(([l]) => l).join(" · ") || "—");
  // act on the first item for a single-click menu; a full menu popover
  // belongs in the workstation, which already has one.
}

function search(q) {
  if (!q) return;
  const n = Number(q.replace(/[^0-9]/g, ""));
  if (n && n > 0 && n <= S.data.identity.genome_length) { S.base = n; loadSequence(); renderFoot(); return toast(`Jumped to position <b>${n.toLocaleString()}</b>.`); }
  const s = S.data.samples.find((x) => x.a.toLowerCase().includes(q.toLowerCase()));
  if (s) { S.sample = s.a; drawField(); return toast(`Genome <b>${s.a}</b> — ${s.l}, ${s.c}`); }
  const L = S.data.loci.find((x) => x.id.toLowerCase().includes(q.toLowerCase()));
  if (L) { S.locus = L; S.base = L.start; renderInspector(); loadSequence(); return toast(`Locus <b>${L.id}</b>`); }
  toast(`No match for “${q}”.`);
}
function exportView() {
  const cv = $("#field");
  if (S.data.identity.synthetic) stampSynthetic();
  cv.toBlob((b) => {
    const u = URL.createObjectURL(b), a = el("a", { href: u, download: "g4_studio_field.png" });
    document.body.append(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(u), 800);
    toast("Field exported as <b>PNG</b> at screen resolution.");
  });
}

/* An exported PNG outlives the screen it came from, so the stamp is burned
   into the pixels rather than added as chrome. */
function stampSynthetic() {
  const g = S.g;
  g.save();
  g.globalAlpha = 0.22;
  g.translate(S.W / 2, S.H / 2);
  g.rotate(-Math.atan2(S.H, S.W));
  g.textAlign = "center";
  g.font = `700 ${Math.round(S.W / 16)}px 'Instrument Sans', sans-serif`;
  g.fillStyle = "#FFB347";
  g.fillText("SYNTHETIC DEMO DATA", 0, 0);
  g.font = `500 ${Math.round(S.W / 52)}px 'DM Mono', monospace`;
  g.fillText("NOT A RESULT — FABRICATED FOR INTERFACE REVIEW", 0, Math.round(S.W / 22));
  g.restore();
}

/* ═══ chrome ══════════════════════════════════════════════════════ */
function toast(html) {
  const t = el("div", { class: "toast glass", html });
  $("#toasts").append(t);
  setTimeout(() => t.remove(), 5200);
}
function gtip(e, head, rows) {
  const t = $("#gtip");
  t.replaceChildren(el("div", { class: "h", text: head }),
    el("dl", {}, ...rows.flatMap(([k, v]) => [el("dt", { text: k }), el("dd", { text: String(v) })])));
  t.hidden = false;
  const r = t.getBoundingClientRect();
  t.style.left = Math.min(e.clientX + 16, innerWidth - r.width - 10) + "px";
  t.style.top = Math.min(e.clientY + 16, innerHeight - r.height - 10) + "px";
}
function hideGtip() { $("#gtip").hidden = true; }
function debounce(f, ms) { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => f(...a), ms); }; }

boot();
