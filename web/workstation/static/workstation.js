/* ═══════════════════════════════════════════════════════════════════
   G4-WATCH workstation — linked-view engine.

   ONE DATASET, MULTIPLE CONNECTED PERSPECTIVES.

   There is a single `sel` object. Every view reads it and nothing else;
   every view writes to it through `brush()`. That is the whole
   architecture — no view knows about any other view, so adding a
   perspective costs one render function and no wiring.

   Selection renders as contrast (Law 2 of the design system): the
   selected set keeps full opacity, everything else drops to --veil. No
   view invents a highlight colour.
   ═══════════════════════════════════════════════════════════════════ */

const $ = (s, r = document) => r.querySelector(s);
const el = (t, a = {}, ...kids) => {
  const n = document.createElement(t);
  for (const [k, v] of Object.entries(a)) {
    if (k === "class") n.className = v;
    else if (k === "text") n.textContent = v;
    else if (k === "html") n.innerHTML = v;
    else if (k.startsWith("on")) n.addEventListener(k.slice(2), v);
    else if (v != null) n.setAttribute(k, v);
  }
  for (const c of kids) if (c) n.append(c);
  return n;
};
const NS = "http://www.w3.org/2000/svg";
const svg = (t, a = {}, ...kids) => {
  const n = document.createElementNS(NS, t);
  for (const [k, v] of Object.entries(a)) if (v != null) n.setAttribute(k, v);
  for (const c of kids) if (c) n.append(c);
  return n;
};
const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const svgText = (attrs, text) => { const n = svg("text", attrs); n.textContent = text; return n; };

/* ── state ───────────────────────────────────────────────────────── */
const D = { data: null, lineageColor: {}, order: [] };

const sel = {
  lineages: new Set(),   // lineage labels
  countries: new Set(),  // country names
  period: null,          // [y0, y1]
  locus: null,           // atlas_id
  sample: null,          // accession
  clade: null,           // {root, n}
};

const HISTORY = [];
let VIEW = "phylogeny";
const OPT = { tipMode: "lineage", showInternal: false, logCounts: false, topCountries: 20 };

/* ── selection resolution ────────────────────────────────────────
   The single source of truth: which samples survive the current brush.
   Every view calls this and nothing else. */
function active() {
  const out = new Set();
  const { samples } = D.data;
  for (let i = 0; i < samples.length; i++) {
    const s = samples[i];
    if (sel.lineages.size && !sel.lineages.has(s.l)) continue;
    if (sel.countries.size && !sel.countries.has(s.c)) continue;
    if (sel.period && (!s.y || s.y < sel.period[0] || s.y > sel.period[1])) continue;
    if (sel.clade && !sel.clade.members.has(i)) continue;
    if (sel.sample && s.a !== sel.sample) continue;
    out.add(i);
  }
  return out;
}
function filtering() {
  return sel.lineages.size || sel.countries.size || sel.period || sel.clade || sel.sample;
}

function brush(mutate, label) {
  mutate();
  if (label) pushHistory(label);
  renderAll();
}

function pushHistory(label) {
  HISTORY.unshift({
    t: new Date(),
    label,
    view: VIEW,
    n: active().size,
    snapshot: JSON.stringify({
      lineages: [...sel.lineages], countries: [...sel.countries],
      period: sel.period, locus: sel.locus, sample: sel.sample,
      clade: sel.clade ? sel.clade.root : null,
    }),
  });
  if (HISTORY.length > 60) HISTORY.pop();
  renderHistory();
}

function restore(entry) {
  const s = JSON.parse(entry.snapshot);
  sel.lineages = new Set(s.lineages);
  sel.countries = new Set(s.countries);
  sel.period = s.period;
  sel.locus = s.locus;
  sel.sample = s.sample;
  sel.clade = s.clade != null ? cladeOf(s.clade) : null;
  renderAll();
}

/* ── boot ────────────────────────────────────────────────────────── */
async function boot() {
  const params = new URLSearchParams(location.search);
  const pathogen = params.get("p") || "fmdv";
  const res = await fetch(`/api/dataset/${pathogen}`);
  if (!res.ok) {
    $("#view").append(el("p", { class: "blank", text: "Dataset unavailable: " + (await res.text()) }));
    return;
  }
  D.data = await res.json();

  D.order = Object.keys(D.data.lineages);
  D.order.forEach((name, i) => { D.lineageColor[name] = css(`--cat-${i % 8}`); });

  renderSpecimen();
  renderSpine();
  renderObservations();
  renderRails();
  renderAll();
  pushHistory("Session opened");

  addEventListener("resize", debounce(() => { renderView(); renderRails(); }, 120));
  $("#btn-params").onclick = () => openConsole("params");
  $("#btn-history").onclick = () => openConsole("history");
  $("#btn-export").onclick = () => openConsole("export");
}

function renderAll() { renderView(); renderRegister(); renderRails(); renderSpine(); }

/* ═══ 1 · SPECIMEN ═══════════════════════════════════════════════ */
function renderSpecimen() {
  const id = D.data.identity;
  $("#rev").textContent = `workstation · atlas v${id.atlas_version}`;

  const fields = [
    ["organism", id.display_name],
    ["reference", id.reference],
    ["genome", `${id.genome_length.toLocaleString()} nt`],
    ["genomes", `${id.n_samples} <em>of ${id.n_raw}</em>`],
    [id.lineage_field + "s", String(id.n_lineages)],
    ["countries", String(id.n_countries)],
    ["period", id.period ? `${id.period[0]}–${id.period[1]}` : "—"],
    ["loci", String(D.data.loci.length)],
  ];
  const host = $("#spec-fields");
  host.replaceChildren();
  for (const [k, v] of fields) {
    host.append(el("div", { class: "spec-field" },
      el("span", { class: "tag", text: k }),
      el("span", { class: "v", html: v }),
    ));
  }

  const gate = D.data.gate;
  const node = $("#spec-gate");
  node.className = "spec-gate " + (gate.permitted ? "open" : "blocked");
  node.replaceChildren(
    el("span", { class: "lamp" }),
    el("span", { class: "txt", text: gate.permitted ? "scoring permitted" : "scoring blocked" }),
    el("span", { class: "tag", text: gate.permission }),
  );
}

/* ═══ 2 · WORKFLOW SPINE ═════════════════════════════════════════
   The pipeline's real shape: phylogeny and variants fork from QC and
   rejoin at the D.H1 gate. Navigation shaped like the process. */
const STAGES = [
  { id: "data",      ord: "0",   name: "Data",      view: "matrix",     col: 0 },
  { id: "qc",        ord: "1",   name: "QC",        view: "temporal",   col: 0 },
  { id: "phylo",     ord: "2",   name: "Phylogeny", view: "phylogeny",  col: -1 },
  { id: "variants",  ord: "3",   name: "Variants",  view: "genome",     col: 1 },
  { id: "gate",      ord: "4",   name: "D.H1",      view: "gate",       col: 0 },
  { id: "model",     ord: "5",   name: "Model",     view: "model",      col: 0 },
  { id: "report",    ord: "6",   name: "Report",    view: "report",     col: 0 },
];

function stageState(stage) {
  const gate = D.data.gate;
  if (stage.id === "gate") return gate.permitted ? "done" : "blocked";
  if (stage.id === "model" || stage.id === "report") return gate.permitted ? "idle" : "blocked";
  if (stage.id === "phylo") return D.data.tree ? "done" : "idle";
  if (stage.id === "qc") return "warn";
  return "done";
}

function renderSpine() {
  const host = $("#spine-nodes");
  host.replaceChildren();
  const marks = [];

  for (const stage of STAGES) {
    const state = stageState(stage);
    const on = stage.view === VIEW;
    const node = el("button", {
      class: `node ${state}${on ? " on" : ""}`,
      title: `${stage.name} — ${state}`,
      onclick: () => { VIEW = stage.view; renderView(); renderSpine(); },
    },
      el("span", { class: "ord", text: stage.ord }),
      el("span", { class: "glyph" }),
      el("span", { class: "name", text: stage.name }),
    );
    node.style.transform = `translateX(${stage.col * 15}px)`;
    host.append(node);
    marks.push({ node, stage });
  }

  requestAnimationFrame(() => {
    const box = $("#spine").getBoundingClientRect();
    const line = $("#spine-svg");
    line.replaceChildren();
    line.setAttribute("height", $("#spine").scrollHeight);
    const pt = marks.map(({ node, stage }) => {
      const r = node.getBoundingClientRect();
      return { x: r.left - box.left + r.width / 2 + stage.col * 0, y: r.top - box.top + 22, stage };
    });
    for (let i = 0; i < pt.length - 1; i++) {
      const a = pt[i], b = pt[i + 1];
      // Fork geometry: right-angle branches, as a cladogram draws them.
      const mid = (a.y + b.y) / 2;
      line.append(svg("path", {
        d: `M${a.x} ${a.y + 12} V${mid} H${b.x} V${b.y - 12}`,
        fill: "none", stroke: css("--edge"), "stroke-width": 1,
      }));
    }
  });
}

/* ═══ 3 · CANVAS VIEWS ═══════════════════════════════════════════ */
const VIEWS = {
  phylogeny: { title: "Phylogeny", render: viewPhylogeny },
  genome:    { title: "Genome architecture", render: viewGenome },
  temporal:  { title: "Temporal distribution", render: viewTemporal },
  matrix:    { title: "Lineage × geography", render: viewMatrix },
  gate:      { title: "D.H1 authorisation", render: viewGate },
  model:     { title: "Surveillance model", render: viewModel },
  report:    { title: "Report", render: viewReport },
};

function renderView() {
  const spec = VIEWS[VIEW] || VIEWS.phylogeny;
  $("#view-title").textContent = spec.title;
  const host = $("#view");
  host.replaceChildren();
  host.className = "view";
  spec.render(host);
}

/* ── 3a · phylogeny (canvas: 1,696 nodes) ───────────────────────── */
function cladeOf(rootIdx) {
  const nodes = D.data.tree.nodes;
  const kids = new Map();
  for (const n of nodes) if (n.p != null) (kids.get(n.p) || kids.set(n.p, []).get(n.p)).push(n.i);
  const members = new Set();
  const stack = [rootIdx];
  while (stack.length) {
    const i = stack.pop();
    const n = nodes[i];
    if (n.s >= 0) members.add(n.s);
    for (const c of kids.get(i) || []) stack.push(c);
  }
  return { root: rootIdx, members, n: members.size };
}

function viewPhylogeny(host) {
  const tree = D.data.tree;
  if (!tree) {
    host.append(el("p", { class: "blank", text: "No rooted tree artifact for this pathogen. Run Stage 2, or supply --rooted_tree." }));
    return;
  }
  $("#view-sub").textContent = `${tree.n_leaves} tips · IQ-TREE rooted · x = divergence from root`;

  const cv = el("canvas");
  host.append(cv);
  const dpr = devicePixelRatio || 1;
  const w = host.clientWidth, h = host.clientHeight;
  cv.width = w * dpr; cv.height = h * dpr;
  const g = cv.getContext("2d");
  g.scale(dpr, dpr);

  const M = { l: 22, r: 130, t: 14, b: 14 };
  const X = (x) => M.l + x * (w - M.l - M.r);
  const Y = (y) => M.t + y * (h - M.t - M.b);

  const chosen = active();
  const dim = filtering();
  const veil = parseFloat(css("--veil"));
  const nodes = tree.nodes;

  const kids = new Map();
  for (const n of nodes) if (n.p != null) { if (!kids.has(n.p)) kids.set(n.p, []); kids.get(n.p).push(n); }

  g.lineWidth = 1;
  // Branches, drawn as a rectangular cladogram.
  for (const n of nodes) {
    if (n.p == null) continue;
    const p = nodes[n.p];
    const lit = n.s < 0 ? true : chosen.has(n.s);
    g.globalAlpha = dim && !lit ? veil : (n.s < 0 ? 0.55 : 0.9);
    g.strokeStyle = n.s >= 0 ? (D.lineageColor[D.data.samples[n.s].l] || css("--ink-3")) : css("--ink-4");
    g.beginPath();
    g.moveTo(X(p.x), Y(n.y));
    g.lineTo(X(n.x), Y(n.y));
    g.stroke();
  }
  for (const [pi, cs] of kids) {
    const p = nodes[pi];
    const ys = cs.map((c) => Y(c.y));
    g.globalAlpha = dim ? 0.4 : 0.55;
    g.strokeStyle = css("--ink-4");
    g.beginPath();
    g.moveTo(X(p.x), Math.min(...ys));
    g.lineTo(X(p.x), Math.max(...ys));
    g.stroke();
  }

  // Tip marks, coloured by lineage.
  for (const n of nodes) {
    if (n.s < 0) continue;
    const lit = chosen.has(n.s);
    g.globalAlpha = dim && !lit ? veil : 1;
    g.fillStyle = D.lineageColor[D.data.samples[n.s].l] || css("--ink-3");
    g.fillRect(X(n.x), Y(n.y) - 1, 3, 2);
  }
  g.globalAlpha = 1;

  // Lineage bar: the phylogeny's own legend, in place.
  const barX = w - M.r + 16;
  for (const n of nodes) {
    if (n.s < 0) continue;
    const lit = chosen.has(n.s);
    g.globalAlpha = dim && !lit ? veil : 1;
    g.fillStyle = D.lineageColor[D.data.samples[n.s].l];
    g.fillRect(barX, Y(n.y) - 1, 10, 2);
  }
  g.globalAlpha = 1;
  g.font = `9px ${css("--mono") || "monospace"}`;
  g.fillStyle = css("--ink-3");
  g.fillText("lineage", barX, M.t - 3);

  host.append(legendStrip());

  cv.style.cursor = "crosshair";
  cv.onmousemove = (ev) => {
    const r = cv.getBoundingClientRect();
    const hit = nearestTip(nodes, ev.clientX - r.left, ev.clientY - r.top, X, Y);
    if (!hit) return hideTip();
    const s = D.data.samples[hit.s];
    showTip(ev, [
      ["accession", s.a], ["lineage", s.l], ["country", s.c],
      ["year", s.y ?? "—"], ["host", s.h],
    ]);
  };
  cv.onmouseleave = hideTip;
  cv.onclick = (ev) => {
    const r = cv.getBoundingClientRect();
    const hit = nearestTip(nodes, ev.clientX - r.left, ev.clientY - r.top, X, Y);
    if (!hit) return;
    const s = D.data.samples[hit.s];
    brush(() => {
      if (ev.shiftKey && hit.p != null) { sel.clade = cladeOf(hit.p); sel.sample = null; }
      else { sel.sample = sel.sample === s.a ? null : s.a; sel.clade = null; }
    }, ev.shiftKey ? `Clade of ${s.a} (${sel.clade ? sel.clade.n : 0} tips)` : `Genome ${s.a}`);
  };
}

function nearestTip(nodes, mx, my, X, Y) {
  let best = null, bd = 64;
  for (const n of nodes) {
    if (n.s < 0) continue;
    const dx = X(n.x) - mx, dy = Y(n.y) - my;
    const d = dx * dx + dy * dy;
    if (d < bd) { bd = d; best = n; }
  }
  return best;
}

/* ── 3b · genome architecture ───────────────────────────────────── */
function viewGenome(host) {
  const id = D.data.identity, loci = D.data.loci;
  $("#view-sub").textContent = `${id.reference} · ${id.genome_length.toLocaleString()} nt · ${loci.length} candidate loci`;

  const w = host.clientWidth, h = host.clientHeight;
  const M = { l: 56, r: 56, t: 40, b: 46 };
  const X = (bp) => M.l + (bp / id.genome_length) * (w - M.l - M.r);
  const axis = M.t + (h - M.t - M.b) * 0.52;

  const root = svg("svg", { viewBox: `0 0 ${w} ${h}`, width: w, height: h });
  const maxScore = Math.max(...loci.map((l) => Math.abs(l.g4hunter || 0)), 1);
  const amp = (h - M.t - M.b) * 0.38;

  // genome band with CDS / UTR partition
  root.append(svg("rect", { x: M.l, y: axis - 13, width: w - M.l - M.r, height: 26, fill: css("--plate-2") }));
  if (id.cds) {
    root.append(svg("rect", { x: M.l, y: axis - 13, width: X(id.cds[0]) - M.l, height: 26, fill: css("--plate-3") }));
    root.append(svg("rect", { x: X(id.cds[1]), y: axis - 13, width: X(id.genome_length) - X(id.cds[1]), height: 26, fill: css("--plate-3") }));
    root.append(svgText({ x: (M.l + X(id.cds[0])) / 2, y: axis + 4, "text-anchor": "middle", fill: css("--ink-3"), "font-family": "DM Mono, monospace", "font-size": 9 }, "5′UTR"));
    root.append(svgText({ x: (X(id.cds[0]) + X(id.cds[1])) / 2, y: axis + 4, "text-anchor": "middle", fill: css("--ink-3"), "font-family": "DM Mono, monospace", "font-size": 10 },
      `polyprotein CDS ${id.cds[0].toLocaleString()}–${id.cds[1].toLocaleString()}`));
  }
  root.append(svg("rect", { x: M.l, y: axis - 13, width: w - M.l - M.r, height: 26, fill: "none", stroke: css("--edge") }));

  for (const locus of loci) {
    const cx = X((locus.start + locus.end) / 2);
    const up = locus.strand === "+";
    const len = (Math.abs(locus.g4hunter || 0) / maxScore) * amp;
    const on = !sel.locus || sel.locus === locus.id;
    const gsel = svg("g", { opacity: on ? 1 : css("--veil"), style: "cursor:pointer" });

    gsel.append(svg("rect", {
      x: cx - 2, y: up ? axis - 13 - len : axis + 13,
      width: 4, height: len, fill: css("--ink"),
    }));
    gsel.append(svg("line", {
      x1: cx, x2: cx,
      y1: up ? axis - 13 - len : axis + 13 + len,
      y2: up ? M.t + 16 : h - M.b - 18,
      stroke: css("--edge"), "stroke-dasharray": "1 3",
    }));
    const label = svg("text", {
      x: cx, y: up ? M.t + 10 : h - M.b - 6, "text-anchor": "middle",
      fill: sel.locus === locus.id ? css("--ink") : css("--ink-2"),
      "font-family": "DM Mono, monospace", "font-size": 10,
    });
    label.textContent = locus.id.replace(/^.*-/, "G4-");
    gsel.append(label);
    const sub = svg("text", {
      x: cx, y: up ? M.t + 22 : h - M.b + 6, "text-anchor": "middle",
      fill: css("--ink-3"), "font-family": "DM Mono, monospace", "font-size": 9,
    });
    sub.textContent = `${locus.start}–${locus.end} · ${locus.g4hunter > 0 ? "+" : ""}${locus.g4hunter}`;
    gsel.append(sub);

    gsel.addEventListener("click", () => brush(() => {
      sel.locus = sel.locus === locus.id ? null : locus.id;
    }, `Locus ${locus.id}`));
    gsel.addEventListener("mousemove", (ev) => showTip(ev, [
      ["locus", locus.id], ["position", `${locus.start}–${locus.end} (${locus.end - locus.start + 1} nt)`],
      ["strand", locus.strand], ["feature", locus.feature],
      ["G4Hunter", locus.g4hunter], ["tools", `${locus.tools} of 3`],
      ["conservation", locus.conservation != null ? locus.conservation + "%" : "—"],
      ["GC flank", locus.gc_flank != null ? locus.gc_flank + "%" : "—"],
      ["confidence", locus.confidence], ["context", locus.context],
    ]));
    gsel.addEventListener("mouseleave", hideTip);
    root.append(gsel);
  }

  // strand annotation + coordinate ticks
  for (const [y, txt] of [[axis - amp - 6, "+ strand"], [axis + amp + 14, "− strand"]]) {
    const t = svg("text", { x: 8, y, fill: css("--ink-4"), "font-family": "DM Mono, monospace", "font-size": 9 });
    t.textContent = txt; root.append(t);
  }
  const step = 1000;
  for (let bp = 0; bp <= id.genome_length; bp += step) {
    root.append(svg("line", { x1: X(bp), x2: X(bp), y1: h - M.b + 14, y2: h - M.b + 18, stroke: css("--ink-4") }));
    const t = svg("text", { x: X(bp), y: h - M.b + 30, "text-anchor": "middle", fill: css("--ink-4"), "font-family": "DM Mono, monospace", "font-size": 9 });
    t.textContent = bp === 0 ? "0" : (bp / 1000) + "k";
    root.append(t);
  }
  root.append(svg("line", { x1: M.l, x2: X(id.genome_length), y1: h - M.b + 14, y2: h - M.b + 14, stroke: css("--edge") }));

  host.append(root);
}

/* ── 3c · temporal ──────────────────────────────────────────────── */
function viewTemporal(host) {
  const years = Object.keys(D.data.years).map(Number).sort((a, b) => a - b);
  $("#view-sub").textContent = `${years.length} sampled years · stacked by ${D.data.identity.lineage_field}`;

  const chosen = active();
  const dim = filtering();
  const counts = new Map();  // year -> lineage -> [total, selected]
  for (let i = 0; i < D.data.samples.length; i++) {
    const s = D.data.samples[i];
    if (!s.y) continue;
    if (!counts.has(s.y)) counts.set(s.y, new Map());
    const row = counts.get(s.y);
    const cur = row.get(s.l) || [0, 0];
    cur[0] += 1;
    if (chosen.has(i)) cur[1] += 1;
    row.set(s.l, cur);
  }

  const w = host.clientWidth, h = host.clientHeight;
  const M = { l: 46, r: 18, t: 18, b: 34 };
  const maxN = Math.max(...[...counts.values()].map((r) => [...r.values()].reduce((a, b) => a + b[0], 0)));
  const scale = (n) => OPT.logCounts ? Math.log10(n + 1) / Math.log10(maxN + 1) : n / maxN;
  const bw = Math.max(2, (w - M.l - M.r) / years.length - 2);
  const X = (i) => M.l + i * ((w - M.l - M.r) / years.length);
  const H = h - M.t - M.b;

  const root = svg("svg", { viewBox: `0 0 ${w} ${h}`, width: w, height: h });
  const veil = css("--veil");

  years.forEach((year, i) => {
    let acc = 0;
    const row = counts.get(year) || new Map();
    for (const lineage of D.order) {
      const cell = row.get(lineage);
      if (!cell) continue;
      const seg = scale(acc + cell[0]) - scale(acc);
      const y = M.t + H - scale(acc + cell[0]) * H;
      const rect = svg("rect", {
        x: X(i), y, width: bw, height: Math.max(seg * H, 0.6),
        fill: D.lineageColor[lineage],
        opacity: dim && cell[1] === 0 ? veil : 1,
        style: "cursor:pointer",
      });
      rect.addEventListener("mousemove", (ev) => showTip(ev, [
        ["year", year], ["lineage", lineage], ["genomes", cell[0]],
        ["in selection", cell[1]],
      ]));
      rect.addEventListener("mouseleave", hideTip);
      rect.addEventListener("click", () => brush(() => {
        sel.period = [year, year];
        sel.lineages = new Set([lineage]);
      }, `${lineage} · ${year}`));
      root.append(rect);
      acc += cell[0];
    }
  });

  // axes
  root.append(svg("line", { x1: M.l, x2: w - M.r, y1: M.t + H, y2: M.t + H, stroke: css("--edge") }));
  years.forEach((year, i) => {
    if (year % 10 !== 0) return;
    const t = svg("text", { x: X(i) + bw / 2, y: h - 12, "text-anchor": "middle", fill: css("--ink-4"), "font-family": "DM Mono, monospace", "font-size": 9 });
    t.textContent = year; root.append(t);
  });
  for (const frac of [0, 0.5, 1]) {
    const y = M.t + H - frac * H;
    root.append(svg("line", { x1: M.l, x2: w - M.r, y1: y, y2: y, stroke: css("--hair") }));
    const t = svg("text", { x: M.l - 8, y: y + 3, "text-anchor": "end", fill: css("--ink-4"), "font-family": "DM Mono, monospace", "font-size": 9 });
    t.textContent = OPT.logCounts ? "" : Math.round(frac * maxN);
    root.append(t);
  }
  host.append(root, legendStrip());
}

/* ── 3d · lineage × geography matrix ────────────────────────────── */
function viewMatrix(host) {
  const chosen = active();
  const dim = filtering();
  const countries = Object.entries(D.data.countries)
    .filter(([c]) => c !== "—").slice(0, OPT.topCountries).map(([c]) => c);
  $("#view-sub").textContent = `${D.order.length} lineages × top ${countries.length} countries · cell = genome count`;

  const grid = new Map();
  let maxCell = 0;
  for (let i = 0; i < D.data.samples.length; i++) {
    const s = D.data.samples[i];
    const key = s.l + " " + s.c;
    const cur = grid.get(key) || [0, 0];
    cur[0] += 1;
    if (chosen.has(i)) cur[1] += 1;
    grid.set(key, cur);
    if (cur[0] > maxCell) maxCell = cur[0];
  }

  host.className = "view scroll";
  const cell = 26, labelW = 108, headH = 74;
  const w = labelW + countries.length * cell + 20;
  const h = headH + D.order.length * cell + 20;
  const root = svg("svg", { viewBox: `0 0 ${w} ${h}`, width: w, height: h, style: `min-width:${w}px` });

  countries.forEach((country, cIdx) => {
    const x = labelW + cIdx * cell;
    const t = svg("text", {
      x: x + cell / 2, y: headH - 8, fill: css("--ink-3"),
      "font-family": "DM Mono, monospace", "font-size": 9.5,
      transform: `rotate(-58 ${x + cell / 2} ${headH - 8})`, "text-anchor": "start",
      style: "cursor:pointer",
    });
    t.textContent = country;
    t.addEventListener("click", () => brush(() => {
      sel.countries.has(country) ? sel.countries.delete(country) : sel.countries.add(country);
    }, `Country ${country}`));
    root.append(t);
  });

  D.order.forEach((lineage, rIdx) => {
    const y = headH + rIdx * cell;
    const label = svg("text", { x: labelW - 10, y: y + cell / 2 + 3, "text-anchor": "end", fill: css("--ink-2"), "font-family": "DM Mono, monospace", "font-size": 11, style: "cursor:pointer" });
    label.textContent = lineage;
    label.addEventListener("click", () => brush(() => {
      sel.lineages.has(lineage) ? sel.lineages.delete(lineage) : sel.lineages.add(lineage);
    }, `Lineage ${lineage}`));
    root.append(label);
    root.append(svg("rect", { x: labelW - 6, y: y + cell / 2 - 4, width: 8, height: 8, fill: D.lineageColor[lineage] }));

    countries.forEach((country, cIdx) => {
      const x = labelW + cIdx * cell;
      const v = grid.get(lineage + " " + country);
      root.append(svg("rect", { x: x + 1, y: y + 1, width: cell - 2, height: cell - 2, fill: css("--plate-2") }));
      if (!v) return;
      const k = Math.sqrt(v[0] / maxCell);
      const r = svg("rect", {
        x: x + 1, y: y + 1, width: cell - 2, height: cell - 2,
        fill: D.lineageColor[lineage], opacity: (dim && v[1] === 0 ? parseFloat(css("--veil")) : 1) * (0.2 + 0.8 * k),
        style: "cursor:pointer",
      });
      r.addEventListener("mousemove", (ev) => showTip(ev, [["lineage", lineage], ["country", country], ["genomes", v[0]], ["in selection", v[1]]]));
      r.addEventListener("mouseleave", hideTip);
      r.addEventListener("click", () => brush(() => {
        sel.lineages = new Set([lineage]); sel.countries = new Set([country]);
      }, `${lineage} × ${country}`));
      root.append(r);
      if (v[0] >= 10) {
        const t = svg("text", { x: x + cell / 2, y: y + cell / 2 + 3, "text-anchor": "middle", fill: css("--ground"), "font-family": "DM Mono, monospace", "font-size": 9, "pointer-events": "none" });
        t.textContent = v[0]; root.append(t);
      }
    });
  });
  host.append(root);
}

/* ── 3e · gate ──────────────────────────────────────────────────── */
function viewGate(host) {
  const gate = D.data.gate, floor = D.data.floor;
  $("#view-sub").textContent = `ledger ${gate.n_ledger_rows} rows · last run ${gate.latest_run ? gate.latest_run.slice(0, 10) : "—"}`;
  host.className = "view scroll";

  const wrap = el("div", { style: "padding:28px 32px; max-width:920px; display:flex; flex-direction:column; gap:26px" });

  wrap.append(el("div", { style: `border:1px solid var(--${gate.permitted ? "ok" : "block"}); padding:18px 20px` },
    el("p", { class: "tag", style: `color:var(--${gate.permitted ? "ok" : "block"})`, text: gate.permission }),
    el("p", { style: `font-size:22px;font-weight:600;margin-top:6px;color:var(--${gate.permitted ? "ok" : "block"})`, text: gate.permitted ? "Scoring permitted" : "Scoring blocked" }),
    el("p", { style: "font-size:12.5px;color:var(--ink-2);line-height:1.6;margin-top:10px;max-width:74ch", text: gate.explanation }),
  ));

  const table = el("div", { style: "display:flex;flex-direction:column" });
  table.append(el("p", { class: "tag", style: "margin-bottom:10px", text: "Appendix C minimum-data floor" }));
  for (const [key, row] of Object.entries(floor)) {
    if (row.value == null) continue;
    const frac = row.unit === "frac";
    const pass = row.value >= row.floor;
    const ratio = Math.min(row.value / (row.floor * 1.6), 1);
    const floorPos = 1 / 1.6;
    table.append(el("div", { style: "display:grid;grid-template-columns:206px 1fr 116px;gap:14px;align-items:center;padding:9px 0;border-bottom:1px solid var(--hair-2)" },
      el("span", { class: "num", style: `font-size:11.5px;color:var(--${pass ? "ink-2" : "block"})`, text: key }),
      el("div", { style: "position:relative;height:10px;background:var(--plate-2)" },
        el("i", { style: `position:absolute;left:0;top:0;bottom:0;width:${ratio * 100}%;background:var(--${pass ? "ok" : "block"})` }),
        el("i", { style: `position:absolute;left:${floorPos * 100}%;top:-3px;bottom:-3px;width:1px;background:var(--ink)` }),
      ),
      el("span", { class: "num", style: `font-size:11px;text-align:right;color:var(--${pass ? "ink-2" : "block"})`,
        text: frac ? `${(row.value * 100).toFixed(1)}% / ${(row.floor * 100).toFixed(0)}%` : `${row.value} / ${row.floor}` }),
    ));
  }
  if (floor.min_sequences_per_lineage.which) {
    table.append(el("p", { class: "hint", style: "margin-top:10px", text: `The binding constraint is lineage ${floor.min_sequences_per_lineage.which}, with ${floor.min_sequences_per_lineage.value} genomes. The vertical mark on each bar is the floor.` }));
  }
  wrap.append(table);
  host.append(wrap);
}

/* ── 3f · model / report (honestly gated) ───────────────────────── */
function viewModel(host) {
  $("#view-sub").textContent = "Stage 5";
  host.append(el("p", { class: "blank" },
    el("span", { style: "display:block;color:var(--block);font-family:var(--mono);font-size:10px;letter-spacing:.16em;text-transform:uppercase;margin-bottom:10px", text: "unwired by design" }),
    el("span", { text: "No surveillance model is fitted and no score is shown, because the D.H1 gate has not returned SUPPORTED. This panel stays empty rather than displaying a provisional number: a score rendered here would be read as a result." }),
  ));
}
function viewReport(host) {
  $("#view-sub").textContent = "Stage 6";
  host.append(el("p", { class: "blank" },
    el("span", { style: "display:block;color:var(--block);font-family:var(--mono);font-size:10px;letter-spacing:.16em;text-transform:uppercase;margin-bottom:10px", text: "scored report blocked" }),
    el("span", { text: "The scored surveillance report is gated behind the same authorisation. The unscored record — Atlas, corpus composition, gate verdict and ledger — can still be exported from the Export console." }),
  ));
}

/* ── legend ─────────────────────────────────────────────────────── */
function legendStrip() {
  const strip = el("div", { class: "legend", style: "position:absolute;left:20px;bottom:10px;background:color-mix(in srgb, var(--ground) 88%, transparent);padding:6px 10px" });
  for (const lineage of D.order) {
    const off = sel.lineages.size && !sel.lineages.has(lineage);
    strip.append(el("button", {
      class: off ? "off" : "",
      onclick: () => brush(() => {
        sel.lineages.has(lineage) ? sel.lineages.delete(lineage) : sel.lineages.add(lineage);
      }, `Lineage ${lineage}`),
    },
      el("i", { style: `background:${D.lineageColor[lineage]}` }),
      el("span", { text: `${lineage} ${D.data.lineages[lineage]}` }),
    ));
  }
  return strip;
}

/* ═══ selection register ═════════════════════════════════════════ */
function renderRegister() {
  const host = $("#register");
  host.replaceChildren();
  const n = active().size;

  const tokens = [];
  for (const lineage of sel.lineages) tokens.push(token(lineage, "lineage", D.lineageColor[lineage], () => sel.lineages.delete(lineage)));
  for (const country of sel.countries) tokens.push(token(country, "country", null, () => sel.countries.delete(country)));
  if (sel.period) tokens.push(token(sel.period[0] === sel.period[1] ? String(sel.period[0]) : sel.period.join("–"), "period", null, () => { sel.period = null; }));
  if (sel.clade) tokens.push(token(`node ${sel.clade.root} · ${sel.clade.n} tips`, "clade", null, () => { sel.clade = null; }));
  if (sel.sample) tokens.push(token(sel.sample, "genome", null, () => { sel.sample = null; }));
  if (sel.locus) tokens.push(token(sel.locus, "locus", null, () => { sel.locus = null; }));

  if (!tokens.length) {
    host.append(el("span", { class: "empty", text: "NO BRUSH — all 848 genomes active. Click any mark to filter every view. Shift-click a tip for its clade." }));
    return;
  }
  host.append(el("span", { class: "tag", text: "brush" }));
  tokens.forEach((t) => host.append(t));
  host.append(el("span", { class: "token", style: "border-color:var(--ink)" }, el("b", { text: `${n} / ${D.data.samples.length} genomes` })));
  host.append(el("button", { class: "k tiny", text: "clear", onclick: () => brush(() => {
    sel.lineages.clear(); sel.countries.clear();
    sel.period = null; sel.clade = null; sel.sample = null; sel.locus = null;
  }, "Brush cleared") }));
}

function token(label, kind, colour, remove) {
  return el("span", { class: "token" },
    colour ? el("i", { class: "swatch", style: `background:${colour}` }) : null,
    el("span", { style: "color:var(--ink-3)", text: kind }),
    el("b", { text: label }),
    el("button", { text: "×", title: "remove", onclick: () => brush(remove, null) }),
  );
}

/* ═══ observations ═══════════════════════════════════════════════ */
function renderObservations() {
  const host = $("#observations");
  const list = D.data.observations;
  $("#obs-count").textContent = `${list.filter((o) => o.severity === "block").length} blocking · ${list.length} total`;
  host.replaceChildren();
  for (const o of list) {
    host.append(el("button", {
      class: "obs " + o.severity,
      onclick: () => focusObservation(o),
    },
      el("span", { class: "sev" }),
      el("div", {},
        el("p", { class: "dom", text: `${o.severity} · ${o.domain}` }),
        el("p", { class: "t", text: o.title }),
        el("p", { class: "d", text: o.detail }),
        el("p", { class: "r", text: o.rule }),
      ),
    ));
  }
}

function focusObservation(o) {
  const f = o.focus || {};
  brush(() => {
    if (f.type === "lineage" && f.value) { sel.lineages = new Set([f.value]); VIEW = "temporal"; }
    else if (f.type === "locus") { sel.locus = f.value; VIEW = "genome"; }
    else if (f.type === "period") { sel.period = f.value; VIEW = "temporal"; }
    else if (f.type === "country") { sel.countries = new Set([f.value]); VIEW = "matrix"; }
    else if (f.type === "region") { VIEW = "genome"; }
    else if (f.type === "gate") { VIEW = "gate"; }
  }, "Observation: " + o.title);
}

/* ═══ coordinate rail ════════════════════════════════════════════ */
function renderRails() {
  railGenome();
  railTime();
}

function railGenome() {
  const node = $("#rail-genome");
  const id = D.data.identity;
  const w = node.clientWidth || 600, h = 62;
  node.setAttribute("viewBox", `0 0 ${w} ${h}`);
  node.replaceChildren();
  const M = 10;
  const X = (bp) => M + (bp / id.genome_length) * (w - 2 * M);

  node.append(svg("rect", { x: M, y: 26, width: w - 2 * M, height: 10, fill: css("--plate-2"), stroke: css("--hair") }));
  if (id.cds) {
    node.append(svg("rect", { x: X(id.cds[0]), y: 26, width: X(id.cds[1]) - X(id.cds[0]), height: 10, fill: css("--plate-3") }));
  }
  for (const locus of D.data.loci) {
    const cx = X((locus.start + locus.end) / 2);
    const on = !sel.locus || sel.locus === locus.id;
    const mark = svg("g", { opacity: on ? 1 : css("--veil"), style: "cursor:pointer" });
    mark.append(svg("rect", { x: cx - 1.5, y: locus.strand === "+" ? 14 : 36, width: 3, height: 12, fill: css("--ink") }));
    const t = svg("text", { x: cx, y: locus.strand === "+" ? 10 : 58, "text-anchor": "middle", fill: css("--ink-3"), "font-family": "DM Mono, monospace", "font-size": 8.5 });
    t.textContent = locus.id.slice(-3);
    mark.append(t);
    mark.addEventListener("click", () => brush(() => { sel.locus = sel.locus === locus.id ? null : locus.id; }, `Locus ${locus.id}`));
    node.append(mark);
  }
  $("#rail-genome-read").textContent = sel.locus
    ? sel.locus
    : `0 – ${id.genome_length.toLocaleString()} nt`;
}

function railTime() {
  const node = $("#rail-time");
  const years = Object.keys(D.data.years).map(Number).sort((a, b) => a - b);
  const w = node.clientWidth || 400, h = 62;
  node.setAttribute("viewBox", `0 0 ${w} ${h}`);
  node.replaceChildren();
  const M = 10;
  const y0 = years[0], y1 = years[years.length - 1];
  const X = (year) => M + ((year - y0) / Math.max(y1 - y0, 1)) * (w - 2 * M);
  const maxN = Math.max(...Object.values(D.data.years));

  for (const year of years) {
    const n = D.data.years[year];
    const bh = (n / maxN) * 34;
    const inSel = !sel.period || (year >= sel.period[0] && year <= sel.period[1]);
    node.append(svg("rect", {
      x: X(year) - 1.5, y: 42 - bh, width: 3, height: bh,
      fill: css("--ink-2"), opacity: inSel ? 1 : css("--veil"),
    }));
  }
  node.append(svg("line", { x1: M, x2: w - M, y1: 42, y2: 42, stroke: css("--edge") }));
  for (const year of [y0, Math.round((y0 + y1) / 2), y1]) {
    const t = svg("text", { x: X(year), y: 55, "text-anchor": "middle", fill: css("--ink-4"), "font-family": "DM Mono, monospace", "font-size": 9 });
    t.textContent = year; node.append(t);
  }
  if (sel.period) {
    node.append(svg("rect", {
      x: X(sel.period[0]) - 2, y: 6, width: Math.max(X(sel.period[1]) - X(sel.period[0]) + 4, 4), height: 36,
      fill: "none", stroke: css("--ink"), "stroke-dasharray": "2 2",
    }));
  }

  // drag to brush a period
  let anchor = null;
  const yearAt = (clientX) => {
    const r = node.getBoundingClientRect();
    const frac = (clientX - r.left - M) / (r.width - 2 * M);
    return Math.round(y0 + frac * (y1 - y0));
  };
  node.onmousedown = (ev) => { anchor = yearAt(ev.clientX); };
  node.onmousemove = (ev) => {
    const year = yearAt(ev.clientX);
    $("#rail-time-read").textContent = anchor == null
      ? `${year}`
      : `${Math.min(anchor, year)} – ${Math.max(anchor, year)}`;
  };
  node.onmouseup = (ev) => {
    if (anchor == null) return;
    const year = yearAt(ev.clientX);
    const lo = Math.max(Math.min(anchor, year), y0), hi = Math.min(Math.max(anchor, year), y1);
    anchor = null;
    brush(() => { sel.period = (lo === hi && sel.period && sel.period[0] === lo) ? null : [lo, hi]; }, `Period ${lo}–${hi}`);
  };
  node.onmouseleave = () => { anchor = null; };
  if (!sel.period) $("#rail-time-read").textContent = `${y0} – ${y1}`;
  else $("#rail-time-read").textContent = `${sel.period[0]} – ${sel.period[1]}`;
}

/* ═══ floating consoles ══════════════════════════════════════════ */
const OPEN = {};

function openConsole(kind) {
  if (OPEN[kind]) { OPEN[kind].remove(); delete OPEN[kind]; return; }
  const body = el("div", { class: "console-body" });
  const panel = el("div", { class: "console" },
    el("div", { class: "console-bar" },
      el("h3", { text: { params: "Parameters", history: "Analysis history", export: "Export" }[kind] }),
      el("button", { class: "x", text: "×", onclick: () => { panel.remove(); delete OPEN[kind]; } }),
    ),
    body,
  );
  const offset = Object.keys(OPEN).length * 24;
  panel.style.right = (28 + offset) + "px";
  panel.style.top = (120 + offset) + "px";
  document.body.append(panel);
  OPEN[kind] = panel;
  dragify(panel, panel.firstChild);

  if (kind === "params") fillParams(body);
  if (kind === "history") { body.id = "hist-body"; renderHistory(); }
  if (kind === "export") fillExport(body);
}

function fillParams(body) {
  body.append(
    ctl("Tip colour", el("select", { onchange: (e) => { OPT.tipMode = e.target.value; renderView(); } },
      el("option", { value: "lineage", text: "lineage" }))),
    ctlSwitch("Logarithmic counts", OPT.logCounts, (v) => { OPT.logCounts = v; renderView(); }),
    ctlRange("Countries in matrix", OPT.topCountries, 5, 57, (v) => { OPT.topCountries = v; renderView(); }),
    el("p", { class: "hint", text: "Parameters affect presentation only. Nothing here changes a pipeline result — the Atlas, the tree and the gate verdict are read from artifacts and are not recomputed by this view." }),
  );
}

function ctl(label, control) {
  return el("div", { class: "ctl" }, el("label", { text: label }), control);
}
function ctlSwitch(label, value, onchange) {
  return el("label", { class: "switch" },
    el("input", { type: "checkbox", checked: value ? "checked" : null, onchange: (e) => onchange(e.target.checked) }),
    el("span", { text: label }));
}
function ctlRange(label, value, min, max, onchange) {
  const out = el("output", { text: String(value) });
  return el("div", { class: "ctl" },
    el("label", { text: label }),
    el("div", { class: "row" },
      el("input", { type: "range", min, max, value, oninput: (e) => { out.textContent = e.target.value; onchange(+e.target.value); } }),
      out));
}

function renderHistory() {
  const body = document.getElementById("hist-body");
  if (!body) return;
  body.replaceChildren();
  const list = el("div", { class: "hist" });
  for (const entry of HISTORY) {
    list.append(el("button", { onclick: () => restore(entry) },
      el("span", { class: "h1", text: entry.label }),
      el("span", { class: "h3", text: entry.t.toLocaleTimeString([], { hour12: false }) }),
      el("span", { class: "h2", text: `${entry.view} · ${entry.n} genomes` }),
    ));
  }
  body.append(list, el("p", { class: "hint", text: "Every brush is recorded. Selecting an entry restores that exact selection state across all views — the analysis is reproducible within the session." }));
}

function fillExport(body) {
  const doExport = (what) => {
    const chosen = [...active()].map((i) => D.data.samples[i]);
    if (what === "svg") {
      const node = $("#view svg");
      if (!node) return alertLine("This view renders to canvas; use PNG.");
      const blob = new Blob([new XMLSerializer().serializeToString(node)], { type: "image/svg+xml" });
      download(blob, `g4watch_${VIEW}.svg`);
    } else if (what === "png") {
      const node = $("#view canvas");
      if (!node) return alertLine("This view renders to SVG; use SVG.");
      node.toBlob((b) => download(b, `g4watch_${VIEW}.png`));
    } else if (what === "tsv") {
      const rows = ["accession\tlineage\tcountry\tyear\thost",
        ...chosen.map((s) => [s.a, s.l, s.c, s.y ?? "", s.h].join("\t"))];
      download(new Blob([rows.join("\n")], { type: "text/tab-separated-values" }), "g4watch_selection.tsv");
    } else if (what === "state") {
      download(new Blob([JSON.stringify({
        pathogen: D.data.identity.pathogen, view: VIEW,
        selection: { lineages: [...sel.lineages], countries: [...sel.countries], period: sel.period, locus: sel.locus, sample: sel.sample },
        n_selected: chosen.length, exported: new Date().toISOString(),
      }, null, 2)], { type: "application/json" }), "g4watch_state.json");
    }
  };
  body.append(
    el("p", { class: "hint", text: `${active().size} genomes are in the current brush.` }),
    el("div", { class: "ctl" }, el("label", { text: "Figure" }),
      el("div", { class: "row" },
        el("button", { class: "k", text: "SVG", onclick: () => doExport("svg") }),
        el("button", { class: "k", text: "PNG", onclick: () => doExport("png") }))),
    el("div", { class: "ctl" }, el("label", { text: "Data" }),
      el("div", { class: "row" },
        el("button", { class: "k", text: "Selection TSV", onclick: () => doExport("tsv") }),
        el("button", { class: "k", text: "State JSON", onclick: () => doExport("state") }))),
    el("p", { class: "hint", text: "Figures export at the vector resolution of the view. The state file records the brush so a colleague can reproduce the same selection." }),
  );
}

function download(blob, name) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = name;
  document.body.append(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
function alertLine(msg) {
  const n = el("div", { class: "tip", style: "left:50%;bottom:130px;transform:translateX(-50%)", text: msg });
  document.body.append(n); setTimeout(() => n.remove(), 2600);
}

function dragify(panel, handle) {
  let dx = 0, dy = 0;
  handle.addEventListener("mousedown", (ev) => {
    if (ev.target.closest("button")) return;
    const r = panel.getBoundingClientRect();
    dx = ev.clientX - r.left; dy = ev.clientY - r.top;
    panel.style.right = "auto";
    const move = (e) => {
      // snap to the 4-unit grid
      panel.style.left = Math.round((e.clientX - dx) / 4) * 4 + "px";
      panel.style.top = Math.round((e.clientY - dy) / 4) * 4 + "px";
    };
    const up = () => { removeEventListener("mousemove", move); removeEventListener("mouseup", up); };
    addEventListener("mousemove", move); addEventListener("mouseup", up);
  });
}

/* ═══ tooltip ════════════════════════════════════════════════════ */
function showTip(ev, rows) {
  const tip = $("#tip");
  tip.replaceChildren();
  for (const [k, v] of rows) {
    tip.append(el("div", {},
      el("span", { class: "k2", text: k + "  " }),
      el("strong", { text: String(v) })));
  }
  tip.hidden = false;
  const pad = 14;
  const r = tip.getBoundingClientRect();
  tip.style.left = Math.min(ev.clientX + pad, innerWidth - r.width - 8) + "px";
  tip.style.top = Math.min(ev.clientY + pad, innerHeight - r.height - 8) + "px";
}
function hideTip() { $("#tip").hidden = true; }

function debounce(fn, ms) {
  let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); };
}

boot();
