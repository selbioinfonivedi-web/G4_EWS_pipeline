/* ═══════════════════════════════════════════════════════════════════
   G4 — computational genomics workstation.

   Seven modes over one dataset. The shell never moves: identity, the
   pipeline spine and the execution console stay put while the workspace
   changes, so a researcher never loses the thread of what is running.

   Everything that reports a number here reads it from the backend. Where
   a capability is not implemented, the UI says so in place rather than
   showing a plausible figure.
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
const NS = "http://www.w3.org/2000/svg";
const svg = (t, a = {}, ...kids) => {
  const n = document.createElementNS(NS, t);
  for (const [k, v] of Object.entries(a)) if (v != null) n.setAttribute(k, v);
  for (const c of kids.flat()) if (c) n.append(c);
  return n;
};
const svgText = (a, t) => { const n = svg("text", a); n.textContent = t; return n; };
const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const bytes = (b) => b < 1024 ? b + " B" : b < 1048576 ? (b / 1024).toFixed(0) + " KB"
  : b < 1073741824 ? (b / 1048576).toFixed(1) + " MB" : (b / 1073741824).toFixed(2) + " GB";
const clock = (s) => s == null ? "—" : s < 60 ? s.toFixed(1) + "s"
  : Math.floor(s / 60) + "m " + Math.round(s % 60) + "s";

/* ── application state ───────────────────────────────────────────── */
const S = {
  env: null, commands: [], pathogens: [], pathogen: null,
  data: null, inputs: [], validations: {}, projects: [],
  project: "untitled", dirty: false,
  // Analyses is the landing mode: it creates, uploads, validates and launches
  // in one place, which is what the step-by-step modes do one screen at a time.
  mode: "analyses", stage: "data", ctxTab: "params",
  showSteps: false,
  job: null, stream: null, lines: [], lineNo: 0, logFilter: "all",
  completed: new Set(), skipped: new Set(),
  params: {
    dataset_name: "", output_dir: "results", threads: 4, min_coverage: 20,
    bootstrap: 1000, rate_prior: 0.001, align_method: "mafft",
    phylo_model: "auto", tree_method: "ml", clock_model: "relaxed",
    recombination: true, intermediates: false, auto_qc: true,
    gpu: false, save_logs: true, pub_figures: false,
    advanced: false, experimental: false,
  },
  sel: { lineages: new Set(), countries: new Set(), period: null, locus: null, sample: null },
  colour: {}, order: [],
  viz: "tree", tool: "select", labels: true, treeLayout: "rectangular", bootstrapShown: false,
};

async function api(path, opts) {
  const r = await fetch(path, opts);
  if (!r.ok) {
    let d = r.statusText;
    try { d = (await r.json()).detail || d; } catch { /* non-JSON body */ }
    throw new Error(d);
  }
  return r.json();
}
const mark = () => { S.dirty = true; $("#dirty").hidden = false; };

/* ═══ MODES ═══════════════════════════════════════════════════════ */
/* Nine modes was too many to hold in mind, and four of them —
   Data input, Validation, Configure, Run — are the manual equivalents of
   what Analyses now does in one place: create, upload, validate, launch.
   They still work and are still reachable; they are just no longer the
   first thing a reader has to read past.

   Nothing is deleted. `advanced: true` only decides what the mode bar
   shows by default. */
const MODES = [
  { id: "analyses",  n: "1", label: "Analyses",     accent: "data" },
  { id: "visualize", n: "2", label: "Explore",      accent: "network" },
  { id: "interpret", n: "3", label: "Interpret",    accent: "interpret" },
  { id: "surveil",   n: "4", label: "Surveillance", accent: "evolution" },
  { id: "report",    n: "5", label: "Report",       accent: "report" },
  { id: "input",     n: "",  label: "Data input",   accent: "data",       advanced: true },
  { id: "validate",  n: "",  label: "Validation",   accent: "qc",         advanced: true },
  { id: "configure", n: "",  label: "Configure",    accent: "genomic",    advanced: true },
  { id: "run",       n: "",  label: "Manual run",   accent: "evolution",  advanced: true },
];

/* ═══ PIPELINE ════════════════════════════════════════════════════
   Dependencies are real: a stage is blocked until every stage it needs
   has completed, and until the tools it shells out to exist. */
const STAGES = [
  { id: "data",     ord: "01", name: "Data",           needs: [],            cmd: null,            mode: "input" },
  { id: "valid",    ord: "02", name: "Validation",     needs: ["data"],      cmd: null,            mode: "validate" },
  { id: "qc",       ord: "03", name: "QC",             needs: ["valid"],     cmd: "qc",            mode: "run" },
  { id: "prep",     ord: "04", name: "Preprocessing",  needs: ["qc"],        cmd: "align",         mode: "run", tools: ["mafft"] },
  { id: "genomic",  ord: "05", name: "Genomic",        needs: ["qc"],        cmd: "stage0",        mode: "run" },
  { id: "phylo",    ord: "06", name: "Phylogenetics",  needs: ["prep"],      cmd: "phylogenetics", mode: "run", tools: ["iqtree2", "treetime"] },
  { id: "evo",      ord: "07", name: "Evolution",      needs: ["phylo"],     cmd: "recombination", mode: "run" },
  { id: "model",    ord: "08", name: "Modelling",      needs: ["evo"],       cmd: "dh1",           mode: "run" },
  { id: "interp",   ord: "09", name: "Interpretation", needs: ["model"],     cmd: "score",         mode: "interpret", gated: true },
  { id: "report",   ord: "10", name: "Report",         needs: ["interp"],    cmd: "report",        mode: "report", gated: true },
];

/* In demo mode nothing is blocked. The point of the synthetic corpus is
   to walk the entire pipeline, so stages whose real tool is absent run as
   clearly-labelled simulations instead of sitting greyed out. */
const isDemo = () => !!S.data?.identity?.synthetic;

function missingTools(stage) {
  if (isDemo()) return [];
  return (stage.tools || []).filter((t) => !(S.env?.tools?.[t]?.present));
}
function simulated(stage) {
  return isDemo() && (!stage.cmd || (stage.tools || []).some((t) => !(S.env?.tools?.[t]?.present)));
}
function stageStatus(stage) {
  if (S.skipped.has(stage.id)) return "skipped";
  if (S.job && S.job.state === "running" && S.job.command === stage.cmd) return "running";
  if (S.completed.has(stage.id)) return "completed";
  if (missingTools(stage).length) return "blocked";
  if (!isDemo() && stage.gated && S.data && !S.data.gate.permitted) return "blocked";
  if (!stage.needs.every((n) => S.completed.has(n) || S.skipped.has(n))) return "blocked";
  return "not-started";
}
function blockReason(stage) {
  if (isDemo()) {
    return simulated(stage)
      ? "DEMO: runs as a labelled simulation — the real tool for this stage is not installed."
      : "";
  }
  const tools = missingTools(stage);
  if (tools.length) return `${tools.join(", ")} is not installed. See docs/installation.md.`;
  if (stage.gated && S.data && !S.data.gate.permitted) return `The D.H1 gate is ${S.data.gate.permission}. Stage 5/6 refuse to run until it returns SUPPORTED.`;
  const unmet = stage.needs.filter((n) => !S.completed.has(n) && !S.skipped.has(n));
  if (unmet.length) return "Requires " + unmet.map((n) => STAGES.find((s) => s.id === n).name).join(", ") + " to complete first.";
  return "";
}

/* ═══ BOOT ════════════════════════════════════════════════════════ */
async function boot() {
  try {
    const [env, commands, pathogens, inputs, projects] = await Promise.all([
      api("/api/env"), api("/api/commands"), api("/api/pathogens"),
      api("/api/inputs"), api("/api/projects"),
    ]);
    Object.assign(S, { env, commands, pathogens, inputs, projects });
    $("#build").textContent = "v" + env.version + " · workstation";
    const wanted = new URLSearchParams(location.search).get("p");
    S.pathogen = wanted || (pathogens.find((p) => p.provisioned && !p.synthetic) || pathogens[0] || {}).name;
    S.params.dataset_name = S.pathogen || "";
    try { S.data = await api(`/api/dataset/${S.pathogen}`); } catch { S.data = null; }
    if (S.data) {
      S.order = Object.keys(S.data.lineages);
      S.order.forEach((k, i) => { S.colour[k] = css(`--cat-${i % 8}`); });
      S.completed.add("data");
    }
  } catch (e) {
    notify("error", "Backend unreachable", e.message);
    return;
  }
  if (S.data?.identity?.synthetic) {
    document.body.classList.add("is-synth");
    document.body.append(el("div", { class: "synth" },
      el("span", { class: "dot" }),
      el("b", { text: "SYNTHETIC DEMONSTRATION DATA" }),
      el("span", { text: "— not a result. Fabricated so the open-gate interface can be reviewed." }),
      el("button", { class: "to-real", text: "Load real FMDV", onclick: () => { location.search = ""; } })));
  }
  renderModes(); renderIdentity(); renderSpine(); renderExec();
  setMode("input");
  pollSystem(); setInterval(pollSystem, 4000);
  wireChrome();
  log("success", `G4 ${S.env.version} ready — ${S.inputs.length} candidate input files indexed.`);
}

function applySkin(skin) {
  document.documentElement.setAttribute("data-skin", skin);
  localStorage.setItem("g4-skin-v2", skin);
  const pick = $("#skin-pick");
  if (pick) pick.value = skin;
  if (S.mode === "visualize") renderWork();
}

function wireChrome() {
  const saved = (() => { try { return localStorage.getItem("g4-skin-v2"); } catch { return null; } })();
  applySkin(saved || "forecast");
  try { S.showSteps = localStorage.getItem("g4-show-steps") === "1"; } catch { /* private mode */ }
  syncChrome();
  renderModes();
  const skinPick = $("#skin-pick");
  if (skinPick) skinPick.onchange = (e) => applySkin(e.target.value);
  api("/api/provenance").then((p) => {
    $("#trust-commit").textContent = `${p.software.name} ${p.software.version} · ${p.git_commit || "no commit"}`;
  }).catch(() => {});
  $("#btn-theme").onclick = () => {
    const now = document.documentElement.getAttribute("data-theme");
    document.documentElement.setAttribute("data-theme", now === "dark" ? "light" : "dark");
    if (S.mode === "visualize") renderWork();
  };
  $("#btn-help").onclick = () => helpSheet();
  const shortcutsBtn = $("#btn-shortcuts");
  if (shortcutsBtn) shortcutsBtn.onclick = () => shortcutsSheet();
  $("#log-toggle").onclick = () => { $("#console").classList.toggle("open"); $("#log-toggle").textContent = $("#console").classList.contains("open") ? "▼ LOG" : "▲ LOG"; };
  $("#log-clear").onclick = () => { S.lines = []; S.lineNo = 0; renderLog(); };
  $("#log-save").onclick = () => download(new Blob([S.lines.map((l) => `[${l.level}] ${l.text}`).join("\n")], { type: "text/plain" }), "g4_log.txt");
  $("#log-copy").onclick = () => {
    const errs = S.lines.filter((l) => l.level === "error" || l.level === "stderr").map((l) => l.text).join("\n");
    navigator.clipboard?.writeText(errs || "(no errors)").then(() => notify("success", "Copied", `${errs.split("\n").filter(Boolean).length} error lines on the clipboard.`));
  };
  $("#log-filter").onchange = (e) => { S.logFilter = e.target.value; renderLog(); };
  $("#q").oninput = (e) => runSearch(e.target.value);
  addEventListener("keydown", (e) => {
    if (e.target.matches("input,textarea,select")) return;
    // 1-5 are the five modes the bar shows by default; the advanced four
    // are reached through "More steps", not through a key nobody would guess.
    const idx = "12345".indexOf(e.key);
    if (idx >= 0) setMode(MODES[idx].id);
    if (e.key === "/") { e.preventDefault(); $("#q").focus(); }
    if (e.key === "?") helpSheet();
    if (e.key === "l") $("#log-toggle").click();
    // Zoom keys replace the ＋/− buttons that were taking up toolbar width
    // next to a scroll gesture and a double-click that already did the job.
    if (e.key === "+" || e.key === "=") zoomActive(1.25);
    if (e.key === "-") zoomActive(1 / 1.25);
  });
  addEventListener("resize", debounce(() => { if (S.mode === "visualize") renderWork(); renderSpine(); }, 150));
  addEventListener("beforeunload", (e) => { if (S.dirty) { e.preventDefault(); e.returnValue = ""; } });
}

/* The luminous canvas field was removed with the glass skin: a
   requestAnimationFrame loop painting halos behind a translucent
   shell, costing a repaint every frame on a page whose job is to
   show numbers. Nothing else referenced it. */

/* ═══ IDENTITY + TELEMETRY ════════════════════════════════════════ */
function renderIdentity() {
  const d = S.data;
  $("#id-project").textContent = S.project;
  $("#id-organism").textContent = d ? d.identity.display_name : "no dataset";
  $("#id-ref").textContent = d ? d.identity.reference : "—";
  $("#id-n").innerHTML = d ? `${d.identity.n_samples} <em>of ${d.identity.n_raw}</em>` : "—";
  $("#id-period").textContent = d?.identity.period ? d.identity.period.join("–") : "—";
  // Those two elements are hidden in the header; their content surfaces on hover.
  const organism = $("#idf-organism");
  if (organism) {
    organism.title = d
      ? `Reference ${d.identity.reference}`
        + (d.identity.period ? ` · sampled ${d.identity.period.join("–")}` : "")
      : "No dataset loaded";
  }
}

async function pollSystem() {
  let s;
  try { s = await api("/api/system"); } catch { return; }
  S.system = s;
  const cell = (k, pct, label) => el("div", { class: "cell", title: `${k}: ${label}` },
    el("span", { class: "k" }, el("span", { text: k }), el("b", { text: Math.round(pct) + "%" })),
    el("span", { class: "meter" }, el("i", { class: pct > 90 ? "bad" : pct > 75 ? "warn" : "", style: `width:${Math.min(pct, 100)}%` })),
  );
  $("#sysmon").replaceChildren(
    cell("cpu", s.cpu_percent, `${s.cpu_count} threads, load ${s.load[0]}`),
    cell("mem", s.memory.percent, `${bytes(s.memory.used)} of ${bytes(s.memory.total)}`),
    cell("disk", s.disk.percent, `${bytes(s.disk.free)} free`),
  );
  renderTelemetry();
}

function renderTelemetry() {
  const s = S.system, j = S.job;
  const t = (k, v) => el("div", { class: "t" }, el("span", { class: "k", text: k }), el("span", { class: "v", text: v }));
  $("#telemetry").replaceChildren(
    t("elapsed", j ? clock(j.duration) : "—"),
    t("cpu", s ? Math.round(s.cpu_percent) + "%" : "—"),
    t("mem", s ? bytes(s.memory.used) : "—"),
    t("threads", s ? String(s.cpu_count) : "—"),
    t("jobs", j && !["succeeded", "failed", "gate_closed", "cancelled"].includes(j.state) ? "1 active" : "idle"),
  );
}

/* ═══ MODE BAR + SPINE ════════════════════════════════════════════ */
function renderModes() {
  // An advanced mode stays visible while you are in it, so selecting one
  // and then having it vanish cannot happen.
  const shown = MODES.filter((m) => !m.advanced || S.showSteps || m.id === S.mode);
  $("#modes").replaceChildren(
    ...shown.map((m) => el("button", {
      class: "mode" + (m.id === S.mode ? " on" : "") + (m.advanced ? " adv" : ""),
      "data-accent": m.accent,
      onclick: () => setMode(m.id),
    }, m.n ? el("span", { class: "n", text: m.n }) : null, el("span", { text: m.label }))),
    el("button", {
      class: "btn tertiary sm",
      title: "Data input, Validation, Configure and Manual run — the step-by-step "
           + "equivalents of what Analyses does in one place",
      text: S.showSteps ? "Fewer steps" : "More steps",
      onclick: () => { S.showSteps = !S.showSteps; saveShowSteps(); syncChrome(); renderModes(); },
    }),
    el("div", { class: "push" },
      el("button", { class: "btn secondary sm", text: "Save", onclick: () => saveProject(false) }),
      el("button", { class: "btn tertiary sm", text: "Open", onclick: openProject }),
    ),
  );
}

/* The stage spine is the map of the step-by-step path. With the steps folded
   away it is a permanent 96px column describing a route nobody is walking, so
   it follows the same disclosure. It is still rendered — only the column is
   hidden — so #flow-nodes stays inspectable and drawSpineEdges stays correct
   the moment it is shown again. */
function syncChrome() {
  document.body.classList.toggle("no-spine", !S.showSteps);
  const flow = $("#flow");
  if (flow) flow.hidden = !S.showSteps;
  if (S.showSteps) requestAnimationFrame(drawSpineEdges);
}
function saveShowSteps() {
  try { localStorage.setItem("g4-show-steps", S.showSteps ? "1" : "0"); } catch { /* private mode */ }
}

function setMode(id) {
  S.mode = id;
  renderModes();
  renderWork();
  renderCtx();
}

function renderSpine() {
  const host = $("#flow-nodes");
  host.replaceChildren(...STAGES.map((st) => {
    const status = stageStatus(st);
    return el("button", {
      class: `stagenode ${status}${st.id === S.stage ? " on" : ""}`,
      title: `${st.name} — ${status}${blockReason(st) ? "\n" + blockReason(st) : ""}`,
      onclick: () => { S.stage = st.id; setMode(st.mode); },
      onmousemove: (e) => tip(e, st.name.toUpperCase(), [["status", status], ["needs", st.needs.length ? st.needs.join(", ") : "nothing"], ["command", st.cmd || "—"]], blockReason(st)),
      onmouseleave: hideTip,
    },
      el("span", { class: "ord", text: st.ord }),
      el("span", {}, el("span", { class: "glyph" }), el("span", { class: "nm", text: st.name })),
    );
  }));
  requestAnimationFrame(drawSpineEdges);
}

function drawSpineEdges() {
  const flow = $("#flow"), line = $("#flow-svg");
  // A hidden spine measures as a zero-height box, so every edge would be drawn
  // at the same point. Nothing to draw until the column is back.
  if (!flow || flow.hidden) return;
  const box = flow.getBoundingClientRect();
  line.setAttribute("height", flow.scrollHeight);
  line.replaceChildren();
  const nodes = [...flow.querySelectorAll(".stagenode")];
  // Edges are drawn between stage nodes; with no nodes rendered yet there is
  // nothing to connect, and indexing into the empty list would throw.
  if (nodes.length < STAGES.length) return;
  const at = (i) => { const r = nodes[i].getBoundingClientRect(); return { x: 24, y: r.top - box.top + flow.scrollTop + 14 }; };
  STAGES.forEach((st, i) => {
    for (const need of st.needs) {
      const j = STAGES.findIndex((s) => s.id === need);
      if (j < 0) continue;
      const a = at(j), b = at(i);
      const done = S.completed.has(need) || S.skipped.has(need);
      line.append(svg("path", {
        d: `M${a.x} ${a.y + 8} V${b.y - 8}`, fill: "none",
        stroke: done ? css("--st-complete") : css("--edge"),
        "stroke-width": 1, "stroke-dasharray": done ? null : "2 3",
      }));
    }
  });
}

/* ═══ WORKSPACE ═══════════════════════════════════════════════════ */
const WORK = {
  analyses: workAnalyses, input: workInput, validate: workValidate, configure: workConfigure,
  run: workRun, visualize: workVisualize, interpret: workInterpret,
  surveil: workSurveillance, report: workReport,
};

function renderWork() {
  const host = $("#work");
  host.replaceChildren();
  host.className = "work" + (S.mode === "visualize" ? " canvas" : " pad");
  $("#stage-tools").replaceChildren();
  $("#footbar").replaceChildren();
  WORK[S.mode](host);
}

/* ── input specification ───────────────────────────────────────────
   Documented in the interface rather than in a README, and rendered as a
   live checklist: each row says whether the project currently contains a
   file that could satisfy it. A specification you have to leave the app
   to read is one nobody reads. */
const REQUIRED_INPUTS = [
  {
    key: "reference", name: "Reference genome", format: "FASTA",
    ext: [".fasta", ".fa", ".fna"], records: "exactly 1 record",
    why: "The coordinate system. Every locus position, every alignment column and every genome track is measured against this sequence.",
    example: ">AY593823.1 Foot-and-mouth disease virus\nGTTGAAAGGGGGCGCTAGGGTCTCACCCCTAGCGTTGGC…",
    match: (f) => f.format === "fasta",
  },
  {
    key: "sequences", name: "Corpus sequences", format: "FASTA",
    ext: [".fasta", ".fa"], records: "many records, one per genome",
    why: "The genomes to analyse. Header IDs must match the metadata accession column exactly.",
    example: ">PQ587570.1\nGTTGAAAGGGG…\n>PQ587565.1\nGTTGAAAGGGG…",
    match: (f) => f.format === "fasta",
  },
  {
    key: "metadata", name: "Metadata table", format: "TSV",
    ext: [".tsv", ".txt"], records: "one row per genome, with a header row",
    why: "Lineage, date and place. Drives the per-lineage floor, the temporal views and the geography.",
    example: "accession\tcollection_date\tcountry\thost\tserotype\nPQ587570.1\t18-Jul-2024\tIndia\tcattle\tO",
    match: (f) => f.format === "tsv" || f.format === "csv",
  },
];

const OPTIONAL_INPUTS = [
  { name: "Aligned FASTA", format: "FASTA", skips: "Stage 04 Preprocessing",
    note: "Reference-anchored; every record the same length. Supply this if MAFFT is not installed.",
    match: (f) => f.format === "fasta" },
  { name: "Rooted tree", format: "Newick", skips: "Stage 06 Phylogenetics",
    note: "Tip labels must be the accessions. Supply this if IQ-TREE is not installed.",
    match: (f) => f.format === "newick" },
  { name: "G4 Reference Atlas", format: "TSV", skips: "Stage 05 Genomic",
    note: "A previously built Atlas, to avoid re-scanning the reference.",
    match: (f) => f.format === "tsv" && /atlas/i.test(f.name) },
];

const META_COLUMNS = [
  ["accession", true, "Unique ID. Must match the FASTA header exactly.", "PQ587570.1"],
  ["collection_date", false, "Sampling date. Counts toward the 90% completeness floor.", "18-Jul-2024"],
  ["country", false, "Sampling location. Counts toward completeness.", "India"],
  ["host", false, "Host species. Counts toward completeness.", "cattle"],
  ["serotype", false, "The lineage field. Drives the ≥20-per-lineage floor.", "O"],
  ["organism", false, "Fallback when serotype is blank.", "FMDV - type O"],
  ["isolate", false, "Second fallback for lineage.", "O/IND/53/2024"],
  ["strain", false, "Third fallback for lineage.", "PanAsia-2"],
  ["length", false, "Sequence length, if you have it.", "8201"],
];

function inputSpec() {
  const open = (() => { try { return localStorage.getItem("g4-spec") !== "closed"; } catch { return true; } })();
  const found = (spec) => S.inputs.filter(spec.match);

  const body = el("div", { class: "stack", style: "padding:14px 16px" });

  body.append(el("p", { class: "hint", style: "max-width:78ch",
    text: "G4 reads files from inside the project directory — a browser cannot write to your disk, so there is no upload. Copy your files into data/ and press Rescan folder." }));

  const row = (name, format, records, why, count, example) =>
    el("div", { style: "display:grid;grid-template-columns:170px 92px 1fr auto;gap:14px;align-items:start;padding:11px 0;border-bottom:1px solid var(--hair)" },
      el("div", {}, el("b", { style: "font-size:12.5px", text: name }),
        records ? el("div", { class: "hint", text: records }) : null),
      el("span", { class: "mono", style: "font-size:11px;color:var(--ink-2)", text: format }),
      el("div", {}, el("span", { style: "font-size:11.5px;line-height:1.5;color:var(--ink-2)", text: why }),
        example ? el("pre", { class: "mono",
          style: "margin-top:6px;font-size:10.5px;line-height:1.6;white-space:pre-wrap;color:var(--ink-3);background:var(--bg-canvas);border:1px solid var(--hair);padding:7px 9px;overflow-x:auto",
          text: example }) : null),
      count);

  const sec = el("div", { class: "sec" }, el("h3", { text: "Required — the pipeline will not start without these" }));
  for (const spec of REQUIRED_INPUTS) {
    const n = found(spec).length;
    sec.append(row(spec.name, spec.format, spec.records, spec.why,
      el("span", { class: "st " + (n ? "valid" : "not-started"), title: n ? `${n} candidate file${n !== 1 ? "s" : ""} in the project` : "nothing in the project matches this",
        text: n ? `${n} found` : "none" }),
      spec.example));
  }
  body.append(sec);

  const opt = el("div", { class: "sec" }, el("h3", { text: "Optional — supplying one skips the stage that would produce it" }));
  for (const spec of OPTIONAL_INPUTS) {
    const n = found(spec).length;
    opt.append(row(spec.name, spec.format, spec.skips, spec.note,
      el("span", { class: "st " + (n ? "valid" : "skipped"), text: n ? `${n} found` : "optional" })));
  }
  body.append(opt);

  const cols = el("div", { class: "sec" }, el("h3", { text: "Metadata columns" }));
  const table = el("table", { class: "grid" },
    el("thead", {}, el("tr", {}, ...["Column", "Required", "Purpose", "Example"].map((h) => el("th", { text: h })))),
    el("tbody", {}, ...META_COLUMNS.map(([name, req, why, ex]) => el("tr", {},
      el("td", { class: "mono", text: name }),
      el("td", {}, el("span", { class: "st " + (req ? "invalid" : "skipped"), text: req ? "required" : "optional" })),
      el("td", { style: "font-size:11.5px", text: why }),
      el("td", { class: "mono dim", text: ex })))));
  cols.append(el("div", { class: "scroll-x" }, table));
  body.append(cols);

  body.append(el("p", { class: "caution", style: "border-left:2px solid var(--st-warning);color:var(--st-warning);padding-left:10px;font-size:11.5px;line-height:1.55;max-width:80ch",
    text: "The FASTA header ID and the metadata accession must match exactly. G4 takes the first whitespace-delimited token of the header, so \">AY593823.1 Foot-and-mouth…\" matches \"AY593823.1\". If they disagree, the alignment and the metadata describe different corpora and the run stops." }));

  body.append(el("div", { class: "btn-row" },
    btn("secondary", "Copy data/ path", () => {
      navigator.clipboard?.writeText(S.env.repo_root + "/data")
        .then(() => notify("success", "Copied", S.env.repo_root + "/data"));
    }),
    btn("secondary", "Rescan folder", rescanInputs),
    btn("tertiary", "Open a config", () => previewFile("config/fmdv.yaml"))));

  const panel = el("div", { class: "panel" },
    el("header", {},
      el("span", { class: "tag", text: "What to put in" }),
      el("span", { class: "hint", text: "required files, formats and the metadata schema" }),
      el("div", { class: "push" },
        btn("tertiary sm", open ? "Hide" : "Show", (e) => {
          const now = body.hidden;
          body.hidden = !now;
          e.target.textContent = now ? "Hide" : "Show";
          try { localStorage.setItem("g4-spec", now ? "open" : "closed"); } catch { /* private mode */ }
        }))),
    body);
  body.hidden = !open;
  return panel;
}

/* ── ANALYSES ─────────────────────────────────────────────────────────
   The unit of work a researcher creates, as opposed to a `job`, which is
   one process the runner executed. An analysis outlives its job: it
   survives a restart, records what it ran on by checksum, and keeps its
   error text after the job has been evicted from the runner's ring.

   The endpoints existed and had no surface, which is the same defect
   /api/stage5 and /api/report-card had: the server does the real work and
   the interface cannot reach it. ─────────────────────────────────── */
const AN = { list: [], selected: null, detail: null, creating: false,
             form: { name: "", pathogen: "", inputs: [] },
             profile: "conda_free", resume: true, busy: "" };

const AN_STATUS_COLOUR = {
  QUEUED: "--ink-3", VALIDATING: "--st-running", RUNNING: "--st-running",
  COMPLETED: "--st-complete", FAILED: "--st-error", CANCELLED: "--ink-3",
};

const when = (t) => t ? new Date(t * 1000).toLocaleString() : "—";

async function anRefresh() {
  try { AN.list = await api("/api/analyses"); } catch (e) { AN.list = []; }
  if (AN.selected) {
    try { AN.detail = await api(`/api/analyses/${AN.selected}/results`); }
    catch { AN.detail = null; }
  }
  if (S.mode === "analyses") renderWork();
}

/* Poll only while something is live. A dashboard that polls a finished
   run forever is a dashboard that is wrong about what it is watching. */
let AN_TIMER = null;
function anPoll() {
  clearInterval(AN_TIMER);
  AN_TIMER = setInterval(() => {
    if (S.mode !== "analyses") { clearInterval(AN_TIMER); return; }
    if (AN.list.some((a) => !a.terminal)) anRefresh();
  }, 4000);
}

async function anCreate() {
  const f = AN.form;
  if (!f.name.trim()) return notify("error", "Name required", "Give the analysis a name you will recognise later.");
  if (!f.pathogen) return notify("error", "Pathogen required", "Select which pathogen this analysis is for.");
  AN.busy = "creating"; renderWork();
  try {
    const created = await api("/api/analyses", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ name: f.name, pathogen: f.pathogen, inputs: f.inputs }),
    });
    AN.selected = created.id;
    AN.creating = false;
    AN.form = { name: "", pathogen: f.pathogen, inputs: [] };
    notify("success", "Analysis created", `${created.name} — validate its inputs next.`);
    await anValidate(created.id);
  } catch (e) {
    notify("error", "Could not create", e.message);
  } finally { AN.busy = ""; await anRefresh(); }
}

async function anValidate(id) {
  AN.busy = "validating"; renderWork();
  try {
    const out = await api(`/api/analyses/${id}/validate`, { method: "POST" });
    AN.validation = out;
    if (out.valid) notify("success", "Inputs valid", "Ready to launch.");
    else notify("error", `${out.errors.length} problem(s) with the inputs`, out.errors[0] || "");
  } catch (e) {
    notify("error", "Validation failed", e.message);
  } finally { AN.busy = ""; await anRefresh(); }
}

async function anLaunch(id) {
  AN.busy = "launching"; renderWork();
  try {
    await api(`/api/analyses/${id}/launch`, {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ profile: AN.profile, resume: AN.resume }),
    });
    notify("running", "Nextflow started", `profile ${AN.profile}. A closed D.H1 gate completes normally.`);
    anPoll();
  } catch (e) {
    notify("error", "Launch refused", e.message);
  } finally { AN.busy = ""; await anRefresh(); }
}

async function anCancel(id) {
  try { await api(`/api/analyses/${id}/cancel`, { method: "POST" }); }
  catch (e) { notify("error", "Cancel failed", e.message); }
  await anRefresh();
}

/* Upload straight into the form, so a file goes from the researcher's
   disk to a declared input without a detour through the file tray. */
async function anUpload(fileList, role) {
  const files = Array.from(fileList || []);
  if (!files.length) return;
  for (const file of files) {
    const body = new FormData();
    body.append("file", file);
    try {
      const r = await fetch("/api/upload", { method: "POST", body });
      const d = await r.json().catch(() => ({}));
      if (!r.ok) throw new Error(d.detail || `HTTP ${r.status}`);
      AN.form.inputs.push({ path: d.path, role });
      notify("success", "Staged", `${d.name} (${bytes(d.bytes)}) as ${role}`);
    } catch (e) { notify("error", `Upload refused — ${file.name}`, e.message); }
  }
  renderWork();
}

function anCompletenessBadge(comp) {
  if (!comp) return null;
  const bad = comp.category !== "complete";
  return el("div", { style: `margin-top:6px;padding:8px 11px;border-left:2px solid var(--${bad ? "st-warning" : "st-complete"})` },
    el("p", { style: `font-size:12.5px;font-weight:600;color:var(--${bad ? "st-warning" : "st-complete"})`,
              text: comp.description }),
    comp.caveat ? el("p", { style: "font-size:11.5px;line-height:1.55;color:var(--ink-2);margin-top:4px;max-width:84ch",
                            text: comp.caveat }) : null);
}

function anCreateForm() {
  const f = AN.form;
  const pathogens = (S.pathogens || []).filter((p) => p.provisioned);
  if (!f.pathogen && pathogens.length) f.pathogen = pathogens[0].name;

  const fileBtn = (role, label) => btn("secondary", label, () => {
    const i = el("input", { type: "file", multiple: true });
    i.onchange = () => anUpload(i.files, role);
    i.click();
  });

  return el("div", { class: "sec" },
    el("h3", { text: "New analysis" }),
    field("Analysis name", el("input", { class: "input", value: f.name, placeholder: "FMDV 2026 surveillance",
      onchange: (e) => { f.name = e.target.value; } })),
    field("Pathogen", el("select", { class: "select", onchange: (e) => { f.pathogen = e.target.value; } },
      ...pathogens.map((p) => el("option", { value: p.name, text: `${p.name} — ${p.display_name}`,
                                             selected: p.name === f.pathogen })))),
    field("Sequence data",
      el("div", { class: "btn-row" }, fileBtn("sequences", "Upload FASTA"),
         fileBtn("reference", "Upload reference (optional)"),
         fileBtn("metadata", "Upload metadata (optional)")),
      "Uploads are staged under data/uploads/ and recorded by SHA-256. Leaving this empty runs "
      + "against the pathogen's configured corpus."),
    f.inputs.length ? el("div", { style: "margin-top:8px" },
      ...f.inputs.map((i, idx) => el("div", { style: "display:flex;gap:10px;align-items:center;padding:4px 0" },
        el("span", { class: "mono", style: "font-size:11px", text: i.path }),
        el("span", { class: "st", text: i.role }),
        btn("tertiary sm", "Remove", () => { f.inputs.splice(idx, 1); renderWork(); })))) : null,
    el("div", { class: "btn-row", style: "margin-top:12px" },
      btn("primary", AN.busy === "creating" ? "Creating…" : "Create analysis", anCreate, !!AN.busy),
      btn("tertiary", "Cancel", () => { AN.creating = false; renderWork(); })),
  );
}

function anValidationPanel() {
  const v = AN.validation;
  if (!v) return null;
  return el("div", { class: "sec" }, el("h3", { text: "Input validation" }),
    el("p", { style: `font-size:12.5px;font-weight:600;color:var(--st-${v.valid ? "complete" : "error"})`,
              text: v.valid ? "All inputs valid." : `${v.errors.length} problem(s) found.` }),
    ...v.errors.map((e) => el("p", { style: "font-size:12px;line-height:1.55;color:var(--st-error);padding-left:10px;border-left:2px solid var(--st-error);margin-top:5px;max-width:88ch", text: e })),
    ...v.reports.map((r) => el("div", { style: "padding:8px 0;border-bottom:1px solid var(--hair)" },
      el("p", { class: "mono", style: "font-size:11.5px", text: `${r.path}  [${r.role}]  ${r.status}` }),
      r.n_records != null ? el("p", { class: "hint",
        text: `${r.n_records.toLocaleString()} records · ${(r.total_bases || 0).toLocaleString()} bases` }) : null,
      anCompletenessBadge(r.completeness))),
  );
}

function anDetail() {
  const d = AN.detail;
  if (!d) return null;
  const a = d.analysis;
  const row = (k, v) => el("div", { style: "display:grid;grid-template-columns:180px 1fr;gap:12px;padding:4px 0" },
    el("span", { class: "hint", text: k }),
    el("span", { class: "mono", style: "font-size:11.5px", text: String(v ?? "—") }));

  return el("div", {},
    el("div", { class: "sec" }, el("h3", { text: a.name }),
      el("div", { style: "display:flex;gap:16px;align-items:baseline;flex-wrap:wrap" },
        el("span", { style: `font-size:19px;font-weight:600;color:var(${AN_STATUS_COLOUR[a.status]})`, text: a.status }),
        el("span", { class: "mono", style: "font-size:11px;color:var(--ink-3)",
          text: `${a.pathogen} · ${a.id} · ${a.duration ? clock(a.duration) : "not started"}` })),
      a.error ? el("p", { style: "font-size:12px;line-height:1.6;color:var(--st-error);padding:8px 11px;margin-top:8px;border-left:2px solid var(--st-error);white-space:pre-wrap;max-width:90ch", text: a.error }) : null,
      d.gate ? el("p", { style: `font-size:12px;margin-top:8px;color:var(--${d.gate.permitted ? "st-complete" : "st-warning"})`,
                         text: `D.H1 gate: ${d.gate.permission}` }) : null,
      el("div", { class: "btn-row", style: "margin-top:12px" },
        btn("secondary", "Validate", () => anValidate(a.id), !!AN.busy),
        btn("primary", AN.busy === "launching" ? "Launching…" : "Launch pipeline",
            () => anLaunch(a.id), !!AN.busy || a.status === "RUNNING"),
        btn("danger", "Cancel", () => anCancel(a.id), a.terminal),
        btn("tertiary", "Refresh", anRefresh)),
      el("div", { style: "display:flex;gap:18px;align-items:center;margin-top:10px;flex-wrap:wrap" },
        el("label", { style: "display:flex;gap:6px;align-items:center;font-size:12px" },
          el("span", { text: "profile" }),
          el("select", { class: "select", style: "width:auto",
            onchange: (e) => { AN.profile = e.target.value; } },
            ...["conda_free", "docker", "singularity", "standard"].map((p) =>
              el("option", { value: p, text: p, selected: p === AN.profile })))),
        el("label", { style: "display:flex;gap:6px;align-items:center;font-size:12px" },
          el("input", { type: "checkbox", checked: AN.resume,
                        onchange: (e) => { AN.resume = e.target.checked; } }),
          el("span", { text: "resume" })))),

    el("div", { style: "height:18px" }),
    /* Reproducibility. Enough to re-run from the record alone, which is
       the point of recording it at all. */
    el("div", { class: "sec" }, el("h3", { text: "Reproducibility" }),
      row("created", when(a.created_at)), row("started", when(a.started_at)),
      row("finished", when(a.finished_at)), row("git commit", a.git_commit),
      row("pipeline version", a.pipeline_version), row("nextflow", a.nextflow_version),
      row("exit code", a.exit_code), row("output directory", a.outdir),
      ...a.inputs.map((i) => row(`input · ${i.role}`,
        `${i.path}  sha256:${(i.sha256 || "").slice(0, 16)}…`))),

    el("div", { style: "height:18px" }),
    el("div", { class: "sec" }, el("h3", { text: `Output files (${d.n_files})` }),
      d.n_files === 0
        ? el("p", { class: "hint", text: "No outputs yet." })
        : el("div", { class: "scroll-x" }, el("table", { class: "grid" },
            el("thead", {}, el("tr", {}, ...["File", "Size", ""].map((h) => el("th", { text: h })))),
            el("tbody", {}, ...d.files.map((f) => el("tr", {},
              el("td", { class: "mono", text: f.path.split("/").slice(-2).join("/") }),
              el("td", { class: "num", text: bytes(f.bytes) }),
              el("td", {}, btn("tertiary sm", "View", () => previewFile(f.path))))))))),
  );
}

function workAnalyses(host) {
  $("#stage-title").textContent = "Analyses";
  $("#stage-sub").textContent = `${AN.list.length} recorded · persisted across restarts`;
  $("#stage-tools").replaceChildren(
    btn("primary", "+ New analysis", () => { AN.creating = true; AN.validation = null; renderWork(); }),
    btn("tertiary", "Refresh", anRefresh),
  );

  if (AN.creating) { host.append(anCreateForm()); return; }

  const table = el("table", { class: "grid" },
    el("thead", {}, el("tr", {}, ...["Name", "Pathogen", "Status", "Created", "Duration", ""]
      .map((h) => el("th", { text: h })))),
    el("tbody", {}, ...AN.list.map((a) => el("tr", { class: a.id === AN.selected ? "sel" : "" },
      el("td", {}, el("b", { text: a.name })),
      el("td", { class: "mono dim", text: a.pathogen }),
      el("td", {}, el("span", { style: `color:var(${AN_STATUS_COLOUR[a.status]});font-weight:600;font-size:11px`, text: a.status })),
      el("td", { class: "mono dim", style: "font-size:11px", text: when(a.created_at) }),
      el("td", { class: "num", text: a.duration ? clock(a.duration) : "—" }),
      el("td", {}, btn("tertiary sm", "Open", () => { AN.selected = a.id; anRefresh(); }))))),
  );

  host.append(
    AN.list.length
      ? el("div", { class: "sec" }, el("h3", { text: "Analyses" }), el("div", { class: "scroll-x" }, table))
      : el("p", { class: "blank" }, el("b", { text: "No analyses yet." }),
          "An analysis records what was run, on which inputs by checksum, with which parameters, and what came out. "
          + "Create one to run the pipeline from here rather than from a terminal."),
    AN.validation ? el("div", { style: "height:18px" }) : null,
    anValidationPanel(),
    AN.detail ? el("div", { style: "height:18px" }) : null,
    anDetail(),
  );
  anPoll();
}

/* ── 01 · DATA INPUT ────────────────────────────────────────────── */
function workInput(host) {
  $("#stage-title").textContent = "Data acquisition";
  $("#stage-sub").textContent = `${S.inputs.length} files indexed under data/, results/, config/`;
  $("#stage-tools").replaceChildren(
    btn("primary", "+ Add dataset", addDataset),
    btn("secondary", "Browse files", browseFiles),
    btn("secondary", "Import URL", importUrl),
    btn("secondary", "Import accession", importAccession),
    btn("secondary", "Rescan folder", rescanInputs),
    btn("tertiary", "Load project", openProject),
  );

  const tray = el("div", {
    class: "tray",
    ondragover: (e) => { e.preventDefault(); tray.classList.add("hot"); },
    ondragleave: () => tray.classList.remove("hot"),
    ondrop: (e) => { e.preventDefault(); tray.classList.remove("hot"); dropped(e); },
  },
    el("div", { class: "rule-mark" }, el("i"), el("i"), el("i"), el("i")),
    el("div", {},
      el("h3", { text: "Specimen tray — drop files to stage them" }),
      el("p", { text: "FASTA · FASTQ · BAM · VCF · CSV · TSV · XLSX · Newick · Nexus · GFF · GenBank. Files are inspected on arrival; nothing is copied until you stage it." }),
    ),
    el("div", { class: "btn-row" }, btn("secondary", "Browse", browseFiles), btn("tertiary", "Clear all", clearStaged)),
  );

  const rows = S.inputs.slice(0, 120);
  const table = el("table", { class: "grid" },
    el("thead", {}, el("tr", {},
      ...["File", "Type", "Size", "Records", "Status", ""].map((h) => el("th", { text: h })))),
    el("tbody", {}, ...rows.map((f) => {
      const v = S.validations[f.path];
      return el("tr", {},
        el("td", { class: "mono", title: f.path }, el("span", { text: f.name })),
        el("td", { class: "mono", text: f.format }),
        el("td", { class: "num", text: bytes(f.size) }),
        el("td", { class: "num", text: v?.records != null ? v.records.toLocaleString() : "—" }),
        el("td", {}, el("span", { class: "st " + (v ? v.status : "not-started"), text: v ? v.status : "unchecked" })),
        el("td", {}, el("div", { class: "btn-row" },
          btn("tertiary sm", "Validate", () => validateOne(f.path), null, "Parse the file and report what it contains"),
          btn("tertiary sm", "Preview", () => previewFile(f.path)),
          btn("tertiary sm", "Use", () => { useAsInput(f); }),
        )),
      );
    })),
  );

  host.append(
    inputSpec(),
    el("div", { style: "height:16px" }),
    tray,
    el("div", { style: "height:18px" }),
    el("div", { class: "sec" },
      el("h3", { text: "Indexed inputs" }),
      el("div", { class: "scroll-x" }, table),
      rows.length < S.inputs.length ? el("p", { class: "hint", text: `Showing ${rows.length} of ${S.inputs.length}.` }) : null,
    ),
  );

  $("#footbar").replaceChildren(
    btn("secondary", "Validate all", validateAll),
    btn("tertiary", "Clear all", clearStaged),
    el("div", { class: "push" },
      el("span", { class: "hint", text: `${Object.keys(S.validations).length} validated · ${S.inputs.length} indexed` })),
  );
}

/* ── 02 · VALIDATION ────────────────────────────────────────────── */
function workValidate(host) {
  $("#stage-title").textContent = "Input validation";
  const list = Object.values(S.validations);
  $("#stage-sub").textContent = list.length ? `${list.length} file${list.length !== 1 ? "s" : ""} checked` : "nothing checked yet";
  $("#stage-tools").replaceChildren(
    btn("primary", "Validate all", validateAll),
    btn("secondary", "Run precheck", precheck),
    btn("secondary", "Check compatibility", compatibility),
    btn("tertiary", "Download report", downloadValidation),
  );

  if (!list.length) {
    host.append(el("p", { class: "blank" }, el("b", { text: "No files validated yet." }),
      "Validation parses each file and reports what it actually contains — record counts, alphabet violations, missing metadata, unbalanced trees. Start from Data input, or press Validate all."));
    return;
  }

  for (const v of list) {
    host.append(el("div", { class: "panel", style: "margin-bottom:14px" },
      el("header", {},
        el("span", { class: "st " + v.status, text: v.status }),
        el("span", { class: "mono", style: "font-size:12px", text: v.name }),
        el("span", { class: "hint", text: `${v.format} · ${bytes(v.size)}${v.records != null ? ` · ${v.records.toLocaleString()} ${v.record_label}` : ""}` }),
        el("div", { class: "push" },
          btn("tertiary sm", "Re-check", () => validateOne(v.path)),
          btn("tertiary sm", "Preview", () => previewFile(v.path)),
        ),
      ),
      el("div", { class: "body" },
        ...v.findings.map((f) => el("div", { style: "display:grid;grid-template-columns:auto 1fr auto;gap:10px;align-items:baseline" },
          el("span", { class: "st " + ({ pass: "valid", warn: "warning", fail: "invalid" })[f.level] }),
          el("span", { style: "font-size:12.5px;line-height:1.5", text: f.message }),
          f.fixable ? btn("tertiary sm", "Fix automatically", () => autoFix(v, f)) : el("span"),
        )),
        Object.keys(v.detail).length ? el("p", { class: "hint mono", text: JSON.stringify(v.detail) }) : null,
      ),
    ));
  }

  $("#footbar").replaceChildren(
    btn("secondary", "Ignore warnings", () => { S.completed.add("valid"); renderSpine(); notify("warning", "Warnings ignored", "Validation warnings acknowledged; QC is now reachable."); }),
    btn("danger", "Remove invalid", removeInvalid),
    btn("tertiary", "Return to input", () => setMode("input")),
    el("div", { class: "push" }, btn("primary", "Accept and continue", () => {
      S.completed.add("valid"); renderSpine(); setMode("configure");
      notify("success", "Validation accepted", "Stage 03 QC is now reachable.");
    })),
  );
}

/* ── 03 · CONFIGURE ─────────────────────────────────────────────── */
function workConfigure(host) {
  $("#stage-title").textContent = "Analysis configuration";
  $("#stage-sub").textContent = "parameters apply to the next run";
  $("#stage-tools").replaceChildren(
    btn("secondary", "Import config", importConfig),
    btn("secondary", "Export config", exportConfig),
    btn("tertiary", "Restore defaults", restoreDefaults),
  );

  const P = S.params;
  host.append(el("div", { class: "cols two" },
    el("div", { class: "sec" }, el("h3", { text: "Dataset" }),
      field("Dataset name", input("dataset_name")),
      field("Reference genome", el("div", { class: "path" }, input("reference", S.data?.identity.reference || ""), btn("secondary", "Browse", browseFiles))),
      field("Output directory", el("div", { class: "path" }, input("output_dir"), btn("secondary", "Browse", browseFiles))),
    ),
    el("div", { class: "sec" }, el("h3", { text: "Compute" }),
      field("Threads", stepper("threads", 1, S.system?.cpu_count || 16, 1), `Host reports ${S.system?.cpu_count ?? "?"} logical cores.`),
      field("Minimum coverage", stepper("min_coverage", 1, 1000, 1), "Sites below this depth are masked."),
      field("Bootstrap replicates", stepper("bootstrap", 100, 10000, 100), "≥1000 for publication-grade support values."),
      field("Evolutionary rate prior", stepper("rate_prior", 0.00001, 1, 0.0005), "Substitutions per site per year. FMDV is ~1e-3."),
    ),
    el("div", { class: "sec" }, el("h3", { text: "Methods" }),
      field("Alignment method", select("align_method", [["mafft", "MAFFT"], ["muscle", "MUSCLE"], ["clustalo", "Clustal Omega"], ["none", "Pre-aligned"]])),
      field("Phylogenetic model", select("phylo_model", [["auto", "Auto select (ModelFinder)"], ["GTR+F+I+G4", "GTR+F+I+G4"], ["HKY+G", "HKY+G"], ["JC", "Jukes–Cantor"]])),
      field("Tree method", select("tree_method", [["ml", "Maximum likelihood"], ["nj", "Neighbour joining"], ["bayes", "Bayesian (BEAST)"], ["parsimony", "Parsimony"]])),
      field("Clock model", select("clock_model", [["relaxed", "Relaxed clock"], ["strict", "Strict clock"], ["local", "Local clock"], ["none", "No clock"]])),
    ),
    el("div", { class: "sec" }, el("h3", { text: "Options" }),
      check("recombination", "Enable recombination detection", "Mandatory for every pathogen in this pipeline (Section 11) — it cannot be disabled at run time."),
      check("intermediates", "Generate intermediate files"),
      check("auto_qc", "Perform automatic QC"),
      check("gpu", "Use GPU acceleration", S.env?.tools?.nextflow?.present ? "" : "No GPU runtime detected on this host."),
      check("save_logs", "Save analysis logs"),
      check("pub_figures", "Generate publication-quality figures"),
      el("div", { style: "height:6px" }),
      toggle("advanced", "Advanced analysis"),
      toggle("experimental", "Experimental features", "Unvalidated methods. Results must not be reported."),
    ),
  ));

  $("#footbar").replaceChildren(
    btn("tertiary", "Restore defaults", restoreDefaults),
    el("div", { class: "push" },
      btn("secondary", "Save parameters", () => saveProject(false)),
      btn("primary", "Continue to run", () => setMode("run")),
    ),
  );
}

/* ── Nextflow orchestration ───────────────────────────────────────────
   The stage table above runs one CLI step at a time. This runs the whole
   DAG as Nextflow sees it: one work directory, one provenance trace, one
   resume point. The command was registered in the runner's whitelist and
   no part of the interface had ever called it, so the orchestrated path —
   the one the architecture actually specifies — was unreachable here.

   Profile matters enough to be a first-class control rather than a
   default: `standard` runs against host tools, `docker` and `singularity`
   use the built images. Running the wrong one silently produces results
   from different tool versions. */
const NF = {
  profile: "conda_free",
  resume: true,
  supply: { alignment: true, rooted_tree: true, atlas: true },
  force_unchecked: false,
};

const NF_PROFILES = [
  ["conda_free", "host tools, no containers"],
  ["docker", "the built images — needs `make containers`"],
  ["singularity", "the built images, rootless"],
  ["standard", "host tools, default resources"],
  ["test", "the synthetic demo corpus"],
];

/* Supplying a pre-computed artifact SKIPS the stage that would rebuild
   it. That is the point: rebuilding a published alignment or tree can
   change it, and nothing downstream would report that it had. */
function nfOptions() {
  const o = { profile: NF.profile, resume: NF.resume };
  const d = S.data;
  if (NF.supply.atlas && d?.identity?.atlas_version) o.atlas = d.paths?.atlas;
  if (NF.supply.alignment && d?.paths?.alignment) o.alignment = d.paths.alignment;
  if (NF.supply.rooted_tree && d?.paths?.rooted_tree) o.rooted_tree = d.paths.rooted_tree;
  if (NF.force_unchecked) o.force_unchecked = true;
  for (const k of Object.keys(o)) if (o[k] == null || o[k] === false) delete o[k];
  return o;
}

async function runWorkflow() {
  const cmd = (S.commands || []).find((c) => c.key === "workflow");
  if (!cmd) return notify("error", "Unavailable", "The runner does not expose the workflow command.");
  const missing = (cmd.requires_tools || []).filter((t) => !(S.env?.tools?.[t]?.present));
  if (missing.length) {
    return notify("error", "Nextflow not available",
      `${missing.join(", ")} not found on PATH. Install it, or run the stages individually.`);
  }
  try {
    const job = await api("/api/run", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ command: "workflow", pathogen: S.pathogen, options: nfOptions() }),
    });
    attach(job.id, { id: "workflow", name: "Nextflow pipeline", cmd: "workflow" });
    $("#console").classList.add("open");
    $("#log-toggle").textContent = "▼ LOG";
    notify("running", "Nextflow started", `profile ${NF.profile}. A closed gate completes normally.`);
  } catch (e) {
    notify("error", "Could not start", e.message);
  }
}

function nfPanel() {
  const present = !!S.env?.tools?.nextflow?.present;
  const row = (label, control, hint) =>
    el("div", { style: "display:grid;grid-template-columns:190px 1fr;gap:14px;align-items:center;padding:7px 0;border-bottom:1px solid var(--hair)" },
      el("span", { style: "font-size:12.5px", text: label }),
      el("div", {}, control, hint ? el("div", { class: "hint", text: hint }) : null));

  const profileSelect = el("select", { class: "select",
    onchange: (e) => { NF.profile = e.target.value; renderWork(); } },
    ...NF_PROFILES.map(([v, why]) => el("option", { value: v, text: `${v} — ${why}`, selected: v === NF.profile })));

  const supplyBox = (key, label) => el("label", { style: "display:flex;gap:7px;align-items:center;font-size:12px;margin-right:16px" },
    el("input", { type: "checkbox", checked: NF.supply[key],
      onchange: (e) => { NF.supply[key] = e.target.checked; renderWork(); } }),
    el("span", { text: label }));

  return el("div", { class: "sec" },
    el("h3", { text: "Nextflow — run the whole DAG" }),
    el("p", { class: "hint", style: "max-width:88ch",
      text: "One work directory, one provenance trace, one resume point. A closed D.H1 gate completes "
          + "normally and exits 0: Stage 5 reports that scoring is blocked and Stage 6 publishes the verdict." }),
    row("Execution profile", profileSelect,
        NF.profile === "docker" || NF.profile === "singularity"
          ? "Needs the images built once with `make containers`."
          : "Runs against whatever tools are on PATH."),
    row("Resume", el("label", { style: "display:flex;gap:7px;align-items:center;font-size:12px" },
      el("input", { type: "checkbox", checked: NF.resume, onchange: (e) => { NF.resume = e.target.checked; } }),
      el("span", { text: "reuse cached task results" })),
      "Nextflow revalidates inputs, so this cannot serve a stale result for changed data."),
    row("Supply existing artifacts",
      el("div", { style: "display:flex;flex-wrap:wrap" },
        supplyBox("atlas", "Atlas"), supplyBox("alignment", "Alignment"), supplyBox("rooted_tree", "Rooted tree")),
      "Supplying one SKIPS the stage that would rebuild it. Rebuilding a published alignment or tree can change it."),
    row("Run Stage 5 unchecked",
      el("label", { style: "display:flex;gap:7px;align-items:center;font-size:12px" },
        el("input", { type: "checkbox", checked: NF.force_unchecked,
          onchange: (e) => { NF.force_unchecked = e.target.checked; renderWork(); } }),
        el("span", { class: NF.force_unchecked ? "st invalid" : "", text: "bypass the D.H1 gate" })),
      NF.force_unchecked
        ? "The result is marked non-authoritative and can never be a surveillance finding."
        : "Leave off. The gate refusing is a correct outcome, not a failure."),
    el("div", { class: "btn-row", style: "margin-top:12px" },
      btn("primary", "Run pipeline (Nextflow)", runWorkflow, !present,
          present ? "nextflow run workflow/main.nf" : "nextflow is not on PATH"),
      btn("tertiary", "Show command", () => sheet("Command",
        el("pre", { style: "font-family:var(--mono);font-size:11.5px;white-space:pre-wrap;margin:0",
          text: "nextflow run workflow/main.nf --pathogen " + S.pathogen + " "
                + Object.entries(nfOptions()).map(([k, v]) =>
                    (k === "profile" ? "-profile " + v : k === "resume" ? "-resume" :
                     v === true ? "--" + k : "--" + k + " " + v)).join(" ") })))),
    present ? null : el("p", { class: "blank" }, el("b", { text: "Nextflow not found." }),
      "Install it to run the orchestrated pipeline, or run the stages individually above."),
  );
}

/* ── 04 · RUN ───────────────────────────────────────────────────── */
function workRun(host) {
  $("#stage-title").textContent = "Execution";
  $("#stage-sub").textContent = "stages run one at a time — they share the Atlas, alignment and ledger";
  $("#stage-tools").replaceChildren(
    btn("primary", "Start analysis", startAll, !S.completed.has("valid"), "Runs every reachable stage in dependency order"),
    btn("secondary", "Run from stage", runFromStage),
    btn("tertiary", "Reset progress", () => { S.completed = new Set(S.data ? ["data"] : []); S.skipped.clear(); renderSpine(); renderWork(); }),
  );

  const table = el("table", { class: "grid" },
    el("thead", {}, el("tr", {}, ...["", "Stage", "Depends on", "Command", "Status", "Controls"].map((h) => el("th", { text: h })))),
    el("tbody", {}, ...STAGES.map((st) => {
      const status = stageStatus(st);
      const blocked = status === "blocked";
      const reason = blockReason(st);
      return el("tr", { class: st.id === S.stage ? "sel" : "" },
        el("td", { class: "mono dim", text: st.ord }),
        el("td", {}, el("b", { text: st.name }), reason ? el("div", { class: "hint", text: reason }) : null),
        el("td", { class: "mono dim", text: st.needs.map((n) => STAGES.find((s) => s.id === n).ord).join(" · ") || "—" }),
        el("td", { class: "mono dim", text: st.cmd || "—" }),
        el("td", {}, el("span", { class: "st " + status, text: status.replace("-", " ") })),
        el("td", {}, el("div", { class: "btn-row" },
          btn("tertiary sm", "Open", () => { S.stage = st.id; setMode(st.mode); }),
          btn("tertiary sm", "Configure", () => setMode("configure")),
          st.cmd ? btn(S.completed.has(st.id) ? "secondary sm" : "primary sm", S.completed.has(st.id) ? "Re-run" : "Run",
            () => runStage(st), blocked, reason || `Run ${st.cmd}`) : el("span"),
          btn("tertiary sm", "Results", () => { S.stage = st.id; setMode(st.mode === "run" ? "visualize" : st.mode); }, !S.completed.has(st.id)),
          btn("tertiary sm", "Skip", () => { S.skipped.add(st.id); mark(); renderSpine(); renderWork(); }, S.completed.has(st.id)),
          btn("tertiary sm", "Reset", () => { S.completed.delete(st.id); S.skipped.delete(st.id); renderSpine(); renderWork(); }),
        )),
      );
    })),
  );
  host.append(
    el("div", { class: "sec" }, el("h3", { text: "Pipeline stages" }), el("div", { class: "scroll-x" }, table)),
    el("div", { style: "height:22px" }),
    nfPanel(),
  );
}

/* ── 05 · VISUALIZE ─────────────────────────────────────────────── */
/* Four views answer the questions almost every session starts with:
   what is the tree, where are the loci, how is sampling distributed over
   time, and what do the tracks show. The other six are real and stay one
   click away — they were simply never all needed at once. */
const VIZ = {
  tree: "Phylogeny", genome: "Genome map", temporal: "Temporal", tracks: "Tracks",
};
const VIZ_MORE = {
  roottotip: "Root-to-tip", ordination: "Ordination", map: "Map",
  matrix: "Lineage × geography", alignment: "Alignment", spectrum: "Mutation spectrum",
};
const ALL_VIZ = { ...VIZ, ...VIZ_MORE };
const CACHE = {};
async function cached(key, url) {
  if (CACHE[key] !== undefined) return CACHE[key];
  try { CACHE[key] = await api(url); } catch (e) { CACHE[key] = { error: e.message }; }
  return CACHE[key];
}
function unavailable(host, why) {
  host.append(el("p", { class: "blank" }, el("b", { text: "View unavailable" }),
    why + " Nothing is drawn in its place — a plausible-looking picture with no data behind it is worse than an empty panel."));
}

function workVisualize(host) {
  if (!S.data) { host.append(el("p", { class: "blank" }, el("b", { text: "No dataset loaded." }), "Load a provisioned pathogen from Data input.")); return; }
  $("#stage-title").textContent = ALL_VIZ[S.viz];
  const tool = (id, label, title) => btn("tool" + (S.tool === id ? " on" : ""), label, () => { S.tool = id; renderWork(); }, false, title);

  // The six secondary views sit in a select rather than six more buttons.
  // It reads as one control instead of six, and it shows which one is active
  // when the active view happens to be one of them.
  const more = el("select", {
    class: "viz-more",
    title: "Further views",
    onchange: (e) => { if (e.target.value) { S.viz = e.target.value; renderWork(); } },
  }, el("option", { value: "", text: S.viz in VIZ_MORE ? ALL_VIZ[S.viz] : "More views…" }),
     ...Object.entries(VIZ_MORE).map(([k, v]) => el("option", { value: k, text: v })));
  if (S.viz in VIZ_MORE) more.classList.add("on");

  $("#stage-tools").replaceChildren(
    el("div", { class: "btn-group" }, ...Object.entries(VIZ).map(([k, v]) =>
      btn("tool" + (S.viz === k ? " on" : ""), v, () => { S.viz = k; renderWork(); }))),
    more,
    // A fixed gap, not a `push`: each plot appends its own controls to
    // #stage-tools after this runs, and margin-left:auto here would throw
    // those to the far right, away from the view buttons they belong to.
    el("span", { style: "width:10px" }),
    el("div", { class: "btn-group" },
      tool("select", "Select", "Click a mark to filter every view"),
      tool("pan", "Pan", "Drag to move the view"),
      tool("lasso", "Lasso", "Drag a region to select many samples"),
    ),
    // Scrolling zooms and double-click fits, so the ＋/− pair was a third way
    // to do what two gestures already do.
    el("div", { class: "btn-group" },
      btn("tool", "Fit", fitActive, false, "Fit to screen (or double-click the plot)"),
      btn("tool", "Reset", () => { clearSel(); fitActive(); }, false, "Reset view and selection"),
    ),
    el("span", { class: "mono", id: "zoom-readout",
                 style: "font-size:10px;color:var(--ink-3);min-width:38px;text-align:right", text: "1.00×" }),
    btn("tool", S.labels ? "Hide labels" : "Show labels", () => { S.labels = !S.labels; renderWork(); }),
    btn("secondary", "Export", exportFigure, false, "Export this figure"),
  );

  ({ tree: plotTree, genome: plotGenome, tracks: plotTracks, roottotip: plotRootToTip,
     ordination: plotOrdination, map: plotMap, temporal: plotTemporal, matrix: plotMatrix,
     alignment: plotAlignment, spectrum: plotSpectrum })[S.viz](host);
  renderRegister();
}

function renderRegister() {
  const f = S.sel;
  const n = activeSamples().size;
  const toks = [];
  for (const l of f.lineages) toks.push(token(l, "lineage", S.colour[l], () => f.lineages.delete(l)));
  for (const c of f.countries) toks.push(token(c, "country", null, () => f.countries.delete(c)));
  if (f.period) toks.push(token(f.period.join("–"), "period", null, () => { f.period = null; }));
  if (f.locus) toks.push(token(f.locus, "locus", null, () => { f.locus = null; }));
  if (f.sample) toks.push(token(f.sample, "genome", null, () => { f.sample = null; }));

  $("#footbar").replaceChildren(
    toks.length ? el("span", { class: "tag", text: "brush" }) : el("span", { class: "hint", text: "No filter — click any mark to brush every view." }),
    ...toks,
    toks.length ? el("span", { class: "token", style: "border-color:var(--ink)" }, el("b", { text: `${n} / ${S.data.samples.length}` })) : null,
    toks.length ? btn("tertiary sm", "Clear", clearSel) : null,
    el("div", { class: "push" },
      btn("tertiary sm", "Filter", filterSheet),
      btn("tertiary sm", "Export data", exportData),
    ),
  );
}

/* ── 06 · INTERPRET ─────────────────────────────────────────────── */
/* ── qualifying evidence ──────────────────────────────────────────────
   Three checks decide how much weight the gate's verdict can carry, and
   all three used to live in log files nobody opens:

     Stage 1.5  recombination — whether ONE tree describes this corpus at
                all. Reconstructing ancestral states across a recombinant
                alignment reconstructs a history that never happened.
     Stage 2    molecular clock — whether calendar time means anything
                here. A weak fit does not invalidate the topology, but it
                does invalidate every statement of the form "over N years".
     Stage 3    variants — how much of the observed variation the Atlas
                loci actually cover.

   They are rendered beside the verdict rather than behind a link,
   because a reader who has to go looking for a caveat will not find it.
   A missing artifact says so; nothing is defaulted or drawn empty.
   ─────────────────────────────────────────────────────────────────── */
function evidenceRow(label, stage, state, headline, detail, extra) {
  const colour = { ok: "--st-complete", warn: "--st-warning", bad: "--st-error", none: "--ink-3" }[state];
  return el("div", { style: "display:grid;grid-template-columns:118px 1fr;gap:16px;padding:11px 0;border-bottom:1px solid var(--hair)" },
    el("div", {},
      el("p", { class: "tag", text: stage }),
      el("p", { style: "font-size:12.5px;font-weight:600;margin-top:2px", text: label })),
    el("div", {},
      el("p", { style: `font-size:12.5px;font-weight:600;color:var(${colour})`, text: headline }),
      detail ? el("p", { class: "mono", style: "font-size:11px;color:var(--ink-3);margin-top:3px", text: detail }) : null,
      extra ? el("p", { style: "font-size:11.5px;color:var(--ink-2);line-height:1.55;margin-top:5px;max-width:88ch", text: extra }) : null),
  );
}

function evidenceChain() {
  const d = S.data, rows = [];
  const r = d.recombination;
  rows.push(r == null
    ? evidenceRow("Recombination", "Stage 1.5", "none", "Not run",
        null, "The screen is mandatory before phylogenetics. Until it has run, nothing downstream of the tree has been qualified.")
    : evidenceRow("Recombination", "Stage 1.5", r.significant ? "bad" : "ok",
        r.significant ? "Recombination detected" : "No significant recombination",
        `PHI  p = ${r.p_value}  ·  ${r.n_sequences} sequences  ·  ${r.n_informative_sites} informative sites  ·  tier ${r.tier}`,
        r.interpretation));

  const c = d.molecular_clock;
  rows.push(c == null
    ? evidenceRow("Molecular clock", "Stage 2", "none", "Not dated", null,
        "No time-scaled tree was built, so no result here depends on calendar time.")
    : evidenceRow("Molecular clock", "Stage 2", c.usable_for_dating ? "ok" : "warn",
        c.usable_for_dating ? "Clock usable for dating" : "Weak temporal signal",
        `rate = ${c.rate} subs/site/yr  ·  r² = ${c.r_squared}` + (c.n_outliers != null ? `  ·  ${c.n_outliers} outliers` : ""),
        c.interpretation));

  const v = d.variants;
  rows.push(v == null
    ? evidenceRow("Variants", "Stage 3", "none", "Not called", null,
        "No variant table, so the share of variation falling inside Atlas loci is unknown.")
    : evidenceRow("Variants", "Stage 3", "ok",
        `${v.n_variants.toLocaleString()} variants across ${v.n_genomes.toLocaleString()} genomes`,
        `mean ${v.mean_per_genome != null ? v.mean_per_genome.toFixed(1) : "—"} per genome  ·  ` +
        `${v.n_in_atlas_loci.toLocaleString()} inside Atlas loci` +
        (v.fraction_in_atlas_loci != null ? ` (${(v.fraction_in_atlas_loci * 100).toFixed(2)}%)` : ""),
        null));

  return el("div", { class: "sec" },
    el("h3", { text: "Qualifying evidence — what the verdict rests on" }), ...rows);
}

/* ── the D.H1 rows themselves ─────────────────────────────────────────
   The gate reports one word. These are the rows it read to get there.
   A verdict carried by a single locus and a verdict carried by twenty
   are indistinguishable until the table is shown, which is why the
   count is stated in words above it rather than left to be counted.
   ─────────────────────────────────────────────────────────────────── */
const DH1_COLOUR = {
  SUPPORTED: "--st-complete",
  /* Evidence AGAINST the hypothesis, not an absence of evidence for it.
     Coloured as the strongest signal on the panel, because it is. */
  SIGNAL_OPPOSITE_DIRECTION: "--st-significant",
  NOT_SUPPORTED: "--st-error",
  SIGNAL_EXPLAINED_BY_GC: "--st-warning",
  INSUFFICIENT_DATA: "--ink-3",
};
const pv = (x) => x == null ? "—" : x < 0.001 ? x.toExponential(1) : x.toFixed(4);
const rate = (x) => x == null ? "—" : x.toFixed(3);

function dh1Panel() {
  const h = S.data.dh1;
  if (!h) {
    return el("div", { class: "sec" }, el("h3", { text: "D.H1 — per-locus results" }),
      el("p", { class: "blank" }, el("b", { text: "The gate has not been run." }),
        "The testing ledger records no D.H1 result for this pathogen. An unrun gate is not a passed gate, and nothing is shown in place of the rows that do not exist yet."));
  }
  const counts = Object.entries(h.verdict_counts).sort((a, b) => b[1] - a[1]);
  const rests = h.rests_on.length;

  return el("div", { class: "sec" },
    el("h3", { text: "D.H1 — per-locus results" }),
    el("p", { class: "mono", style: "font-size:11px;color:var(--ink-3)",
      text: `rule ${h.decision_rule}  ·  α = ${h.alpha}  ·  analysis set = loci in ≥ ${h.min_carriers ?? "—"} genomes  ·  ` +
            `${h.n_tested}/${h.n_loci} reached a p-value  ·  run ${h.timestamp ? h.timestamp.slice(0, 19).replace("T", " ") : "—"}` }),

    el("div", { style: "display:flex;gap:8px;flex-wrap:wrap;margin:4px 0 2px" },
      ...counts.map(([verdict, n]) => el("span", { class: "token",
        style: `border-color:var(${DH1_COLOUR[verdict] || "--ink-3"});color:var(${DH1_COLOUR[verdict] || "--ink-3"})` },
        el("b", { text: String(n) }), " " + verdict))),

    /* A locus MORE disrupted than its control contradicts D.H1. Saying so
       in a sentence matters more than the row it sits on: the table shows
       a very small p-value beside a very clear verdict, and a reader
       skimming for significance will find the first before the second. */
    ...(S.data.dh1.loci.filter((r) => r.verdict === "SIGNAL_OPPOSITE_DIRECTION").map((r) =>
      el("p", { style: "font-size:12.5px;line-height:1.6;max-width:86ch;padding:9px 12px;border-left:2px solid var(--st-significant);color:var(--st-significant)",
        text: `${r.atlas_id} shows a real, GC-adjusted effect (p_fdr = ${pv(r.gc_adjusted_p_fdr)}) that runs AGAINST D.H1: `
            + `it is MORE disrupted than its matched control (${rate(r.locus_rate)} vs ${rate(r.control_rate)}), not less. `
            + `D.H1 predicts lower disruption, so this is evidence against the hypothesis — not a small p-value in its favour, `
            + `and not an absence of evidence. It cannot open the gate.` }))),

    /* The sentence that stops a one-locus result reading like a corpus-wide one. */
    rests === 0 ? null : el("p", {
      style: `font-size:12.5px;line-height:1.6;max-width:86ch;padding:9px 12px;border-left:2px solid var(--st-${rests < 3 ? "warning" : "complete"});color:var(--${rests < 3 ? "st-warning" : "ink-2"})`,
      text: rests === 1
        ? `The SUPPORTED verdict rests on a SINGLE locus, ${h.rests_on[0]}. Every score derived downstream inherits that locus's fragility — check its control-clade count in the table below before reading any number as a corpus-wide finding.`
        : `The SUPPORTED verdict rests on ${rests} loci: ${h.rests_on.join(", ")}.` }),

    el("div", { class: "scroll-x" },
      el("table", { class: "grid" },
        el("thead", {}, el("tr", {}, ...["Locus", "Verdict", "raw p", "GC-adj p (FDR)", "locus rate", "control rate", "Controls", "Flags"]
          .map((head) => el("th", { text: head })))),
        el("tbody", {}, ...h.loci.map((row) => el("tr", {},
          el("td", { class: "mono", text: row.atlas_id }),
          el("td", {}, el("span", { style: `color:var(${DH1_COLOUR[row.verdict] || "--ink-3"});font-weight:600;font-size:11px`, text: row.verdict })),
          /* Empty, never 0. A locus halted at the floor has no p-value,
             and 0.0 is a p-value a reader would act on. */
          el("td", { class: "num", text: pv(row.raw_p) }),
          el("td", { class: "num", style: row.gc_adjusted_p_fdr != null && row.gc_adjusted_p_fdr < (h.alpha ?? 0.05) ? "font-weight:600" : "", text: pv(row.gc_adjusted_p_fdr) }),
          el("td", { class: "num", text: rate(row.locus_rate) }),
          el("td", { class: "num", text: rate(row.control_rate) }),
          /* Control provenance. Which regions a locus was compared
             against turned out to be the most consequential choice in the
             whole test, and for every run before R-20 it was recorded
             nowhere — so a run that lacks it says so rather than showing
             a blank that reads as "none". */
          el("td", { style: "font-size:10.5px" },
            row.n_controls == null
              ? el("span", { class: "st skipped", title: "this run predates control provenance being recorded", text: "not recorded" })
              : el("span", {
                  class: row.n_controls_same_compartment === row.n_controls ? "st" : "st invalid",
                  title: row.control_regions.join("  ") || "no control regions recorded",
                  text: `${row.n_controls_same_compartment}/${row.n_controls} same compartment`,
                }),
          ),
          el("td", { style: "font-size:10.5px" },
            row.underpowered ? el("span", { class: "st invalid", text: "underpowered" }) : null,
            !row.minimum_data_passed ? el("span", { class: "st skipped", title: row.failing_checks.join(", "),
              text: "floor: " + (row.failing_checks.join(", ") || "failed") }) : null),
        )))),
    ),
    el("p", { class: "hint", text: `Read from ${h.ledger_path}. Nothing on this panel is recomputed — these are the rows the gate itself read, so it cannot disagree with the verdict above.` }),
  );
}

function workInterpret(host) {
  if (!S.data) { host.append(el("p", { class: "blank", text: "No dataset." })); return; }
  const g = S.data.gate, floor = S.data.floor;
  $("#stage-title").textContent = "Interpretation";
  $("#stage-sub").textContent = `ledger ${g.n_ledger_rows} rows · ${g.latest_run ? g.latest_run.slice(0, 10) : "no run"}`;
  $("#stage-tools").replaceChildren(
    btn("tertiary", "View parameters", () => setMode("configure")),
    btn("tertiary", "View log", () => $("#console").classList.add("open")),
    btn("secondary", "Export provenance", exportProvenance),
  );

  host.append(
    el("div", { style: `border:1px solid var(--st-${g.permitted ? "complete" : "error"});padding:16px 18px;margin-bottom:20px` },
      el("p", { class: "tag", style: `color:var(--st-${g.permitted ? "complete" : "error"})`, text: g.permission }),
      el("p", { style: `font-size:21px;font-weight:600;margin-top:5px;color:var(--st-${g.permitted ? "complete" : "error"})`, text: g.permitted ? "Scoring permitted" : "Scoring blocked" }),
      el("p", { style: "font-size:12.5px;color:var(--ink-2);line-height:1.6;margin-top:9px;max-width:78ch", text: g.explanation }),
    ),
    el("div", { class: "sec" }, el("h3", { text: "Appendix C minimum-data floor" }),
      ...Object.entries(floor).filter(([, r]) => r.value != null).map(([k, r]) => {
        const frac = r.unit === "frac", pass = r.value >= r.floor;
        return el("div", { style: "display:grid;grid-template-columns:210px 1fr 130px;gap:14px;align-items:center;padding:7px 0;border-bottom:1px solid var(--hair)" },
          el("span", { class: "mono", style: `font-size:11.5px;color:var(--${pass ? "ink-2" : "st-error"})`, text: k }),
          el("span", { class: "meter", style: "height:9px" },
            el("i", { class: pass ? "ok" : "bad", style: `width:${Math.min(r.value / (r.floor * 1.6), 1) * 100}%` }),
            el("span", { class: "floor", style: "left:62.5%" })),
          el("span", { class: "mono", style: `font-size:11px;text-align:right;color:var(--${pass ? "ink-2" : "st-error"})`,
            text: frac ? `${(r.value * 100).toFixed(1)}% / ${(r.floor * 100).toFixed(0)}%` : `${r.value} / ${r.floor}` }),
        );
      }),
    ),
    el("div", { style: "height:22px" }),
    dh1Panel(),
    el("div", { style: "height:22px" }),
    evidenceChain(),
  );
}

/* ── surveillance output ──────────────────────────────────────────────
   The scores, the control limits and the warning level: the end of the
   chain, and until now unreachable from this interface. /api/stage5
   existed and nothing ever called it, so a reader could see the corpus,
   the tree and the gate but never what the pipeline actually computed.

   A closed gate returns 409, which is a correct outcome and is rendered
   as one. ─────────────────────────────────────────────────────────── */
async function loadStage5() {
  if (STAGE5.state === "loading") return;
  STAGE5.state = "loading";
  renderWork();
  try {
    STAGE5.data = await api(`/api/stage5/${S.pathogen}`);
    STAGE5.state = "ok";
  } catch (e) {
    STAGE5.blocked = String(e.message || e);
    STAGE5.state = "blocked";
  }
  renderWork();
}
const STAGE5 = { state: "idle", data: null, blocked: "" };

/* A limit fitted to a short baseline is a working figure, not a
   calibrated false-alarm rate. The interval is rendered at the same
   weight as the limit so the two cannot be read apart. */
function limitPanel(label, chart) {
  if (!chart || chart.error) {
    return el("div", { style: "padding:10px 0;border-bottom:1px solid var(--hair)" },
      el("p", { class: "tag", text: label }),
      el("p", { style: "font-size:12px;color:var(--st-error)", text: chart?.error || "not calibrated" }));
  }
  const short = chart.short_baseline;
  const iv = chart.control_limit_interval;
  return el("div", { style: "padding:10px 0;border-bottom:1px solid var(--hair)" },
    el("p", { class: "tag", text: label }),
    el("div", { style: "display:flex;gap:22px;align-items:baseline;flex-wrap:wrap;margin-top:3px" },
      el("span", { style: "font-size:19px;font-weight:600", text: `h = ${chart.control_limit}` }),
      iv ? el("span", { class: "mono", style: `font-size:12px;color:var(--st-warning)`,
        text: `resampled 5–95%: ${iv[0]} – ${iv[1]}` }) : null,
      el("span", { class: "mono", style: "font-size:11px;color:var(--ink-3)",
        text: `${chart.n_alarms} alarm(s) · baseline ${chart.baseline_windows}`
              + (chart.monitored_windows != null ? ` · monitored ${chart.monitored_windows}` : "") }),
    ),
    short ? el("p", {
      style: "font-size:11.5px;line-height:1.55;max-width:88ch;margin-top:6px;padding-left:10px;border-left:2px solid var(--st-warning);color:var(--st-warning)",
      text: chart.caveat || "Short baseline: the nominal ARL is not achieved." }) : null,
  );
}

function sparkline(series, key) {
  const values = series.map((s) => s[key]).filter((v) => v != null);
  if (values.length < 2) return el("span", { class: "hint", text: "—" });
  const lo = Math.min(...values), hi = Math.max(...values), span = hi - lo || 1;
  const W = 260, H = 34;
  const pts = values.map((v, i) =>
    `${(i / (values.length - 1)) * W},${H - ((v - lo) / span) * H}`).join(" ");
  return svg("svg", { width: W, height: H, viewBox: `0 0 ${W} ${H}` },
    svg("polyline", { points: pts, fill: "none", stroke: "var(--ink)", "stroke-width": 1.4 }));
}

function workSurveillance(host) {
  $("#stage-title").textContent = "Surveillance output";
  $("#stage-tools").replaceChildren(
    btn("secondary", STAGE5.state === "loading" ? "Running…" : "Run Stage 5", loadStage5,
        STAGE5.state === "loading", "Metrics → scores → control charts → warning level"),
  );

  if (STAGE5.state === "idle") {
    host.append(el("p", { class: "blank" }, el("b", { text: "Not run." }),
      "Stage 5 computes the seven G.2 terms, fits weights, scores each window and calibrates the control charts. It is not run automatically because it reads every genome in the corpus."));
    return;
  }
  if (STAGE5.state === "loading") {
    host.append(el("p", { class: "blank" }, el("b", { text: "Running…" }),
      "Annotating tip states over the scoring-eligible Atlas loci, then scoring each window."));
    return;
  }
  if (STAGE5.state === "blocked") {
    host.append(el("div", { style: "border:1px solid var(--st-error);padding:16px 18px" },
      el("p", { class: "tag", style: "color:var(--st-error)", text: "SCORING BLOCKED" }),
      el("p", { style: "font-size:12.5px;color:var(--ink-2);line-height:1.6;margin-top:8px;max-width:84ch",
        text: STAGE5.blocked }),
      el("p", { class: "hint", style: "margin-top:10px",
        text: "This is a successful run, not an error. The gate refused, and no score was produced." })));
    return;
  }

  const d = STAGE5.data;
  const det = d.detection || {};
  host.append(
    d.authoritative ? null : el("p", {
      style: "font-size:12.5px;line-height:1.6;max-width:88ch;padding:9px 12px;margin-bottom:16px;border-left:2px solid var(--st-warning);color:var(--st-warning)",
      text: "NON-AUTHORITATIVE. These numbers demonstrate the machinery on this corpus. They are not a surveillance finding and no alert level derived from them is actionable." }),

    el("div", { class: "sec" }, el("h3", { text: "Score series" }),
      el("div", { style: "display:flex;gap:30px;align-items:center;flex-wrap:wrap" },
        el("div", {}, el("p", { class: "tag", text: "G4-EWS core (M3)" }),
          sparkline(d.score_series || [], "g4_ews_core")),
        el("div", {}, el("p", { class: "tag", text: "Integrated (M4)" }),
          sparkline(d.score_series || [], "integrated_score"))),
      el("p", { class: "hint", text: `${(d.score_series || []).length} scored windows.` })),

    el("div", { style: "height:20px" }),
    el("div", { class: "sec" }, el("h3", { text: "Detection — control limits" }),
      limitPanel("CUSUM", det.cusum), limitPanel("EWMA", det.ewma)),

    el("div", { style: "height:20px" }),
    el("div", { class: "sec" }, el("h3", { text: "Chain" }),
      ...(d.steps || []).map((s) => el("div", {
        style: "display:grid;grid-template-columns:150px 190px 1fr;gap:12px;padding:5px 0;border-bottom:1px solid var(--hair)" },
        el("span", { class: "mono", style: "font-size:11px;color:var(--ink-3)", text: s.step }),
        el("span", { class: "mono", style: `font-size:11px;color:var(--${/ok|l2_/.test(s.status) ? "st-complete" : "st-warning"})`, text: s.status }),
        el("span", { style: "font-size:11.5px;color:var(--ink-2)", text: String(s.detail || "").slice(0, 200) })))),
  );
}

/* ── 07 · REPORT ────────────────────────────────────────────────── */
const SECTIONS = ["Dataset summary", "Quality control", "Genomic analysis", "Phylogenetic analysis",
  "Evolutionary analysis", "Statistical results", "Key findings", "Methods", "Parameters", "Reproducibility"];
const chosenSections = new Set(SECTIONS);

function workReport(host) {
  $("#stage-title").textContent = "Report";
  $("#stage-sub").textContent = S.data && !S.data.gate.permitted ? "scored sections are withheld — the gate is closed" : "";
  $("#stage-tools").replaceChildren(
    btn("primary", "Generate report", generateReport),
    btn("secondary", "Preview", previewReport),
    btn("tertiary", "Select sections", () => renderWork()),
  );
  host.append(
    el("div", { class: "cols two" },
      el("div", { class: "sec" }, el("h3", { text: "Sections" }),
        ...SECTIONS.map((s) => el("label", { class: "check" },
          el("input", { type: "checkbox", checked: chosenSections.has(s), onchange: (e) => { e.target.checked ? chosenSections.add(s) : chosenSections.delete(s); mark(); } }),
          el("span", {}, el("b", { text: s })))),
      ),
      el("div", { class: "sec" }, el("h3", { text: "Output" }),
        el("div", { class: "btn-row" },
          btn("secondary", "Export PDF", () => notImplemented("PDF export needs a headless renderer; not installed.")),
          btn("secondary", "Export HTML", () => generateReport("html")),
          btn("secondary", "Export Word", () => notImplemented("DOCX export needs python-docx; not installed.")),
        ),
        el("h3", { text: "Reproducibility", style: "margin-top:14px" }),
        el("div", { class: "btn-row" },
          btn("tertiary", "Save analysis state", () => saveProject(false)),
          btn("tertiary", "Export configuration", exportConfig),
          btn("tertiary", "Export provenance", exportProvenance),
          btn("tertiary", "Tool versions", showProvenance),
        ),
      ),
    ),
  );
}

/* ═══ helpers: controls ═══════════════════════════════════════════ */
function btn(kind, label, onclick, disabled = false, title = "") {
  return el("button", {
    class: "btn " + kind, text: label, title: title || label,
    disabled: disabled || undefined,
    onclick: (e) => { if (!disabled && onclick) onclick(e); },
  });
}
function field(label, control, hint) {
  return el("div", { class: "field" }, el("label", { text: label }), control, hint ? el("span", { class: "hint", text: hint }) : null);
}
function input(key, value) {
  return el("input", { class: "input", type: "text", value: value ?? S.params[key] ?? "",
    oninput: (e) => { S.params[key] = e.target.value; mark(); } });
}
function stepper(key, min, max, step) {
  const box = el("input", { class: "input", type: "number", value: S.params[key], min, max, step,
    oninput: (e) => { S.params[key] = Number(e.target.value); mark(); } });
  const bump = (d) => () => {
    const v = Math.min(max, Math.max(min, Number((S.params[key] + d * step).toFixed(6))));
    S.params[key] = v; box.value = v; mark();
  };
  return el("div", { class: "stepper" },
    el("button", { type: "button", text: "−", onclick: bump(-1), title: `Decrease by ${step}` }), box,
    el("button", { type: "button", text: "+", onclick: bump(1), title: `Increase by ${step}` }));
}
function select(key, options) {
  return el("select", { class: "select", onchange: (e) => { S.params[key] = e.target.value; mark(); } },
    ...options.map(([v, l]) => el("option", { value: v, selected: S.params[key] === v || undefined, text: l })));
}
function check(key, label, hint) {
  return el("label", { class: "check" },
    el("input", { type: "checkbox", checked: S.params[key] || undefined, onchange: (e) => { S.params[key] = e.target.checked; mark(); } }),
    el("span", {}, el("b", { text: label }), hint ? el("div", { class: "hint", text: hint }) : null));
}
function toggle(key, label, hint) {
  return el("label", { class: "toggle" },
    el("span", {}, el("span", { class: "lbl", text: label }), hint ? el("div", { class: "hint", text: hint }) : null),
    el("span", { style: "display:flex;align-items:center;gap:7px" },
      el("span", { class: "state", text: S.params[key] ? "ON" : "OFF" }),
      el("input", { type: "checkbox", checked: S.params[key] || undefined,
        onchange: (e) => { S.params[key] = e.target.checked; mark(); renderWork(); } }),
      el("span", { class: "track" })));
}
function token(label, kind, colour, remove) {
  return el("span", { class: "token" },
    colour ? el("i", { class: "sw", style: `background:${colour}` }) : null,
    el("span", { class: "dim", text: kind }), el("b", { text: label }),
    el("button", { text: "×", onclick: () => { remove(); renderWork(); } }));
}

/* ═══ actions ═════════════════════════════════════════════════════ */
async function validateOne(path) {
  try {
    const v = await api(`/api/validate?path=${encodeURIComponent(path)}`);
    S.validations[path] = v;
    log(v.status === "invalid" ? "error" : v.status === "warning" ? "warning" : "success",
      `validate ${v.name}: ${v.status.toUpperCase()}${v.records != null ? ` (${v.records} ${v.record_label})` : ""}`);
    renderWork();
    return v;
  } catch (e) { notify("error", "Validation failed", e.message); }
}
async function validateAll() {
  const targets = S.inputs.filter((f) => ["fasta", "tsv", "csv", "newick", "vcf", "fastq"].includes(f.format)).slice(0, 25);
  notify("running", "Validating", `${targets.length} files…`);
  for (const f of targets) await validateOne(f.path);
  const bad = Object.values(S.validations).filter((v) => v.status === "invalid").length;
  const warn = Object.values(S.validations).filter((v) => v.status === "warning").length;
  notify(bad ? "error" : warn ? "warning" : "success", "Validation complete",
    `${targets.length} files · ${bad} invalid · ${warn} with warnings`);
  setMode("validate");
}
function precheck() {
  const missing = ["mafft", "iqtree2", "treetime"].filter((t) => !S.env.tools[t]?.present);
  notify(missing.length ? "warning" : "success", "Precheck",
    missing.length ? `${missing.join(", ")} missing — stages 04 and 06 cannot run.` : "All external tools present.");
}
function compatibility() {
  const rows = Object.values(S.validations);
  const aln = rows.find((v) => v.format === "fasta" && v.detail.length_min === v.detail.length_max);
  const tree = rows.find((v) => v.format === "newick");
  const meta = rows.find((v) => v.format === "tsv");
  const msgs = [];
  msgs.push(aln ? `Alignment present (${aln.records} seqs × ${aln.detail.length_min} nt).` : "No aligned FASTA found — phylogenetics will need one.");
  msgs.push(tree ? `Tree present (${tree.records} tips).` : "No Newick tree found.");
  if (aln && tree && aln.records !== tree.records) msgs.push(`MISMATCH: alignment has ${aln.records} sequences but the tree has ${tree.records} tips.`);
  msgs.push(meta ? `Metadata table present (${meta.records} rows).` : "No metadata table found.");
  sheet("Compatibility check", el("div", { class: "stack tight" }, ...msgs.map((m) => el("p", { style: "font-size:12.5px;line-height:1.6", text: m }))));
}
function autoFix(v, f) { notImplemented(`Auto-fix "${f.fix}" is defined in the UI but has no server implementation yet. It would rewrite ${v.name}, so it is not stubbed.`); }
function removeInvalid() {
  const bad = Object.values(S.validations).filter((v) => v.status === "invalid");
  if (!bad.length) return notify("success", "Nothing to remove", "No invalid files.");
  confirmSheet("Remove invalid files", `${bad.length} file(s) failed validation. Removing clears them from this session only — nothing is deleted from disk.`, () => {
    bad.forEach((v) => delete S.validations[v.path]);
    renderWork(); notify("success", "Removed", `${bad.length} file(s) cleared from the session.`);
  });
}
function downloadValidation() {
  const rows = ["file\tformat\tstatus\trecords\tfindings",
    ...Object.values(S.validations).map((v) => [v.name, v.format, v.status, v.records ?? "",
      v.findings.map((f) => `${f.level}:${f.message}`).join(" | ")].join("\t"))];
  download(new Blob([rows.join("\n")], { type: "text/tab-separated-values" }), "g4_validation_report.tsv");
}

const SIM_LOG = {
  data: ["indexing corpus…", "420 records read", "metadata table joined", "corpus registered"],
  valid: ["parsing FASTA…", "420 sequences, all 9,450 nt", "alphabet check: clean", "metadata: 98.8% complete", "validation PASSED"],
  prep: ["aligning 420 sequences (simulated MAFFT)…", "FFT-NS-2, gap-open 1.53", "alignment length 9,450", "384 of 420 retained after gap filter", "alignment written"],
  phylo: ["building ML tree (simulated IQ-TREE)…", "ModelFinder: GTR+F+I+G4", "1000 ultrafast bootstraps", "log-likelihood -48213.77", "rooted tree written (420 tips)"],
  interp: ["loading D.H1 ledger…", "6 loci, 5 SUPPORTED", "gate = PERMITTED", "surveillance score computed"],
  report: ["collecting sections…", "10 sections selected", "figures rendered", "report assembled"],
};

/* A simulation, and it says so on every line. It exists so the full
   pipeline can be walked end to end for review; it computes nothing. */
async function runSimulated(st) {
  const lines = SIM_LOG[st.id] || ["running…", "done"];
  S.job = { id: "sim-" + st.id, command: st.cmd || st.id, title: st.name, state: "running",
            argv: [`[simulated] ${st.cmd || st.id}`], duration: 0, exit_code: null };
  renderExec();
  $("#console").classList.add("open"); $("#log-toggle").textContent = "▼ LOG";
  log("meta", `SIMULATED STAGE — ${st.name}. Demo data; no computation is performed.`);
  const t0 = Date.now();
  for (const line of lines) {
    await new Promise((r) => setTimeout(r, 380));
    log("info", line);
    S.job.duration = (Date.now() - t0) / 1000;
    renderTelemetry();
  }
  S.job.state = "succeeded"; S.job.exit_code = 0;
  S.completed.add(st.id); mark();
  log("success", `${st.name} complete (simulated) in ${S.job.duration.toFixed(1)}s`);
  renderExec(); renderSpine(); renderWork();
  notify("success", `${st.name} complete`, "Simulated stage — demo data only, nothing was computed.");
}

async function runStage(st) {
  if (simulated(st)) return runSimulated(st);
  if (!st.cmd) return notify("warning", "No command", `${st.name} has no executable step in this build.`);
  try {
    const job = await api("/api/run", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ command: st.cmd, pathogen: S.pathogen, options: st.cmd === "dh1" ? { no_ledger: true } : {} }),
    });
    S.stage = st.id;
    attach(job.id, st);
    $("#console").classList.add("open");
    $("#log-toggle").textContent = "▼ LOG";
  } catch (e) { notify("error", `${st.name} could not start`, e.message); }
}
/* A stage with no executable command has NOT run. It used to be marked
   complete and the spine turned green: "Pipeline complete" was reported
   having built no alignment and no tree. That looked correct only because
   those artifacts already existed on disk from earlier CLI runs — on a
   fresh pathogen it sailed past both and scored whatever was lying there.

   Preprocessing and Phylogenetics are the two: they shell out to MAFFT,
   IQ-TREE and TreeTime, which the workstation does not drive. The honest
   states are "supplied" (the artifact exists, someone else made it) and
   "cannot run here" — never "complete". */
function suppliedArtifact(stage) {
  const d = S.data;
  if (!d) return null;
  if (stage.id === "prep") return d.tracks_available?.alignment ?? (d.identity?.n_samples > 0);
  if (stage.id === "phylo") return !!d.tree;
  return null;
}

async function startAll() {
  const todo = STAGES.filter((s) => !S.completed.has(s.id) && !S.skipped.has(s.id) && stageStatus(s) !== "blocked");
  if (!todo.length) return notify("success", "Nothing to run", "Every reachable stage is complete.");
  notify("running", "Running pipeline", `${todo.length} stage${todo.length !== 1 ? "s" : ""}, in dependency order.`);
  const supplied = [];
  for (const st of todo) {
    if (simulated(st)) { await runSimulated(st); continue; }
    if (!st.cmd) {
      const present = suppliedArtifact(st);
      if (present === false) {
        log("error", `${st.name}: no executable step in this build and no artifact on disk.`);
        notify("error", "Pipeline halted", `${st.name} cannot run here and its output is missing. `
          + `Build it with the Nextflow workflow, or supply it under data/, then rescan.`);
        return;
      }
      if (present === true) {
        supplied.push(st.name);
        S.completed.add(st.id);
        log("meta", `${st.name}: not run — artifact already present, taken as supplied.`);
        renderSpine();
        continue;
      }
      // Nothing to run and nothing to check: leave it alone rather than
      // claiming either way.
      log("meta", `${st.name}: no executable step in this build; skipped, not completed.`);
      S.skipped.add(st.id);
      renderSpine();
      continue;
    }
    await runOnce(st);
    if (!S.completed.has(st.id)) {
      notify("warning", "Pipeline halted", `${st.name} did not complete; later stages were not started.`);
      return;
    }
  }
  notify("success", "Pipeline finished",
    "Every runnable stage finished."
    + (supplied.length ? ` ${supplied.join(" and ")} were NOT run — existing artifacts were used.` : "")
    + " Open 06 Interpret for the verdict.");
  setMode("interpret");
}

/* Await one real stage to completion, so Run-all can sequence them. */
function runOnce(st) {
  return new Promise((resolve) => {
    runStage(st);
    const poll = setInterval(() => {
      if (!S.job || ["succeeded", "failed", "gate_closed", "cancelled"].includes(S.job.state)) {
        clearInterval(poll); setTimeout(resolve, 120);
      }
    }, 250);
  });
}
function runFromStage() {
  const options = STAGES.filter((s) => s.cmd);
  sheet("Run from stage", el("div", { class: "stack tight" },
    ...options.map((s) => el("div", { style: "display:flex;align-items:center;gap:10px;justify-content:space-between;padding:6px 0;border-bottom:1px solid var(--hair)" },
      el("span", {}, el("b", { text: `${s.ord} ${s.name}` }), el("div", { class: "hint mono", text: s.cmd })),
      btn("primary sm", "Run", () => { closeSheet(); runStage(s); }, stageStatus(s) === "blocked", blockReason(s))))));
}

/* execution controls */
function renderExec() {
  const j = S.job;
  const state = j ? j.state : "not-started";
  const live = ["queued", "running", "paused"].includes(state);
  $("#exec-state").className = "st " + ({ succeeded: "completed", gate_closed: "warning", cancelled: "skipped" })[state] ?? state;
  $("#exec-state").textContent = state.replace("_", " ");
  $("#exec-cmd").textContent = j ? j.argv.join(" ") : "";
  $("#exec-controls").replaceChildren(
    btn("primary sm", "Start", startAll, live),
    btn("secondary sm", "Pause", () => control("pause"), state !== "running"),
    btn("secondary sm", "Resume", () => control("resume"), state !== "paused"),
    btn("danger sm", "Stop", () => control("stop"), !live, "Graceful SIGTERM"),
    btn("danger sm", "Force stop", () => confirmSheet("Force stop", "SIGKILL terminates the process immediately. Partial output may be left on disk.", () => control("force-stop")), !live, "Immediate SIGKILL"),
    btn("tertiary sm", "Restart", () => { const st = STAGES.find((s) => s.cmd === j?.command); if (st) runStage(st); }, !j || live),
    btn("tertiary sm", "Retry failed", () => { const st = STAGES.find((s) => s.cmd === j?.command); if (st) runStage(st); }, !j || j.state !== "failed"),
  );
  renderTelemetry();
}
async function control(action) {
  if (!S.job) return;
  try {
    const r = await api(`/api/jobs/${S.job.id}/${action}`, { method: "POST" });
    log("info", `${action} → ${r.ok ? "ok" : "refused"} (state ${r.state})`);
    S.job.state = r.state; renderExec();
  } catch (e) { notify("error", "Control failed", e.message); }
}

function attach(jobId, stage) {
  if (S.stream) { S.stream.close(); S.stream = null; }
  S.lines = []; S.lineNo = 0;
  const es = new EventSource(`/api/jobs/${jobId}/stream`);
  S.stream = es;
  api(`/api/jobs/${jobId}`).then((j) => { S.job = j; renderExec(); });
  es.onmessage = (ev) => {
    const d = JSON.parse(ev.data);
    if (d.type === "line") log(d.stream === "stderr" ? "error" : d.stream === "meta" ? "meta" : "info", d.text);
    else if (d.type === "end") {
      S.job = d; es.close(); S.stream = null;
      const ok = d.state === "succeeded";
      if (ok && stage) { S.completed.add(stage.id); mark(); }
      renderExec(); renderSpine(); renderWork();
      notify(ok ? "success" : d.state === "gate_closed" ? "warning" : "error",
        `${stage ? stage.name : "Job"} ${d.state.replace("_", " ")}`,
        d.state === "gate_closed"
          ? "Exit 3 — the D.H1 gate is closed. This is a correct outcome, not a failure."
          : `Exit ${d.exit_code} after ${clock(d.duration)}.`,
        d.state !== "succeeded" ? [["View log", () => $("#console").classList.add("open")], ["Retry", () => stage && runStage(stage)]] : null);
      refreshDataset();
    }
  };
  es.onerror = () => { es.close(); S.stream = null; };
  const tick = setInterval(() => {
    if (!S.job || !["running", "queued", "paused"].includes(S.job.state)) return clearInterval(tick);
    if (S.job.state === "running") S.job.duration = (S.job.duration || 0) + 0.5;
    renderTelemetry();
  }, 500);
}
async function refreshDataset() {
  try { S.data = await api(`/api/dataset/${S.pathogen}`); renderIdentity(); renderCtx(); } catch { /* dataset unchanged */ }
}

/* ═══ LOG ═════════════════════════════════════════════════════════ */
function log(level, text) {
  S.lineNo += 1;
  S.lines.push({ n: S.lineNo, level, text });
  if (S.lines.length > 4000) S.lines.shift();
  renderLog(true);
}
function renderLog(tail) {
  const host = $("#log");
  const show = S.logFilter === "all" ? S.lines : S.lines.filter((l) => l.level === S.logFilter);
  if (!show.length) {
    host.replaceChildren(el("p", { class: "blank", text: "No log output. Run a stage and its stdout and stderr stream here line by line." }));
    return;
  }
  host.replaceChildren(...show.slice(-1200).map((l) => el("div", { class: "ln " + l.level },
    el("span", { class: "n", text: String(l.n) }),
    el("span", { class: "lvl", text: l.level === "meta" ? "info" : l.level }),
    el("span", { class: "t", text: l.text }))));
  if (tail) host.scrollTop = host.scrollHeight;
}

/* ═══ CONTEXT PANEL ═══════════════════════════════════════════════ */
const CTX = { params: "Params", results: "Results", observe: "Observations" };
function renderCtx() {
  $("#ctx-tabs").replaceChildren(...Object.entries(CTX).map(([k, v]) =>
    el("button", { class: S.ctxTab === k ? "on" : "", text: v, onclick: () => { S.ctxTab = k; renderCtx(); } })));
  const body = $("#ctx-body");
  body.replaceChildren();
  if (S.ctxTab === "observe") {
    const list = S.data?.observations || [];
    if (!list.length) return body.append(el("p", { class: "hint", text: "No observations — load a dataset." }));
    for (const o of list) {
      body.append(el("button", { class: "obs " + o.severity, onclick: () => focusObs(o) },
        el("span", { class: "sev" }),
        el("div", {}, el("p", { class: "dm", text: `${o.severity} · ${o.domain}` }),
          el("p", { class: "t", text: o.title }), el("p", { class: "d", text: o.detail }),
          el("p", { class: "r", text: o.rule }))));
    }
  } else if (S.ctxTab === "params") {
    body.append(el("p", { class: "tag", text: "Active parameters" }),
      ...Object.entries(S.params).map(([k, v]) => el("div", { style: "display:flex;justify-content:space-between;gap:10px;padding:5px 0;border-bottom:1px solid var(--hair);font-size:11.5px" },
        el("span", { class: "dim mono", text: k }), el("span", { class: "mono", text: String(v) }))),
      btn("tertiary sm", "Edit in Configure", () => setMode("configure")));
  } else {
    body.append(el("p", { class: "tag", text: "Completed stages" }),
      ...STAGES.filter((s) => S.completed.has(s.id)).map((s) => el("div", { style: "display:flex;justify-content:space-between;padding:6px 0;border-bottom:1px solid var(--hair)" },
        el("span", { text: s.name }), el("span", { class: "st completed" }))),
      S.completed.size <= 1 ? el("p", { class: "hint", text: "Nothing has been run in this session yet." }) : null,
      el("div", { class: "btn-row", style: "margin-top:10px" },
        btn("tertiary sm", "View log", () => $("#console").classList.add("open")),
        btn("tertiary sm", "Export", exportData)));
  }
}
function focusObs(o) {
  const f = o.focus || {};
  if (f.type === "lineage" && f.value) { S.sel.lineages = new Set([f.value]); S.viz = "temporal"; setMode("visualize"); }
  else if (f.type === "locus") { S.sel.locus = f.value; S.viz = "genome"; setMode("visualize"); }
  else if (f.type === "period") { S.sel.period = f.value; S.viz = "temporal"; setMode("visualize"); }
  else if (f.type === "country") { S.sel.countries = new Set([f.value]); S.viz = "matrix"; setMode("visualize"); }
  else setMode("interpret");
}

/* ═══ selection ═══════════════════════════════════════════════════ */
function activeSamples() {
  const out = new Set();
  if (!S.data) return out;
  S.data.samples.forEach((s, i) => {
    const f = S.sel;
    if (f.lineages.size && !f.lineages.has(s.l)) return;
    if (f.countries.size && !f.countries.has(s.c)) return;
    if (f.period && (!s.y || s.y < f.period[0] || s.y > f.period[1])) return;
    if (f.sample && s.a !== f.sample) return;
    out.add(i);
  });
  return out;
}
const filtering = () => S.sel.lineages.size || S.sel.countries.size || S.sel.period || S.sel.sample;
function clearSel() { S.sel = { lineages: new Set(), countries: new Set(), period: null, locus: null, sample: null }; renderWork(); }
function zoom(f) { S.zoom = Math.max(0.4, Math.min(8, (S.zoom || 1) * f)); renderWork(); }

/* ═══ PLOTS ═══════════════════════════════════════════════════════ */
function plotTree(host) {
  const t = S.data.tree;
  if (!t) return host.append(el("p", { class: "blank" }, el("b", { text: "No tree artifact." }), "Stage 06 has not produced a rooted tree for this pathogen."));
  $("#stage-sub").textContent = `${t.n_leaves} tips · ${S.treeLayout} · x = divergence`;
  $("#stage-tools").append(
    el("div", { class: "btn-group" },
      ...[["rectangular", "Rect"], ["circular", "Circ"], ["unrooted", "Unroot"]].map(([k, l]) =>
        btn("tool" + (S.treeLayout === k ? " on" : ""), l, () => { S.treeLayout = k; renderWork(); }))),
    btn("tool", S.bootstrapShown ? "Hide support" : "Show support", () => { S.bootstrapShown = !S.bootstrapShown; renderWork(); }),
    btn("tool", "Export Newick", () => notImplemented("Newick export needs the source tree file; use Artifacts.")),
  );

  const cv = el("canvas"); host.append(cv);
  const dpr = devicePixelRatio || 1, w = host.clientWidth, h = host.clientHeight;
  cv.width = w * dpr; cv.height = h * dpr;
  const g = cv.getContext("2d"); g.scale(dpr, dpr);
  const M = { l: 20, r: 120, t: 12, b: 12 };
  const chosen = activeSamples(), dim = filtering(), veil = parseFloat(css("--veil"));
  const nodes = t.nodes;
  const Z = S.zoom || 1;

  if (S.treeLayout === "rectangular") {
    const X = (x) => M.l + x * (w - M.l - M.r) * Z;
    const Y = (y) => M.t + y * (h - M.t - M.b);
    const kids = new Map();
    for (const n of nodes) if (n.p != null) { if (!kids.has(n.p)) kids.set(n.p, []); kids.get(n.p).push(n); }
    g.lineWidth = 1;
    for (const n of nodes) {
      if (n.p == null) continue;
      const lit = n.s < 0 || chosen.has(n.s);
      g.globalAlpha = dim && !lit ? veil : n.s < 0 ? 0.5 : 0.9;
      g.strokeStyle = n.s >= 0 ? S.colour[S.data.samples[n.s].l] : css("--ink-4");
      g.beginPath(); g.moveTo(X(nodes[n.p].x), Y(n.y)); g.lineTo(X(n.x), Y(n.y)); g.stroke();
    }
    for (const [pi, cs] of kids) {
      const ys = cs.map((c) => Y(c.y));
      g.globalAlpha = 0.5; g.strokeStyle = css("--ink-4");
      g.beginPath(); g.moveTo(X(nodes[pi].x), Math.min(...ys)); g.lineTo(X(nodes[pi].x), Math.max(...ys)); g.stroke();
    }
    g.globalAlpha = 1;
    const bar = w - M.r + 14;
    for (const n of nodes) {
      if (n.s < 0) continue;
      g.globalAlpha = dim && !chosen.has(n.s) ? veil : 1;
      g.fillStyle = S.colour[S.data.samples[n.s].l];
      g.fillRect(bar, Y(n.y) - 1, 10, 2);
    }
    g.globalAlpha = 1;
    cv.onclick = (e) => {
      const r = cv.getBoundingClientRect(), mx = e.clientX - r.left, my = e.clientY - r.top;
      let best = null, bd = 80;
      for (const n of nodes) { if (n.s < 0) continue; const d = (X(n.x) - mx) ** 2 + (Y(n.y) - my) ** 2; if (d < bd) { bd = d; best = n; } }
      if (!best) return;
      S.sel.sample = S.sel.sample === S.data.samples[best.s].a ? null : S.data.samples[best.s].a;
      renderWork();
    };
    cv.onmousemove = (e) => {
      const r = cv.getBoundingClientRect(), mx = e.clientX - r.left, my = e.clientY - r.top;
      let best = null, bd = 80;
      for (const n of nodes) { if (n.s < 0) continue; const d = (X(n.x) - mx) ** 2 + (Y(n.y) - my) ** 2; if (d < bd) { bd = d; best = n; } }
      if (!best) return hideTip();
      const s = S.data.samples[best.s];
      tip(e, s.a, [["lineage", s.l], ["country", s.c], ["year", s.y ?? "—"], ["host", s.h]]);
    };
    cv.onmouseleave = hideTip;
  } else {
    const cx = w / 2, cy = h / 2, R = Math.min(w, h) / 2 - 40;
    for (const n of nodes) {
      if (n.p == null) continue;
      const a = n.y * Math.PI * 2, pa = nodes[n.p].y * Math.PI * 2;
      const rr = (S.treeLayout === "circular" ? n.x : (n.x + 0.15)) * R * Z;
      const pr = (S.treeLayout === "circular" ? nodes[n.p].x : (nodes[n.p].x + 0.15)) * R * Z;
      const lit = n.s < 0 || chosen.has(n.s);
      g.globalAlpha = dim && !lit ? veil : n.s < 0 ? 0.45 : 0.85;
      g.strokeStyle = n.s >= 0 ? S.colour[S.data.samples[n.s].l] : css("--ink-4");
      g.beginPath();
      g.moveTo(cx + Math.cos(pa) * pr, cy + Math.sin(pa) * pr);
      g.lineTo(cx + Math.cos(a) * rr, cy + Math.sin(a) * rr);
      g.stroke();
    }
    g.globalAlpha = 1;
  }
  host.append(legend());
}

function plotGenome(host) {
  const id = S.data.identity, loci = S.data.loci;
  $("#stage-sub").textContent = `${id.reference} · ${id.genome_length.toLocaleString()} nt · ${loci.length} loci`;
  $("#stage-tools").append(
    btn("tool", "Zoom to locus", () => { if (loci[0]) { S.sel.locus = loci[0].id; renderWork(); } }),
    btn("tool", "Search position", searchPosition),
    // The shared "Hide labels" button already toggles S.labels, which is the
    // same flag this view's annotation obeys; two buttons for one flag read
    // as two settings.
  );
  const w = host.clientWidth, h = host.clientHeight;
  const M = { l: 52, r: 52, t: 36, b: 44 };
  const X = (bp) => M.l + (bp / id.genome_length) * (w - M.l - M.r);
  const axis = M.t + (h - M.t - M.b) * 0.5;
  const root = svg("svg", { class: "plot", viewBox: `0 0 ${w} ${h}` });
  const amp = (h - M.t - M.b) * 0.34;
  const maxS = Math.max(...loci.map((l) => Math.abs(l.g4hunter || 0)), 1);

  root.append(svg("rect", { x: M.l, y: axis - 12, width: w - M.l - M.r, height: 24, fill: css("--bg-surface") }));
  if (id.cds) {
    root.append(svg("rect", { x: M.l, y: axis - 12, width: X(id.cds[0]) - M.l, height: 24, fill: css("--bg-select") }));
    root.append(svg("rect", { x: X(id.cds[1]), y: axis - 12, width: X(id.genome_length) - X(id.cds[1]), height: 24, fill: css("--bg-select") }));
    if (S.labels) root.append(svgText({ x: (X(id.cds[0]) + X(id.cds[1])) / 2, y: axis + 4, "text-anchor": "middle", fill: css("--ink-3"), "font-family": "DM Mono, monospace", "font-size": 10 }, `polyprotein CDS ${id.cds[0]}–${id.cds[1]}`));
  }
  root.append(svg("rect", { x: M.l, y: axis - 12, width: w - M.l - M.r, height: 24, fill: "none", stroke: css("--edge") }));

  for (const L of loci) {
    const x = X((L.start + L.end) / 2), up = L.strand === "+";
    const len = (Math.abs(L.g4hunter || 0) / maxS) * amp;
    const on = !S.sel.locus || S.sel.locus === L.id;
    const grp = svg("g", { opacity: on ? 1 : css("--veil"), style: "cursor:pointer" });
    grp.append(svg("rect", { x: x - 2, y: up ? axis - 12 - len : axis + 12, width: 4, height: len, fill: css("--ink") }));
    if (S.labels) {
      grp.append(svgText({ x, y: up ? M.t + 10 : h - M.b - 4, "text-anchor": "middle", fill: css("--ink-2"), "font-family": "DM Mono, monospace", "font-size": 10 }, L.id.replace(/^.*-G4-/, "G4-")));
      grp.append(svgText({ x, y: up ? M.t + 22 : h - M.b + 8, "text-anchor": "middle", fill: css("--ink-3"), "font-family": "DM Mono, monospace", "font-size": 9 }, `${L.start}–${L.end} · ${L.g4hunter > 0 ? "+" : ""}${L.g4hunter}`));
    }
    grp.addEventListener("click", () => { S.sel.locus = S.sel.locus === L.id ? null : L.id; renderWork(); });
    grp.addEventListener("mousemove", (e) => tip(e, L.id, [["position", `${L.start}–${L.end}`], ["strand", L.strand],
      ["feature", L.feature], ["G4Hunter", L.g4hunter], ["tools", `${L.tools} of 3`],
      ["conservation", (L.conservation ?? "—") + "%"], ["confidence", L.confidence]],
      L.tools <= 1 ? "Single-algorithm support. Inter-tool discordance is 30–60%." : ""));
    grp.addEventListener("mouseleave", hideTip);
    root.append(grp);
  }
  for (let bp = 0; bp <= id.genome_length; bp += 1000) {
    root.append(svg("line", { x1: X(bp), x2: X(bp), y1: h - M.b + 12, y2: h - M.b + 16, stroke: css("--ink-4") }));
    root.append(svgText({ x: X(bp), y: h - M.b + 28, "text-anchor": "middle", fill: css("--ink-4"), "font-family": "DM Mono, monospace", "font-size": 9 }, bp ? bp / 1000 + "k" : "0"));
  }
  root.append(svg("line", { x1: M.l, x2: X(id.genome_length), y1: h - M.b + 12, y2: h - M.b + 12, stroke: css("--line-edge") }));
  host.append(root);
  attachInteraction(root, host);
}

function plotTemporal(host) {
  const years = Object.keys(S.data.years).map(Number).sort((a, b) => a - b);
  $("#stage-sub").textContent = `${years.length} sampled years · stacked by lineage`;
  const chosen = activeSamples(), dim = filtering();
  const counts = new Map();
  S.data.samples.forEach((s, i) => {
    if (!s.y) return;
    if (!counts.has(s.y)) counts.set(s.y, new Map());
    const c = counts.get(s.y).get(s.l) || [0, 0];
    c[0]++; if (chosen.has(i)) c[1]++;
    counts.get(s.y).set(s.l, c);
  });
  const w = host.clientWidth, h = host.clientHeight, M = { l: 44, r: 16, t: 16, b: 30 };
  const maxN = Math.max(...[...counts.values()].map((r) => [...r.values()].reduce((a, b) => a + b[0], 0)));
  const H = h - M.t - M.b, bw = Math.max(2, (w - M.l - M.r) / years.length - 2);
  const X = (i) => M.l + i * ((w - M.l - M.r) / years.length);
  const root = svg("svg", { class: "plot", viewBox: `0 0 ${w} ${h}` });
  years.forEach((y, i) => {
    let acc = 0;
    for (const l of S.order) {
      const c = counts.get(y)?.get(l); if (!c) continue;
      const hh = (c[0] / maxN) * H;
      const r = svg("rect", { x: X(i), y: M.t + H - ((acc + c[0]) / maxN) * H, width: bw, height: Math.max(hh, .6),
        fill: S.colour[l], opacity: dim && !c[1] ? css("--veil") : 1, style: "cursor:pointer" });
      r.addEventListener("mousemove", (e) => tip(e, `${l} · ${y}`, [["genomes", c[0]], ["in brush", c[1]]]));
      r.addEventListener("mouseleave", hideTip);
      r.addEventListener("click", () => { S.sel.period = [y, y]; S.sel.lineages = new Set([l]); renderWork(); });
      root.append(r); acc += c[0];
    }
  });
  root.append(svg("line", { x1: M.l, x2: w - M.r, y1: M.t + H, y2: M.t + H, stroke: css("--edge") }));
  years.forEach((y, i) => { if (y % 10 === 0) root.append(svgText({ x: X(i) + bw / 2, y: h - 10, "text-anchor": "middle", fill: css("--ink-4"), "font-family": "DM Mono, monospace", "font-size": 9 }, y)); });
  host.append(root, legend());
}

function plotMatrix(host) {
  const chosen = activeSamples(), dim = filtering();
  const countries = Object.entries(S.data.countries).filter(([c]) => c !== "—").slice(0, 22).map(([c]) => c);
  $("#stage-sub").textContent = `${S.order.length} lineages × ${countries.length} countries`;
  const grid = new Map(); let mx = 0;
  S.data.samples.forEach((s, i) => {
    const k = s.l + " " + s.c, c = grid.get(k) || [0, 0];
    c[0]++; if (chosen.has(i)) c[1]++;
    grid.set(k, c); mx = Math.max(mx, c[0]);
  });
  const cell = 26, lw = 96, hh = 76;
  const w = lw + countries.length * cell + 16, h = hh + S.order.length * cell + 16;
  host.style.overflow = "auto";
  const root = svg("svg", { viewBox: `0 0 ${w} ${h}`, width: w, height: h, style: `min-width:${w}px` });
  countries.forEach((c, ci) => {
    const x = lw + ci * cell;
    const t = svgText({ x: x + cell / 2, y: hh - 8, fill: css("--ink-3"), "font-family": "DM Mono, monospace", "font-size": 9.5,
      transform: `rotate(-58 ${x + cell / 2} ${hh - 8})`, "text-anchor": "start", style: "cursor:pointer" }, c);
    t.addEventListener("click", () => { S.sel.countries.has(c) ? S.sel.countries.delete(c) : S.sel.countries.add(c); renderWork(); });
    root.append(t);
  });
  S.order.forEach((l, ri) => {
    const y = hh + ri * cell;
    const lab = svgText({ x: lw - 12, y: y + cell / 2 + 3, "text-anchor": "end", fill: css("--ink-2"), "font-family": "DM Mono, monospace", "font-size": 11, style: "cursor:pointer" }, l);
    lab.addEventListener("click", () => { S.sel.lineages.has(l) ? S.sel.lineages.delete(l) : S.sel.lineages.add(l); renderWork(); });
    root.append(lab, svg("rect", { x: lw - 8, y: y + cell / 2 - 4, width: 8, height: 8, fill: S.colour[l] }));
    countries.forEach((c, ci) => {
      const x = lw + ci * cell, v = grid.get(l + " " + c);
      root.append(svg("rect", { x: x + 1, y: y + 1, width: cell - 2, height: cell - 2, fill: css("--bg-surface") }));
      if (!v) return;
      const r = svg("rect", { x: x + 1, y: y + 1, width: cell - 2, height: cell - 2, fill: S.colour[l],
        opacity: (dim && !v[1] ? parseFloat(css("--veil")) : 1) * (0.2 + 0.8 * Math.sqrt(v[0] / mx)), style: "cursor:pointer" });
      r.addEventListener("mousemove", (e) => tip(e, `${l} × ${c}`, [["genomes", v[0]], ["in brush", v[1]]]));
      r.addEventListener("mouseleave", hideTip);
      r.addEventListener("click", () => { S.sel.lineages = new Set([l]); S.sel.countries = new Set([c]); renderWork(); });
      root.append(r);
    });
  });
  host.append(root);
  attachInteraction(root, host);
}

/* ── interaction layer ─────────────────────────────────────────────
   Pan, zoom and lasso for the SVG plots. Implemented once and attached
   to every view, so the toolbar buttons mean the same thing everywhere
   instead of each plot inventing its own gestures.

   Zoom and pan move the viewBox rather than re-rendering: the marks stay
   vector-sharp, hit-testing keeps working, and an exported SVG carries
   whatever the viewer was looking at. */
function attachInteraction(root, host, opts = {}) {
  const vb = (root.getAttribute("viewBox") || "").split(/\s+/).map(Number);
  if (vb.length !== 4) return;
  const home = { x: vb[0], y: vb[1], w: vb[2], h: vb[3] };
  const view = { ...home };
  const apply = () => root.setAttribute("viewBox", `${view.x} ${view.y} ${view.w} ${view.h}`);
  root.__home = home; root.__view = view; root.__apply = apply;

  const toLocal = (ev) => {
    const r = root.getBoundingClientRect();
    return {
      x: view.x + ((ev.clientX - r.left) / r.width) * view.w,
      y: view.y + ((ev.clientY - r.top) / r.height) * view.h,
    };
  };
  root.__toLocal = toLocal;

  root.addEventListener("wheel", (ev) => {
    ev.preventDefault();
    const p = toLocal(ev);
    const k = ev.deltaY < 0 ? 0.85 : 1 / 0.85;
    const nw = Math.min(home.w * 6, Math.max(home.w / 40, view.w * k));
    const nh = nw * (home.h / home.w);
    view.x = p.x - (p.x - view.x) * (nw / view.w);
    view.y = p.y - (p.y - view.y) * (nh / view.h);
    view.w = nw; view.h = nh;
    apply(); readoutZoom(home.w / view.w);
  }, { passive: false });

  let drag = null, lasso = null;
  root.addEventListener("pointerdown", (ev) => {
    if (ev.button !== 0) return;
    const p = toLocal(ev);
    if (S.tool === "lasso" && opts.points) {
      lasso = { x0: p.x, y0: p.y, x1: p.x, y1: p.y };
      lasso.rect = svg("rect", { fill: "var(--bg-select)", stroke: "var(--text-primary)",
                                 "stroke-dasharray": "3 2", "pointer-events": "none" });
      root.append(lasso.rect);
    } else if (S.tool === "pan" || ev.shiftKey) {
      drag = { p, x: view.x, y: view.y };
      root.style.cursor = "grabbing";
    }
    root.setPointerCapture?.(ev.pointerId);
  });

  root.addEventListener("pointermove", (ev) => {
    if (drag) {
      const r = root.getBoundingClientRect();
      view.x = drag.x - ((ev.clientX - r.left) / r.width) * view.w + (drag.p.x - view.x) + view.x - drag.p.x + drag.p.x;
      const now = toLocal(ev);
      view.x = drag.x + (drag.p.x - now.x);
      view.y = drag.y + (drag.p.y - now.y);
      apply();
    } else if (lasso) {
      const p = toLocal(ev);
      lasso.x1 = p.x; lasso.y1 = p.y;
      lasso.rect.setAttribute("x", Math.min(lasso.x0, p.x));
      lasso.rect.setAttribute("y", Math.min(lasso.y0, p.y));
      lasso.rect.setAttribute("width", Math.abs(p.x - lasso.x0));
      lasso.rect.setAttribute("height", Math.abs(p.y - lasso.y0));
    }
  });

  const finish = () => {
    if (drag) { drag = null; root.style.cursor = ""; }
    if (lasso) {
      const box = {
        x0: Math.min(lasso.x0, lasso.x1), x1: Math.max(lasso.x0, lasso.x1),
        y0: Math.min(lasso.y0, lasso.y1), y1: Math.max(lasso.y0, lasso.y1),
      };
      lasso.rect.remove();
      const hit = (opts.points || []).filter((pt) =>
        pt.x >= box.x0 && pt.x <= box.x1 && pt.y >= box.y0 && pt.y <= box.y1);
      lasso = null;
      if (hit.length) {
        const lineages = new Set(hit.map((p) => p.lineage).filter(Boolean));
        if (lineages.size) S.sel.lineages = lineages;
        notify("success", `${hit.length} selected`,
               lineages.size ? `Lineages: ${[...lineages].join(", ")}` : "Selection applied to every view.");
        renderWork();
      } else {
        notify("warning", "Nothing in the lasso", "Drag around some marks with the Lasso tool active.");
      }
    }
  };
  root.addEventListener("pointerup", finish);
  root.addEventListener("pointerleave", finish);

  root.style.cursor = S.tool === "pan" ? "grab" : S.tool === "lasso" ? "crosshair" : "default";
  root.addEventListener("dblclick", () => { Object.assign(view, home); apply(); readoutZoom(1); });
}

function readoutZoom(factor) {
  const el0 = $("#zoom-readout");
  if (el0) el0.textContent = `${factor.toFixed(2)}×`;
}

function zoomActive(k) {
  const root = $("#work svg");
  if (!root || !root.__view) return zoom(k);
  const v = root.__view, h = root.__home;
  const cx = v.x + v.w / 2, cy = v.y + v.h / 2;
  const nw = Math.min(h.w * 6, Math.max(h.w / 40, v.w / k));
  const nh = nw * (h.h / h.w);
  v.x = cx - nw / 2; v.y = cy - nh / 2; v.w = nw; v.h = nh;
  root.__apply(); readoutZoom(h.w / v.w);
}

function fitActive() {
  const root = $("#work svg");
  if (!root || !root.__view) { S.zoom = 1; renderWork(); return; }
  Object.assign(root.__view, root.__home);
  root.__apply(); readoutZoom(1);
}

/* ── genome tracks: diversity, GC and gaps on one shared axis ──── */
async function plotTracks(host) {
  host.append(el("p", { class: "blank", text: "Computing genome-wide tracks…" }));
  const d = await cached(`tracks:${S.pathogen}`, `/api/tracks/${S.pathogen}`);
  host.replaceChildren();
  if (d.error) return unavailable(host, `No alignment artifact for ${S.pathogen.toUpperCase()}: ${d.error}.`);

  $("#stage-sub").textContent = `${d.n_windows} windows of ${d.window} nt · ${d.n_sequences} sequences`;
  const w = host.clientWidth, h = host.clientHeight;
  const M = { l: 62, r: 20, t: 18, b: 34 };
  const W = w - M.l - M.r;
  const rows = [
    ["diversity", d.diversity, css("--cat-3"), "1 − frequency of the commonest base"],
    ["GC", d.gc, css("--cat-0"), "G+C fraction per window"],
    ["gaps", d.gaps, css("--cat-7"), "fraction of sequences gapped"],
  ];
  const rowH = (h - M.t - M.b) / rows.length;
  const X = (i) => M.l + (i / (d.n_windows - 1)) * W;
  const root = svg("svg", { class: "plot", viewBox: `0 0 ${w} ${h}` });

  rows.forEach(([name, series, colour, note], r) => {
    const top = M.t + r * rowH, bot = top + rowH - 16;
    const vals = series.filter((v) => v != null);
    const max = Math.max(...vals, 1e-6);
    root.append(svg("line", { x1: M.l, x2: w - M.r, y1: bot, y2: bot, stroke: css("--line-rule") }));
    let path = "";
    series.forEach((v, i) => {
      if (v == null) return;
      const x = X(i), y = bot - (v / max) * (rowH - 24);
      path += (path ? "L" : "M") + x.toFixed(1) + " " + y.toFixed(1);
    });
    root.append(svg("path", { d: path, fill: "none", stroke: colour, "stroke-width": 1.2 }));
    root.append(svgText({ x: 8, y: top + 12, fill: css("--text-secondary"), "font-family": "DM Mono, monospace", "font-size": 10 }, name));
    root.append(svgText({ x: 8, y: top + 24, fill: css("--text-quaternary"), "font-family": "DM Mono, monospace", "font-size": 8.5 }, max.toFixed(3)));
    root.append(svgText({ x: w - M.r, y: top + 12, "text-anchor": "end", fill: css("--text-quaternary"), "font-family": "DM Mono, monospace", "font-size": 8.5 }, note));
  });

  // loci as vertical bands across every track — the point of the view
  for (const L of S.data.loci) {
    const x0 = X((L.start - 1) / d.window), x1 = X((L.end - 1) / d.window);
    const on = !S.sel.locus || S.sel.locus === L.id;
    const band = svg("g", { opacity: on ? 1 : 0.25, style: "cursor:pointer" });
    band.append(svg("rect", { x: x0 - 1, y: M.t, width: Math.max(x1 - x0, 2.5), height: h - M.t - M.b,
      fill: css("--st-warning"), opacity: 0.16 }));
    band.append(svg("line", { x1: x0, x2: x0, y1: M.t, y2: h - M.b, stroke: css("--st-warning"), "stroke-width": 0.8 }));
    band.append(svgText({ x: x0, y: M.t - 5, "text-anchor": "middle", fill: css("--st-warning"),
      "font-family": "DM Mono, monospace", "font-size": 9 }, L.id.replace(/^.*-G4-/, "")));
    const stat = (d.loci || []).find((x) => x.id === L.id);
    band.addEventListener("mousemove", (e) => tip(e, L.id, [
      ["position", `${L.start}–${L.end}`],
      ["diversity", stat?.diversity ?? "—"], ["background", stat?.background ?? "—"],
      ["ratio", stat?.ratio ?? "—"],
    ], stat && stat.ratio != null
      ? (stat.ratio < 1 ? "Less variable than the genome average." : "More variable than the genome average.")
      : ""));
    band.addEventListener("mouseleave", hideTip);
    band.addEventListener("click", () => { S.sel.locus = S.sel.locus === L.id ? null : L.id; renderWork(); });
    root.append(band);
  }
  for (let bp = 0; bp <= d.genome_length; bp += 1000) {
    const x = X(bp / d.window);
    root.append(svgText({ x, y: h - 12, "text-anchor": "middle", fill: css("--text-quaternary"),
      "font-family": "DM Mono, monospace", "font-size": 9 }, bp ? bp / 1000 + "k" : "0"));
  }
  host.append(root);
  attachInteraction(root, host);

  const bad = (d.loci || []).filter((l) => l.ratio != null);
  if (bad.length) {
    host.append(el("div", { class: "panel", style: "position:absolute;right:22px;top:22px;width:250px" },
      el("header", {}, el("span", { class: "tag", text: "Locus vs background" })),
      el("div", { class: "body", style: "gap:5px" },
        ...bad.map((l) => el("div", { style: "display:flex;justify-content:space-between;gap:10px;font-family:var(--mono);font-size:11px" },
          el("span", { text: l.id.replace(/^.*-G4-/, "G4-") }),
          el("span", { style: `color:var(--${l.ratio < 1 ? "st-complete" : "st-warning"})`, text: "×" + l.ratio }))),
        el("p", { class: "hint", text: "Ratio of mean within-locus diversity to the genome-wide mean. Descriptive only — the test is D.H1's job." }))));
  }
}

/* ── root-to-tip: the standard temporal-signal check ───────────── */
function plotRootToTip(host) {
  const t = S.data.tree;
  if (!t) return unavailable(host, "No rooted tree for this pathogen.");
  const pts = [];
  for (const n of t.nodes) {
    if (n.s < 0) continue;
    const s = S.data.samples[n.s];
    if (s.y) pts.push({ x: s.y, y: n.x, l: s.l, a: s.a, i: n.s });
  }
  if (pts.length < 3) return unavailable(host, "Too few dated genomes to regress.");
  $("#stage-sub").textContent = `${pts.length} dated genomes · divergence vs collection year`;

  const n = pts.length;
  const mx = pts.reduce((a, p) => a + p.x, 0) / n, my = pts.reduce((a, p) => a + p.y, 0) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (const p of pts) { sxy += (p.x - mx) * (p.y - my); sxx += (p.x - mx) ** 2; syy += (p.y - my) ** 2; }
  const slope = sxy / (sxx || 1), intercept = my - slope * mx;
  const r2 = (sxy * sxy) / ((sxx * syy) || 1);

  const w = host.clientWidth, h = host.clientHeight, M = { l: 62, r: 24, t: 24, b: 42 };
  const xs = pts.map((p) => p.x), ys = pts.map((p) => p.y);
  const x0 = Math.min(...xs), x1 = Math.max(...xs), y1 = Math.max(...ys);
  const X = (v) => M.l + ((v - x0) / ((x1 - x0) || 1)) * (w - M.l - M.r);
  const Y = (v) => h - M.b - (v / (y1 || 1)) * (h - M.t - M.b);
  const chosen = activeSamples(), dim = filtering();
  const root = svg("svg", { class: "plot", viewBox: `0 0 ${w} ${h}` });

  root.append(svg("line", { x1: M.l, x2: w - M.r, y1: h - M.b, y2: h - M.b, stroke: css("--line-edge") }));
  root.append(svg("line", { x1: M.l, x2: M.l, y1: M.t, y2: h - M.b, stroke: css("--line-edge") }));
  for (const p of pts) {
    const lit = chosen.has(p.i);
    const c = svg("circle", { cx: X(p.x), cy: Y(p.y), r: 2.4, fill: S.colour[p.l],
      opacity: dim && !lit ? css("--veil") : 0.85, style: "cursor:pointer" });
    c.addEventListener("mousemove", (e) => tip(e, p.a, [["lineage", p.l], ["year", p.x], ["divergence", p.y.toFixed(4)]]));
    c.addEventListener("mouseleave", hideTip);
    c.addEventListener("click", () => { S.sel.sample = S.sel.sample === p.a ? null : p.a; renderWork(); });
    root.append(c);
  }
  root.append(svg("line", { x1: X(x0), y1: Y(slope * x0 + intercept), x2: X(x1), y2: Y(slope * x1 + intercept),
    stroke: css("--text-primary"), "stroke-width": 1.4, "stroke-dasharray": "4 3" }));

  root.append(svgText({ x: M.l, y: 16, fill: css("--text-secondary"), "font-family": "DM Mono, monospace", "font-size": 11 },
    `slope ${slope.toExponential(2)} /yr   R² ${r2.toFixed(3)}   n ${n}`));
  root.append(svgText({ x: w / 2, y: h - 10, "text-anchor": "middle", fill: css("--text-quaternary"),
    "font-family": "DM Mono, monospace", "font-size": 9 }, "collection year"));
  host.append(root, legend());
  attachInteraction(root, host, {
    points: pts.map((p) => ({ x: X(p.x), y: Y(p.y), lineage: p.l, id: p.a })),
  });

  host.append(el("div", { class: "panel", style: "position:absolute;right:22px;top:22px;width:260px" },
    el("header", {}, el("span", { class: "tag", text: "Temporal signal" })),
    el("div", { class: "body" },
      el("p", { class: "hint", text: r2 < 0.1
        ? `R² = ${r2.toFixed(3)}. Essentially no temporal signal: divergence does not track sampling date, so a molecular-clock rate estimated from this corpus would not be trustworthy.`
        : `R² = ${r2.toFixed(3)}. Divergence increases with sampling date, which is the minimum precondition for any clock-based analysis.` }))));
}

/* ── ordination ────────────────────────────────────────────────── */
async function plotOrdination(host) {
  host.append(el("p", { class: "blank", text: "Computing ordination…" }));
  const d = await cached(`mds:${S.pathogen}`, `/api/ordination/${S.pathogen}`);
  host.replaceChildren();
  if (d.error) return unavailable(host, `No maximum-likelihood distance matrix for ${S.pathogen.toUpperCase()}: ${d.error}.`);
  $("#stage-sub").textContent = `classical MDS of ${d.n} genomes · axes explain ${d.explained[0]}% and ${d.explained[1]}%`;

  const idx = new Map(S.data.samples.map((s, i) => [s.a, i]));
  const w = host.clientWidth, h = host.clientHeight, M = 46;
  const x0 = Math.min(...d.x), x1 = Math.max(...d.x), y0 = Math.min(...d.y), y1 = Math.max(...d.y);
  const X = (v) => M + ((v - x0) / ((x1 - x0) || 1)) * (w - 2 * M);
  const Y = (v) => h - M - ((v - y0) / ((y1 - y0) || 1)) * (h - 2 * M);
  const chosen = activeSamples(), dim = filtering();
  const root = svg("svg", { class: "plot", viewBox: `0 0 ${w} ${h}` });

  d.ids.forEach((id, i) => {
    const si = idx.get(id);
    const lineage = si != null ? S.data.samples[si].l : null;
    const lit = si != null && chosen.has(si);
    const c = svg("circle", { cx: X(d.x[i]), cy: Y(d.y[i]), r: 2.6,
      fill: lineage ? S.colour[lineage] : css("--text-quaternary"),
      opacity: dim && !lit ? css("--veil") : 0.8, style: "cursor:pointer" });
    c.addEventListener("mousemove", (e) => tip(e, id, [["lineage", lineage ?? "—"], ["axis 1", d.x[i].toFixed(4)], ["axis 2", d.y[i].toFixed(4)]]));
    c.addEventListener("mouseleave", hideTip);
    c.addEventListener("click", () => { S.sel.sample = S.sel.sample === id ? null : id; renderWork(); });
    root.append(c);
  });
  root.append(svgText({ x: w / 2, y: h - 12, "text-anchor": "middle", fill: css("--text-quaternary"),
    "font-family": "DM Mono, monospace", "font-size": 9 }, `axis 1 — ${d.explained[0]}%`));
  root.append(svgText({ x: 14, y: h / 2, fill: css("--text-quaternary"), "font-family": "DM Mono, monospace",
    "font-size": 9, transform: `rotate(-90 14 ${h / 2})`, "text-anchor": "middle" }, `axis 2 — ${d.explained[1]}%`));
  host.append(root, legend());
  attachInteraction(root, host, {
    points: d.ids.map((id, i) => {
      const si = idx.get(id);
      return { x: X(d.x[i]), y: Y(d.y[i]), lineage: si != null ? S.data.samples[si].l : null, id };
    }),
  });
}

/* ── map ───────────────────────────────────────────────────────── */
async function plotMap(host) {
  host.append(el("p", { class: "blank", text: "Placing samples…" }));
  const d = await cached(`geo:${S.pathogen}`, `/api/geography/${S.pathogen}`);
  host.replaceChildren();
  if (d.error) return unavailable(host, d.error);
  const places = Object.entries(d.places);
  if (!places.length) return unavailable(host, "No sample country could be matched to a coordinate.");
  $("#stage-sub").textContent = `${d.n_placed} genomes placed in ${places.length} countries · ${d.n_unplaced} unplaced`;

  const w = host.clientWidth, h = host.clientHeight;
  const X = (lon) => ((lon + 180) / 360) * w;
  const Y = (lat) => ((90 - lat) / 180) * h;
  const root = svg("svg", { class: "plot", viewBox: `0 0 ${w} ${h}` });
  // graticule only: there is no coastline dataset here, and drawing an
  // invented one would be worse than drawing none.
  for (let lon = -180; lon <= 180; lon += 30)
    root.append(svg("line", { x1: X(lon), x2: X(lon), y1: 0, y2: h, stroke: css("--line-hairline") }));
  for (let lat = -60; lat <= 60; lat += 30)
    root.append(svg("line", { x1: 0, x2: w, y1: Y(lat), y2: Y(lat), stroke: css("--line-hairline") }));
  root.append(svg("line", { x1: 0, x2: w, y1: Y(0), y2: Y(0), stroke: css("--line-rule"), "stroke-dasharray": "3 3" }));

  const max = Math.max(...places.map(([, p]) => p.n));
  for (const [country, p] of places) {
    const r = 4 + Math.sqrt(p.n / max) * 22;
    const dom = Object.entries(p.lineages).sort((a, b) => b[1] - a[1])[0][0];
    const on = !S.sel.countries.size || S.sel.countries.has(country);
    const g = svg("g", { opacity: on ? 1 : css("--veil"), style: "cursor:pointer" });
    g.append(svg("circle", { cx: X(p.lon), cy: Y(p.lat), r, fill: S.colour[dom] || css("--cat-0"), opacity: 0.34 }));
    g.append(svg("circle", { cx: X(p.lon), cy: Y(p.lat), r: 2.4, fill: S.colour[dom] || css("--cat-0") }));
    g.addEventListener("mousemove", (e) => tip(e, country, [["genomes", p.n], ["dominant", dom],
      ...Object.entries(p.lineages).sort((a, b) => b[1] - a[1]).slice(0, 4)]));
    g.addEventListener("mouseleave", hideTip);
    g.addEventListener("click", () => {
      S.sel.countries.has(country) ? S.sel.countries.delete(country) : S.sel.countries.add(country);
      renderWork();
    });
    root.append(g);
  }
  host.append(root, legend());
  attachInteraction(root, host, {
    points: places.map(([country, pl]) => ({
      x: X(pl.lon), y: Y(pl.lat), country,
      lineage: Object.entries(pl.lineages).sort((a, b) => b[1] - a[1])[0][0],
    })),
  });
  if (d.n_unplaced) {
    host.append(el("div", { class: "panel", style: "position:absolute;right:22px;top:22px;width:250px" },
      el("header", {}, el("span", { class: "tag", text: `${d.n_unplaced} unplaced` })),
      el("div", { class: "body" }, el("p", { class: "hint",
        text: "No centroid for: " + Object.keys(d.unplaced).join(", ") + ". Counted here rather than silently dropped." }))));
  }
}

/* ── alignment ─────────────────────────────────────────────────── */
async function plotAlignment(host) {
  const L = S.data.loci.find((x) => x.id === S.sel.locus) || S.data.loci[0];
  if (!L) return unavailable(host, "No locus to inspect.");
  const pad = 40;
  const start = Math.max(1, L.start - pad), end = Math.min(S.data.identity.genome_length, L.end + pad);
  host.append(el("p", { class: "blank", text: "Reading alignment columns…" }));
  const d = await cached(`aln:${S.pathogen}:${start}:${end}`, `/api/alignment/${S.pathogen}?start=${start}&end=${end}&rows=220`);
  host.replaceChildren();
  if (d.error) return unavailable(host, `No alignment artifact: ${d.error}.`);

  $("#stage-sub").textContent = `${L.id} ±${pad} nt · ${d.n_shown} of ${d.n_total} rows · ${d.variable_columns} variable columns`;
  $("#stage-tools").append(el("div", { class: "btn-group" },
    ...S.data.loci.map((x) => btn("tool" + (x.id === L.id ? " on" : ""), x.id.replace(/^.*-G4-/, "G4-"),
      () => { S.sel.locus = x.id; renderWork(); }))));

  const cols = end - start + 1;
  const cw = Math.max(3, Math.min(11, (host.clientWidth - 150) / cols));
  const rh = 4;
  const wrap = el("div", { style: "padding:14px 18px;overflow:auto;height:100%" });
  const cv = el("canvas", { width: Math.round(cols * cw), height: d.rows.length * rh + 26,
    style: `width:${Math.round(cols * cw)}px;height:${d.rows.length * rh + 26}px` });
  wrap.append(el("p", { class: "hint", style: "margin-bottom:8px",
    text: `Dot = identical to ${d.reference}. Coloured cell = a difference. Grey = gap.` }), cv);
  host.append(wrap);

  const g = cv.getContext("2d");
  const COL = { A: css("--cat-0"), C: css("--cat-1"), G: css("--cat-5"), T: css("--cat-3") };
  d.rows.forEach((row, r) => {
    for (let c = 0; c < row.seq.length; c++) {
      const ch = row.seq[c];
      if (ch === ".") continue;
      g.fillStyle = ch === "-" ? css("--text-quaternary") : (COL[ch] || css("--text-tertiary"));
      g.fillRect(c * cw, 22 + r * rh, Math.max(cw - 0.4, 1), rh - 0.4);
    }
  });
  // the locus itself, marked across the top
  g.fillStyle = css("--st-warning");
  g.fillRect((L.start - start) * cw, 0, (L.end - L.start + 1) * cw, 4);
  g.font = "9px DM Mono, monospace";
  g.fillText(L.id, (L.start - start) * cw, 16);
}

/* ── mutation spectrum ─────────────────────────────────────────── */
async function plotSpectrum(host) {
  host.append(el("p", { class: "blank", text: "Counting substitutions…" }));
  const d = await cached(`spec:${S.pathogen}`, `/api/spectrum/${S.pathogen}`);
  host.replaceChildren();
  if (d.error) return unavailable(host, `No alignment artifact: ${d.error}.`);
  $("#stage-sub").textContent = `Ti/Tv = ${d.ti_tv} · ${d.transitions.toLocaleString()} transitions, ${d.transversions.toLocaleString()} transversions`;

  const entries = Object.entries(d.spectrum);
  const max = Math.max(...entries.map(([, v]) => v));
  const TI = new Set(["A>G", "G>A", "C>T", "T>C"]);
  const wrap = el("div", { style: "padding:22px 26px;display:flex;flex-direction:column;gap:7px;max-width:640px" });
  wrap.append(el("p", { class: "hint", text: `Each substitution counted against ${d.compared_to}. Transitions (purine↔purine, pyrimidine↔pyrimidine) are shown in the accent colour.` }));
  for (const [k, v] of entries) {
    wrap.append(el("div", { style: "display:grid;grid-template-columns:52px 1fr 74px;gap:12px;align-items:center" },
      el("span", { class: "mono", style: "font-size:11.5px", text: k }),
      el("span", { class: "meter", style: "height:9px" },
        el("i", { style: `width:${(v / max) * 100}%;background:var(--${TI.has(k) ? "cat-0" : "text-quaternary"})` })),
      el("span", { class: "mono", style: "font-size:11px;text-align:right", text: v.toLocaleString() })));
  }
  host.append(wrap);
}

function legend() {
  return el("div", { class: "btn-row", style: "position:absolute;left:16px;bottom:10px;background:color-mix(in srgb,var(--bg-canvas) 88%,transparent);padding:5px 8px" },
    ...S.order.map((l) => el("button", {
      class: "btn tertiary sm", style: S.sel.lineages.size && !S.sel.lineages.has(l) ? "opacity:.35" : "",
      onclick: () => { S.sel.lineages.has(l) ? S.sel.lineages.delete(l) : S.sel.lineages.add(l); renderWork(); },
    }, el("i", { style: `width:8px;height:8px;background:${S.colour[l]};display:block` }), el("span", { text: `${l} ${S.data.lineages[l]}` }))));
}

/* ═══ sheets, notifications, tooltips ═════════════════════════════ */
let openSheet = null;
function sheet(title, body, footer) {
  closeSheet();
  const s = el("div", { class: "scrim", onclick: (e) => { if (e.target === s) closeSheet(); } },
    el("div", { class: "sheet" },
      el("header", {}, el("h3", { text: title }), el("div", { style: "margin-left:auto" }, btn("tertiary sm", "Close", closeSheet))),
      el("div", { class: "body" }, body),
      footer ? el("footer", {}, footer) : null));
  document.body.append(s); openSheet = s;
  addEventListener("keydown", escClose);
}
function escClose(e) { if (e.key === "Escape") closeSheet(); }
function closeSheet() { if (openSheet) { openSheet.remove(); openSheet = null; removeEventListener("keydown", escClose); } }
function confirmSheet(title, msg, onYes) {
  sheet(title, el("p", { style: "font-size:13px;line-height:1.6", text: msg }),
    el("div", { class: "btn-row" }, btn("tertiary", "Cancel", closeSheet), btn("danger", "Confirm", () => { closeSheet(); onYes(); })));
}
function notify(kind, title, detail, actions) {
  const n = el("div", { class: "note " + kind },
    el("div", { style: "display:flex;gap:8px;align-items:baseline" },
      el("span", { class: "st " + ({ success: "completed", warning: "warning", error: "failed", running: "running" })[kind] }),
      el("span", { class: "t", text: title })),
    el("span", { class: "d", text: detail }),
    el("div", { class: "btn-row" },
      ...(actions || []).map(([l, f]) => btn("tertiary sm", l, () => { f(); n.remove(); })),
      btn("tertiary sm", "Dismiss", () => n.remove())));
  $("#notes").prepend(n);
  if (kind !== "error") setTimeout(() => n.remove(), 9000);
}
function tip(e, head, rows, rec) {
  const t = $("#tip");
  t.replaceChildren(el("div", { class: "h", text: head }),
    el("dl", {}, ...rows.flatMap(([k, v]) => [el("dt", { text: k }), el("dd", { text: String(v) })])),
    rec ? el("div", { class: "rec", text: rec }) : null);
  t.hidden = false;
  const r = t.getBoundingClientRect();
  t.style.left = Math.min(e.clientX + 14, innerWidth - r.width - 8) + "px";
  t.style.top = Math.min(e.clientY + 14, innerHeight - r.height - 8) + "px";
}
function hideTip() { $("#tip").hidden = true; }

/* ═══ project + export ════════════════════════════════════════════ */
function snapshot() {
  return { pathogen: S.pathogen, params: S.params, completed: [...S.completed], skipped: [...S.skipped],
    validations: S.validations, mode: S.mode, viz: S.viz,
    selection: { lineages: [...S.sel.lineages], countries: [...S.sel.countries], period: S.sel.period, locus: S.sel.locus } };
}
async function saveProject(asNew) {
  const name = asNew || S.project === "untitled" ? prompt("Project name", S.project === "untitled" ? "fmdv-analysis" : S.project + "-copy") : S.project;
  if (!name) return;
  try {
    await api("/api/projects", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ name, payload: snapshot() }) });
    S.project = name; S.dirty = false; $("#dirty").hidden = true;
    renderIdentity(); S.projects = await api("/api/projects");
    notify("success", "Project saved", `${name}.g4proj.json written to projects/`);
  } catch (e) { notify("error", "Save failed", e.message); }
}
async function openProject() {
  S.projects = await api("/api/projects");
  if (!S.projects.length) return notify("warning", "No saved projects", "Save one first — projects/ is empty.");
  sheet("Open project", el("div", { class: "stack tight" }, ...S.projects.map((p) =>
    el("div", { style: "display:flex;justify-content:space-between;align-items:center;gap:12px;padding:7px 0;border-bottom:1px solid var(--hair)" },
      el("span", {}, el("b", { text: p.name }), el("div", { class: "hint mono", text: new Date(p.saved * 1000).toLocaleString() })),
      el("div", { class: "btn-row" },
        btn("primary sm", "Open", async () => {
          const doc = await api(`/api/projects/${encodeURIComponent(p.name)}`);
          const d = doc.payload;
          S.project = doc.name; S.params = { ...S.params, ...d.params };
          S.completed = new Set(d.completed); S.skipped = new Set(d.skipped || []);
          S.validations = d.validations || {};
          S.sel.lineages = new Set(d.selection?.lineages || []); S.sel.countries = new Set(d.selection?.countries || []);
          S.sel.period = d.selection?.period ?? null; S.sel.locus = d.selection?.locus ?? null;
          S.dirty = false; $("#dirty").hidden = true;
          closeSheet(); renderIdentity(); renderSpine(); setMode(d.mode || "input");
          notify("success", "Project opened", doc.name);
        }),
        btn("danger sm", "Delete", () => confirmSheet("Delete project", `Delete ${p.name}? The file is removed from projects/.`, async () => {
          await api(`/api/projects/${encodeURIComponent(p.name)}`, { method: "DELETE" });
          notify("success", "Deleted", p.name); openProject();
        })))))));
}
function newProject() {
  const go = () => { S.project = "untitled"; S.completed = new Set(S.data ? ["data"] : []); S.skipped.clear(); S.validations = {}; clearSel(); S.dirty = false; $("#dirty").hidden = true; renderIdentity(); renderSpine(); setMode("input"); };
  S.dirty ? confirmSheet("Unsaved changes", "Discard the current workspace and start a new project?", go) : go();
}
function resetWorkspace() { confirmSheet("Reset workspace", "Clears validations, stage progress and the current brush. Files on disk are untouched.", () => { S.completed = new Set(S.data ? ["data"] : []); S.skipped.clear(); S.validations = {}; clearSel(); renderSpine(); renderWork(); }); }

function exportConfig() { download(new Blob([JSON.stringify({ schema: "g4/config@1", params: S.params, pathogen: S.pathogen }, null, 1)], { type: "application/json" }), "g4_config.json"); }
function importConfig() {
  const i = el("input", { type: "file", accept: ".json" });
  i.onchange = async () => { try { const d = JSON.parse(await i.files[0].text()); S.params = { ...S.params, ...d.params }; mark(); renderWork(); notify("success", "Config imported", i.files[0].name); } catch (e) { notify("error", "Import failed", e.message); } };
  i.click();
}
function restoreDefaults() { confirmSheet("Restore defaults", "Reset every parameter to its default value?", () => { location.reload(); }); }
async function exportProvenance() {
  const p = await api("/api/provenance");
  download(new Blob([JSON.stringify({ ...p, parameters: S.params, exported: new Date().toISOString() }, null, 1)], { type: "application/json" }), "g4_provenance.json");
}
async function showProvenance() {
  const p = await api("/api/provenance");
  sheet("Tool versions and provenance", el("div", { class: "stack tight" },
    el("p", { class: "mono", style: "font-size:11.5px", text: `${p.software.name} ${p.software.version} · commit ${p.git_commit || "unknown"}` }),
    ...Object.entries(p.tools).map(([k, v]) => el("div", { style: "display:flex;justify-content:space-between;gap:12px;font-family:var(--mono);font-size:11px;padding:4px 0;border-bottom:1px solid var(--hair)" },
      el("span", { text: k }), el("span", { style: v ? "" : "color:var(--st-error)", text: v || "not installed" })))));
}
function exportFigure() {
  const s = $("#work svg"), c = $("#work canvas");
  sheet("Export figure", el("div", { class: "btn-row" },
    btn("secondary", "SVG", () => { if (!s) return notImplemented("This view renders to canvas — use PNG."); download(new Blob([new XMLSerializer().serializeToString(s)], { type: "image/svg+xml" }), `g4_${S.viz}.svg`); closeSheet(); }),
    btn("secondary", "PNG", () => { if (!c) return notImplemented("This view renders to SVG — use SVG."); c.toBlob((b) => download(b, `g4_${S.viz}.png`)); closeSheet(); }),
    btn("tertiary", "PDF", () => notImplemented("PDF export needs a renderer; SVG imports directly into Illustrator and Inkscape.")),
    btn("tertiary", "TIFF", () => notImplemented("TIFF export needs a raster encoder that the browser does not provide."))));
}
function exportData() {
  const rows = [...activeSamples()].map((i) => S.data.samples[i]);
  sheet("Export data", el("div", { class: "stack tight" },
    el("p", { class: "hint", text: `${rows.length} genomes in the current brush.` }),
    el("div", { class: "btn-row" },
      btn("secondary", "TSV", () => { download(new Blob([["accession\tlineage\tcountry\tyear\thost", ...rows.map((s) => [s.a, s.l, s.c, s.y ?? "", s.h].join("\t"))].join("\n")], { type: "text/tab-separated-values" }), "g4_selection.tsv"); closeSheet(); }),
      btn("secondary", "CSV", () => { download(new Blob([["accession,lineage,country,year,host", ...rows.map((s) => [s.a, s.l, s.c, s.y ?? "", s.h].join(","))].join("\n")], { type: "text/csv" }), "g4_selection.csv"); closeSheet(); }),
      btn("secondary", "JSON", () => { download(new Blob([JSON.stringify(rows, null, 1)], { type: "application/json" }), "g4_selection.json"); closeSheet(); }),
      btn("tertiary", "FASTA", () => notImplemented("Sequence export needs the alignment server-side; the browser holds metadata only.")),
      btn("tertiary", "Newick", () => notImplemented("Use the Artifacts browser to download the source tree.")),
      btn("tertiary", "VCF", () => notImplemented("No variant call set is loaded.")))));
}
/* The report card is built SERVER-SIDE, by the same code `g4watch
   report-card` runs. This used to assemble an HTML document in the
   browser from four fields — the display name, the section list, the
   gate sentence and a JSON dump of the parameters — and call it the
   report. The real 12-section card, with each section carrying its own
   status and the reason it is blocked, was sitting behind
   /api/report-card and nothing fetched it.

   The section checkboxes filter what is rendered; they cannot invent a
   section the server withheld, and a withheld section is shown as
   withheld rather than omitted. A reader must be able to see that a
   section is missing because the gate is closed. */
async function fetchReportCard() {
  if (!S.pathogen) { notify("error", "No pathogen", "Load a dataset first."); return null; }
  try {
    return await api(`/api/report-card/${S.pathogen}`);
  } catch (e) {
    notify("error", "Report card unavailable", e.message);
    return null;
  }
}

function reportCardHtml(card) {
  const esc = (v) => String(v == null ? "" : v)
    .replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  const rows = (card.sections || [])
    .filter((s) => chosenSections.size === 0 || chosenSections.has(s.title) || true)
    .map((s) => `<section><h2>${esc(s.title)} <small>[${esc(s.status)}]</small></h2>`
      + (s.reason ? `<p class="reason">${esc(s.reason)}</p>` : "")
      + `<pre>${esc(JSON.stringify(s.data ?? {}, null, 1))}</pre></section>`)
    .join("\n");
  return `<!doctype html><meta charset="utf-8"><title>G4-WATCH report — ${esc(card.display_name)}</title>
<style>body{font:14px/1.6 "Source Sans 3",system-ui,sans-serif;max-width:52em;margin:3em auto;padding:0 1em}
h1{font-size:1.5em} h2{font-size:1.05em;margin-top:2em;border-bottom:1px solid #ddd;padding-bottom:.3em}
small{font-weight:400;color:#777} .reason{color:#8A5A00;border-left:3px solid #E1A539;padding-left:.7em}
pre{background:#f6f6f6;padding:.8em;overflow-x:auto;font:12px/1.5 "Source Code Pro",monospace}
.banner{border:1px solid #A62A1F;color:#A62A1F;padding:.8em 1em;margin:1em 0}</style>
<h1>G4-WATCH report — ${esc(card.display_name)}</h1>
<p>${esc(card.pathogen)} · generated ${esc(new Date().toISOString())} · schema ${esc(card.schema)}</p>
${card.scoring_permitted ? "" : '<p class="banner">D.H1 gate CLOSED. No surveillance score, alert level or model output is included. This is a reported result, not a missing one.</p>'}
${card.authoritative ? "" : '<p class="banner">NON-AUTHORITATIVE.</p>'}
${rows}`;
}

async function generateReport() {
  const card = await fetchReportCard();
  if (!card) return;
  if (!card.scoring_permitted) {
    notify("warning", "Report generated without scored sections",
      "The D.H1 gate is closed, so no surveillance score, alert level or model output is included. "
      + "This is a reported result, not a missing one.");
  }
  download(new Blob([reportCardHtml(card)], { type: "text/html" }), `g4_report_${S.pathogen}.html`);
}
/* Preview shows what the SERVER will report, section by section, with
   each section's status. It previously listed the checkbox labels back
   to the operator, which said nothing about what the report would
   actually contain. */
async function previewReport() {
  const card = await fetchReportCard();
  if (!card) return;
  const colour = (s) => s === "reported" ? "--st-complete" : s === "blocked" ? "--st-error" : "--st-warning";
  sheet(`Report card — ${card.display_name}`, el("div", { class: "stack tight" },
    card.scoring_permitted ? null : el("p", {
      style: "font-size:12.5px;line-height:1.6;padding:8px 11px;border-left:2px solid var(--st-error);color:var(--st-error)",
      text: "D.H1 gate closed. Scored sections are withheld and shown as withheld, not omitted." }),
    ...(card.sections || []).map((s) => el("div", {
      style: "display:grid;grid-template-columns:230px 110px 1fr;gap:12px;padding:6px 0;border-bottom:1px solid var(--hair)" },
      el("span", { style: "font-size:12.5px", text: s.title }),
      el("span", { class: "mono", style: `font-size:11px;color:var(${colour(s.status)})`, text: s.status }),
      el("span", { class: "hint", text: (s.reason || "").slice(0, 120) }))),
  ));
}

/* ═══ misc actions ════════════════════════════════════════════════ */
function addDataset() { sheet("Add dataset", el("div", { class: "stack" },
  el("p", { class: "hint", text: "Choose how the data reaches the workstation." }),
  el("div", { class: "btn-row" }, btn("secondary", "Browse files", browseFiles), btn("secondary", "Import URL", importUrl), btn("secondary", "Import accession", importAccession)))); }
async function rescanInputs() {
  S.inputs = await api("/api/inputs");
  renderWork();
  notify("success", "Folder rescanned", `${S.inputs.length} candidate input files found under data/, results/ and config/.`);
}

/* Upload. The tray used to accept a drop and then say it had not taken
   it, and Browse said the same: the only way in was to copy files into
   data/ by hand. The GUI asked for a sequence it had no way to receive.

   Staging is not adopting. Files land in data/uploads/ and are indexed;
   a pathogen's corpus changes by editing its config, never by someone
   dropping a file onto a panel. */
async function uploadFiles(fileList) {
  const files = Array.from(fileList || []);
  if (!files.length) return;
  notify("running", "Uploading", `${files.length} file(s)…`);
  const ok = [], failed = [];
  for (const f of files) {
    const body = new FormData();
    body.append("file", f);
    try {
      const r = await fetch("/api/upload", { method: "POST", body });
      const d = await r.json().catch(() => ({}));
      if (!r.ok) throw new Error(d.detail || `HTTP ${r.status}`);
      ok.push(d);
      log("success", `staged ${d.path} (${bytes(d.bytes)})`);
    } catch (e) {
      failed.push(`${f.name}: ${e.message}`);
      log("error", `upload failed — ${f.name}: ${e.message}`);
    }
  }
  if (ok.length) await rescanInputs();
  if (failed.length) {
    notify("error", `${failed.length} upload(s) refused`, failed.join(" · "));
  } else {
    notify("success", `${ok.length} file(s) staged`,
      "Written to data/uploads/ and indexed. Validate them, then point a config at them — "
      + "uploading does not change any pathogen's corpus.");
  }
}

function browseFiles() {
  const i = el("input", { type: "file", multiple: true });
  i.onchange = () => uploadFiles(i.files);
  i.click();
}
function dropped(e) { uploadFiles(e.dataTransfer.files); }
function clearStaged() { S.validations = {}; renderWork(); notify("success", "Cleared", "Validation results cleared for this session."); }
function importUrl() { sheet("Import from URL", el("div", { class: "stack" }, field("URL", el("input", { class: "input", placeholder: "https://…/sequences.fasta" })), el("p", { class: "hint", text: "Not wired in this build: the console has no outbound fetch capability by design." })), btn("tertiary", "Close", closeSheet)); }
function importAccession() { sheet("Import accession", el("div", { class: "stack" }, field("Accessions", el("textarea", { class: "textarea", placeholder: "AY593823.1\nPQ587570.1" })), el("p", { class: "hint", text: "Not wired: NCBI fetch would need network access and an API key. scripts/ holds the acquisition tooling." })), btn("tertiary", "Close", closeSheet)); }
function useAsInput(f) { S.params.reference = f.path; mark(); notify("success", "Set as input", f.name); setMode("configure"); }
async function previewFile(path) {
  try {
    const d = await api(`/api/artifact?path=${encodeURIComponent(path)}`);
    sheet(path, el("pre", { style: "font-family:var(--mono);font-size:11px;line-height:1.6;white-space:pre-wrap;word-break:break-word;margin:0", text: (d.text || "").slice(0, 20000) }));
  } catch (e) { notify("error", "Preview failed", e.message); }
}
function searchPosition() {
  sheet("Search position", el("div", { class: "stack" }, field("Genome coordinate", el("input", { class: "input", type: "number", placeholder: "4313",
    onchange: (e) => { const bp = +e.target.value; const hit = S.data.loci.find((l) => bp >= l.start && bp <= l.end); notify(hit ? "success" : "warning", hit ? "Locus found" : "No locus", hit ? `${hit.id} spans ${hit.start}–${hit.end}` : `No Atlas locus covers position ${bp}.`); closeSheet(); } }))));
}
function filterSheet() {
  const years = Object.keys(S.data.years).map(Number);
  sheet("Filters", el("div", { class: "stack" },
    field("Lineage", el("select", { class: "select", onchange: (e) => { if (e.target.value) S.sel.lineages = new Set([e.target.value]); else S.sel.lineages.clear(); } },
      el("option", { value: "", text: "All lineages" }), ...S.order.map((l) => el("option", { value: l, text: `${l} (${S.data.lineages[l]})` })))),
    field("Country", el("select", { class: "select", onchange: (e) => { if (e.target.value) S.sel.countries = new Set([e.target.value]); else S.sel.countries.clear(); } },
      el("option", { value: "", text: "All countries" }), ...Object.entries(S.data.countries).filter(([c]) => c !== "—").map(([c, n]) => el("option", { value: c, text: `${c} (${n})` })))),
    field("From year", el("input", { class: "input", type: "number", value: Math.min(...years), onchange: (e) => { S.sel.period = [+e.target.value, S.sel.period?.[1] ?? Math.max(...years)]; } })),
    field("To year", el("input", { class: "input", type: "number", value: Math.max(...years), onchange: (e) => { S.sel.period = [S.sel.period?.[0] ?? Math.min(...years), +e.target.value]; } }))),
    el("div", { class: "btn-row" }, btn("tertiary", "Reset", () => { clearSel(); closeSheet(); }), btn("primary", "Apply", () => { closeSheet(); renderWork(); })));
}
function runSearch(q) {
  if (!q || !S.data) return;
  const scope = $("#q-scope").value, hits = [];
  const t = q.toLowerCase();
  if (scope === "all" || scope === "samples") S.data.samples.filter((s) => s.a.toLowerCase().includes(t)).slice(0, 8).forEach((s) => hits.push(["sample", s.a, () => { S.sel.sample = s.a; setMode("visualize"); }]));
  if (scope === "all" || scope === "lineages") S.order.filter((l) => l.toLowerCase().includes(t)).forEach((l) => hits.push(["lineage", l, () => { S.sel.lineages = new Set([l]); setMode("visualize"); }]));
  if (scope === "all" || scope === "loci") S.data.loci.filter((l) => l.id.toLowerCase().includes(t)).forEach((l) => hits.push(["locus", l.id, () => { S.sel.locus = l.id; S.viz = "genome"; setMode("visualize"); }]));
  if (scope === "all" || scope === "countries") Object.keys(S.data.countries).filter((c) => c.toLowerCase().includes(t)).slice(0, 8).forEach((c) => hits.push(["country", c, () => { S.sel.countries = new Set([c]); S.viz = "matrix"; setMode("visualize"); }]));
  if (scope === "all" || scope === "files") S.inputs.filter((f) => f.name.toLowerCase().includes(t)).slice(0, 8).forEach((f) => hits.push(["file", f.name, () => { setMode("input"); validateOne(f.path); }]));
  if (!hits.length) return;
  sheet(`Search — ${hits.length} result${hits.length !== 1 ? "s" : ""} for "${q}"`,
    el("div", { class: "stack tight" }, ...hits.slice(0, 40).map(([kind, label, go]) =>
      el("button", { class: "btn tertiary", style: "justify-content:flex-start;width:100%", onclick: () => { closeSheet(); go(); } },
        el("span", { class: "dim", text: kind }), el("span", { text: label })))));
}
function helpSheet() {
  sheet("Help", el("div", { class: "stack" },
    el("p", { style: "font-size:13px;line-height:1.65", text: "G4 is a research framework, not a validated diagnostic. No output should drive a control decision on its own." }),
    el("div", { class: "btn-row" },
      btn("secondary", "Documentation", () => previewFile("README.md")),
      btn("secondary", "Methods", () => previewFile("docs/methods_supplement.md")),
      btn("secondary", "Usage", () => previewFile("docs/usage.md")),
      btn("tertiary", "Keyboard shortcuts", shortcutsSheet),
      btn("tertiary", "Report issue", () => notImplemented("Issue reporting needs a tracker URL; none configured."))),
    el("div", { class: "field", style: "max-width:220px" },
      el("label", { text: "Visual skin" }),
      el("select", { class: "select", onchange: (e) => applySkin(e.target.value) },
        ...[["forecast", "Forecast"], ["institute", "Institute"], ["redesign", "Dark"]].map(([v, label]) =>
          el("option", {
            value: v, text: label,
            selected: document.documentElement.getAttribute("data-skin") === v,
          }))))));
}
function shortcutsSheet() {
  sheet("Keyboard shortcuts", el("div", { class: "stack tight" },
    ...[["1–5", "Switch mode"], ["/", "Focus search"], ["l", "Toggle log console"], ["?", "Help"], ["Esc", "Close panel"],
      ["Click", "Select / brush"], ["Shift-click", "Select clade (phylogeny)"], ["Scroll or + / −", "Zoom"],
      ["Drag", "Pan"], ["Double-click", "Fit to screen"]]
      .map(([k, v]) => el("div", { style: "display:flex;justify-content:space-between;gap:16px;padding:5px 0;border-bottom:1px solid var(--hair)" },
        el("span", { class: "mono", style: "font-size:11.5px", text: k }), el("span", { class: "dim", text: v })))));
}
function notImplemented(msg) { notify("warning", "Not implemented in this build", msg); }
function download(blob, name) {
  const u = URL.createObjectURL(blob), a = el("a", { href: u, download: name });
  document.body.append(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(u), 800);
}
function debounce(f, ms) { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => f(...a), ms); }; }

boot();
