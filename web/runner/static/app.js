/* G4-WATCH operator console — front end.
 *
 * No framework and no build step: the console has to work from a plain
 * `uvicorn` process on a machine that may have no network access.
 */

const $ = (sel) => document.querySelector(sel);
const el = (tag, attrs = {}, ...kids) => {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === "class") node.className = v;
    else if (k === "text") node.textContent = v;
    else if (k.startsWith("on")) node.addEventListener(k.slice(2), v);
    else if (v !== null && v !== undefined) node.setAttribute(k, v);
  }
  for (const kid of kids) if (kid) node.append(kid);
  return node;
};

const state = {
  commands: [],
  pathogens: [],
  tools: {},
  pathogen: null,
  jobId: null,
  stream: null,
  jobs: [],
  ticker: null,
  lineNo: 0,
};

async function api(path, options) {
  const res = await fetch(path, options);
  if (!res.ok) {
    let detail = res.statusText;
    try { detail = (await res.json()).detail || detail; } catch { /* non-JSON error body */ }
    throw new Error(detail);
  }
  return res.json();
}

function alert_(message) {
  const node = el("div", { class: "alert", text: message });
  document.body.append(node);
  setTimeout(() => node.remove(), 5500);
}

// ── boot ─────────────────────────────────────────────────────────
async function boot() {
  try {
    const [env, commands, pathogens] = await Promise.all([
      api("/api/env"), api("/api/commands"), api("/api/pathogens"),
    ]);
    state.tools = env.tools;
    state.commands = commands;
    state.pathogens = pathogens;
    $("#version").textContent = "operator console · v" + env.version;
    renderLamps();
    renderPathogens();
    renderProtocol();
    await Promise.all([refreshGate(), refreshJobs(), refreshArtifacts()]);
  } catch (err) {
    alert_("Console API unreachable: " + err.message);
  }
}

function renderLamps() {
  const strip = $("#tools");
  strip.replaceChildren();
  for (const [name, info] of Object.entries(state.tools)) {
    strip.append(el("span", {
      class: "lamp " + (info.present ? "on" : "off"),
      title: info.present ? `${info.path}\n${info.stage}` : `Not installed — required for ${info.stage}`,
    },
      el("i", {}),
      el("span", { text: name }),
    ));
  }
}

function renderPathogens() {
  const select = $("#pathogen");
  select.replaceChildren();
  for (const p of state.pathogens) {
    select.append(el("option", {
      value: p.name,
      text: p.name.toUpperCase() + (p.provisioned ? "" : " · not provisioned"),
    }));
  }
  const provisioned = state.pathogens.find((p) => p.provisioned);
  state.pathogen = provisioned ? provisioned.name : (state.pathogens[0] || {}).name;
  select.value = state.pathogen;
  select.addEventListener("change", () => {
    state.pathogen = select.value;
    refreshGate();
    renderProtocol();
  });
}

// ── gate stamp ───────────────────────────────────────────────────
async function refreshGate() {
  const host = $("#gate");
  if (!state.pathogen) { host.replaceChildren(); return; }
  try {
    const gate = await api(`/api/gate/${state.pathogen}`);
    const blocked = !gate.scoring_permitted;
    host.replaceChildren(el("div", { class: "stamp " + (blocked ? "blocked" : "open") },
      el("p", { class: "stamp-top", text: state.pathogen.toUpperCase() + " · stage 5/6 authorisation" }),
      el("div", { class: "stamp-body" },
        el("p", { class: "stamp-verdict", text: blocked ? "Scoring blocked" : "Scoring permitted" }),
        el("p", { class: "stamp-code", text: gate.permission }),
        el("p", { class: "stamp-note", text: gate.explanation }),
        gate.failing_checks.length
          ? el("p", { class: "stamp-checks" },
              el("b", { text: "failing — " }),
              el("span", { text: gate.failing_checks.join(", ") }))
          : null,
      ),
    ));
  } catch (err) {
    host.replaceChildren(el("div", { class: "stamp unknown" },
      el("p", { class: "stamp-top", text: "gate" }),
      el("div", { class: "stamp-body" }, el("p", { class: "stamp-note", text: "Unavailable: " + err.message })),
    ));
  }
}

// ── protocol ─────────────────────────────────────────────────────
function toolMissing(command) {
  return command.requires_tools.filter((t) => !(state.tools[t] || {}).present);
}

/* Stage commands carry a real ordinal; utilities do not. The margin
   numeral is information, so only the sequence gets one. */
function marker(command) {
  const match = /^Stage\s+([\d.]+)/.exec(command.stage);
  if (match) return { text: match[1], util: false };
  return { text: "▫", util: true };
}

function renderProtocol() {
  const host = $("#stages");
  host.replaceChildren();

  const phases = new Map();
  for (const command of state.commands) {
    if (!phases.has(command.stage)) phases.set(command.stage, []);
    phases.get(command.stage).push(command);
  }
  for (const [phase, commands] of phases) {
    const block = el("div", { class: "phase" }, el("span", { class: "lab", text: phase }));
    for (const command of commands) block.append(step(command));
    host.append(block);
  }
}

function step(command) {
  const missing = toolMissing(command);
  const values = {};
  const mark = marker(command);

  const row = el("div", { class: "step" + (missing.length ? " off" : "") });

  const name = el("button", { class: "step-name", type: "button", onclick: () => row.classList.toggle("open") },
    el("span", { text: command.title }));
  if (command.gate_aware || command.mutates) {
    const flags = el("span", { class: "flags" });
    if (command.gate_aware) flags.append(el("span", { class: "flag gate", text: "gated" }));
    if (command.mutates) flags.append(el("span", { class: "flag write", text: "writes" }));
    name.append(flags);
  }

  const runKey = el("button", {
    class: "key primary", type: "button", text: "Run",
    disabled: missing.length ? "disabled" : null,
    title: missing.length ? `${missing.join(", ")} not installed` : "",
    onclick: () => launch(command, values),
  });

  row.append(
    el("span", { class: "step-no" + (mark.util ? " util" : ""), text: mark.text }),
    el("div", { class: "step-main" }, name, el("p", { class: "step-desc", text: command.summary })),
    runKey,
  );

  if (command.options.length || missing.length || command.danger_note) {
    const inner = el("div", { class: "params-inner" });

    if (missing.length) {
      inner.append(el("p", { class: "caution", text: `${missing.join(", ")} is not installed. See docs/installation.md.` }));
    }
    if (command.danger_note) {
      inner.append(el("p", { class: "caution", text: command.danger_note }));
    }

    for (const option of command.options) {
      const id = `${command.key}-${option.name}`;
      if (option.type === "flag") {
        if (option.default === true) values[option.name] = true;
        inner.append(el("div", { class: "param check" },
          el("input", {
            type: "checkbox", id,
            checked: option.default === true ? "checked" : null,
            onchange: (ev) => { values[option.name] = ev.target.checked; },
          }),
          el("div", {},
            el("label", { for: id, text: option.label }),
            option.help ? el("span", { class: "hint", text: option.help }) : null,
          ),
        ));
      } else {
        inner.append(el("div", { class: "param" },
          el("label", { for: id, text: option.label }),
          el("input", {
            type: "text", id, placeholder: option.placeholder || "",
            oninput: (ev) => {
              const v = ev.target.value.trim();
              if (v === "") delete values[option.name]; else values[option.name] = v;
            },
          }),
          option.help ? el("span", { class: "hint", text: option.help }) : null,
        ));
      }
    }
    row.append(el("div", { class: "params" }, inner));
  }
  return row;
}

// ── launching ────────────────────────────────────────────────────
async function launch(command, values) {
  const payload = {
    command: command.key,
    pathogen: command.needs_pathogen ? state.pathogen : null,
    options: { ...values },
  };

  if (command.mutates) {
    const writesLedger = command.key === "dh1" && values.no_ledger !== true;
    const note = writesLedger
      ? "This appends permanent rows to the append-only study ledger."
      : command.danger_note || "This writes pipeline artifacts.";
    if (!(await confirmRun(command, payload, note))) return;
  }

  try {
    const job = await api("/api/run", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    await refreshJobs();
    selectJob(job.id);
  } catch (err) {
    alert_(err.message);
  }
}

function confirmRun(command, payload, note) {
  return new Promise((resolve) => {
    const preview = ["g4watch", command.key, payload.pathogen ? `--pathogen ${payload.pathogen}` : "",
      ...Object.entries(payload.options).map(([k, v]) => (v === true ? `--${k}` : v === false ? "" : `--${k} ${v}`)),
    ].filter(Boolean).join(" ");

    const close = (answer) => { scrim.remove(); document.removeEventListener("keydown", onKey); resolve(answer); };
    const onKey = (e) => { if (e.key === "Escape") close(false); };

    const scrim = el("div", { class: "scrim", onclick: (e) => { if (e.target === scrim) close(false); } },
      el("div", { class: "sheet" },
        el("header", {}, el("h3", { text: "Confirm — " + command.title })),
        el("div", { class: "content" },
          el("p", { class: "prose" },
            el("span", { text: note }),
            el("code", { class: "argv", text: preview }),
          ),
        ),
        el("footer", {},
          el("button", { class: "key", type: "button", text: "Cancel", onclick: () => close(false) }),
          el("button", { class: "key primary", type: "button", text: "Run it", onclick: () => close(true) }),
        ),
      ),
    );
    document.body.append(scrim);
    document.addEventListener("keydown", onKey);
  });
}

// ── job selection and streaming ──────────────────────────────────
async function selectJob(jobId) {
  if (state.stream) { state.stream.close(); state.stream = null; }
  stopTicker();
  state.jobId = jobId;
  state.lineNo = 0;
  renderJobs();

  $("#console").replaceChildren();

  let job;
  try { job = await api(`/api/jobs/${jobId}`); }
  catch (err) { alert_(err.message); return; }

  setHeader(job);
  for (const line of job.lines) appendLine(line);

  if (["queued", "running"].includes(job.state)) {
    startTicker(job);
    const stream = new EventSource(`/api/jobs/${jobId}/stream`);
    state.stream = stream;
    const seen = new Set(job.lines.map((l) => l.seq));
    stream.onmessage = (event) => {
      const data = JSON.parse(event.data);
      if (data.type === "line") {
        if (seen.has(data.seq)) return;
        seen.add(data.seq);
        appendLine(data);
      } else if (data.type === "end") {
        stopTicker();
        setHeader(data);
        stream.close();
        state.stream = null;
        refreshJobs(); refreshArtifacts(); refreshGate();
      }
    };
    stream.onerror = () => { stream.close(); state.stream = null; stopTicker(); };
  }
}

const MARK = { stderr: "!", meta: "›", stdout: "" };

function appendLine(line) {
  const out = $("#console");
  state.lineNo += 1;
  out.append(el("div", { class: "pl " + line.stream },
    el("span", { class: "n", text: String(state.lineNo) }),
    el("span", { class: "m", text: MARK[line.stream] ?? "" }),
    el("span", { class: "t", text: line.text }),
  ));
  if ($("#follow").checked) out.scrollTop = out.scrollHeight;
}

function setHeader(job) {
  const annun = $("#job-state");
  annun.className = "annun " + job.state;
  annun.textContent = job.state.replace("_", " ");

  $("#job-cmd").textContent = job.argv.join(" ");
  $("#job-cmd").title = job.argv.join(" ");
  $("#exit").textContent = (job.exit_code === null || job.exit_code === undefined) ? "—" : String(job.exit_code);
  $("#elapsed").textContent = job.duration ? job.duration.toFixed(1) + " s" : "—";

  const cancel = $("#cancel");
  cancel.disabled = !["queued", "running"].includes(job.state);
  cancel.onclick = () => api(`/api/jobs/${job.id}/cancel`, { method: "POST" }).catch((e) => alert_(e.message));

  $("#job-meta").textContent = job.pathogen ? job.pathogen.toUpperCase() : "";

  document.querySelector(".verdict-band")?.remove();
  if (job.state === "gate_closed") {
    $("#console").after(el("div", { class: "verdict-band" },
      el("b", { text: "Exit 3 is the correct outcome. " }),
      el("span", { text: "The D.H1 gate is closed, so Stage 5/6 declined to run. This is not a pipeline failure — it means the evidence does not yet permit scoring." }),
    ));
  }
}

function startTicker(job) {
  const t0 = (job.started || job.created) * 1000;
  state.ticker = setInterval(() => {
    $("#elapsed").textContent = ((Date.now() - t0) / 1000).toFixed(1) + " s";
  }, 100);
}
function stopTicker() {
  if (state.ticker) { clearInterval(state.ticker); state.ticker = null; }
}

// ── run log ──────────────────────────────────────────────────────
async function refreshJobs() {
  try { state.jobs = await api("/api/jobs"); } catch { return; }
  renderJobs();
}

function renderJobs() {
  const host = $("#jobs");
  host.replaceChildren();
  if (!state.jobs.length) {
    host.append(el("p", { class: "blank", text: "Nothing run in this session." }));
    return;
  }
  for (const job of state.jobs) {
    const when = new Date(job.created * 1000).toLocaleTimeString([], { hour12: false });
    const bits = [when];
    if (job.pathogen) bits.push(job.pathogen);
    if (job.duration) bits.push(job.duration.toFixed(1) + "s");
    host.append(el("div", {
      class: "run" + (job.id === state.jobId ? " sel" : ""),
      onclick: () => selectJob(job.id),
    },
      el("span", { class: "rn", text: job.title }),
      el("span", { class: "rs " + job.state, text: job.state.replace("_", " ") }),
      el("span", { class: "rm", text: bits.join("  ·  ") }),
    ));
  }
}

// ── artifacts ────────────────────────────────────────────────────
async function refreshArtifacts() {
  const host = $("#artifacts");
  let files;
  try { files = await api("/api/artifacts"); } catch { return; }
  host.replaceChildren();
  if (!files.length) {
    host.append(el("p", { class: "blank", text: "No outputs yet. Files written by a run appear here." }));
    return;
  }
  for (const file of files) {
    host.append(el("div", { class: "file", onclick: () => openArtifact(file.path) },
      el("span", { class: "fp", title: file.path, text: file.path }),
      el("span", { class: "fs", text: humanSize(file.size) }),
    ));
  }
}

function humanSize(bytes) {
  if (bytes < 1024) return bytes + " B";
  if (bytes < 1024 * 1024) return (bytes / 1024).toFixed(0) + " KB";
  return (bytes / 1024 / 1024).toFixed(1) + " MB";
}

async function openArtifact(path) {
  let data;
  try { data = await api(`/api/artifact?path=${encodeURIComponent(path)}`); }
  catch (err) { alert_(err.message); return; }

  const scrim = el("div", { class: "scrim", onclick: (e) => { if (e.target === scrim) scrim.remove(); } },
    el("div", { class: "sheet" },
      el("header", {},
        el("h3", { text: path }),
        el("button", { class: "key sm", type: "button", text: "Close", onclick: () => scrim.remove() }),
      ),
      el("div", { class: "content mono", text: data.text }),
    ),
  );
  document.body.append(scrim);
}

$("#refresh-gate").addEventListener("click", refreshGate);
$("#refresh-artifacts").addEventListener("click", refreshArtifacts);

boot();
