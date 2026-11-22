const palette = {
  core: "#327b78",
  accessory: "#c48d56",
  rare: "#8177bd",
  singleton: "#cc716a",
};
const labels = {
  core: "Core",
  accessory: "Accessory",
  rare: "Rare",
  singleton: "Singleton",
};
let data;
let treeData, statsData, similarityData;
let datasetId =
  new URLSearchParams(location.search).get("dataset") || "demo-lp-v1";
let selectedGenomeIds =
  new URLSearchParams(location.search)
    .get("strains")
    ?.split(",")
    .filter(Boolean) || [];
const $ = (selector, root = document) => root.querySelector(selector);
const $$ = (selector, root = document) => [...root.querySelectorAll(selector)];
const fmt = (n) => new Intl.NumberFormat("en-US").format(n);
const escapeHtml = (value) =>
  String(value).replace(
    /[&<>"']/g,
    (character) =>
      ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[
        character
      ],
  );

async function init() {
  try {
    const catalogResponse = await fetch("/api/datasets");
    const catalog = await catalogResponse.json();
    $("#dataset-select").innerHTML = catalog.datasets
      .map(
        (d) =>
          `<option value="${escapeHtml(d.id)}">${d.mode === "real" ? "● " : ""}${d.mode === "real" ? "Real · " : "Demo · "}${escapeHtml(d.name)}</option>`,
      )
      .join("");
    $("#dataset-select").value = datasetId;
    const query = `?dataset=${encodeURIComponent(datasetId)}`;
    const responses = await Promise.all(
      ["/api/overview", "/api/tree", "/api/statistics", "/api/similarity"].map(
        (url) => fetch(url + query),
      ),
    );
    if (responses.some((response) => !response.ok))
      throw new Error("Atlas API request failed");
    [data, treeData, statsData, similarityData] = await Promise.all(
      responses.map((response) => response.json()),
    );
    const validIds = new Set(data.strains.map((strain) => strain.id));
    selectedGenomeIds = selectedGenomeIds.filter((id) => validIds.has(id));
    render();
    bind();
  } catch {
    $(".content").innerHTML =
      '<div class="error-state"><h1>Atlas could not connect</h1><p>Start the local API with <code>make dev</code>, then reload this page.</p></div>';
  }
}
function render() {
  $("#crumb-current").textContent = "Overview";
  updateDatasetContext();
  const real = data.dataset.mode === "real";
  $("#metrics").innerHTML = [
    metric(
      "COHORT SIZE",
      fmt(data.strains.length),
      real ? "validated RefSeq genomes" : "strains in demo cohort",
      "◉",
      "teal",
    ),
    metric(
      "GENE FAMILIES",
      fmt(data.families.length),
      real ? "Panaroo gene families" : "illustrative family rows",
      "⌘",
      "violet",
    ),
    metric(
      "CORE FAMILIES",
      fmt(data.counts.core),
      `${((data.counts.core / data.families.length) * 100).toFixed(1)}% of ${real ? "analyzed families" : "demo families"}`,
      "◎",
      "blue",
    ),
    metric(
      "FAMILY COPIES",
      fmt(data.strains.reduce((a, s) => a + s.genes, 0)),
      "sum across strains",
      "⌁",
      "amber",
    ),
  ].join("");
  $("#family-total").textContent = data.families.length;
  $("#family-count-tag").textContent = `${fmt(data.families.length)} families`;
  $("#strain-count").textContent = data.strains.length;
  drawDonut();
  drawLandscape();
  drawCircular();
  renderStrains();
  renderFamilies();
  renderMatrix();
  renderTree();
  renderStatistics();
  renderCompareOptions();
}
function metric(label, value, note, icon, tone) {
  return `<article class="metric"><div class="metric-top"><span>${label}</span><i class="metric-icon ${tone}">${icon}</i></div><strong>${value}</strong><small>${note}</small><div class="metric-rule ${tone}"></div></article>`;
}
function drawDonut() {
  const keys = ["core", "accessory", "rare", "singleton"],
    total = data.families.length;
  let start = -Math.PI / 2;
  const arcs = keys.map((k) => {
    const v = data.counts[k],
      end = start + (2 * Math.PI * v) / total,
      path = arcPath(110, 110, 78, start, end);
    const part = `<path d="${path}" fill="none" stroke="${palette[k]}" stroke-width="20" stroke-linecap="butt"><title>${labels[k]}: ${v} families</title></path>`;
    start = end;
    return part;
  });
  $("#donut").innerHTML =
    `<circle cx="110" cy="110" r="78" fill="none" stroke="var(--line-soft)" stroke-width="20"/>${arcs.join("")}`;
  $("#legend").innerHTML = keys
    .map(
      (k) =>
        `<div class="legend-row"><span class="legend-name"><i class="legend-dot" style="--dot:${palette[k]}"></i>${labels[k]}</span><strong>${data.counts[k]}</strong><small>${((100 * data.counts[k]) / total).toFixed(0)}%</small></div>`,
    )
    .join("");
}
function arcPath(cx, cy, r, a0, a1) {
  const x0 = cx + r * Math.cos(a0),
    y0 = cy + r * Math.sin(a0),
    x1 = cx + r * Math.cos(a1),
    y1 = cy + r * Math.sin(a1);
  return `M ${x0} ${y0} A ${r} ${r} 0 ${a1 - a0 > Math.PI ? 1 : 0} 1 ${x1} ${y1}`;
}
function drawLandscape() {
  if (!data.families.every((f) => Number.isFinite(f.position))) {
    $("#landscape-chart").innerHTML =
      '<div class="unavailable-chart"><strong>Coordinate track unavailable</strong><span>This dataset has no normalized family-to-locus coordinate map.</span></div>';
    return;
  }
  const width = 690,
    height = 170,
    left = 34,
    right = 12,
    base = 145,
    scale = (width - left - right) / data.genome_length;
  const guides = [0, 1, 2, 3, 4]
    .map(
      (m) =>
        `<line x1="${left + m * 1e6 * scale}" y1="16" x2="${left + m * 1e6 * scale}" y2="${base}" class="chart-guide"/>`,
    )
    .join("");
  const pts = data.families
    .map((f) => {
      const x = left + f.position * scale,
        y = base - f.prevalence * 105,
        r = f.group === "core" ? 2.8 : 2.1;
      return `<circle cx="${x.toFixed(1)}" cy="${y.toFixed(1)}" r="${r}" fill="${palette[f.group]}" opacity=".85"><title>${escapeHtml(f.id)} · ${labels[f.group]} · ${(f.prevalence * 100).toFixed(0)}% prevalence</title></circle>`;
    })
    .join("");
  $("#landscape-chart").innerHTML =
    `<svg viewBox="0 0 ${width} ${height}" role="img" aria-label="Synthetic family prevalence positions">${guides}<line x1="${left}" x2="${width - right}" y1="${base}" y2="${base}" class="chart-axis"/><path d="M${left},${base - 105} H${width - right}" class="chart-threshold"/>${pts}<text x="${left - 6}" y="${base - 101}" text-anchor="end" class="axis-label">100%</text></svg>`;
}
function drawCircular() {
  if (!data.families.every((f) => Number.isFinite(f.position))) {
    $("#circular-genome").innerHTML =
      '<div class="unavailable-chart"><strong>No circular coordinate view for this analysis</strong><span>Family presence does not provide chromosome positions. Import gene loci and contig-aware coordinates to enable this track.</span></div>';
    $("#genome-legend").innerHTML = "";
    $("#genome-track-kicker").textContent = "COORDINATE DATA UNAVAILABLE";
    $("#genome-footnote").textContent =
      "The current family matrix has no normalized family-to-locus coordinate map; no genome markers are drawn.";
    return;
  }
  const cx = 210,
    cy = 210,
    r = 124,
    rin = 89;
  let out = `<svg viewBox="0 0 420 420" role="img" aria-label="Circular illustrative family coordinate map"><circle cx="${cx}" cy="${cy}" r="${r}" fill="none" stroke="var(--line-soft)" stroke-width="1"/><circle cx="${cx}" cy="${cy}" r="${rin}" fill="none" stroke="var(--line-soft)" stroke-width="1"/>`;
  for (let i = 0; i < 12; i++) {
    const a = (i * Math.PI) / 6 - Math.PI / 2;
    const x1 = cx + rin * Math.cos(a),
      y1 = cy + rin * Math.sin(a),
      x2 = cx + (r + 7) * Math.cos(a),
      y2 = cy + (r + 7) * Math.sin(a);
    out += `<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" class="chart-guide"/>`;
  }
  for (const f of data.families) {
    const a = (2 * Math.PI * f.position) / data.genome_length - Math.PI / 2;
    const x1 = cx + rin * Math.cos(a),
      y1 = cy + rin * Math.sin(a),
      x2 = cx + r * Math.cos(a),
      y2 = cy + r * Math.sin(a);
    out += `<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="${palette[f.group]}" stroke-width="${f.group === "core" ? 2.5 : 1.6}" opacity=".85"><title>${escapeHtml(f.id)} · ${labels[f.group]}</title></line>`;
  }
  out += `<text x="${cx}" y="${cy - 5}" text-anchor="middle" class="circle-center">4.45 Mb</text><text x="${cx}" y="${cy + 15}" text-anchor="middle" class="circle-sub">SYNTHETIC AXIS</text></svg>`;
  $("#circular-genome").innerHTML = out;
  $("#genome-legend").innerHTML = Object.keys(labels)
    .map(
      (k) =>
        `<span><i class="legend-dot" style="--dot:${palette[k]}"></i>${labels[k]}</span>`,
    )
    .join("");
}
function renderStrains() {
  const real = data.dataset.mode === "real";
  $("#view-strains .version").textContent = `${data.strains.length} GENOMES`;
  $("#view-strains .lede").textContent = real
    ? "Compare family counts, assembly accessions and source records across the analyzed RefSeq cohort."
    : "Compare family counts and provenance fields across the illustrative cohort.";
  $("#view-strains .callout strong").textContent = real
    ? "RefSeq assembly accessions"
    : "Not accession records";
  $("#view-strains .callout span:not(.callout-icon)").textContent = real
    ? "Assembly accessions and NCBI source records are shown with the analyzed genomes. Selection exclusions are recorded in the local manifest."
    : "Strain identifiers and source labels are synthetic placeholders, not NCBI accessions or isolate metadata.";
  $("#strain-strip").innerHTML = data.strains
    .slice(0, 6)
    .map(
      (s, i) =>
        `<button class="strain-card" data-strain="${escapeHtml(s.id)}"><span class="strain-index">${String(i + 1).padStart(2, "0")}</span><span><strong>${escapeHtml(s.id)}</strong><small>${escapeHtml(s.source)} label</small></span><b>${s.families}<small> fam.</small></b></button>`,
    )
    .join("");
  $("#strain-table").innerHTML = data.strains
    .map(
      (s) =>
        `<tr><td class="mono primary-cell">${escapeHtml(s.id)}</td><td>${escapeHtml(s.name)}</td><td><span class="source-pill">${escapeHtml(s.source)}</span></td><td>${s.families}</td><td>${s.genes}</td><td class="mono">${(s.length / 1e6).toFixed(2)} Mb</td></tr>`,
    )
    .join("");
}
function renderMatrix() {
  const filter = $("#class-filter")?.value || "all",
    query = ($("#family-search")?.value || "").toLowerCase();
  const fs = data.families.filter(
    (f) =>
      (filter === "all" || f.group === filter) &&
      f.id.toLowerCase().includes(query),
  );
  const genomes = selectedGenomeIds.length
    ? data.strains.filter((s) => selectedGenomeIds.includes(s.id))
    : data.strains;
  const head = `<div class="matrix-row matrix-header" style="--n:${genomes.length}"><div class="matrix-sticky">FAMILY <span>↓</span></div>${genomes.map((s) => `<div class="strain-col" title="${escapeHtml(s.id)}">${escapeHtml(s.id.slice(-2))}</div>`).join("")}</div>`;
  const rows = fs
    .map(
      (f) =>
        `<button class="matrix-row family-row ${f.id === selectedFamily ? "selected" : ""}" style="--n:${genomes.length}" data-family="${escapeHtml(f.id)}" aria-label="Select ${escapeHtml(f.id)}"><div class="matrix-sticky"><b class="class-tick" style="--dot:${palette[f.group]}"></b><strong>${escapeHtml(f.id)}</strong><small>${labels[f.group]}</small></div>${genomes
          .map((s) => {
            const c = f.copies[s.id] || 0;
            return `<span class="matrix-cell ${c ? "is-present" : ""} ${c > 1 ? "is-multi" : ""}" title="${escapeHtml(s.id)}: ${c ? `${c} copy${c > 1 ? "ies" : ""}` : "absent"}">${c > 1 ? c : ""}</span>`;
          })
          .join("")}</button>`,
    )
    .join("");
  $("#matrix-wrap").innerHTML =
    `<div class="matrix-grid" style="--n:${genomes.length}">${head}${rows || '<div class="no-results">No families match this search.</div>'}</div>`;
  $("#matrix-count").textContent =
    `Showing ${fs.length} of ${data.families.length} families · ${genomes.length} of ${data.strains.length} strains`;
  $$(".family-row").forEach((el) =>
    el.addEventListener("click", () => selectFamily(el.dataset.family)),
  );
}
let selectedFamily = null;
function selectFamily(id) {
  selectedFamily = id;
  const f = data.families.find((x) => x.id === id);
  if (!f) return;
  renderMatrix();
  const present = f.members.length;
  $("#family-detail").innerHTML =
    `<div class="detail-main"><span class="section-kicker">SELECTED FAMILY</span><strong>${escapeHtml(f.id)}</strong><span class="detail-class" style="--dot:${palette[f.group]}">${labels[f.group]}</span></div><div class="detail-stat"><span>PREVALENCE</span><strong>${(100 * f.prevalence).toFixed(1)}%</strong></div><div class="detail-stat"><span>GENOMES</span><strong>${present} <small>/ ${data.strains.length}</small></strong></div><div class="detail-members"><span>MEMBERS</span><div>${f.members.map((id) => `<span>${escapeHtml(id)}</span>`).join("") || "<small>None</small>"}</div></div>`;
}
function renderFamilies() {
  const real = data.dataset.mode === "real";
  $("#view-families .version").textContent =
    `${fmt(data.families.length)} FAMILIES`;
  $("#view-families .lede").textContent =
    `Filter the ${real ? "Panaroo-derived" : "demo"} family catalogue and inspect prevalence, copy number, and genome membership.`;
  const filter = $("#family-class-filter")?.value || "all",
    query = ($("#catalog-search")?.value || "").toLowerCase();
  const fs = data.families.filter(
    (f) =>
      (filter === "all" || f.group === filter) &&
      f.id.toLowerCase().includes(query),
  );
  $("#family-table").innerHTML = fs
    .map(
      (f) =>
        `<tr><td><button class="table-family" data-family="${escapeHtml(f.id)}">${escapeHtml(f.id)}</button></td><td><span class="detail-class" style="--dot:${palette[f.group]}">${labels[f.group]}</span></td><td>${f.members.length} / ${data.strains.length}</td><td>${(100 * f.prevalence).toFixed(1)}%</td><td>${Object.values(f.copies).reduce((a, b) => a + b, 0)}</td></tr>`,
    )
    .join("");
  $$(".table-family").forEach((b) =>
    b.addEventListener("click", () => {
      selectFamily(b.dataset.family);
      openView("matrix");
    }),
  );
}
function renderTree() {
  const width = 880,
    top = 28,
    row = 31,
    height = top * 2 + (treeData.root.members.length - 1) * row,
    maxH = Math.max(treeData.root.height, 1e-9),
    x0 = 42,
    x1 = width - 182;
  const yById = new Map(data.strains.map((s, i) => [s.id, top + i * row]));
  let svg = `<svg viewBox="0 0 ${width} ${height}" role="img" aria-label="UPGMA clustering by gene-content Jaccard distance"><line x1="${x0}" y1="${top - 15}" x2="${x1}" y2="${top - 15}" class="tree-axis"/><text x="${x0}" y="${top - 20}" class="tree-scale">${maxH.toFixed(3)} distance</text><text x="${x1}" y="${top - 20}" text-anchor="end" class="tree-scale">0 distance</text>`;
  const drawn = [];
  const visit = (node) => {
    if (node.name) {
      const y = yById.get(node.name),
        x = x1,
        selected = selectedGenomeIds.includes(node.name);
      svg += `<g class="tree-leaf ${selected ? "selected" : ""}" data-genome="${escapeHtml(node.name)}" tabindex="0" role="button" aria-label="Toggle ${escapeHtml(node.name)} selection"><circle cx="${x}" cy="${y}" r="4"/><text x="${x + 11}" y="${y + 4}">${escapeHtml(node.name)}</text></g>`;
      return { x, y, min: y, max: y };
    }
    const x = x0 + (1 - node.height / maxH) * (x1 - x0);
    const children = node.children.map((c) => visit(c.node));
    const y = children.reduce((sum, c) => sum + c.y, 0) / children.length;
    svg += `<path d="M${x},${children[0].y} V${children.at(-1).y}" class="tree-branch"/>`;
    for (const c of children)
      svg += `<path d="M${x},${c.y} H${c.x}" class="tree-branch"/>`;
    return { x, y, min: children[0].min, max: children.at(-1).max };
  };
  visit(treeData.root);
  svg += "</svg>";
  $("#tree-canvas").innerHTML = svg;
  $$(".tree-leaf").forEach((el) => {
    const toggle = () => toggleGenome(el.dataset.genome);
    el.addEventListener("click", toggle);
    el.addEventListener("keydown", (e) => {
      if (e.key === "Enter" || e.key === " ") {
        e.preventDefault();
        toggle();
      }
    });
  });
  $("#tree-selection-label").textContent = selectedGenomeIds.length
    ? `${selectedGenomeIds.length} selected: ${selectedGenomeIds.join(", ")}`
    : "No strains selected";
}
function toggleGenome(id) {
  selectedGenomeIds = selectedGenomeIds.includes(id)
    ? selectedGenomeIds.filter((x) => x !== id)
    : [...selectedGenomeIds, id];
  syncSelection();
  renderTree();
  renderCompareOptions();
  renderStatistics();
}
function syncSelection() {
  const url = new URL(location.href);
  selectedGenomeIds.length
    ? url.searchParams.set("strains", selectedGenomeIds.join(","))
    : url.searchParams.delete("strains");
  history.replaceState(null, "", url);
}
async function renderStatistics() {
  let current = statsData;
  if (
    selectedGenomeIds.length &&
    selectedGenomeIds.length < data.strains.length
  ) {
    const query = new URLSearchParams({
      dataset: datasetId,
      genomes: selectedGenomeIds.join(","),
    });
    const response = await fetch(`/api/statistics?${query}`);
    if (!response.ok)
      throw new Error("Could not calculate statistics for selected genomes");
    current = await response.json();
  }
  const section = $("#view-statistics .section-kicker");
  if (section)
    section.textContent = `PANGENOME STATISTICS · ${current.selection_mode.toUpperCase()}`;
  const frequency = current.family_frequency;
  drawBars(
    "#frequency-chart",
    frequency.map((x) => ({ label: x.genomes_present, value: x.families })),
    { color: "#327b78", xLabel: "Genomes present", yLabel: "Families" },
  );
  drawLines(
    "#accumulation-chart",
    current.accumulation.map((x) => ({
      x: x.genomes_sampled,
      mean: x.pangenome_mean,
      min: x.pangenome_min,
      max: x.pangenome_max,
    })),
    { color: "#327b78", xLabel: "Genomes sampled", yLabel: "Families" },
  );
  drawLines(
    "#core-chart",
    current.accumulation.map((x) => ({
      x: x.genomes_sampled,
      mean: x.core_mean,
      min: x.core_min,
      max: x.core_max,
    })),
    { color: "#8177bd", xLabel: "Genomes sampled", yLabel: "Families" },
  );
  const names = similarityData.genomes;
  $("#similarity-heatmap").innerHTML =
    `<table class="similarity-table"><thead><tr><th></th>${names.map((id) => `<th>${escapeHtml(id.slice(-2))}</th>`).join("")}</tr></thead><tbody>${similarityData.matrix.map((values, i) => `<tr><th>${escapeHtml(names[i])}</th>${values.map((v, j) => `<td style="--similarity:${v}" title="${escapeHtml(names[i])} × ${escapeHtml(names[j])}: ${(v * 100).toFixed(1)}%"><span>${v.toFixed(2)}</span></td>`).join("")}</tr>`).join("")}</tbody></table>`;
}
function drawBars(selector, values, { color, xLabel, yLabel }) {
  const width = 560,
    height = 250,
    l = 42,
    r = 12,
    t = 18,
    b = 35,
    max = Math.max(...values.map((d) => d.value), 1),
    innerW = width - l - r,
    innerH = height - t - b,
    gap = 3,
    bw = Math.max(3, innerW / values.length - gap);
  const bars = values
    .map((d, i) => {
      const x = l + (i * innerW) / values.length,
        y = t + innerH - (d.value / max) * innerH;
      return `<rect x="${x}" y="${y}" width="${bw}" height="${(d.value / max) * innerH}" rx="2" fill="${color}"><title>${d.value} families present in ${d.label} genomes</title></rect>`;
    })
    .join("");
  $(selector).innerHTML =
    `<svg viewBox="0 0 ${width} ${height}" role="img" aria-label="${yLabel} by ${xLabel}"><line x1="${l}" y1="${t + innerH}" x2="${width - r}" y2="${t + innerH}" class="chart-axis"/>${bars}<text x="${l}" y="${height - 9}" class="axis-label">1</text><text x="${width - r}" y="${height - 9}" text-anchor="end" class="axis-label">${values.length} ${xLabel.toLowerCase()}</text></svg>`;
}
function drawLines(selector, values, { color, xLabel, yLabel }) {
  const width = 560,
    height = 250,
    l = 42,
    r = 12,
    t = 18,
    b = 35,
    innerW = width - l - r,
    innerH = height - t - b,
    max = Math.max(...values.map((d) => d.max), 1),
    min = 0;
  const px = (i) => l + (i * innerW) / Math.max(values.length - 1, 1),
    py = (v) => t + innerH - ((v - min) / (max - min || 1)) * innerH;
  const area = `M${values.map((d, i) => `${px(i)},${py(d.min)}`).join(" L")} L${values
      .slice()
      .reverse()
      .map((d, i) => `${px(values.length - 1 - i)},${py(d.max)}`)
      .join(" L")} Z`,
    line = values
      .map((d, i) => `${i ? "L" : "M"}${px(i)},${py(d.mean)}`)
      .join(" ");
  $(selector).innerHTML =
    `<svg viewBox="0 0 ${width} ${height}" role="img" aria-label="${yLabel} by ${xLabel}, randomized order mean and observed range"><path d="${area}" fill="${color}" opacity=".14"/><path d="${line}" fill="none" stroke="${color}" stroke-width="2"/><line x1="${l}" y1="${t + innerH}" x2="${width - r}" y2="${t + innerH}" class="chart-axis"/><text x="${l}" y="${height - 9}" class="axis-label">1 genome</text><text x="${width - r}" y="${height - 9}" text-anchor="end" class="axis-label">${values.length} genomes</text><text x="${l - 5}" y="${t + 4}" text-anchor="end" class="axis-label">${max.toFixed(0)}</text></svg>`;
}
function renderCompareOptions() {
  $("#compare-options").innerHTML = data.strains
    .map(
      (s) =>
        `<label class="compare-option"><input type="checkbox" value="${escapeHtml(s.id)}" ${selectedGenomeIds.includes(s.id) ? "checked" : ""}><span><strong>${escapeHtml(s.id)}</strong><small>${s.families} families</small></span></label>`,
    )
    .join("");
  $$(".compare-option input").forEach((input) =>
    input.addEventListener("change", () => {
      selectedGenomeIds = $$(".compare-option input:checked").map(
        (x) => x.value,
      );
      syncSelection();
      renderTree();
      renderStatistics();
    }),
  );
}
async function runCompare() {
  const ids = $$(".compare-option input:checked").map((x) => x.value);
  $("#compare-error").textContent = "";
  try {
    const response = await fetch("/api/compare", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ genome_ids: ids, dataset_id: datasetId }),
    });
    const result = await response.json();
    if (!response.ok) throw new Error(result.error || "Comparison failed");
    renderComparison(result);
  } catch (error) {
    $("#compare-error").textContent = error.message;
    $("#comparison-results").innerHTML = "";
  }
}
function renderComparison(result) {
  const unit = data.dataset.mode === "real" ? "genomes" : "demo genomes";
  const unique = Object.entries(result.unique_families)
    .map(
      ([id, families]) =>
        `<div class="compare-result-row"><strong>${escapeHtml(id)}</strong><span>${families.length}</span><small>families unique within selection</small></div>`,
    )
    .join("");
  const pairs = result.genomes
    .map(
      (a) =>
        `<tr><th>${escapeHtml(a)}</th>${result.genomes.map((b) => `<td>${(100 * result.pairwise_similarity[a][b]).toFixed(1)}%</td>`).join("")}</tr>`,
    )
    .join("");
  $("#comparison-results").innerHTML =
    `<div class="compare-summary"><article class="metric"><span class="section-kicker">SHARED BY ALL SELECTED</span><strong>${result.shared_family_count}</strong><small>gene families</small></article><article class="metric"><span class="section-kicker">FAMILY UNION</span><strong>${result.union_family_count}</strong><small>distinct families across selection</small></article><article class="metric"><span class="section-kicker">STRAINS</span><strong>${result.genomes.length}</strong><small>selected ${unit}</small></article></div><div class="comparison-grid"><article class="panel"><div class="panel-head"><div><span class="section-kicker">STRAIN-SPECIFIC WITHIN SELECTION</span><h3>Families unique to each</h3></div></div>${unique}</article><article class="panel"><div class="panel-head"><div><span class="section-kicker">PAIRWISE</span><h3>Jaccard similarity</h3></div></div><div class="table-scroll"><table><tbody>${pairs}</tbody></table></div></article></div><p class="comparison-note">${escapeHtml(result.metric)}</p>`;
}
function openView(name) {
  $$(".page-view").forEach((v) =>
    v.classList.toggle("shown", v.id === `view-${name}`),
  );
  $$(".nav-item").forEach((b) =>
    b.classList.toggle("active", b.dataset.view === name),
  );
  const title = {
    overview: "Overview",
    matrix: "Presence / absence",
    strains: "Strain explorer",
    genome: "Genome landscape",
    synteny: "Synteny explorer",
    families: "Gene families",
    tree: "Gene-content tree",
    statistics: "Statistics",
    compare: "Compare strains",
    methods: "Methods & provenance",
    sources: "Data sources",
  }[name];
  $("#crumb-current").textContent = title;
  window.scrollTo({ top: 0, behavior: "smooth" });
}
function updateDatasetContext() {
  const real = data.dataset.mode === "real";
  $("#dataset-warning-title").textContent = real
    ? "Real RefSeq / Panaroo analysis"
    : "Illustrative dataset";
  $("#dataset-warning-body").textContent = real
    ? `${data.strains.length} RefSeq assemblies processed by Panaroo. Results are source-based; sequence phylogeny, synteny and normalized locus coordinates remain unavailable.`
    : "All values below come from a deterministic synthetic fixture. They do not describe real Lactiplantibacillus plantarum isolates.";
  $("#matrix-warning-title").textContent = real
    ? "Panaroo-derived family matrix"
    : "Synthetic values";
  $("#matrix-warning-body").textContent = real
    ? "Family membership and per-genome copy counts were imported from Panaroo gene_presence_absence.csv."
    : "This matrix is generated from a fixed seed and is not biological evidence.";
  $("#view-matrix .lede").textContent =
    `A family-by-strain view of the ${real ? "Panaroo-derived RefSeq cohort" : "demo cohort"}. Cells encode copy number; select a family to inspect its distribution.`;
  $("#family-distribution-note").textContent = real
    ? `Classes are calculated from Panaroo family membership across ${data.strains.length} validated genomes.`
    : "Classification is calculated from the selected demo cohort.";
  $("#landscape-note").textContent = real
    ? "The source GFF3 files are preserved, but family membership has not been mapped to a normalized coordinate axis."
    : "Synthetic coordinates for interface demonstration.";
  $("#landscape-axis").hidden = real;
  $("#footer-data-note").textContent = real
    ? "RefSeq and Panaroo output · see methods and provenance"
    : "All demo values clearly labeled synthetic";
  $(".profile small").textContent = real
    ? "Local RefSeq analysis"
    : "Single-user demo";
  $(".section-heading h2 i").textContent = data.dataset.name;
  $("#view-tree .tree-panel .section-kicker").textContent = real
    ? `${data.strains.length} VALIDATED REFSEQ GENOMES`
    : "DEMO COHORT";
  $("#view-tree .tree-panel .panel-foot").textContent = real
    ? "Select leaves to filter the matrix and seed a strain comparison. Distances derive from binary family presence in the Panaroo matrix; this is not a sequence phylogeny."
    : "Select leaves to filter the matrix and seed a strain comparison. Tree labels and distances derive from the synthetic fixture.";
  $("#view-statistics .version").textContent = real
    ? "REFSEQ · PANAROO"
    : "DEMO CALCULATIONS";
  $("#view-compare .lede").textContent =
    `Compare shared and strain-specific family membership within this ${real ? "RefSeq cohort" : "demo cohort"}.`;
  $("#view-compare .compare-controls > div:first-child p").textContent =
    "Choose two to eight genomes. Tree selections are carried into this view.";
  $(".dataset-select").title = real
    ? "NCBI RefSeq data processed by Panaroo"
    : "Deterministic synthetic demonstration dataset";
  const query = `?dataset=${encodeURIComponent(datasetId)}`;
  $('a[href^="/api/export/families.tsv"]').href =
    `/api/export/families.tsv${query}`;
  $('a[href^="/api/export/tree.newick"]').href =
    `/api/export/tree.newick${query}`;
  $("#view-tree .lede").textContent = real
    ? "Hierarchical clustering of binary family-presence profiles from the selected RefSeq/Panaroo analysis. This is a gene-content dendrogram, not a sequence-based evolutionary phylogeny."
    : "Hierarchical clustering of binary family-presence profiles. This is a gene-content dendrogram, not a sequence-based evolutionary phylogeny.";
  const cards = $$(".method-card");
  if (real) {
    $(".methodology h3").textContent = "Interpretation boundary";
    $(".methodology p").textContent =
      "This convenience cohort is selected from public RefSeq assemblies and is not a curated biological study set. Panaroo family membership supports gene-content comparisons; this build does not infer sequence phylogeny, synteny, or rearrangements.";
  }
  if (real && cards.length >= 4) {
    cards[0].querySelector("h3").textContent = "NCBI RefSeq inputs and Panaroo";
    cards[0].querySelector("p").textContent =
      `${data.strains.length} validated RefSeq assemblies were processed with Panaroo ${data.analysis?.engine_version || "(version in manifest)"} (NCBI Datasets ${data.analysis?.datasets_cli_version || "version in manifest"}). Source accessions, validation exclusions, input checksums, parameters and logs are retained in the local analysis artifact directory.`;
    cards[0].querySelector(".status-pill").textContent =
      "Real assemblies processed";
    cards[1].querySelector("h3").textContent =
      "Family membership and prevalence";
    cards[1].querySelector("p").textContent =
      "Membership and copy counts come from Panaroo gene_presence_absence.csv. Prevalence is distinct genomes with the family divided by analyzed genomes. Presence-based Jaccard and UPGMA are derived from that matrix.";
    cards[1].querySelector(".status-pill").textContent =
      "Derived from recorded Panaroo output";
    cards[2].querySelector("h3").textContent = "Unavailable sequence analyses";
    cards[2].querySelector("p").textContent =
      "This analysis did not run a core-genome inference program, protein alignment, or whole-genome alignment. Sequence phylogeny, sequence identity and synteny are unavailable; the gene-content dendrogram is not a phylogeny.";
    cards[2].querySelector(".status-pill").textContent = "Not computed";
    cards[3].querySelector("h3").textContent = "Reproducibility artifacts";
    cards[3].querySelector("p").textContent =
      "The input selection, NCBI Datasets version, Panaroo version, selected parameters, file checksums, logs and normalized family matrix are stored locally.";
    $("#synteny-state-title").textContent =
      "Gene order has not been reconciled";
    $("#synteny-state-body").textContent =
      "Source GFF3 annotations are preserved, but this run does not currently reconcile Panaroo member IDs to contig-aware loci or compute aligned syntenic blocks. No rearrangement calls are available.";
    $("#ncbi-source-status").textContent =
      `Used · Datasets ${data.analysis?.datasets_cli_version || "recorded in manifest"}`;
    $("#ncbi-source-status").classList.remove("muted");
    $("#panaroo-source-status").textContent =
      `Used · ${data.analysis?.engine_version || "version in manifest"}`;
    $("#panaroo-source-status").classList.remove("muted");
  } else if (!real && cards.length >= 4) {
    $(".methodology h3").textContent = "Illustrations are not observations";
    $(".methodology p").textContent =
      "The demo is useful for exploring interaction patterns and the data model. Its synthetic membership patterns, coordinates and source labels must not be cited or interpreted as findings about L. plantarum.";
    cards[0].querySelector("h3").textContent =
      "Deterministic synthetic fixture";
    cards[0].querySelector("p").textContent =
      "Dataset ID demo-lp-v1, seed 20261122. Strain names, family membership, copy counts and coordinate offsets are generated for UI demonstration. The species label provides context only; records are not real isolates.";
    cards[0].querySelector(".status-pill").textContent =
      "No real genomic data loaded";
    cards[1].querySelector("h3").textContent = "Family prevalence";
    cards[1].querySelector("p").textContent =
      "Prevalence is the number of distinct genomes with a family divided by the number of genomes in the selected cohort. Core means present in all genomes; accessory uses a 15% frequency floor; singleton means one genome.";
    cards[2].querySelector("h3").textContent = "Real analysis pipeline";
    cards[2].querySelector("p").textContent =
      "NCBI acquisition, Panaroo clustering, sequence-based phylogenetic inference and sequence-based synteny are not available in this demo. No evolutionary tree or biological relationship is claimed.";
    cards[2].querySelector(".status-pill").textContent =
      "Unavailable for demo data";
    cards[3].querySelector("h3").textContent = "API and export";
    cards[3].querySelector("p").textContent =
      "The local API exposes dataset summary, matrix, family, similarity and statistics views. TSV/CSV exports preserve integer copy counts with zero for absent families; prevalence is a fraction in [0, 1].";
    $("#synteny-state-title").textContent =
      "This dataset has no gene-order evidence";
    $("#synteny-state-body").textContent =
      "The demo fixture contains family membership and illustrative markers, but no contigs, genes, strands, or validated alignments. It cannot support synteny blocks or rearrangement calls.";
    $("#panaroo-source-status").textContent = "Executable not installed";
    $("#panaroo-source-status").classList.add("muted");
  }
}
function bind() {
  $$(".nav-item").forEach((b) =>
    b.addEventListener("click", () => openView(b.dataset.view)),
  );
  $$("[data-open]").forEach((b) =>
    b.addEventListener("click", () => openView(b.dataset.open)),
  );
  $("#class-filter").addEventListener("change", renderMatrix);
  $("#family-search").addEventListener("input", renderMatrix);
  $("#family-class-filter").addEventListener("change", renderFamilies);
  $("#catalog-search").addEventListener("input", renderFamilies);
  $("#clear-tree-selection").addEventListener("click", () => {
    selectedGenomeIds = [];
    syncSelection();
    renderTree();
    renderCompareOptions();
  });
  $("#run-compare").addEventListener("click", runCompare);
  $("#dataset-select").addEventListener("change", (event) => {
    const url = new URL(location.href);
    url.searchParams.set("dataset", event.target.value);
    url.searchParams.delete("strains");
    location.href = url;
  });
  $("#theme-toggle").addEventListener("click", () => {
    document.documentElement.classList.toggle("dark");
    localStorage.setItem(
      "gca-dark",
      document.documentElement.classList.contains("dark") ? "1" : "0",
    );
  });
  if (localStorage.getItem("gca-dark") === "1")
    document.documentElement.classList.add("dark");
}
init();
