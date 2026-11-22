"""Small dependency-free API for the deterministic Gene Content Atlas demo."""
from __future__ import annotations

import csv
import io
import json
import math
import os
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from analysis import compare_genomes, pangenome_statistics, similarity_matrix, upgma

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "demo.json"
REAL_DATA = ROOT / "data" / "real" / "analysis" / "pangenome.json"
WEB = ROOT / "apps" / "web"


def demo_data() -> dict:
    """Generate a reproducible illustrative matrix and coordinate model."""
    strains = [
        {"id": f"LP-{i:02d}", "name": f"L. plantarum demo {i:02d}", "source": ["dairy", "plant", "fermentation"][i % 3]}
        for i in range(1, 13)
    ]
    families = []
    for i in range(1, 81):
        if i <= 38:
            members = list(range(12))
            kind = "core"
        elif i <= 52:
            members = [j for j in range(12) if (j + i) % 5 != 0]
            kind = "accessory"
        elif i <= 68:
            members = [j for j in range(12) if (j * 3 + i) % 4 == 0]
            kind = "rare"
        else:
            members = [(i * 7) % 12]
            kind = "singleton"
        # Add reproducible copy-number variation while retaining one family row per strain.
        copies = {str(j): (2 if ((i * 11 + j * 7) % 29 == 0) else 1) for j in members}
        families.append({
            "id": f"GF{i:04d}", "label": f"Family {i:03d}", "members": [strains[j]["id"] for j in members],
            "copies": {strains[int(j)]["id"]: n for j, n in copies.items()}, "group": kind,
            # Distribute illustrative markers across the full axis instead of
            # implying a localized genomic cluster through a poor demo layout.
            "position": ((i - 1) * 4_450_000 // 80) + ((i * 7919) % 20_000),
        })
    return {
        "dataset": {"id": "demo-lp-v1", "name": "Lactiplantibacillus plantarum", "label": "ILLUSTRATIVE DEMO · SYNTHETIC DATA", "mode": "demo", "seed": 20261122, "source": "Deterministic local fixture; no biological conclusions are represented."},
        "strains": strains,
        "families": families,
        "genome_length": 4_450_000,
        "methods": {"prevalence": "Distinct genomes containing a family / genomes in the selected dataset.", "distance": "1 − Jaccard similarity of family presence sets; copy number is not used.", "tree": "UPGMA clustering of the demo presence/absence Jaccard distances; illustrative only."},
    }


def load_data(dataset_id: str = "demo-lp-v1") -> dict:
    if dataset_id == "real-lp-refseq-panaroo":
        if not REAL_DATA.is_file():
            raise KeyError(dataset_id)
        return json.loads(REAL_DATA.read_text())
    if dataset_id != "demo-lp-v1":
        raise KeyError(dataset_id)
    if not DATA.exists():
        DATA.parent.mkdir(parents=True, exist_ok=True)
        DATA.write_text(json.dumps(demo_data(), indent=2) + "\n")
    return json.loads(DATA.read_text())


def derive(d: dict) -> dict:
    n = len(d["strains"])
    if n < 1:
        raise ValueError("Dataset has no genomes")
    counts = {k: 0 for k in ("core", "accessory", "rare", "singleton")}
    partitions = {k: 0 for k in ("strict_core", "soft_core", "shell", "cloud", "singleton")}
    soft_core_floor = math.ceil(.95 * n)
    shell_floor = math.ceil(.15 * n)
    for f in d["families"]:
        present = len(f["members"])
        f["prevalence"] = present / n
        if present == n:
            partition = "strict_core"
        elif present >= soft_core_floor:
            partition = "soft_core"
        elif present == 1:
            partition = "singleton"
        elif present >= shell_floor:
            partition = "shell"
        else:
            partition = "cloud"
        partitions[partition] += 1
        f["group"] = "core" if partition in ("strict_core", "soft_core") else "singleton" if partition == "singleton" else "rare" if partition == "cloud" else "accessory"
        counts[f["group"]] += 1
    genome_stats = []
    for strain in d["strains"]:
        present = [f for f in d["families"] if strain["id"] in f["members"]]
        genes = sum(f["copies"][strain["id"]] for f in present)
        genome_stats.append({**strain, "families": len(present), "genes": genes, "length": strain.get("length") or d.get("genome_length")})
    genome_length = d.get("genome_length") or max((s.get("length") or 0 for s in genome_stats), default=0)
    return {"dataset": d["dataset"], "strains": genome_stats, "families": d["families"], "counts": counts, "partitions": partitions, "genome_length": genome_length, "methods": d["methods"], "analysis": d.get("analysis")}


class Handler(BaseHTTPRequestHandler):
    def respond(self, status: int, body: bytes, content_type: str) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("X-Frame-Options", "DENY")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Content-Security-Policy", "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'self'; object-src 'none'; base-uri 'self'; form-action 'self'; frame-ancestors 'none'")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)
        path = parsed.path
        dataset_id = parse_qs(parsed.query).get("dataset", ["demo-lp-v1"])[0]
        if path in {"/api/overview", "/api/matrix", "/api/strains", "/api/families", "/api/similarity", "/api/tree", "/api/statistics", "/api/export/families.tsv", "/api/export/families.csv", "/api/export/tree.newick"}:
            try:
                load_data(dataset_id)
            except KeyError:
                self.respond(404, json.dumps({"error": f"Dataset unavailable: {dataset_id}"}).encode(), "application/json; charset=utf-8")
                return
        if path == "/api/health":
            self.respond(200, b'{"status":"ok","dataset":"demo-lp-v1"}', "application/json")
        elif path == "/api/datasets":
            datasets = [{"id": "demo-lp-v1", "name": "Lactiplantibacillus plantarum", "mode": "demo", "label": "Illustrative synthetic demonstration"}]
            if REAL_DATA.is_file():
                real = json.loads(REAL_DATA.read_text())["dataset"]
                datasets.append({"id": real["id"], "name": real["name"], "mode": "real", "label": real["label"]})
            self.respond(200, json.dumps({"datasets": datasets}).encode(), "application/json; charset=utf-8")
        elif path in ("/api/overview", "/api/matrix", "/api/strains", "/api/families"):
            try:
                payload = derive(load_data(dataset_id))
            except KeyError:
                self.respond(404, json.dumps({"error": f"Dataset unavailable: {dataset_id}"}).encode(), "application/json; charset=utf-8")
                return
            if path.endswith("/strains"):
                payload = {"dataset": payload["dataset"], "strains": payload["strains"]}
            elif path.endswith("/families"):
                payload = {"dataset": payload["dataset"], "families": payload["families"]}
            self.respond(200, json.dumps(payload).encode(), "application/json; charset=utf-8")
        elif path == "/api/similarity":
            try:
                payload = similarity_matrix(derive(load_data(dataset_id)))
            except KeyError:
                self.respond(404, b'{"error":"Dataset unavailable"}', "application/json; charset=utf-8")
                return
            self.respond(200, json.dumps(payload).encode(), "application/json; charset=utf-8")
        elif path == "/api/tree":
            try:
                payload = upgma(derive(load_data(dataset_id)))
            except KeyError:
                self.respond(404, b'{"error":"Dataset unavailable"}', "application/json; charset=utf-8")
                return
            self.respond(200, json.dumps(payload).encode(), "application/json; charset=utf-8")
        elif path == "/api/statistics":
            try:
                selected = parse_qs(parsed.query).get("genomes", [None])[0]
                genome_ids = selected.split(",") if selected else None
                payload = pangenome_statistics(derive(load_data(dataset_id)), genome_ids=genome_ids)
            except KeyError:
                self.respond(404, b'{"error":"Dataset unavailable"}', "application/json; charset=utf-8")
                return
            except ValueError as error:
                self.respond(400, json.dumps({"error": str(error)}).encode(), "application/json; charset=utf-8")
                return
            self.respond(200, json.dumps(payload).encode(), "application/json; charset=utf-8")
        elif path == "/api/export/families.tsv":
            d = derive(load_data(dataset_id))
            out = io.StringIO(newline="")
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            writer.writerow(["family_id", "class", "genomes_present", "prevalence", *[s["id"] for s in d["strains"]]])
            for f in d["families"]:
                writer.writerow([f["id"], f["group"], len(f["members"]), f"{f['prevalence']:.6f}", *[f["copies"].get(s["id"], 0) for s in d["strains"]]])
            self.respond(200, out.getvalue().encode(), "text/tab-separated-values; charset=utf-8")
        elif path == "/api/export/families.csv":
            d = derive(load_data(dataset_id))
            out = io.StringIO(newline="")
            writer = csv.writer(out)
            writer.writerow(["family_id", "class", "genomes_present", "prevalence", *[s["id"] for s in d["strains"]]])
            for f in d["families"]:
                writer.writerow([f["id"], f["group"], len(f["members"]), f"{f['prevalence']:.6f}", *[f["copies"].get(s["id"], 0) for s in d["strains"]]])
            self.respond(200, out.getvalue().encode(), "text/csv; charset=utf-8")
        elif path == "/api/export/tree.newick":
            self.respond(200, (upgma(derive(load_data(dataset_id)))["newick"] + "\n").encode(), "text/plain; charset=utf-8")
        else:
            target = WEB / ("index.html" if path == "/" else path.lstrip("/"))
            if not target.resolve().is_relative_to(WEB.resolve()) or not target.is_file():
                self.respond(404, b"Not found", "text/plain; charset=utf-8")
                return
            kind = "text/html" if target.suffix == ".html" else "text/css" if target.suffix == ".css" else "text/javascript"
            self.respond(200, target.read_bytes(), f"{kind}; charset=utf-8")

    def do_POST(self) -> None:  # noqa: N802
        path = urlparse(self.path).path
        if path != "/api/compare":
            self.respond(404, b'{"error":"Not found"}', "application/json; charset=utf-8")
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if length < 2 or length > 16_384:
                raise ValueError("Request must contain a small JSON genome selection")
            body = json.loads(self.rfile.read(length))
            genome_ids = body.get("genome_ids") if isinstance(body, dict) else None
            if not isinstance(genome_ids, list) or not all(isinstance(item, str) for item in genome_ids):
                raise ValueError("genome_ids must be a list of genome ID strings")
            dataset_id = body.get("dataset_id", "demo-lp-v1")
            try:
                dataset = derive(load_data(dataset_id))
            except KeyError:
                self.respond(404, json.dumps({"error": f"Dataset unavailable: {dataset_id}"}).encode(), "application/json; charset=utf-8")
                return
            result = compare_genomes(dataset, genome_ids)
            self.respond(200, json.dumps(result).encode(), "application/json; charset=utf-8")
        except (ValueError, json.JSONDecodeError) as error:
            response = json.dumps({"error": str(error)}).encode()
            self.respond(400, response, "application/json; charset=utf-8")

    def log_message(self, fmt: str, *args) -> None:
        print(f"[{self.log_date_time_string()}] {fmt % args}")


if __name__ == "__main__":
    load_data()
    host = os.environ.get("HOST", "127.0.0.1")
    port = int(os.environ.get("PORT", "8000"))
    print(f"Gene Content Atlas local API: http://{host}:{port}")
    ThreadingHTTPServer((host, port), Handler).serve_forever()
