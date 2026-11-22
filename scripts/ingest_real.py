"""Discover a RefSeq cohort and acquire its official NCBI Datasets package."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

ROOT = Path(__file__).resolve().parents[1]
REAL = ROOT / "data" / "real"
SPECIES = "Lactiplantibacillus plantarum"


def walk(value: Any) -> Iterator[tuple[str, Any]]:
    if isinstance(value, dict):
        for key, item in value.items():
            yield key.casefold().replace("-", "_").replace(" ", "_"), item
            yield from walk(item)
    elif isinstance(value, list):
        for item in value:
            yield from walk(item)


def find_value(value: Any, *keys: str) -> Any:
    targets = {key.casefold().replace("-", "_").replace(" ", "_") for key in keys}
    return next((item for key, item in walk(value) if key in targets and item not in (None, "")), None)


def report_accession(report: dict[str, Any]) -> str | None:
    direct = report.get("accession")
    if isinstance(direct, str) and re.fullmatch(r"GCF_\d+\.\d+", direct):
        return direct
    return next((value for key, value in walk(report) if key in {"accession", "assembly_accession"} and isinstance(value, str) and re.fullmatch(r"GCF_\d+\.\d+", value)), None)


def parse_json_lines(text: str) -> list[dict[str, Any]]:
    entries = []
    for line_number, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        try:
            value = json.loads(line)
        except json.JSONDecodeError as error:
            raise ValueError(f"NCBI Datasets returned invalid JSON on line {line_number}") from error
        reports = value.get("reports") if isinstance(value, dict) else None
        if isinstance(reports, list):
            entries.extend(row for row in reports if isinstance(row, dict))
        elif isinstance(value, dict):
            entries.append(value)
    return entries


def checked_run(args: list[str], timeout: int) -> subprocess.CompletedProcess[str]:
    try:
        result = subprocess.run(args, text=True, capture_output=True, timeout=timeout, check=False)
    except FileNotFoundError as error:
        raise RuntimeError(f"Required executable not found: {args[0]}. Install the official NCBI Datasets CLI first.") from error
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(f"Command timed out after {timeout}s: {args[0]}") from error
    if result.returncode:
        detail = (result.stderr or result.stdout).strip()[-4000:]
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(args)}\n{detail}")
    return result


def normalized_record(report: dict[str, Any]) -> dict[str, Any] | None:
    accession = report_accession(report)
    if not accession:
        return None
    organism = find_value(report, "organism_name", "scientific_name") or SPECIES
    assembly_level = find_value(report, "assembly_level") or "unknown"
    strain = find_value(report, "strain", "isolate") or ""
    seq_len = find_value(report, "total_sequence_length", "sequence_length")
    n_contigs = find_value(report, "number_of_contigs", "contig_count")
    n50 = find_value(report, "contig_n50", "scaffold_n50")
    seq_len = int(seq_len) if str(seq_len).isdigit() else None
    n_contigs = int(n_contigs) if str(n_contigs).isdigit() else None
    n50 = int(n50) if str(n50).isdigit() else None
    level = str(assembly_level).casefold()
    return {
        "accession": accession,
        "assembly_version": accession.rsplit(".", 1)[-1],
        "organism": str(organism),
        "strain": str(strain),
        "assembly_level": str(assembly_level),
        "genome_length": seq_len,
        "contig_count": n_contigs,
        "contig_n50": n50,
        "source": "NCBI RefSeq",
        "source_url": f"https://www.ncbi.nlm.nih.gov/datasets/genome/{accession}/",
        "original_metadata": report,
        "_level_rank": {"complete genome": 0, "chromosome": 1, "scaffold": 2, "contig": 3}.get(level, 4),
    }


def select_records(records: list[dict[str, Any]], count: int) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    normalized = [record for raw in records if (record := normalized_record(raw))]
    normalized.sort(key=lambda r: (r["_level_rank"], r["contig_count"] if r["contig_count"] is not None else 10**9, -(r["contig_n50"] or 0), r["accession"]))
    selected = []
    exclusions = []
    seen_strains = set()
    for record in normalized:
        strain = record["strain"].strip().casefold()
        isolate_source = str(find_value(record["original_metadata"], "isolation_source") or "").casefold()
        risk_fields = " ".join(str(v).casefold() for k, v in walk(record["original_metadata"]) if "host" in k or "isolation_source" in k or "sample" in k)
        if any(token in risk_fields for token in ("clinical", "hospital", "patient", "blood culture")):
            exclusions.append({"accession": record["accession"], "reason": "Metadata indicates a clinical or patient-associated source; excluded for manual review."})
            continue
        dedupe_key = strain or record["accession"]
        if dedupe_key in seen_strains:
            exclusions.append({"accession": record["accession"], "reason": f"Duplicate strain label: {strain}"})
            continue
        if not record["genome_length"] or not record["contig_count"]:
            exclusions.append({"accession": record["accession"], "reason": "Assembly size or contig count missing from metadata."})
            continue
        seen_strains.add(dedupe_key)
        record["selection_rank"] = len(selected) + 1
        record["selection_reason"] = "Ranked by RefSeq assembly level, lower contig count, higher N50, then accession."
        record.pop("_level_rank", None)
        selected.append(record)
        if len(selected) >= count:
            break
    if not selected:
        raise RuntimeError("NCBI returned no usable RefSeq assemblies for the configured species and filters.")
    return selected, exclusions


def safe_extract(package: Path, destination: Path, max_uncompressed: int = 30_000_000_000) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(package) as archive:
        total = sum(info.file_size for info in archive.infolist())
        if total > max_uncompressed:
            raise RuntimeError(f"NCBI package expands to {total:,} bytes; limit is {max_uncompressed:,} bytes")
        for info in archive.infolist():
            target = (destination / info.filename).resolve()
            if not target.is_relative_to(destination.resolve()):
                raise RuntimeError(f"Unsafe path in NCBI package: {info.filename}")
            if info.is_dir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(info) as source, target.open("wb") as output:
                shutil.copyfileobj(source, output)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count", type=int, default=24, help="Target number of unique strains (1–32)")
    parser.add_argument("--discover-only", action="store_true", help="Write the selection manifest without downloading sequence files")
    parser.add_argument("--timeout", type=int, default=3600, help="Maximum seconds for each NCBI Datasets command")
    args = parser.parse_args()
    if not 1 <= args.count <= 32:
        parser.error("--count must be between 1 and 32")
    version = checked_run(["datasets", "--version"], timeout=30).stdout.strip()
    summary = checked_run([
        "datasets", "summary", "genome", "taxon", SPECIES,
        "--assembly-source", "RefSeq", "--annotated",
        "--assembly-level", "complete,chromosome,scaffold",
        "--exclude-atypical", "--exclude-multi-isolate",
        "--as-json-lines", "--limit", "all",
    ], timeout=args.timeout)
    raw_reports = parse_json_lines(summary.stdout)
    selected, exclusions = select_records(raw_reports, args.count)
    manifest = {
        "dataset_id": "real-lp-refseq",
        "species": SPECIES,
        "mode": "real",
        "source": "NCBI RefSeq",
        "source_release": None,
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
        "datasets_cli_version": version,
        "selection_criteria": ["RefSeq", "annotated", "complete/chromosome/scaffold", "exclude atypical", "exclude multi-isolate", "unique strain label", "assembly-level then contig-count/N50 ranking"],
        "selected_count": len(selected),
        "candidate_count": len(raw_reports),
        "selected": selected,
        "exclusions": exclusions,
        "selection_review": "Review the accession list and metadata before treating this convenience cohort as a curated biological study set.",
    }
    REAL.mkdir(parents=True, exist_ok=True)
    manifest_path = REAL / "selection_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Selected {len(selected)} unique RefSeq assemblies from {len(raw_reports)} candidate reports.")
    print(f"Selection manifest: {manifest_path}")
    if args.discover_only:
        return 0
    accessions = REAL / "selected_accessions.txt"
    accessions.write_text("\n".join(item["accession"] for item in selected) + "\n")
    package = REAL / "ncbi_dataset.zip"
    checked_run([
        "datasets", "download", "genome", "accession", "--inputfile", str(accessions),
        "--include", "genome,gff3,protein,cds,seq-report", "--filename", str(package),
        "--no-progressbar",
    ], timeout=args.timeout)
    digest = file_sha256(package)
    extract_dir = REAL / "package"
    safe_extract(package, extract_dir)
    manifest["package"] = {"path": str(package.relative_to(ROOT)), "sha256": digest, "extracted_to": str(extract_dir.relative_to(ROOT))}
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Downloaded and checksum-recorded NCBI package ({package.stat().st_size:,} bytes).")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, ValueError, zipfile.BadZipFile) as error:
        print(f"Real-data acquisition failed: {error}", file=sys.stderr)
        raise SystemExit(1)
