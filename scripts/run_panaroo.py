"""Validate acquired RefSeq inputs, run Panaroo, and normalize its results."""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REAL = ROOT / "data" / "real"
sys.path.insert(0, str(ROOT / "services" / "api"))
from ingest import genome_metadata, parse_fasta, parse_gff3, parse_panaroo_csv  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def find_inputs(package_dir: Path, accessions: list[str], input_dir: Path) -> tuple[list[tuple[str, Path, Path]], list[dict], list[dict[str, str]]]:
    source_root = package_dir / "ncbi_dataset" / "data"
    inputs = []
    genome_meta = []
    exclusions = []
    input_dir.mkdir(parents=True, exist_ok=True)
    for accession in accessions:
        folder = source_root / accession
        if not folder.is_dir():
            raise RuntimeError(f"Missing downloaded assembly directory: {folder}")
        fastas = [p for p in folder.rglob("*_genomic.fna") if "cds_from_genomic" not in p.name]
        gffs = list(folder.rglob("*.gff")) + list(folder.rglob("*.gff3"))
        if len(fastas) != 1 or len(gffs) != 1:
            raise RuntimeError(f"Expected exactly one genomic FASTA and one GFF3 for {accession}; found {len(fastas)} FASTA and {len(gffs)} GFF files")
        original_fasta, original_gff = fastas[0], gffs[0]
        # Stable, accession-based filenames are also the sample IDs Panaroo writes.
        fasta = input_dir / f"{accession}.fna"
        gff = input_dir / f"{accession}.gff"
        for target, source in ((fasta, original_fasta), (gff, original_gff)):
            if target.exists() or target.is_symlink():
                target.unlink()
            target.symlink_to(source.resolve())
        try:
            with fasta.open() as handle:
                sequences = parse_fasta(handle)
            with gff.open() as handle:
                genes = parse_gff3(handle, {key: len(value) for key, value in sequences.items()})
        except ValueError as error:
            exclusions.append({"accession": accession, "reason": f"FASTA/GFF3 validation failed: {error}"})
            fasta.unlink(missing_ok=True)
            gff.unlink(missing_ok=True)
            continue
        cds_count = len({(gene["contig"], gene["id"]) for gene in genes})
        if not cds_count:
            exclusions.append({"accession": accession, "reason": "No CDS features parsed from the NCBI GFF3 annotation."})
            fasta.unlink(missing_ok=True)
            gff.unlink(missing_ok=True)
            continue
        metadata = genome_metadata(original_fasta, sequences)
        metadata.update({"accession": accession, "fasta_path": str(original_fasta.relative_to(ROOT)), "gff_path": str(original_gff.relative_to(ROOT)), "fasta_sha256": sha256(original_fasta), "gff_sha256": sha256(original_gff), "cds_count": cds_count})
        inputs.append((accession, gff, fasta))
        genome_meta.append(metadata)
    return inputs, genome_meta, exclusions


def main() -> int:
    manifest_path = REAL / "selection_manifest.json"
    package_dir = REAL / "package"
    if not manifest_path.is_file() or not package_dir.is_dir():
        raise RuntimeError("Run `make ingest-real` first to acquire the RefSeq package.")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("dataset_id") != "real-lp-refseq" or not manifest.get("selected"):
        raise RuntimeError("Selection manifest is missing or does not describe the expected RefSeq dataset.")
    # The first invocation in a fresh conda environment can spend substantial
    # time importing compiled scientific dependencies before printing a version.
    version_result = subprocess.run(["panaroo", "--version"], text=True, capture_output=True, timeout=120)
    if version_result.returncode:
        raise RuntimeError(f"Panaroo version check failed: {version_result.stderr.strip()}")
    panaroo_version = (version_result.stdout or version_result.stderr).strip()
    run_dir = REAL / "analysis"
    input_dir = run_dir / "inputs"
    outputs = run_dir / "panaroo"
    if outputs.exists():
        shutil.rmtree(outputs)
    inputs, genome_metadata_records, exclusions = find_inputs(package_dir, [row["accession"] for row in manifest["selected"]], input_dir)
    if len(inputs) < 2:
        raise RuntimeError(f"Only {len(inputs)} assemblies passed FASTA/GFF3 validation; at least two are required. Exclusions: {exclusions}")
    manifest["analysis_included_accessions"] = [accession for accession, _, _ in inputs]
    manifest["analysis_exclusions"] = exclusions
    manifest["analysis_included_count"] = len(inputs)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    input_list = run_dir / "panaroo_inputs.txt"
    input_list.parent.mkdir(parents=True, exist_ok=True)
    input_list.write_text("\n".join(f"{gff}\t{fasta}" for _, gff, fasta in inputs) + "\n")
    thread_count = max(1, min(4, os.cpu_count() or 1))
    command = ["panaroo", "-i", str(input_list), "-o", str(outputs), "--clean-mode", "strict", "--remove-invalid-genes", "-t", str(thread_count)]
    started = datetime.now(timezone.utc).isoformat()
    result = subprocess.run(command, cwd=run_dir, text=True, capture_output=True, timeout=int(os.environ.get("PANAROO_TIMEOUT_SECONDS", "86400")), check=False)
    (run_dir / "panaroo.stdout.log").write_text(result.stdout)
    (run_dir / "panaroo.stderr.log").write_text(result.stderr)
    if result.returncode:
        raise RuntimeError(f"Panaroo failed with exit code {result.returncode}; see {run_dir / 'panaroo.stderr.log'}")
    matrix_path = outputs / "gene_presence_absence.csv"
    if not matrix_path.is_file():
        raise RuntimeError("Panaroo exited successfully but gene_presence_absence.csv was not produced")
    aliases = {accession: accession for accession, _, _ in inputs}
    with matrix_path.open(newline="") as handle:
        parsed = parse_panaroo_csv(handle, aliases)
    genome_by_id = {row["accession"]: row for row in manifest["selected"]}
    strains = []
    for meta in genome_metadata_records:
        source = genome_by_id[meta["accession"]]
        strains.append({"id": meta["accession"], "name": source.get("strain") or source.get("organism") or meta["accession"], "source": "NCBI RefSeq", "length": meta["genome_length"], "gc_fraction": meta["gc_fraction"], "contig_count": meta["contig_count"], "cds_count": meta["cds_count"], "sequence_sha256": meta["sequence_sha256"], "assembly_accession": meta["accession"], "assembly_level": source.get("assembly_level"), "source_url": source["source_url"]})
    families = parsed["families"]
    for family in families:
        prevalence = len(family["members"]) / len(strains)
        family["prevalence"] = prevalence
        family["group"] = "core" if prevalence == 1 else "accessory" if prevalence >= .15 else "singleton" if len(family["members"]) == 1 else "rare"
    artifact = {
        "dataset": {"id": "real-lp-refseq-panaroo", "name": "Lactiplantibacillus plantarum", "label": "REAL DATA · NCBI REFSEQ + PANAROO", "source": "NCBI RefSeq; Panaroo", "mode": "real"},
        "strains": strains,
        "families": families,
        "methods": {"family_assignment": "Panaroo gene_presence_absence.csv", "prevalence": "Distinct genomes containing a family / genomes in the analysis.", "classification": "Frequency thresholds: core=all, accessory=at least 15% but below all, rare=2+ below 15%, singleton=1 genome.", "coordinates": "Source GFF3 retained in the downloaded data package; not represented as a single pangenome coordinate axis."},
        "analysis": {"started_at": started, "completed_at": datetime.now(timezone.utc).isoformat(), "datasets_cli_version": manifest.get("datasets_cli_version"), "engine": "Panaroo", "engine_version": panaroo_version, "command": command, "thread_count": thread_count, "selection_manifest_sha256": sha256(manifest_path), "input_manifest_sha256": sha256(input_list), "matrix_sha256": sha256(matrix_path), "excluded_assemblies": exclusions, "warnings": parsed["warnings"], "artifacts": {"matrix": str(matrix_path.relative_to(ROOT)), "graph": str((outputs / "final_graph.gml").relative_to(ROOT)) if (outputs / "final_graph.gml").exists() else None, "core_alignment": str((outputs / "core_gene_alignment.aln").relative_to(ROOT)) if (outputs / "core_gene_alignment.aln").exists() else None}},
    }
    artifact_path = run_dir / "pangenome.json"
    artifact_path.write_text(json.dumps(artifact, indent=2) + "\n")
    print(f"Imported {len(families)} families from Panaroo for {len(strains)} RefSeq genomes.")
    if exclusions:
        print(f"Excluded {len(exclusions)} assemblies after FASTA/GFF3 validation; see the selection manifest.")
    print(f"Normalized result: {artifact_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, ValueError, OSError, subprocess.TimeoutExpired) as error:
        print(f"Panaroo analysis failed: {error}", file=sys.stderr)
        raise SystemExit(1)
