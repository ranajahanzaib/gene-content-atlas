"""Small, dependency-free parsers used to validate downloaded bacterial data."""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from urllib.parse import unquote
from typing import Any, TextIO

PANAROO_METADATA = {
    "Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
    "Avg sequences per isolate", "Genome Fragment", "Order within Fragment",
    "Accessory Fragment", "Accessory Order with Fragment", "QC",
}


def parse_fasta(handle: TextIO) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    identifier: str | None = None
    for line_number, raw in enumerate(handle, 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            identifier = line[1:].split(None, 1)[0]
            if not identifier:
                raise ValueError(f"Empty FASTA identifier at line {line_number}")
            if identifier in records:
                raise ValueError(f"Duplicate FASTA identifier: {identifier}")
            records[identifier] = []
        elif identifier is None:
            raise ValueError(f"Sequence encountered before a FASTA header at line {line_number}")
        else:
            sequence = "".join(line.split()).upper()
            if any(base not in "ACGTRYSWKMBDHVN.-" for base in sequence):
                raise ValueError(f"Invalid FASTA character at line {line_number}")
            records[identifier].append(sequence)
    if not records:
        raise ValueError("FASTA file contains no records")
    return {name: "".join(parts) for name, parts in records.items()}


def gff_attributes(raw: str) -> dict[str, str]:
    result = {}
    for item in raw.split(";"):
        if not item:
            continue
        key, separator, value = item.partition("=")
        if not separator:
            continue
        result[unquote(key)] = unquote(value)
    return result


def parse_gff3(handle: TextIO, contig_lengths: dict[str, int]) -> list[dict[str, Any]]:
    genes: list[dict[str, Any]] = []
    feature_parts: dict[tuple[str, str], tuple[str, set[tuple[int, int]]]] = {}
    for line_number, raw in enumerate(handle, 1):
        line = raw.rstrip("\r\n")
        if line == "##FASTA":
            break
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            raise ValueError(f"GFF3 line {line_number} must contain 9 tab-separated fields")
        seqid, source, feature_type, start_text, end_text, score, strand, phase, attributes_text = fields
        if feature_type != "CDS":
            continue
        if seqid not in contig_lengths:
            raise ValueError(f"GFF3 line {line_number} references unknown contig {seqid}")
        try:
            start_1, end_1 = int(start_text), int(end_text)
        except ValueError as error:
            raise ValueError(f"Invalid GFF3 coordinates on line {line_number}") from error
        if start_1 < 1 or end_1 < start_1 or end_1 > contig_lengths[seqid]:
            raise ValueError(f"GFF3 coordinates out of contig bounds on line {line_number}")
        if strand not in ("+", "-", ".", "?"):
            raise ValueError(f"Invalid GFF3 strand on line {line_number}")
        attributes = gff_attributes(attributes_text)
        feature_id = attributes.get("ID") or attributes.get("locus_tag")
        if not feature_id:
            raise ValueError(f"CDS feature has no ID or locus_tag on line {line_number}")
        unique = (seqid, feature_id)
        parent = attributes.get("Parent", "")
        interval = (start_1, end_1)
        prior = feature_parts.get(unique)
        if prior:
            prior_parent, prior_intervals = prior
            # GFF3 allows a discontinuous CDS to use the same ID on multiple
            # rows (for example, a programmed frameshift). Keep every segment
            # while rejecting duplicate intervals or inconsistent parents.
            if prior_parent != parent or interval in prior_intervals:
                raise ValueError(f"Inconsistent repeated CDS identifier {feature_id} on {seqid}")
            prior_intervals.add(interval)
        else:
            feature_parts[unique] = (parent, {interval})
        # GFF3 is 1-based inclusive; internal coordinates are 0-based half-open.
        genes.append({"id": feature_id, "contig": seqid, "start": start_1 - 1, "end": end_1, "strand": strand, "phase": phase, "source": source, "attributes": attributes})
    return genes


def genome_metadata(fasta_path: Path, sequences: dict[str, str]) -> dict[str, Any]:
    digest = hashlib.sha256()
    with fasta_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    sequence = "".join(sequences.values())
    gc = sum(base in "GC" for base in sequence)
    canonical = sum(base in "ACGT" for base in sequence)
    return {
        "sequence_sha256": digest.hexdigest(),
        "genome_length": len(sequence),
        "contig_count": len(sequences),
        "gc_fraction": gc / canonical if canonical else None,
        "contigs": [{"id": key, "length": len(value)} for key, value in sequences.items()],
    }


def parse_panaroo_csv(handle: TextIO, genome_aliases: dict[str, str] | None = None) -> dict[str, Any]:
    """Read Panaroo's Roary-style family table and retain copy counts.

    Membership cells contain semicolon-separated sequence identifiers. A blank
    cell is absent; it is never treated as unknown.
    """
    reader = csv.DictReader(handle)
    if not reader.fieldnames:
        raise ValueError("Panaroo family CSV has no header")
    fields = [name.strip() for name in reader.fieldnames]
    genome_columns = [name for name in fields if name not in PANAROO_METADATA]
    if not genome_columns:
        raise ValueError("No genome columns found in Panaroo gene_presence_absence.csv")
    aliases = genome_aliases or {}
    families = []
    warnings = []
    for row_number, raw in enumerate(reader, 2):
        row = {(key or "").strip(): (value or "").strip() for key, value in raw.items()}
        family_id = row.get("Gene")
        if not family_id:
            raise ValueError(f"Missing family identifier on CSV row {row_number}")
        copies = {}
        members = []
        genes = {}
        for genome_column in genome_columns:
            cell = row.get(genome_column, "").strip()
            # Roary-compatible Panaroo output separates paralog identifiers with semicolons.
            ids = [value.strip() for value in cell.split(";") if value.strip() and value.strip() != "-"]
            normalized = aliases.get(genome_column, genome_column)
            if ids:
                if normalized in copies:
                    raise ValueError(f"Genome alias collision for {normalized}")
                members.append(normalized)
                copies[normalized] = len(ids)
                genes[normalized] = ids
        # Panaroo 1.8.1 omits Roary's summary metadata columns by default.
        declared = row.get("No. isolates") or ""
        if declared.isdigit() and int(declared) != len(members):
            warnings.append(f"{family_id}: declared isolate count {declared} differs from parsed membership {len(members)}")
        families.append({"id": family_id, "annotation": row.get("Annotation", ""), "members": members, "copies": copies, "gene_ids": genes})
    if not families:
        raise ValueError("Panaroo family CSV has no data rows")
    return {"genome_columns": [aliases.get(name, name) for name in genome_columns], "families": families, "warnings": warnings}
