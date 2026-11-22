# Data model and coordinates

## Demo and real datasets

Every family identifier is scoped to its dataset and analysis. `members` is the set of distinct genome IDs containing the family. `copies` maps each member genome to the number of distinct gene IDs Panaroo reported for that row. Absence is represented by a missing mapping and exported as integer zero; it is not the same as unknown.

The demo dataset has 12 synthetic genome labels and 80 reproducibly generated family rows. It includes no genes, proteins, contigs, or real accessions. The real artifact stores accession-version genome IDs, the Panaroo family table, sequence metadata, source URLs, and an analysis provenance object.

## Coordinates

GFF3 input uses one-based inclusive start/end coordinates. The parser stores `start = start_1 - 1` and `end = end_1`, yielding zero-based half-open coordinates. Contig identifiers and strand are retained. CDS features referring to unknown contigs or positions outside the corresponding FASTA record are rejected. Contig-relative coordinates must not be compared as though they share a cross-strain coordinate system.

## Provenance

The acquisition selection manifest stores NCBI Datasets CLI version, retrieval time, accession versions, source URLs, original assembly reports, selected and excluded rows, and the downloaded package checksum. The Panaroo result stores executable version, exact command, input manifest checksum, output matrix checksum, logs, per-input FASTA/GFF checksums, and any matrix consistency warnings.
