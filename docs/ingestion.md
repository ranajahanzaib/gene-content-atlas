# RefSeq ingestion and Panaroo

## Acquisition

`make ingest-real` calls the supported `datasets` CLI (version and syntax are checked at runtime) for annotated RefSeq assemblies of _Lactiplantibacillus plantarum_ at complete, chromosome, or scaffold level. A candidate manifest records the selection logic and exclusions before downloading genome, GFF3, protein, CDS, and sequence-report files. The current ranked selection is a starting cohort, not a hand-curated biological study; review it before analysis.

NCBI publishes a genome package ZIP and per-assembly reports. The package is checksum recorded and extracted with path traversal and expanded-size checks. Original NCBI FASTA and GFF3 files are kept under `data/real/package`.

## Validation and analysis

`make analyze` uses the installed Panaroo executable. For each selected accession it expects exactly one genomic nucleotide FASTA and one GFF/GFF3. FASTA IDs are unique; GFF3 CDS coordinates must resolve to the same contig identifiers and remain within sequence bounds. GFF3 coordinates normalize to zero-based half-open intervals for internal use. Invalid assemblies are excluded with explicit reasons in the selection manifest; the source annotations are never rewritten. A stable tab-separated list of valid GFF3/FASTA input pairs is recorded and passed with argument arrays, never through shell interpolation.

Panaroo runs with strict cleaning, invalid-gene filtering, and a bounded thread count. Its version, command, stdout/stderr, input checksums, exit status, result matrix checksum, and family warnings are retained. A non-zero exit or missing expected matrix stops the import. The parser handles the Roary-compatible metadata columns and semicolon-separated gene IDs in sample cells. Per-genome copy counts are the number of IDs in that cell; blank cells are absent.

The artifact powers the family matrix and gene-content analyses. Panaroo member IDs are not reconciled to contig-aware GFF loci, so the application does not draw synteny blocks. No core-alignment inference is run.

## Reproducibility and limits

The selection and analysis artifacts are local files under `data/real` and should be backed up or checksummed separately. The workflow has been run locally with NCBI Datasets CLI 18.41.1 and Panaroo 1.8.1. This build has no database transactions, background worker, cancellation flow, upload handling, GTDB mapping, or standardized reannotation. NCBI and Panaroo command syntax can evolve; the commands in the scripts fail visibly and version strings are captured for each run.
