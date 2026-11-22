# Troubleshooting

## Demo page does not open

Run `make doctor`, then `make seed-demo` and `make dev`. If port 8000 is occupied, use `PORT=8765 make dev`. Stop the server with Ctrl+C.

## `make ingest-real` says `datasets` is missing

Install the official NCBI Datasets CLI for your OS from the official installation guide, add it to `PATH`, and check the resolved version with `make doctor`. Acquisition also needs internet access to NCBI.

## `make analyze` says `panaroo` is missing

Install Panaroo and its supported dependencies using a maintained bioinformatics distribution. Check `panaroo --version` and `make doctor` before running. The demo remains available independently.

## Panaroo input validation fails

The analysis requires one genome FASTA and one annotation GFF per selected assembly with matching sequence identifiers, valid CDS coordinates, and parseable features. Assemblies that fail this validation are excluded, with accession and reason recorded under `analysis_exclusions` in `data/real/selection_manifest.json` and in the normalized artifact. Inspect the named accession in that manifest and the local input package before interpreting the cohort.

## Real analysis does not show in the selector

The normalized `data/real/analysis/pangenome.json` artifact is written only after successful Panaroo matrix validation. Check exit status and `panaroo.stderr.log`; there is no automatic fallback from an incomplete real dataset to the demo.
