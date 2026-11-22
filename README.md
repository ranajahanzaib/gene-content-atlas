# Gene Content Atlas

**Explore the genomic diversity within a species.**

Gene Content Atlas is a local comparative genomics workspace. Its deterministic demo is ready to explore immediately. Optional real-data commands discover annotated RefSeq assemblies and run Panaroo when the required external tools and network access are available. The demo is always labeled synthetic and is never substituted for a missing real dataset.

## Quick start

Requirements: Python 3.10 or newer. The demo API and web client use only the Python standard library and browser JavaScript.

```sh
make setup
make seed-demo
make dev
```

Open <http://127.0.0.1:8000>. If that port is occupied, start with `PORT=8765 make dev`. Stop the server with Ctrl+C. Optional environment variables are listed in `.env.example`; export them in the shell before running Make targets.

The workspace includes a family-class overview, searchable copy-count matrix, strain and family tables, synthetic coordinate illustrations, gene-content UPGMA clustering, randomized pangenome statistics, pairwise strain comparisons, CSV/TSV/Newick exports, methods, and source pages. The circular coordinates in the demo are illustrative and do not describe real chromosomes.

## Real RefSeq workflow

Install the official NCBI Datasets CLI and Panaroo using their supported distribution channels. Check resolved executable versions with `make doctor`; the acquisition manifest and analysis artifact record the versions used. NCBI Datasets currently publishes macOS universal CLI binaries and documents both CLI syntax and downloadable genome package contents. Panaroo accepts GFF-based input lists and writes a Roary-style `gene_presence_absence.csv` matrix plus other optional artifacts. See the official [NCBI install guide](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/), [NCBI genome download reference](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/), [Panaroo quickstart](https://github.com/gtonkinhill/panaroo/blob/master/docs/gettingstarted/quickstart.md), and [Panaroo output guide](https://github.com/gtonkinhill/panaroo/blob/master/docs/gettingstarted/output.md).

Then run:

```sh
make ingest-real    # Discover, rank, manifest and download RefSeq assemblies
make analyze        # Validate local FASTA/GFF3, run Panaroo, normalize its matrix
```

`make ingest-real` selects up to 24 unique strain labels from annotated RefSeq complete/chromosome/scaffold candidates, excludes atypical and multi-isolate assemblies, flags clearly clinical/patient-associated metadata for exclusion, and ranks by assembly level, contig count and N50. Review `data/real/selection_manifest.json`: this automated convenience cohort is not a curated biological study set. The original accession versions, metadata, retrieval timestamp, CLI version, package checksum and selection/exclusion reasons are retained. Downloads are stored under `data/real/package/`.

`make analyze` requires the downloaded RefSeq package plus the `panaroo` executable. It validates FASTA records and CDS coordinate bounds, records and excludes any selected assembly with invalid annotations, writes a stable per-accession GFF/FASTA input list, runs strict-mode Panaroo with invalid-gene filtering, captures logs and tool versions, verifies expected output files, parses family membership and multicopy cells, and writes `data/real/analysis/pangenome.json`. Panaroo's graph and source GFFs remain in the local analysis directory. The normalized result appears in the dataset selector only after the artifact exists.

The UI distinguishes the Jaccard/UPGMA gene-content tree from sequence phylogeny. It does not claim sequence identity, syntenic blocks, rearrangements, or a circular family coordinate axis without supporting computed data. GTDB classification and sequence-based phylogeny are not provided in this build.

## Local commands

```text
make setup
make dev
make doctor
make seed-demo
make ingest-real
make analyze
make test
make lint
make build
```

`make doctor` reports installed runtime and bioinformatics tool versions. Missing external executables are reported as unavailable; no fake analysis is run in their place. A failed acquisition or analysis keeps diagnostics and exits non-zero.

## API

The local API serves the web app and exposes:

| Route                                     | Purpose                                                       |
| ----------------------------------------- | ------------------------------------------------------------- |
| `/api/health`                             | Local service status                                          |
| `/api/datasets`                           | Demo dataset and any available real Panaroo artifact          |
| `/api/overview?dataset=...`               | Dataset summary, family classes, strain summaries             |
| `/api/matrix?dataset=...`                 | Family membership and exact per-genome copy counts            |
| `/api/strains?dataset=...`                | Genome catalogue                                              |
| `/api/families?dataset=...`               | Family catalogue                                              |
| `/api/tree?dataset=...`                   | UPGMA tree from binary Jaccard distances                      |
| `/api/statistics?dataset=...&genomes=...` | Family frequency, randomized accumulation/core curves         |
| `/api/similarity?dataset=...`             | Pairwise Jaccard similarity matrix                            |
| `POST /api/compare`                       | Shared, unique and pairwise families for 2–8 selected genomes |
| `/api/export/families.tsv?dataset=...`    | Tab-delimited family matrix                                   |
| `/api/export/families.csv?dataset=...`    | CSV family matrix                                             |
| `/api/export/tree.newick?dataset=...`     | Gene-content dendrogram in Newick syntax                      |

The API is intentionally small and local-only. It does not currently provide FastAPI-generated OpenAPI, pagination, PostgreSQL, authentication, or background analysis jobs.

## Calculation notes

- Family prevalence is distinct genomes containing the family divided by genomes in the selected set.
- Strict core = present in every genome; soft core = present in at least 95% but not every genome; shell = at least 15% below soft-core; cloud = at least two genomes below 15%; singleton = exactly one genome. The primary visualization groups strict+soft as core and shell as accessory; cloud and singleton stay separate.
- Pairwise similarity is Jaccard on binary family presence; copy number does not affect distance.
- The gene-content tree is UPGMA on `1 − Jaccard similarity`. It is not a sequence-based evolutionary phylogeny.
- Accumulation and core decay use 100 seeded random genome-order permutations. The displayed min/max is the observed range, not a confidence interval.
- GFF3 source coordinates are normalized to zero-based, half-open intervals in the parser.

## Architecture and current limits

- `apps/web/` is the responsive JavaScript/SVG client.
- `services/api/app.py` is the standard-library local HTTP API.
- `services/api/analysis.py` holds pure Jaccard, UPGMA, comparison and permutation calculations.
- `services/api/ingest.py` validates FASTA/GFF3 and parses Panaroo matrix output.
- `scripts/ingest_real.py` uses NCBI Datasets for RefSeq discovery and acquisition.
- `scripts/run_panaroo.py` validates downloaded inputs, runs Panaroo and writes a provenance-bearing normalized artifact.
- `services/api/tests/` covers fixture calculations, input parsing, matrix import, selection logic, and package path safety.

The local real-data workflow was run on 24 RefSeq candidates: four assemblies with out-of-bounds GFF3 CDS coordinates were explicitly excluded, and Panaroo 1.8.1 imported 5,488 families for the remaining 20 genomes. NCBI Datasets CLI 18.41.1 produced the selection and package. The resulting matrix is available in the app and API. This convenience cohort is not a curated biological study set. Do not cite the synthetic demo as biological evidence.

This is a local, single-user exploration build. It has no accounts, database, background jobs, upload workflow, or hosted deployment. GTDB classification, sequence-based phylogeny, synteny, rearrangement analysis, and pangenome graph visualization are unavailable. The analysis was verified locally; hosted and production behavior has not been verified.

Focused notes are in [`docs/`](docs/). The live source pages link to NCBI, Panaroo, and GTDB materials.

## License

MIT. See [LICENSE](LICENSE).
