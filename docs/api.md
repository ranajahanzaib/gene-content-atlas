# Local API

The Python standard-library HTTP server hosts the UI and JSON API at the same origin. Query `dataset=demo-lp-v1` selects the synthetic fixture; `dataset=real-lp-refseq-panaroo` is exposed only when its normalized Panaroo artifact exists. Unsupported IDs return 404 and never fall back to demo data.

`POST /api/compare` accepts JSON `{ "dataset_id": "demo-lp-v1", "genome_ids": ["LP-01", "LP-02"] }`. The genome list must contain 2–8 unique IDs in the selected dataset. Bad selections return 400; unavailable datasets return 404.

`GET /api/statistics?dataset=...&genomes=LP-01,LP-02` calculates for that subset. Omitting `genomes` calculates for the full dataset.

The API is local-only and intentionally dependency-free. It currently has no generated OpenAPI contract, server-side pagination, authentication, rate limiting, database, or job queue. Do not expose it on an untrusted network.
