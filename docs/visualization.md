# Visualization behavior

The web client separates the compact demo API data from SVG/HTML presentation functions in `apps/web/app.js`. Matrix rows encode absence, one copy, and multiple copies with both color and accessible labels. The selected genome list is serialized in the URL and updates the matrix, clustering highlights, comparison form, and selected-subset statistics.

The gene-content tree is a rectangular UPGMA dendrogram calculated from binary family presence. Its label and Newick export explicitly identify it as gene-content clustering. It does not stand in for a core-genome phylogeny. The pairwise heatmap shows exact Jaccard similarities. Randomized accumulation and core-decay charts use seeded permutations and observed min/max ranges.

The circular demo diagram is marked synthetic. Real pangenome artifacts omit family coordinates until locus reconciliation is implemented; the UI shows a coordinate-unavailable state instead of reusing demo positions. Synteny similarly remains unavailable until a validated contig-aware gene-order analysis exists.

The app respects reduced-motion preference and includes light/dark styles. The current renderer has only been checked with the local 20-genome cohort. Larger matrices have not been performance validated, and virtualization or tiled rendering is not available.
