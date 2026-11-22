# Scientific definitions

## Frequency classes

For `N` selected genomes and `k` genomes containing a family, prevalence is `k / N`. The exclusive partition applies in order: strict core (`k = N`); soft core (`k < N` and `k >= ceil(0.95N)`); singleton (`k = 1`); shell (`k >= ceil(0.15N)`); cloud (remaining `k >= 2`). This ordering resolves overlaps in small cohorts. The main chart combines strict and soft core, combines shell as accessory, and displays cloud and singleton separately. PPanGGOLiN model-derived persistent/shell/cloud categories are not implemented and must not be equated with these frequency bins.

## Similarity and tree

Each genome is represented by its set of present family IDs; copy number is ignored. Jaccard similarity is `|A ∩ B| / |A ∪ B|`; Jaccard distance is `1 - similarity`. If both sets are empty, similarity is defined as 1. UPGMA joins the closest clusters and updates cluster distances using size-weighted arithmetic means. The resulting branches summarize gene-content dissimilarity, not evolutionary time. It is labeled as a gene-content clustering tree.

## Accumulation and core decay

For each deterministic pseudo-random genome ordering, the pangenome accumulation curve counts the union of family sets as genomes are added; core decay counts their intersection. The plotted center is the mean across 100 permutations. The shaded lower and upper values are the observed minimum and maximum across permutations, not an inferential confidence interval. Selecting a genome subset recomputes the curves over that subset.

## Missing analyses

Presence alone does not establish sequence identity, orthology, synteny, or phylogeny. Circular positions in the demo are solely a visualization fixture. Real Panaroo membership supports frequency, Jaccard, UPGMA, and membership-based statistics. The application does not reconcile source GFF loci to Panaroo families and does not calculate alignment-based synteny or sequence phylogeny.
