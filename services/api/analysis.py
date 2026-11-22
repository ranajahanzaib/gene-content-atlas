"""Transparent calculations over family-by-genome presence data.

These routines intentionally do not infer sequence phylogeny, synteny, or
sequence conservation. Their inputs are family membership and copy counts only.
"""
from __future__ import annotations

import random
from collections import Counter
from typing import Any


def presence_sets(dataset: dict[str, Any]) -> dict[str, set[str]]:
    genomes = {strain["id"]: set() for strain in dataset["strains"]}
    for family in dataset["families"]:
        for genome_id in family["members"]:
            if genome_id not in genomes:
                raise ValueError(f"Family {family['id']} references unknown genome {genome_id}")
            genomes[genome_id].add(family["id"])
    return genomes


def jaccard_similarity(left: set[str], right: set[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 1.0


def similarity_matrix(dataset: dict[str, Any]) -> dict[str, Any]:
    sets = presence_sets(dataset)
    ids = list(sets)
    matrix = [[jaccard_similarity(sets[a], sets[b]) for b in ids] for a in ids]
    return {"genomes": ids, "metric": "Jaccard similarity of binary family presence", "matrix": matrix}


def upgma(dataset: dict[str, Any]) -> dict[str, Any]:
    """Build deterministic UPGMA from Jaccard distances, preserving labels."""
    sets = presence_sets(dataset)
    if not sets:
        raise ValueError("At least one genome is required")
    ids = sorted(sets)
    clusters: dict[tuple[str, ...], dict[str, Any]] = {
        (genome_id,): {"name": genome_id, "members": (genome_id,), "height": 0.0}
        for genome_id in ids
    }
    distances: dict[tuple[tuple[str, ...], tuple[str, ...]], float] = {}

    def key(a: tuple[str, ...], b: tuple[str, ...]):
        return (a, b) if a < b else (b, a)

    for i, left in enumerate(ids):
        for right in ids[i + 1:]:
            distances[key((left,), (right,))] = 1.0 - jaccard_similarity(sets[left], sets[right])

    while len(clusters) > 1:
        candidates = [(dist, a, b) for (a, b), dist in distances.items() if a in clusters and b in clusters]
        distance, a, b = min(candidates, key=lambda row: (row[0], row[1], row[2]))
        ca, cb = clusters.pop(a), clusters.pop(b)
        members = tuple(sorted(a + b))
        height = distance / 2.0
        node = {
            "name": None,
            "members": list(members),
            "height": height,
            "distance": distance,
            "children": [
                {"node": ca, "branch_length": max(0.0, height - ca["height"])},
                {"node": cb, "branch_length": max(0.0, height - cb["height"])},
            ],
        }
        for other in list(clusters):
            da = distances[key(a, other)]
            db = distances[key(b, other)]
            distances[key(members, other)] = (len(a) * da + len(b) * db) / (len(a) + len(b))
        clusters[members] = node
        # Old distances are no longer needed and would only grow quadratic memory.
        for old_key in list(distances):
            if a in old_key or b in old_key:
                distances.pop(old_key, None)
    root = next(iter(clusters.values()))
    return {"method": "UPGMA on binary gene-content Jaccard distance", "interpretation": "Gene-content clustering only; not a sequence-based evolutionary phylogeny.", "root": root, "newick": to_newick(root) + ";"}


def to_newick(node: dict[str, Any]) -> str:
    if node.get("name") is not None:
        return node["name"].replace("'", "''")
    pieces = []
    for child in node["children"]:
        pieces.append(f"{to_newick(child['node'])}:{child['branch_length']:.8f}")
    return f"({','.join(pieces)})"


def pangenome_statistics(dataset: dict[str, Any], genome_ids: list[str] | None = None, permutations: int = 100, seed: int = 20261122) -> dict[str, Any]:
    available = presence_sets(dataset)
    genome_ids = list(available) if genome_ids is None else genome_ids
    if not genome_ids:
        raise ValueError("Select at least one genome")
    if len(set(genome_ids)) != len(genome_ids):
        raise ValueError("Genome selections must be unique")
    unknown = [genome_id for genome_id in genome_ids if genome_id not in available]
    if unknown:
        raise ValueError(f"Unknown genome ID: {unknown[0]}")
    sets = {genome_id: available[genome_id] for genome_id in genome_ids}
    genome_ids = list(sets)
    if not genome_ids:
        raise ValueError("At least one genome is required")
    selected_ids = set(genome_ids)
    frequencies = Counter(sum(member in selected_ids for member in f["members"]) for f in dataset["families"] if any(member in selected_ids for member in f["members"]))
    rng = random.Random(seed)
    accumulation = [[] for _ in genome_ids]
    core_decay = [[] for _ in genome_ids]
    for _ in range(max(1, permutations)):
        order = genome_ids.copy()
        rng.shuffle(order)
        union: set[str] = set()
        intersection: set[str] | None = None
        for index, genome_id in enumerate(order):
            current = sets[genome_id]
            union.update(current)
            intersection = current.copy() if intersection is None else intersection & current
            accumulation[index].append(len(union))
            core_decay[index].append(len(intersection))
    curve = []
    for index in range(len(genome_ids)):
        values = accumulation[index]
        core_values = core_decay[index]
        curve.append({
            "genomes_sampled": index + 1,
            "pangenome_mean": sum(values) / len(values),
            "pangenome_min": min(values),
            "pangenome_max": max(values),
            "core_mean": sum(core_values) / len(core_values),
            "core_min": min(core_values),
            "core_max": max(core_values),
        })
    return {
        "genomes": genome_ids,
        "selection_mode": "selected genomes" if len(genome_ids) != len(available) else "full dataset",
        "family_frequency": [{"genomes_present": k, "families": frequencies.get(k, 0)} for k in range(1, len(genome_ids) + 1)],
        "accumulation": curve,
        "permutations": max(1, permutations),
        "seed": seed,
        "variation_note": "Min/max are the observed range across randomized genome orders, not a confidence interval.",
        "gene_copies_per_genome": [{"genome_id": genome_id, "copies": sum(f["copies"].get(genome_id, 0) for f in dataset["families"])} for genome_id in genome_ids],
    }


def compare_genomes(dataset: dict[str, Any], genome_ids: list[str]) -> dict[str, Any]:
    available = presence_sets(dataset)
    if len(genome_ids) < 2 or len(genome_ids) > 8:
        raise ValueError("Select between 2 and 8 genomes")
    if len(set(genome_ids)) != len(genome_ids):
        raise ValueError("Genome selections must be unique")
    unknown = [genome_id for genome_id in genome_ids if genome_id not in available]
    if unknown:
        raise ValueError(f"Unknown genome ID: {unknown[0]}")
    sets = [available[genome_id] for genome_id in genome_ids]
    shared = set.intersection(*sets)
    union = set.union(*sets)
    return {
        "genomes": genome_ids,
        "shared_family_count": len(shared),
        "union_family_count": len(union),
        "unique_families": {genome_id: sorted(available[genome_id] - set.union(*(s for i, s in enumerate(sets) if i != index))) for index, genome_id in enumerate(genome_ids)},
        "pairwise_similarity": {a: {b: jaccard_similarity(available[a], available[b]) for b in genome_ids} for a in genome_ids},
        "family_ids": sorted(union),
        "metric": "Binary gene-family presence; Jaccard similarity. Copy counts are not used in distances.",
    }
