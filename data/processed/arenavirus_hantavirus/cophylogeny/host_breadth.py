#!/usr/bin/env python3
"""
Compute host phylogenetic breadth metrics per virus from a host tree and virus-host associations.

Inputs:
- Host tree (Newick, branch lengths required; rooted or unrooted)
- Associations TSV with columns: virus_id, host_tip_label
- Optional family map TSV for family-level metrics

CLI:
python host_breadth.py \
  --tree host_tree.newick \
  --assoc associations.tsv \
  --out metrics.csv \
  [--family-map host_to_family.tsv] [--n-perm 10000] [--seed 1]
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

try:
    import dendropy
except ImportError as exc:  # pragma: no cover
    raise SystemExit("Missing dependency 'dendropy'. Install with: pip install dendropy") from exc


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(description="Host phylogenetic breadth metrics per virus")
    p.add_argument("--tree", required=True, help="Host tree in Newick format with branch lengths")
    p.add_argument("--assoc", required=True, help="Associations TSV with virus_id and host_tip_label")
    p.add_argument("--out", required=True, help="Output metrics CSV path")
    p.add_argument("--family-map", default=None, help="Optional TSV mapping host to family")
    p.add_argument("--n-perm", type=int, default=10000, help="Permutations for PD null model (default: 10000)")
    p.add_argument("--seed", type=int, default=1, help="Random seed (default: 1)")
    return p.parse_args()


def load_tree(path: Path) -> dendropy.Tree:
    """Load host tree from Newick and midpoint-root if unrooted."""
    if not path.exists():
        raise FileNotFoundError(f"Tree file not found: {path}")

    tree = dendropy.Tree.get(path=str(path), schema="newick", preserve_underscores=True)
    if not tree.is_rooted:
        tree.reroot_at_midpoint(update_bipartitions=False)
    return tree


def load_assoc(path: Path) -> pd.DataFrame:
    """Load associations TSV and validate required columns."""
    if not path.exists():
        raise FileNotFoundError(f"Associations file not found: {path}")

    df = pd.read_csv(path, sep="\t", dtype=str)
    required = {"virus_id", "host_tip_label"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Associations TSV missing columns: {sorted(missing)}")

    df = df[["virus_id", "host_tip_label"]].copy()
    df["virus_id"] = df["virus_id"].astype(str).str.strip()
    df["host_tip_label"] = df["host_tip_label"].astype(str).str.strip()
    df = df[(df["virus_id"] != "") & (df["host_tip_label"] != "")]
    df = df.drop_duplicates()
    return df


def tip_label_to_node(tree: dendropy.Tree) -> Dict[str, dendropy.Node]:
    """Map tip labels to tip nodes, raising on duplicate labels."""
    mapping: Dict[str, dendropy.Node] = {}
    for leaf in tree.leaf_node_iter():
        if leaf.taxon is None or leaf.taxon.label is None:
            raise ValueError("Tree contains a tip with missing label")
        label = leaf.taxon.label
        if label in mapping:
            raise ValueError(f"Duplicate tip label in tree: {label}")
        mapping[label] = leaf
    return mapping


def resolve_host_labels(
    labels: Iterable[str],
    tree_tips: Set[str],
) -> Tuple[Dict[str, str], List[str], int]:
    """
    Resolve host labels against tree tips.

    Resolution order:
    1) exact match
    2) replace spaces with underscores
    3) replace underscores with spaces

    Returns:
        mapping from original label -> resolved tree label
        list of missing labels
        number of labels resolved via non-exact normalization
    """
    out: Dict[str, str] = {}
    missing: List[str] = []
    non_exact = 0

    for h in labels:
        if h in tree_tips:
            out[h] = h
            continue

        c1 = h.replace(" ", "_")
        if c1 in tree_tips:
            out[h] = c1
            non_exact += 1
            continue

        c2 = h.replace("_", " ")
        if c2 in tree_tips:
            out[h] = c2
            non_exact += 1
            continue

        missing.append(h)

    return out, missing, non_exact


def node_depths(tree: dendropy.Tree) -> Dict[dendropy.Node, float]:
    """Compute root-to-node depths using branch lengths."""
    depths = {tree.seed_node: 0.0}
    for node in tree.preorder_node_iter():
        if node is tree.seed_node:
            continue
        parent = node.parent_node
        if parent is None:
            depths[node] = 0.0
            continue
        edge_len = 0.0 if node.edge.length is None else float(node.edge.length)
        depths[node] = depths[parent] + edge_len
    return depths


def total_tree_pd(tree: dendropy.Tree) -> float:
    """Sum branch lengths for all edges in the tree (Faith's PD over all tips)."""
    total = 0.0
    for node in tree.preorder_node_iter():
        if node is tree.seed_node:
            continue
        total += 0.0 if node.edge.length is None else float(node.edge.length)
    return float(total)


def mrca_two_nodes(
    a: dendropy.Node,
    b: dendropy.Node,
    depths: Dict[dendropy.Node, float],
) -> dendropy.Node:
    """Find MRCA of two nodes by lifting deeper node upward."""
    na, nb = a, b
    da, db = depths[na], depths[nb]

    while da > db and na.parent_node is not None:
        na = na.parent_node
        da = depths[na]
    while db > da and nb.parent_node is not None:
        nb = nb.parent_node
        db = depths[nb]

    while na is not nb:
        if na.parent_node is None or nb.parent_node is None:
            return na if na.parent_node is None else nb
        na = na.parent_node
        nb = nb.parent_node
    return na


def mrca_many_nodes(nodes: Sequence[dendropy.Node], depths: Dict[dendropy.Node, float]) -> Optional[dendropy.Node]:
    """Find MRCA for a sequence of nodes."""
    if len(nodes) == 0:
        return None
    m = nodes[0]
    for n in nodes[1:]:
        m = mrca_two_nodes(m, n, depths)
    return m


def faith_pd_from_nodes(
    nodes: Sequence[dendropy.Node],
    depths: Dict[dendropy.Node, float],
) -> float:
    """
    Faith's PD of minimal subtree spanning nodes.

    PD is computed as the sum of unique branch lengths in the union of
    paths from each node to the MRCA of the set.
    """
    n = len(nodes)
    if n == 0:
        return np.nan
    if n == 1:
        return 0.0

    m = mrca_many_nodes(nodes, depths)
    if m is None:
        return np.nan

    seen_edges: Set[int] = set()
    pd_sum = 0.0

    for node in nodes:
        cur = node
        while cur is not m:
            edge = cur.edge
            edge_id = id(edge)
            if edge_id not in seen_edges:
                seen_edges.add(edge_id)
                pd_sum += 0.0 if edge.length is None else float(edge.length)
            if cur.parent_node is None:
                break
            cur = cur.parent_node

    return float(pd_sum)


def node_distance(
    a: dendropy.Node,
    b: dendropy.Node,
    depths: Dict[dendropy.Node, float],
) -> float:
    """Patristic distance between any two nodes."""
    m = mrca_two_nodes(a, b, depths)
    return float(depths[a] + depths[b] - 2.0 * depths[m])


def build_host_distance_matrix(
    tree: dendropy.Tree,
    host_labels: Sequence[str],
) -> Tuple[np.ndarray, Dict[str, int]]:
    """Precompute host patristic matrix for provided host labels."""
    pdm = tree.phylogenetic_distance_matrix()
    tax_by_label = {tx.label: tx for tx in tree.taxon_namespace}

    n = len(host_labels)
    mat = np.zeros((n, n), dtype=float)
    for i in range(n):
        tx_i = tax_by_label[host_labels[i]]
        for j in range(i + 1, n):
            tx_j = tax_by_label[host_labels[j]]
            d = float(pdm.distance(tx_i, tx_j))
            mat[i, j] = d
            mat[j, i] = d

    idx = {h: i for i, h in enumerate(host_labels)}
    return mat, idx


def mpd_mntd_from_indices(indices: Sequence[int], dist_mat: np.ndarray) -> Tuple[float, float]:
    """Compute MPD and MNTD from a submatrix index list."""
    k = len(indices)
    if k < 2:
        return np.nan, np.nan

    idx = np.asarray(indices, dtype=int)
    sub = dist_mat[np.ix_(idx, idx)]

    upper = np.triu_indices(k, 1)
    mpd = float(np.mean(sub[upper]))

    near = sub.copy()
    np.fill_diagonal(near, np.inf)
    mntd = float(np.mean(np.min(near, axis=1)))
    return mpd, mntd


def load_family_map(path: Path, tree_hosts: Set[str]) -> Dict[str, str]:
    """Load optional host->family map from TSV with flexible column names."""
    if not path.exists():
        raise FileNotFoundError(f"Family map file not found: {path}")

    df = pd.read_csv(path, sep="\t", dtype=str)
    cols = {c.lower(): c for c in df.columns}

    host_col = None
    for c in ["host_tip_label", "host_id", "host", "tip", "tip_label"]:
        if c in cols:
            host_col = cols[c]
            break

    fam_col = None
    for c in ["host_family", "family", "fam"]:
        if c in cols:
            fam_col = cols[c]
            break

    if host_col is None or fam_col is None:
        if df.shape[1] >= 2:
            host_col, fam_col = df.columns[0], df.columns[1]
        else:
            raise ValueError(
                "Family map must contain host and family columns (e.g., host_tip_label, host_family)."
            )

    out: Dict[str, str] = {}
    resolver, _, _ = resolve_host_labels(df[host_col].astype(str).tolist(), tree_hosts)
    for row in df[[host_col, fam_col]].itertuples(index=False):
        raw_h = str(row[0]).strip()
        f = str(row[1]).strip()
        h = resolver.get(raw_h, "")
        if h and f:
            out[h] = f
    return out


def summarize_pd_null(obs: float, null_vals: np.ndarray) -> Tuple[float, float, float, float, float]:
    """Return null mean, sd, z, upper-tail p, two-sided p for PD."""
    arr = np.asarray(null_vals, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0 or not np.isfinite(obs):
        return np.nan, np.nan, np.nan, np.nan, np.nan

    mu = float(np.mean(arr))
    sd = float(np.std(arr, ddof=1)) if arr.size > 1 else np.nan
    z = float((obs - mu) / sd) if np.isfinite(sd) and sd > 0 else np.nan

    p_high = float((np.sum(arr >= obs) + 1) / (arr.size + 1))
    p_low = float((np.sum(arr <= obs) + 1) / (arr.size + 1))
    p_two = float(min(1.0, 2.0 * min(p_high, p_low)))
    return mu, sd, z, p_high, p_two


def methods_paragraph(n_perm: int, seed: int) -> str:
    """Manuscript-ready Methods paragraph for PD and null model."""
    return (
        "Host phylogenetic breadth was quantified per virus on the mammal host phylogeny using branch-length-aware "
        "metrics. Faith's phylogenetic diversity (PD) was calculated as the sum of branch lengths in the minimal subtree "
        "spanning all host tips associated with each virus; mean pairwise distance (MPD) and mean nearest-taxon distance "
        "(MNTD) were also computed from patristic distances among associated hosts. To account for differences in host "
        "richness, observed PD was normalized by the total PD of the full host tree and by the PD of all hosts observed in "
        "the association dataset. A null model was implemented by sampling random host sets without replacement from the "
        "observed host pool, matching each virus-specific host count, with {n_perm} permutations per host-count class. "
        "Empirical significance of PD was summarized with z-scores and permutation p-values, and analyses were made "
        "reproducible with a fixed random seed ({seed})."
    ).format(n_perm=n_perm, seed=seed)


def main() -> None:
    args = parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)

    tree = load_tree(Path(args.tree))
    assoc = load_assoc(Path(args.assoc))

    tip_to_node = tip_label_to_node(tree)
    tree_tips = set(tip_to_node.keys())

    assoc_hosts_raw = sorted(set(assoc["host_tip_label"]))
    host_resolver, missing_hosts, n_non_exact = resolve_host_labels(assoc_hosts_raw, tree_tips)

    missing_path = out_path.with_suffix(".missing_hosts.tsv")
    if missing_hosts:
        pd.DataFrame({"missing_host_tip_label": missing_hosts}).to_csv(missing_path, sep="\t", index=False)
        raise ValueError(
            f"{len(missing_hosts)} host_tip_label values from associations were not found in tree tips. "
            f"Details: {missing_path}"
        )

    # Harmonize association host labels to tree tip labels.
    assoc["host_tip_label"] = assoc["host_tip_label"].map(host_resolver)
    assoc = assoc.dropna(subset=["host_tip_label"]).copy()

    if n_non_exact > 0:
        print(
            f"INFO: Resolved {n_non_exact} host labels via space/underscore normalization.", file=sys.stderr
        )

    # Restrict to hosts observed in associations (background B for null model).
    bg_hosts = sorted(set(assoc["host_tip_label"]))
    host_dist, host_idx = build_host_distance_matrix(tree, bg_hosts)

    bg_nodes = [tip_to_node[h] for h in bg_hosts]
    bg_host_count = len(bg_hosts)

    depths = node_depths(tree)
    pd_tree_total = total_tree_pd(tree)
    pd_assoc_total = faith_pd_from_nodes(bg_nodes, depths)

    # Virus -> unique host tips
    virus_to_hosts: Dict[str, Set[str]] = defaultdict(set)
    for row in assoc.itertuples(index=False):
        virus_to_hosts[row.virus_id].add(row.host_tip_label)

    # Optional family mapping and family representative nodes (global family MRCAs).
    host_to_family: Dict[str, str] = {}
    family_to_rep_node: Dict[str, dendropy.Node] = {}
    if args.family_map is not None:
        host_to_family = load_family_map(Path(args.family_map), tree_tips)
        fam_to_hosts: Dict[str, List[str]] = defaultdict(list)
        for h, fam in host_to_family.items():
            fam_to_hosts[fam].append(h)
        for fam, fam_hosts in fam_to_hosts.items():
            fam_nodes = [tip_to_node[h] for h in fam_hosts if h in tip_to_node]
            if len(fam_nodes) == 0:
                continue
            rep = mrca_many_nodes(fam_nodes, depths)
            if rep is not None:
                family_to_rep_node[fam] = rep

    # Null PD distributions by unique n_hosts (reuse for speed).
    unique_k = sorted({len(hs) for hs in virus_to_hosts.values()})
    null_pd_by_k: Dict[int, np.ndarray] = {}

    for k in unique_k:
        if k <= 0:
            null_pd_by_k[k] = np.array([np.nan], dtype=float)
            continue
        if k == 1:
            null_pd_by_k[k] = np.zeros(args.n_perm, dtype=float)
            continue
        if k == bg_host_count:
            null_pd_by_k[k] = np.full(args.n_perm, pd_assoc_total, dtype=float)
            continue

        vals = np.empty(args.n_perm, dtype=float)
        for i in range(args.n_perm):
            sample_ix = rng.choice(bg_host_count, size=k, replace=False)
            sample_nodes = [bg_nodes[j] for j in sample_ix]
            vals[i] = faith_pd_from_nodes(sample_nodes, depths)
        null_pd_by_k[k] = vals

    results = []

    for virus in sorted(virus_to_hosts.keys()):
        hosts = sorted(virus_to_hosts[virus])
        n_hosts = len(hosts)
        host_nodes = [tip_to_node[h] for h in hosts]

        pd_species = faith_pd_from_nodes(host_nodes, depths) if n_hosts >= 1 else np.nan

        if n_hosts >= 2:
            ix = [host_idx[h] for h in hosts]
            mpd_species, mntd_species = mpd_mntd_from_indices(ix, host_dist)
        else:
            mpd_species, mntd_species = np.nan, np.nan

        norm_pd_tree = (pd_species / pd_tree_total) if np.isfinite(pd_species) and pd_tree_total > 0 else np.nan
        norm_pd_assoc = (pd_species / pd_assoc_total) if np.isfinite(pd_species) and pd_assoc_total > 0 else np.nan

        null_vals = null_pd_by_k.get(n_hosts, np.array([np.nan], dtype=float))
        null_mean, null_sd, null_z, null_p_high, null_p_two = summarize_pd_null(pd_species, null_vals)

        # Family-level metrics (MRCA representative per family)
        n_families = np.nan
        pd_family = np.nan
        mpd_family = np.nan
        mntd_family = np.nan

        if host_to_family:
            fams = sorted({host_to_family[h] for h in hosts if h in host_to_family and host_to_family[h] in family_to_rep_node})
            if fams:
                n_families = float(len(fams))
                reps = [family_to_rep_node[f] for f in fams]

                pd_family = faith_pd_from_nodes(reps, depths)

                if len(reps) >= 2:
                    dvals = []
                    nearest = []
                    for i in range(len(reps)):
                        ri = reps[i]
                        nn = np.inf
                        for j in range(len(reps)):
                            if i == j:
                                continue
                            rj = reps[j]
                            d = node_distance(ri, rj, depths)
                            dvals.append(d) if j > i else None
                            if d < nn:
                                nn = d
                        nearest.append(nn)
                    mpd_family = float(np.mean(dvals)) if dvals else np.nan
                    mntd_family = float(np.mean(nearest)) if nearest else np.nan

        results.append(
            {
                "virus_id": virus,
                "n_hosts": n_hosts,
                "Faith_PD_species": pd_species,
                "MPD_species": mpd_species,
                "MNTD_species": mntd_species,
                "Normalized_PD_species_tree": norm_pd_tree,
                "Normalized_PD_species_assoc": norm_pd_assoc,
                "n_families": n_families,
                "Faith_PD_family": pd_family,
                "MPD_family": mpd_family,
                "MNTD_family": mntd_family,
                "null_PD_mean": null_mean,
                "null_PD_sd": null_sd,
                "null_PD_z": null_z,
                "null_PD_p_high": null_p_high,
                "null_PD_p_two_sided": null_p_two,
                "n_perm": args.n_perm,
                "seed": args.seed,
            }
        )

    out_df = pd.DataFrame(results)
    out_df.to_csv(out_path, index=False)

    methods_path = out_path.with_suffix(".methods.txt")
    methods_path.write_text(methods_paragraph(args.n_perm, args.seed) + "\n", encoding="utf-8")

    metadata_path = out_path.with_suffix(".metadata.tsv")
    pd.DataFrame(
        [
            {
                "tree": str(Path(args.tree).resolve()),
                "assoc": str(Path(args.assoc).resolve()),
                "out": str(out_path.resolve()),
                "n_tips_tree": len(tree_tips),
                "n_hosts_assoc": len(bg_hosts),
                "n_viruses": len(virus_to_hosts),
                "pd_tree_total": pd_tree_total,
                "pd_assoc_total": pd_assoc_total,
                "n_perm": args.n_perm,
                "seed": args.seed,
                "family_map_used": bool(args.family_map),
            }
        ]
    ).to_csv(metadata_path, sep="\t", index=False)

    print(f"Wrote metrics: {out_path}")
    print(f"Wrote methods: {methods_path}")
    print(f"Wrote metadata: {metadata_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
