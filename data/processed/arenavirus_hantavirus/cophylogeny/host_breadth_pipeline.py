#!/usr/bin/env python3
"""
Quantify per-virus host phylogenetic breadth and test phylogenetic signal on virus trees.

Requirements implemented:
- Host patristic distances from a Newick host tree (branch lengths required)
- Per-virus observed metrics: MPD, MNTD, Faith's PD, MRCA depth, minimal clade cover
- Null models:
    * Null0: uniform host sampling
    * Null2: degree-weighted host sampling
    * Null3: degree-preserving bipartite rewiring via double-edge swaps
- SES and empirical p-values for MPD/MNTD
- Pagel's lambda (BM likelihood) for SES traits on the virus tree
- Reproducible outputs, plots, and methods snippet

Python: 3.11+
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    import dendropy
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "Missing dependency: dendropy. Install with `pip install dendropy` and rerun."
    ) from exc


@dataclass
class NullStats:
    mean: float
    sd: float
    ses: float
    p_two_sided: float
    p_high: float
    n_null: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Host breadth pipeline: observed metrics, null models, SES, and Pagel's lambda."
    )
    parser.add_argument("--host_tree", required=True, help="Host tree Newick path with branch lengths.")
    parser.add_argument("--virus_tree", required=True, help="Virus tree Newick path with branch lengths.")
    parser.add_argument("--assoc", required=True, help="Associations TSV path.")
    parser.add_argument("--outdir", required=True, help="Output directory.")

    parser.add_argument("--n_reps", type=int, default=999, help="Null replicates (default: 999).")
    parser.add_argument("--alpha", type=float, default=1.0, help="Degree exponent for Null2 (default: 1.0).")
    parser.add_argument("--q", type=float, default=0.7, help="Depth fraction cutoff for minimal clade cover (default: 0.7).")
    parser.add_argument("--seed", type=int, default=1, help="Random seed (default: 1).")

    parser.add_argument(
        "--background_mode",
        choices=["observed", "file"],
        default="observed",
        help="Background hosts for Null0/Null2: observed or from file (default: observed).",
    )
    parser.add_argument(
        "--background_hosts",
        default=None,
        help="Optional file with one host_id per line (used when --background_mode file).",
    )

    parser.add_argument(
        "--swaps_per_replicate",
        type=int,
        default=None,
        help="Accepted swaps per Null3 replicate (default: 20*|E|).",
    )
    parser.add_argument(
        "--lambda_trait_null",
        choices=["null0", "null2", "null3"],
        default="null3",
        help="Which SES null model to use for Pagel's lambda traits (default: null3).",
    )
    return parser.parse_args()


def load_tree_newick(path: Path, label: str) -> dendropy.Tree:
    if not path.exists():
        raise FileNotFoundError(f"{label} tree not found: {path}")
    tree = dendropy.Tree.get(path=str(path), schema="newick", preserve_underscores=True)
    if not tree.is_rooted:
        tree.reroot_at_midpoint(update_bipartitions=False)
    return tree


def tip_labels(tree: dendropy.Tree) -> List[str]:
    labels = []
    for leaf in tree.leaf_node_iter():
        if leaf.taxon is None or leaf.taxon.label is None:
            raise ValueError("Tree has a tip without label.")
        labels.append(leaf.taxon.label)
    return labels


def node_depths(tree: dendropy.Tree) -> Dict[dendropy.Node, float]:
    depths: Dict[dendropy.Node, float] = {tree.seed_node: 0.0}
    for node in tree.preorder_node_iter():
        if node is tree.seed_node:
            continue
        parent = node.parent_node
        if parent is None:
            depths[node] = 0.0
            continue
        elen = 0.0 if node.edge.length is None else float(node.edge.length)
        depths[node] = depths[parent] + elen
    return depths


def tree_height(tree: dendropy.Tree, depths: Dict[dendropy.Node, float]) -> float:
    heights = [depths[n] for n in tree.leaf_node_iter()]
    if not heights:
        raise ValueError("Tree has no tips.")
    return float(max(heights))


def build_tip_node_map(tree: dendropy.Tree) -> Dict[str, dendropy.Node]:
    out: Dict[str, dendropy.Node] = {}
    for leaf in tree.leaf_node_iter():
        label = leaf.taxon.label
        if label in out:
            raise ValueError(f"Duplicate tip label in tree: {label}")
        out[label] = leaf
    return out


def compute_host_distance_matrix(
    tree: dendropy.Tree, host_labels: Sequence[str]
) -> Tuple[np.ndarray, Dict[str, int]]:
    """Compute full patristic matrix for host tips listed in host_labels."""
    pdm = tree.phylogenetic_distance_matrix()

    taxon_by_label = {}
    for tx in tree.taxon_namespace:
        taxon_by_label[tx.label] = tx

    n = len(host_labels)
    mat = np.zeros((n, n), dtype=float)
    for i in range(n):
        tx_i = taxon_by_label[host_labels[i]]
        for j in range(i + 1, n):
            tx_j = taxon_by_label[host_labels[j]]
            d = float(pdm.distance(tx_i, tx_j))
            mat[i, j] = d
            mat[j, i] = d

    idx = {h: i for i, h in enumerate(host_labels)}
    return mat, idx


def _pick_column(lower_to_original: Dict[str, str], aliases: Sequence[str]) -> Optional[str]:
    for a in aliases:
        if a in lower_to_original:
            return lower_to_original[a]
    return None


def resolve_labels(values: Iterable[str], labels: Set[str]) -> Tuple[Dict[str, str], List[str], int]:
    """
    Resolve IDs to tree labels allowing underscore/space normalization.

    Resolution order:
    1) exact
    2) replace spaces with underscores
    3) replace underscores with spaces
    """
    mapping: Dict[str, str] = {}
    missing: List[str] = []
    n_normalized = 0

    for v in values:
        if v in labels:
            mapping[v] = v
            continue

        c1 = v.replace(" ", "_")
        if c1 in labels:
            mapping[v] = c1
            n_normalized += 1
            continue

        c2 = v.replace("_", " ")
        if c2 in labels:
            mapping[v] = c2
            n_normalized += 1
            continue

        missing.append(v)

    return mapping, sorted(missing), n_normalized


def parse_associations(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Associations file not found: {path}")

    # Auto-detect delimiter so comma-separated residual tables also work.
    df = pd.read_csv(path, sep=None, engine="python", dtype=str)
    cols = [str(c).strip() for c in df.columns]
    lower_to_original = {c.lower(): c for c in cols}

    virus_col = _pick_column(lower_to_original, ["virus_id", "virus", "parasite", "pathogen"])
    host_col = _pick_column(lower_to_original, ["host_id", "host_tip_label", "host", "tip"])

    if virus_col is None or host_col is None:
        raise ValueError(
            "Associations file must contain virus/host columns. "
            "Accepted names include virus_id|virus and host_id|host_tip_label|host. "
            f"Found columns: {cols}"
        )

    evidence_col = _pick_column(
        lower_to_original, ["evidence_weight", "weight", "evidence", "w"]
    )
    family_col = _pick_column(
        lower_to_original, ["host_family", "family", "fam", "hostfam"]
    )

    df = df.rename(columns={virus_col: "virus_id", host_col: "host_id"}).copy()
    if evidence_col is not None:
        df = df.rename(columns={evidence_col: "evidence_weight"})
    if family_col is not None:
        df = df.rename(columns={family_col: "host_family"})

    if "evidence_weight" not in df.columns:
        df["evidence_weight"] = ""
    if "host_family" not in df.columns:
        df["host_family"] = ""

    df["virus_id"] = df["virus_id"].astype(str).str.strip()
    df["host_id"] = df["host_id"].astype(str).str.strip()
    df["host_family"] = df["host_family"].fillna("").astype(str).str.strip()

    w = pd.to_numeric(df["evidence_weight"], errors="coerce")
    df["evidence_weight"] = w.fillna(1.0).astype(float)

    df = df[(df["virus_id"] != "") & (df["host_id"] != "")].copy()

    # Collapse duplicate virus-host rows; retain strongest evidence and first non-empty family label.
    def first_non_empty(vals: pd.Series) -> str:
        for v in vals:
            vv = str(v).strip()
            if vv:
                return vv
        return ""

    collapsed = (
        df.groupby(["virus_id", "host_id"], as_index=False)
        .agg(
            evidence_weight=("evidence_weight", "max"),
            host_family=("host_family", first_non_empty),
        )
        .copy()
    )
    return collapsed


def validate_association_ids(
    assoc: pd.DataFrame,
    host_labels: Set[str],
    virus_labels: Set[str],
    outdir: Path,
) -> Tuple[pd.DataFrame, int, int]:
    host_map, unmatched_hosts, n_host_normalized = resolve_labels(set(assoc["host_id"]), host_labels)
    virus_map, unmatched_viruses, n_virus_normalized = resolve_labels(set(assoc["virus_id"]), virus_labels)

    report = outdir / "unmatched_ids.tsv"
    rows = []
    rows.extend([{"type": "host", "id": h} for h in unmatched_hosts])
    rows.extend([{"type": "virus", "id": v} for v in unmatched_viruses])
    pd.DataFrame(rows).to_csv(report, sep="\t", index=False)

    if unmatched_hosts or unmatched_viruses:
        msg = [
            "Association IDs do not fully match tree tip labels.",
            f"Unmatched hosts: {len(unmatched_hosts)}",
            f"Unmatched viruses: {len(unmatched_viruses)}",
            f"Details written to: {report}",
        ]
        if unmatched_hosts:
            msg.append("Example unmatched host IDs: " + ", ".join(unmatched_hosts[:10]))
        if unmatched_viruses:
            msg.append("Example unmatched virus IDs: " + ", ".join(unmatched_viruses[:10]))
        raise ValueError("\n".join(msg))

    out = assoc.copy()
    out["host_id"] = out["host_id"].map(host_map)
    out["virus_id"] = out["virus_id"].map(virus_map)
    return out, n_host_normalized, n_virus_normalized


def host_background(
    assoc: pd.DataFrame,
    host_tip_labels: Set[str],
    mode: str,
    host_file: Optional[Path],
) -> List[str]:
    observed_hosts = sorted(set(assoc["host_id"]))
    if mode == "observed":
        return observed_hosts

    if host_file is None:
        raise ValueError("--background_mode file requires --background_hosts")
    if not host_file.exists():
        raise FileNotFoundError(f"Background host list not found: {host_file}")

    chosen = []
    with open(host_file, "r", encoding="utf-8") as handle:
        for line in handle:
            h = line.strip()
            if h:
                chosen.append(h)

    chosen = sorted(set(chosen))
    chosen_map, missing, _ = resolve_labels(chosen, host_tip_labels)
    if missing:
        raise ValueError(
            "Background host list contains IDs missing from host tree. "
            f"Example: {', '.join(missing[:10])}"
        )
    resolved = sorted(set(chosen_map.values()))
    if len(resolved) == 0:
        raise ValueError("Background host list is empty.")
    return resolved


def compute_mpd_mntd(indices: Sequence[int], dist_mat: np.ndarray) -> Tuple[float, float]:
    k = len(indices)
    if k < 2:
        return np.nan, np.nan

    idx = np.asarray(indices, dtype=int)
    sub = dist_mat[np.ix_(idx, idx)]

    iu = np.triu_indices(k, 1)
    mpd = float(np.mean(sub[iu]))

    tmp = sub.copy()
    np.fill_diagonal(tmp, np.inf)
    nearest = np.min(tmp, axis=1)
    mntd = float(np.mean(nearest))
    return mpd, mntd


def faith_pd(hosts: Sequence[str], tip_to_node: Dict[str, dendropy.Node], host_tree: dendropy.Tree) -> float:
    if len(hosts) == 0:
        return np.nan
    if len(hosts) == 1:
        return 0.0

    mrca = host_tree.mrca(taxon_labels=list(hosts))
    if mrca is None:
        return np.nan

    seen_edges: Set[int] = set()
    total = 0.0

    for h in hosts:
        node = tip_to_node[h]
        while node is not mrca:
            edge = node.edge
            eid = id(edge)
            if eid not in seen_edges:
                seen_edges.add(eid)
                total += 0.0 if edge.length is None else float(edge.length)
            parent = node.parent_node
            if parent is None:
                break
            node = parent

    return float(total)


def mrca_depth_norm(
    hosts: Sequence[str],
    host_tree: dendropy.Tree,
    depths: Dict[dendropy.Node, float],
    total_height: float,
) -> float:
    if len(hosts) == 0 or total_height <= 0:
        return np.nan
    mrca = host_tree.mrca(taxon_labels=list(hosts))
    if mrca is None:
        return np.nan
    return float(depths[mrca] / total_height)


def clade_cut_assignment(
    host_tree: dendropy.Tree,
    tip_to_node: Dict[str, dendropy.Node],
    depths: Dict[dendropy.Node, float],
    q: float,
) -> Dict[str, dendropy.Node]:
    if not (0.0 <= q <= 1.0):
        raise ValueError("q must be within [0, 1].")

    total_height = tree_height(host_tree, depths)
    threshold = q * total_height
    assign: Dict[str, dendropy.Node] = {}

    for label, node in tip_to_node.items():
        cur = node
        chosen = host_tree.seed_node
        while cur.parent_node is not None:
            d_cur = depths[cur]
            d_parent = depths[cur.parent_node]
            if d_cur >= threshold and d_parent < threshold:
                chosen = cur
                break
            cur = cur.parent_node
        assign[label] = chosen

    return assign


def minimal_clade_cover(hosts: Sequence[str], assignment: Dict[str, dendropy.Node]) -> float:
    if len(hosts) == 0:
        return np.nan
    ids = {id(assignment[h]) for h in hosts}
    return float(len(ids))


def family_summary(
    hosts: Sequence[str],
    host_to_family: Dict[str, str],
    host_index: Dict[str, int],
    host_dist: np.ndarray,
) -> Tuple[float, str, float]:
    fams = []
    for h in hosts:
        fam = host_to_family.get(h, "")
        fam = fam.strip()
        if fam:
            fams.append(fam)

    if not fams:
        return np.nan, "", np.nan

    counts = Counter(fams)
    n_unique = float(len(counts))
    fam_counts = ";".join(f"{k}:{v}" for k, v in sorted(counts.items()))

    # Between-family MPD among host pairs with non-empty family labels.
    labeled_hosts = [h for h in hosts if host_to_family.get(h, "").strip()]
    if len(labeled_hosts) < 2:
        return n_unique, fam_counts, np.nan

    vals = []
    for i in range(len(labeled_hosts)):
        hi = labeled_hosts[i]
        fi = host_to_family[hi].strip()
        ii = host_index[hi]
        for j in range(i + 1, len(labeled_hosts)):
            hj = labeled_hosts[j]
            fj = host_to_family[hj].strip()
            if fi == fj:
                continue
            jj = host_index[hj]
            vals.append(host_dist[ii, jj])

    between = float(np.mean(vals)) if vals else np.nan
    return n_unique, fam_counts, between


def empirical_stats(obs: float, null_vals: Sequence[float]) -> NullStats:
    arr = np.asarray(null_vals, dtype=float)
    arr = arr[np.isfinite(arr)]
    if not np.isfinite(obs) or arr.size == 0:
        return NullStats(np.nan, np.nan, np.nan, np.nan, np.nan, int(arr.size))

    mean = float(np.mean(arr))
    sd = float(np.std(arr, ddof=1)) if arr.size > 1 else np.nan
    ses = float((obs - mean) / sd) if np.isfinite(sd) and sd > 0 else np.nan

    p_two = float((np.sum(np.abs(arr - mean) >= abs(obs - mean)) + 1) / (arr.size + 1))
    p_high = float((np.sum(arr >= obs) + 1) / (arr.size + 1))
    return NullStats(mean, sd, ses, p_two, p_high, int(arr.size))


def simulate_null_sampling(
    virus_to_k: Dict[str, int],
    background_hosts: List[str],
    host_index: Dict[str, int],
    host_dist: np.ndarray,
    rng: np.random.Generator,
    n_reps: int,
    probs: Optional[np.ndarray] = None,
) -> Tuple[Dict[str, List[float]], Dict[str, List[float]]]:
    mpd_null = {v: [] for v in virus_to_k}
    mntd_null = {v: [] for v in virus_to_k}

    bg_indices = np.array([host_index[h] for h in background_hosts], dtype=int)
    n_bg = len(bg_indices)

    for v, k in virus_to_k.items():
        if k < 2 or k > n_bg:
            mpd_null[v] = [np.nan] * n_reps
            mntd_null[v] = [np.nan] * n_reps
            continue

        if probs is not None:
            p = probs.copy()
            if not np.isclose(np.sum(p), 1.0):
                s = np.sum(p)
                if s <= 0:
                    mpd_null[v] = [np.nan] * n_reps
                    mntd_null[v] = [np.nan] * n_reps
                    continue
                p = p / s

            # Need at least k hosts with non-zero probability for sampling without replacement.
            if np.sum(p > 0) < k:
                mpd_null[v] = [np.nan] * n_reps
                mntd_null[v] = [np.nan] * n_reps
                continue
        else:
            p = None

        vals_mpd = []
        vals_mntd = []
        for _ in range(n_reps):
            chosen = rng.choice(bg_indices, size=k, replace=False, p=p)
            mpd, mntd = compute_mpd_mntd(chosen, host_dist)
            vals_mpd.append(mpd)
            vals_mntd.append(mntd)

        mpd_null[v] = vals_mpd
        mntd_null[v] = vals_mntd

    return mpd_null, mntd_null


def run_degree_preserving_rewiring(
    virus_to_hosts: Dict[str, Set[str]],
    host_index: Dict[str, int],
    host_dist: np.ndarray,
    n_reps: int,
    rng: np.random.Generator,
    swaps_per_replicate: int,
) -> Tuple[Dict[str, List[float]], Dict[str, List[float]], float, int, int]:
    """
    Null3: degree-preserving bipartite rewiring via double-edge swaps.

    Returns:
        mpd_null, mntd_null, acceptance_rate, total_accepted, total_attempted
    """
    viruses = sorted(virus_to_hosts.keys())
    virus_to_ix = {v: i for i, v in enumerate(viruses)}

    all_hosts = sorted({h for hs in virus_to_hosts.values() for h in hs})
    host_to_ix = {h: i for i, h in enumerate(all_hosts)}

    base_edges = []
    for v in viruses:
        vi = virus_to_ix[v]
        for h in sorted(virus_to_hosts[v]):
            hi = host_to_ix[h]
            base_edges.append((vi, hi))

    if not base_edges:
        raise ValueError("No edges in association graph for Null3 rewiring.")

    E = len(base_edges)
    mpd_null = {v: [] for v in viruses}
    mntd_null = {v: [] for v in viruses}

    total_attempted = 0
    total_accepted = 0

    for _rep in range(n_reps):
        edges = list(base_edges)
        edge_set = set(edges)

        accepted = 0
        attempted = 0
        max_attempts = max(1000, swaps_per_replicate * 100)

        while accepted < swaps_per_replicate and attempted < max_attempts:
            i, j = rng.integers(0, E, size=2)
            if i == j:
                continue
            (v1, h1) = edges[i]
            (v2, h2) = edges[j]
            if v1 == v2 or h1 == h2:
                attempted += 1
                continue

            new1 = (v1, h2)
            new2 = (v2, h1)

            attempted += 1

            if new1 == new2:
                continue
            if new1 in edge_set or new2 in edge_set:
                continue

            edge_set.remove((v1, h1))
            edge_set.remove((v2, h2))
            edge_set.add(new1)
            edge_set.add(new2)

            edges[i] = new1
            edges[j] = new2
            accepted += 1

        total_attempted += attempted
        total_accepted += accepted

        v_hosts_int = [set() for _ in viruses]
        for vi, hi in edges:
            v_hosts_int[vi].add(hi)

        for v in viruses:
            vi = virus_to_ix[v]
            hosts_ix = [host_index[all_hosts[h]] for h in v_hosts_int[vi]]
            mpd, mntd = compute_mpd_mntd(hosts_ix, host_dist)
            mpd_null[v].append(mpd)
            mntd_null[v].append(mntd)

    acceptance_rate = float(total_accepted / total_attempted) if total_attempted > 0 else np.nan
    return mpd_null, mntd_null, acceptance_rate, total_accepted, total_attempted


def bm_covariance(
    tree: dendropy.Tree,
    tip_order: Sequence[str],
    tip_to_node: Dict[str, dendropy.Node],
    depths: Dict[dendropy.Node, float],
) -> np.ndarray:
    n = len(tip_order)
    C = np.zeros((n, n), dtype=float)

    for i, ti in enumerate(tip_order):
        ni = tip_to_node[ti]
        C[i, i] = depths[ni]

    for i in range(n):
        ti = tip_order[i]
        for j in range(i + 1, n):
            tj = tip_order[j]
            m = tree.mrca(taxon_labels=[ti, tj])
            if m is None:
                raise ValueError(f"Could not find MRCA for virus tips: {ti}, {tj}")
            cov = depths[m]
            C[i, j] = cov
            C[j, i] = cov
    return C


def loglik_pagel_lambda(y: np.ndarray, C: np.ndarray, lam: float, jitter: float = 1e-8) -> float:
    n = y.shape[0]
    if n < 2:
        return np.nan

    V = lam * C
    np.fill_diagonal(V, np.diag(C))

    # Numerical stabilization.
    scale = float(np.mean(np.diag(V))) if np.mean(np.diag(V)) > 0 else 1.0
    V = V + np.eye(n) * (jitter * scale)

    try:
        sign, logdet = np.linalg.slogdet(V)
        if sign <= 0:
            return np.nan
        Vinv = np.linalg.inv(V)
    except np.linalg.LinAlgError:
        return np.nan

    one = np.ones(n)
    denom = float(one @ Vinv @ one)
    if denom <= 0:
        return np.nan
    mu = float((one @ Vinv @ y) / denom)

    resid = y - mu
    quad = float(resid @ Vinv @ resid)
    ll = -0.5 * (n * math.log(2.0 * math.pi) + logdet + quad)
    return ll


def estimate_lambda(y: np.ndarray, C: np.ndarray) -> Dict[str, float]:
    out = {
        "lambda_hat": np.nan,
        "ll_hat": np.nan,
        "ci_low": np.nan,
        "ci_high": np.nan,
        "lr_vs_0": np.nan,
        "p_vs_0": np.nan,
        "lr_vs_1": np.nan,
        "p_vs_1": np.nan,
        "n_taxa": float(len(y)),
    }

    if len(y) < 3 or np.nanstd(y) < 1e-12:
        return out

    coarse = np.linspace(0.0, 1.0, 201)
    ll_coarse = np.array([loglik_pagel_lambda(y, C, float(l)) for l in coarse])
    if np.all(~np.isfinite(ll_coarse)):
        return out

    best_idx = int(np.nanargmax(ll_coarse))
    best_l = float(coarse[best_idx])

    lo = max(0.0, best_l - 0.05)
    hi = min(1.0, best_l + 0.05)
    fine_local = np.linspace(lo, hi, 401)
    ll_local = np.array([loglik_pagel_lambda(y, C, float(l)) for l in fine_local])

    if np.all(~np.isfinite(ll_local)):
        l_hat = best_l
        ll_hat = float(ll_coarse[best_idx])
    else:
        i2 = int(np.nanargmax(ll_local))
        l_hat = float(fine_local[i2])
        ll_hat = float(ll_local[i2])

    # Profile likelihood CI (chi-square df=1 threshold 3.84).
    grid = np.linspace(0.0, 1.0, 1001)
    ll_grid = np.array([loglik_pagel_lambda(y, C, float(l)) for l in grid])
    finite = np.isfinite(ll_grid)
    if np.any(finite):
        ll_max = float(np.nanmax(ll_grid))
        cutoff = ll_max - 0.5 * 3.841458820694124
        valid = grid[(ll_grid >= cutoff) & finite]
        ci_low = float(np.min(valid)) if valid.size else np.nan
        ci_high = float(np.max(valid)) if valid.size else np.nan
    else:
        ci_low = np.nan
        ci_high = np.nan

    ll0 = loglik_pagel_lambda(y, C, 0.0)
    ll1 = loglik_pagel_lambda(y, C, 1.0)

    # df=1 chi-square survival: p = erfc(sqrt(x/2))
    lr0 = float(2.0 * (ll_hat - ll0)) if np.isfinite(ll0) else np.nan
    p0 = float(math.erfc(math.sqrt(max(lr0, 0.0) / 2.0))) if np.isfinite(lr0) else np.nan

    lr1 = float(2.0 * (ll_hat - ll1)) if np.isfinite(ll1) else np.nan
    p1 = float(math.erfc(math.sqrt(max(lr1, 0.0) / 2.0))) if np.isfinite(lr1) else np.nan

    out.update(
        {
            "lambda_hat": l_hat,
            "ll_hat": ll_hat,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "lr_vs_0": lr0,
            "p_vs_0": p0,
            "lr_vs_1": lr1,
            "p_vs_1": p1,
        }
    )
    return out


def write_methods_snippet(
    outpath: Path,
    n_reps: int,
    alpha: float,
    q: float,
    seed: int,
    background_mode: str,
    swaps_per_replicate: int,
    lambda_trait_null: str,
) -> None:
    text = (
        "Host phylogenetic breadth was quantified for each virus from the set of unique associated hosts "
        "using host patristic distances on a branch-length mammal phylogeny. We computed mean pairwise "
        "distance (MPD), mean nearest taxon distance (MNTD), Faith's phylogenetic diversity (PD), normalized "
        "MRCA depth, and a minimal clade-cover statistic obtained by cutting the host tree at q={q:.3f} of total "
        "root height. Viruses with fewer than two hosts were retained but MPD/MNTD were set to NA. "
        "Null expectations for MPD and MNTD were generated with {n_reps} replicates under three models: "
        "(Null0) uniform sampling of k hosts from background B ({background_mode}); "
        "(Null2) degree-weighted sampling without replacement with probabilities proportional to host degree^alpha "
        "(alpha={alpha:.3f}); and "
        "(Null3) degree-preserving bipartite rewiring by double-edge swaps with {swaps_per_replicate} accepted swaps "
        "per replicate. Standardized effect sizes were computed as SES=(observed-mean(null))/sd(null), with empirical "
        "two-sided and upper-tail p-values from null ranks. Phylogenetic signal of breadth was evaluated on the virus "
        "tree by maximum-likelihood Pagel's lambda under Brownian-motion covariance, using traits SES_MPD and SES_MNTD "
        "from {lambda_trait_null}. Lambda confidence intervals were obtained by profile likelihood (df=1), with LR tests "
        "against lambda=0 and lambda=1. All analyses used random seed={seed}."
    ).format(
        q=q,
        n_reps=n_reps,
        background_mode=background_mode,
        alpha=alpha,
        swaps_per_replicate=swaps_per_replicate,
        lambda_trait_null=lambda_trait_null,
        seed=seed,
    )
    outpath.write_text(text + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    np.random.seed(args.seed)
    rng = np.random.default_rng(args.seed)

    host_tree = load_tree_newick(Path(args.host_tree), label="Host")
    virus_tree = load_tree_newick(Path(args.virus_tree), label="Virus")

    host_tips = tip_labels(host_tree)
    virus_tips = tip_labels(virus_tree)
    host_tip_set = set(host_tips)
    virus_tip_set = set(virus_tips)

    assoc = parse_associations(Path(args.assoc))
    assoc, n_host_norm, n_virus_norm = validate_association_ids(
        assoc, host_tip_set, virus_tip_set, outdir
    )
    if n_host_norm > 0 or n_virus_norm > 0:
        print(
            f"INFO: normalized labels (host={n_host_norm}, virus={n_virus_norm}) "
            "using underscore/space harmonization.",
            file=sys.stderr,
        )

    # Restrict associations to virus tips currently in virus tree (all should already match).
    assoc = assoc[assoc["virus_id"].isin(virus_tip_set) & assoc["host_id"].isin(host_tip_set)].copy()

    # Host family map from associations.
    host_to_family = {}
    if "host_family" in assoc.columns:
        fam_df = assoc[["host_id", "host_family"]].copy()
        fam_df["host_family"] = fam_df["host_family"].fillna("").astype(str)
        fam_df = fam_df[fam_df["host_family"].str.strip() != ""]
        if not fam_df.empty:
            fam_mode = (
                fam_df.groupby("host_id")["host_family"]
                .agg(lambda s: s.value_counts().index[0])
                .to_dict()
            )
            host_to_family.update(fam_mode)

    # Virus -> hosts and weighted richness.
    virus_to_hosts: Dict[str, Set[str]] = defaultdict(set)
    virus_to_weighted_k: Dict[str, float] = defaultdict(float)

    for row in assoc.itertuples(index=False):
        virus_to_hosts[row.virus_id].add(row.host_id)

    weights_per_vh = assoc.groupby(["virus_id", "host_id"], as_index=False)["evidence_weight"].max()
    for row in weights_per_vh.itertuples(index=False):
        virus_to_weighted_k[row.virus_id] += float(row.evidence_weight)

    # Ensure all virus tips are present in table, even if no associations.
    for v in virus_tips:
        virus_to_hosts.setdefault(v, set())
        virus_to_weighted_k.setdefault(v, 0.0)

    # Precompute host patristic distances.
    host_tip_to_node = build_tip_node_map(host_tree)
    host_depths = node_depths(host_tree)
    host_height = tree_height(host_tree, host_depths)
    host_dist, host_index = compute_host_distance_matrix(host_tree, host_tips)

    # Precompute clade assignment at q for minimal clade cover.
    clade_assign = clade_cut_assignment(host_tree, host_tip_to_node, host_depths, q=args.q)

    # Observed metrics per virus.
    rows = []
    virus_to_k = {}
    for v in virus_tips:
        hosts = sorted(virus_to_hosts[v])
        k = len(hosts)
        virus_to_k[v] = k

        if k >= 2:
            idx = [host_index[h] for h in hosts]
            mpd_obs, mntd_obs = compute_mpd_mntd(idx, host_dist)
        else:
            mpd_obs, mntd_obs = np.nan, np.nan

        pd_obs = faith_pd(hosts, host_tip_to_node, host_tree) if k >= 1 else np.nan
        mrca_d = mrca_depth_norm(hosts, host_tree, host_depths, host_height) if k >= 1 else np.nan
        mcc = minimal_clade_cover(hosts, clade_assign) if k >= 1 else np.nan

        n_fam, fam_counts, between_fam_mpd = family_summary(hosts, host_to_family, host_index, host_dist)

        rows.append(
            {
                "virus_id": v,
                "k_v": k,
                "k_weighted": float(virus_to_weighted_k.get(v, 0.0)),
                "mpd_obs": mpd_obs,
                "mntd_obs": mntd_obs,
                "pd_obs": pd_obs,
                "mrca_depth_obs": mrca_d,
                "minimal_clade_cover": mcc,
                "n_host_families": n_fam,
                "host_family_counts": fam_counts,
                "between_family_mpd": between_fam_mpd,
            }
        )

    breadth_df = pd.DataFrame(rows)

    # Null backgrounds and host degrees for Null2.
    B = host_background(
        assoc=assoc,
        host_tip_labels=host_tip_set,
        mode=args.background_mode,
        host_file=Path(args.background_hosts) if args.background_hosts else None,
    )

    host_degree = assoc.groupby("host_id")["virus_id"].nunique().to_dict()
    b_degrees = np.array([float(host_degree.get(h, 0.0)) for h in B], dtype=float)
    weighted = np.power(b_degrees, args.alpha)
    if np.sum(weighted) > 0:
        weighted = weighted / np.sum(weighted)
    else:
        weighted = np.zeros_like(weighted)

    # Null0 / Null2.
    null0_mpd, null0_mntd = simulate_null_sampling(
        virus_to_k=virus_to_k,
        background_hosts=B,
        host_index=host_index,
        host_dist=host_dist,
        rng=rng,
        n_reps=args.n_reps,
        probs=None,
    )
    null2_mpd, null2_mntd = simulate_null_sampling(
        virus_to_k=virus_to_k,
        background_hosts=B,
        host_index=host_index,
        host_dist=host_dist,
        rng=rng,
        n_reps=args.n_reps,
        probs=weighted,
    )

    # Null3 rewiring.
    E = int(sum(len(hs) for hs in virus_to_hosts.values()))
    swaps_per_rep = args.swaps_per_replicate if args.swaps_per_replicate is not None else max(1, 20 * E)

    null3_mpd, null3_mntd, acceptance_rate, swaps_accepted, swaps_attempted = run_degree_preserving_rewiring(
        virus_to_hosts=virus_to_hosts,
        host_index=host_index,
        host_dist=host_dist,
        n_reps=args.n_reps,
        rng=rng,
        swaps_per_replicate=swaps_per_rep,
    )

    # Attach SES and p-values for each null model.
    def add_null_columns(df: pd.DataFrame, null_mpd: Dict[str, List[float]], null_mntd: Dict[str, List[float]], tag: str) -> None:
        out_mpd = []
        out_mntd = []
        for row in df.itertuples(index=False):
            v = row.virus_id
            st_mpd = empirical_stats(row.mpd_obs, null_mpd.get(v, []))
            st_mntd = empirical_stats(row.mntd_obs, null_mntd.get(v, []))
            out_mpd.append(st_mpd)
            out_mntd.append(st_mntd)

        df[f"mpd_null_mean_{tag}"] = [x.mean for x in out_mpd]
        df[f"mpd_null_sd_{tag}"] = [x.sd for x in out_mpd]
        df[f"mpd_ses_{tag}"] = [x.ses for x in out_mpd]
        df[f"mpd_p_two_sided_{tag}"] = [x.p_two_sided for x in out_mpd]
        df[f"mpd_p_high_{tag}"] = [x.p_high for x in out_mpd]
        df[f"mpd_n_null_{tag}"] = [x.n_null for x in out_mpd]

        df[f"mntd_null_mean_{tag}"] = [x.mean for x in out_mntd]
        df[f"mntd_null_sd_{tag}"] = [x.sd for x in out_mntd]
        df[f"mntd_ses_{tag}"] = [x.ses for x in out_mntd]
        df[f"mntd_p_two_sided_{tag}"] = [x.p_two_sided for x in out_mntd]
        df[f"mntd_p_high_{tag}"] = [x.p_high for x in out_mntd]
        df[f"mntd_n_null_{tag}"] = [x.n_null for x in out_mntd]

    add_null_columns(breadth_df, null0_mpd, null0_mntd, "null0")
    add_null_columns(breadth_df, null2_mpd, null2_mntd, "null2")
    add_null_columns(breadth_df, null3_mpd, null3_mntd, "null3")

    # Add run metadata columns.
    breadth_df["seed"] = args.seed
    breadth_df["n_reps"] = args.n_reps
    breadth_df["alpha"] = args.alpha
    breadth_df["q"] = args.q
    breadth_df["background_mode"] = args.background_mode
    breadth_df["null3_swaps_per_replicate"] = swaps_per_rep
    breadth_df["null3_swap_acceptance_rate"] = acceptance_rate
    breadth_df["null3_swaps_accepted_total"] = swaps_accepted
    breadth_df["null3_swaps_attempted_total"] = swaps_attempted

    breadth_out = outdir / "virus_host_breadth.tsv"
    breadth_df.to_csv(breadth_out, sep="\t", index=False)

    # Lambda summary using chosen SES columns.
    ses_mpd_col = f"mpd_ses_{args.lambda_trait_null}"
    ses_mntd_col = f"mntd_ses_{args.lambda_trait_null}"

    virus_tip_to_node = build_tip_node_map(virus_tree)
    virus_depths = node_depths(virus_tree)

    lambda_rows = []
    for trait_col, trait_name in [
        (ses_mpd_col, "SES_MPD"),
        (ses_mntd_col, "SES_MNTD"),
        ("mrca_depth_obs", "MRCA_depth_obs"),
    ]:
        sub = breadth_df[["virus_id", trait_col]].copy()
        sub = sub[np.isfinite(sub[trait_col])]
        sub = sub[sub["virus_id"].isin(virus_tips)]

        if len(sub) < 3:
            lambda_rows.append(
                {
                    "trait": trait_name,
                    "source_column": trait_col,
                    "n_taxa": len(sub),
                    "lambda_hat": np.nan,
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "lr_vs_0": np.nan,
                    "p_vs_0": np.nan,
                    "lr_vs_1": np.nan,
                    "p_vs_1": np.nan,
                }
            )
            continue

        tip_order = sorted(sub["virus_id"].tolist())
        y = sub.set_index("virus_id").loc[tip_order, trait_col].to_numpy(dtype=float)

        C = bm_covariance(virus_tree, tip_order, virus_tip_to_node, virus_depths)
        est = estimate_lambda(y, C)

        lambda_rows.append(
            {
                "trait": trait_name,
                "source_column": trait_col,
                "n_taxa": int(est["n_taxa"]),
                "lambda_hat": est["lambda_hat"],
                "ci_low": est["ci_low"],
                "ci_high": est["ci_high"],
                "lr_vs_0": est["lr_vs_0"],
                "p_vs_0": est["p_vs_0"],
                "lr_vs_1": est["lr_vs_1"],
                "p_vs_1": est["p_vs_1"],
            }
        )

    lambda_df = pd.DataFrame(lambda_rows)
    lambda_out = outdir / "lambda_summary.tsv"
    lambda_df.to_csv(lambda_out, sep="\t", index=False)

    # Plots.
    plot_df = breadth_df.copy()

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    v1 = plot_df[ses_mpd_col].to_numpy(dtype=float)
    v1 = v1[np.isfinite(v1)]
    v2 = plot_df[ses_mntd_col].to_numpy(dtype=float)
    v2 = v2[np.isfinite(v2)]

    axes[0].hist(v1, bins=30, color="#4c72b0", alpha=0.8, edgecolor="white")
    axes[0].set_title("SES MPD")
    axes[0].set_xlabel(ses_mpd_col)
    axes[0].set_ylabel("Count")

    axes[1].hist(v2, bins=30, color="#dd8452", alpha=0.8, edgecolor="white")
    axes[1].set_title("SES MNTD")
    axes[1].set_xlabel(ses_mntd_col)
    axes[1].set_ylabel("Count")

    fig.savefig(outdir / "ses_histograms.png", dpi=300)
    plt.close(fig)

    fig2, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)
    x = plot_df["k_v"].to_numpy(dtype=float)
    y = plot_df[ses_mpd_col].to_numpy(dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    ax.scatter(x[mask], y[mask], s=20, alpha=0.8, color="#4c72b0")
    ax.set_xlabel("k_v (number of associated hosts)")
    ax.set_ylabel(ses_mpd_col)
    ax.set_title("SES MPD vs host count")
    fig2.savefig(outdir / "ses_mpd_vs_kv.png", dpi=300)
    plt.close(fig2)

    # Methods snippet.
    write_methods_snippet(
        outpath=outdir / "methods_snippet.txt",
        n_reps=args.n_reps,
        alpha=args.alpha,
        q=args.q,
        seed=args.seed,
        background_mode=args.background_mode,
        swaps_per_replicate=swaps_per_rep,
        lambda_trait_null=args.lambda_trait_null,
    )

    # Run metadata.
    metadata = {
        "host_tree": str(Path(args.host_tree).resolve()),
        "virus_tree": str(Path(args.virus_tree).resolve()),
        "assoc": str(Path(args.assoc).resolve()),
        "outdir": str(outdir.resolve()),
        "n_reps": args.n_reps,
        "alpha": args.alpha,
        "q": args.q,
        "seed": args.seed,
        "background_mode": args.background_mode,
        "background_hosts": str(Path(args.background_hosts).resolve()) if args.background_hosts else "",
        "swaps_per_replicate": swaps_per_rep,
        "null3_acceptance_rate": acceptance_rate,
        "null3_swaps_accepted_total": swaps_accepted,
        "null3_swaps_attempted_total": swaps_attempted,
        "lambda_trait_null": args.lambda_trait_null,
        "n_hosts_in_host_tree": len(host_tips),
        "n_viruses_in_virus_tree": len(virus_tips),
        "n_assoc_edges_unique": len(assoc),
    }
    pd.DataFrame([metadata]).to_csv(outdir / "run_metadata.tsv", sep="\t", index=False)

    print(f"Wrote: {breadth_out}")
    print(f"Wrote: {lambda_out}")
    print(f"Wrote: {outdir / 'ses_histograms.png'}")
    print(f"Wrote: {outdir / 'ses_mpd_vs_kv.png'}")
    print(f"Wrote: {outdir / 'methods_snippet.txt'}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
