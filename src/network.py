"""Build co-occurrence networks from extracted entities.

Two-pass approach:

1. **Aggregation.** Group entities by (pmid). For each abstract, list
   the unique normalised entities present. Filter to entities meeting
   a minimum global frequency (`min_count`) — rare entities produce
   noisy edges with no statistical support.
2. **Edge weighting.** For each pair of entities co-occurring in an
   abstract, increment the edge weight. Optionally compute **PMI**
   (pointwise mutual information) to surface associations that are
   stronger than chance — raw co-occurrence is dominated by hub
   entities (e.g. "cancer" edges with everything).

We provide both raw count and PMI views; the reader chooses what
matters for the question they're asking.
"""

from __future__ import annotations

import json
import math
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from typing import Iterable

import networkx as nx
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROCESSED_DIR = ROOT / "data" / "processed"


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def load_entities(path: Path | None = None) -> pd.DataFrame:
    path = path or (PROCESSED_DIR / "entities.jsonl")
    if not path.exists():
        raise FileNotFoundError(
            f"No entities at {path}. Run: python -m src.ner"
        )
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
    return pd.DataFrame(rows)


def entity_doc_frequency(entities: pd.DataFrame) -> pd.DataFrame:
    """For each (norm, type), count distinct PMIDs it appears in (DF).

    We use document frequency rather than raw mention frequency because
    one paper mentioning EGFR fifteen times shouldn't outweigh fifteen
    papers mentioning EGFR once.
    """
    df = (entities
          .drop_duplicates(["pmid", "norm", "type"])
          .groupby(["norm", "type"], as_index=False)
          .agg(doc_freq=("pmid", "size"))
          .sort_values("doc_freq", ascending=False))
    return df


def filter_entities(
    entities: pd.DataFrame,
    min_count: int = 5,
    keep_types: list[str] | None = None,
) -> pd.DataFrame:
    """Filter to entities with DF >= `min_count` and (optionally) a type subset."""
    df_freq = entity_doc_frequency(entities)
    keep = df_freq.loc[df_freq["doc_freq"] >= min_count, ["norm", "type"]]
    if keep_types is not None:
        keep = keep.loc[keep["type"].isin(keep_types)]
    keep_set = set(map(tuple, keep[["norm", "type"]].itertuples(index=False, name=None)))

    mask = entities.apply(
        lambda r: (r["norm"], r["type"]) in keep_set, axis=1
    )
    return entities.loc[mask].copy()


# ---------------------------------------------------------------------------
# Co-occurrence
# ---------------------------------------------------------------------------

def cooccurrence_edges(
    entities: pd.DataFrame,
    n_docs_total: int | None = None,
) -> pd.DataFrame:
    """Build edge list from a filtered entities frame.

    Returns a DataFrame with columns:
        a, b, type_a, type_b, weight, pmi

    `weight` is the document co-occurrence count.
    `pmi` is PMI based on global document frequencies, treating
    each abstract as a "document".
    """
    # Keys = (norm, type) tuples — we can have an entity appear with
    # different types from different models, treat them as distinct.
    by_pmid: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for _, r in entities[["pmid", "norm", "type"]].iterrows():
        by_pmid[r["pmid"]].add((r["norm"], r["type"]))

    n_docs = n_docs_total or len(by_pmid)

    # Per-entity document frequency
    df_counts: Counter = Counter()
    for ents in by_pmid.values():
        for e in ents:
            df_counts[e] += 1

    # Pair counts
    pair_counts: Counter = Counter()
    for ents in by_pmid.values():
        # Sort to make pair tuples canonical
        sorted_ents = sorted(ents)
        for a, b in combinations(sorted_ents, 2):
            pair_counts[(a, b)] += 1

    rows = []
    for (a, b), w in pair_counts.items():
        # PMI = log2(P(a, b) / (P(a) * P(b)))
        # using P(x) = df(x) / n_docs and P(a, b) = w / n_docs
        p_a = df_counts[a] / n_docs
        p_b = df_counts[b] / n_docs
        p_ab = w / n_docs
        pmi = math.log2(p_ab / (p_a * p_b)) if p_a * p_b > 0 else 0.0
        rows.append({
            "a": a[0], "type_a": a[1],
            "b": b[0], "type_b": b[1],
            "weight": w, "pmi": pmi,
        })

    return (pd.DataFrame(rows)
            .sort_values("weight", ascending=False)
            .reset_index(drop=True))


def build_graph(
    edges: pd.DataFrame,
    min_weight: int = 3,
    weight_col: str = "weight",
) -> nx.Graph:
    """Construct a NetworkX graph from a filtered edge list."""
    sub = edges.loc[edges["weight"] >= min_weight].copy()
    G = nx.Graph()
    for _, r in sub.iterrows():
        for node, ntype in [(r["a"], r["type_a"]), (r["b"], r["type_b"])]:
            if node not in G:
                G.add_node(node, type=ntype)
        G.add_edge(r["a"], r["b"],
                   weight=int(r["weight"]),
                   pmi=float(r["pmi"]))
    return G


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

def centrality_table(G: nx.Graph) -> pd.DataFrame:
    """Per-node degree, weighted degree, betweenness — for ranking hubs."""
    deg = dict(G.degree())
    wdeg = dict(G.degree(weight="weight"))
    btw = nx.betweenness_centrality(G, weight=None, normalized=True)
    rows = []
    for n, attrs in G.nodes(data=True):
        rows.append({
            "node": n,
            "type": attrs.get("type"),
            "degree": deg[n],
            "weighted_degree": wdeg[n],
            "betweenness": btw[n],
        })
    return (pd.DataFrame(rows)
            .sort_values("weighted_degree", ascending=False)
            .reset_index(drop=True))


def communities(G: nx.Graph, resolution: float = 1.0) -> dict[str, int]:
    """Greedy modularity communities — returns {node: community_id}."""
    # Louvain is in networkx as of 3.0+ via community.louvain_communities,
    # but greedy_modularity_communities is always available and
    # adequate for a few hundred nodes.
    comms = nx.community.greedy_modularity_communities(
        G, weight="weight", resolution=resolution
    )
    out = {}
    for cid, comm in enumerate(comms):
        for node in comm:
            out[node] = cid
    return out
