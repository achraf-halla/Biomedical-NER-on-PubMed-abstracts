"""Hand-curated normalisation rules for scispaCy entity output.

The raw scispaCy output has three recurring problems on this corpus:

1. **Fragment artifacts.** scispaCy sometimes emits "tyrosine" as a
   standalone CHEMICAL when the actual phrase is "tyrosine kinase
   inhibitor". This produces a single artificial hub that absorbs most
   of the network's centrality.
2. **Acronym/expansion duplicates.** "NSCLC" and "non-small cell lung
   cancer" are the same disease but appear as separate nodes with
   thousands of edges between them.
3. **Generic descriptors classified as entities.** "toxicity", "death",
   "antitumor", "cellular", etc. — not entities, just adjectives or
   weakly-defined nouns that scispaCy's models surface anyway.

This module addresses (2) and (3) and partially (1). It is **not** a
replacement for proper UMLS linking via `scispacy.linking` — it is a
pragmatic patch that materially improves the network without that
dependency.

Limitations called out honestly in the README:
- Coverage is limited to the high-frequency cases I observed.
- Synonyms across disease subtypes (e.g. "DLBCL" vs "diffuse large
  B-cell lymphoma" vs "B-cell lymphoma") are merged conservatively;
  some legitimate distinctions are lost.
- Drug-class names (e.g. "TKIs", "EGFR-TKIs") are merged to a single
  "tyrosine kinase inhibitor" node, which is appropriate for class-
  level analysis but loses specificity.
"""

from __future__ import annotations

import pandas as pd

# ---------------------------------------------------------------------------
# Stop-list (extends the one in src/ner.py for the post-extraction pass)
# ---------------------------------------------------------------------------

EXTENDED_STOPLIST = {
    # Generic descriptors
    "toxicity", "toxicities", "death", "deaths", "antitumor", "anticancer",
    "malignancies", "malignancy",

    # Cell-related fragments
    "cellular", "line", "lines", "cell line", "cell lines", "single-cell",
    "tumor cell", "tumor cells", "cancer cell", "cancer cells", "cell",
    "immune cell", "tumour cells", "tumour cell",

    # Trial/clinical artifacts that scispaCy surfaces as entities
    "orr", "ici", "icis", "second-line", "first-line", "third-line",
    "frontline", "front-line",

    # Treatment / regimen non-entities
    "treatment", "treatments", "regimen", "regimens", "therapy", "therapies",
    "chemotherapy", "monotherapy", "combination therapy", "maintenance",

    # Generic class words
    "compound", "compounds", "drug", "drugs", "agent", "agents", "molecule",
    "molecules", "inhibitor", "inhibitors", "antibody", "antibodies",

    # Outcome / measurement words
    "response", "responses", "survival", "outcome", "outcomes", "efficacy",
    "expression", "growth",

    # Anatomical / organism non-entities surfaced by ANATOMY mapping
    "organ", "organs",

    # The big offender — "tyrosine" alone is almost always a fragment
    # of "tyrosine kinase" / "tyrosine kinase inhibitor"
    "tyrosine",

    # Other fragment artifacts
    "kinase", "kinases", "receptor", "receptors", "growth factor",
    "epidermal growth factor", "growth factor receptor",
}

# ---------------------------------------------------------------------------
# Alias map: noisy scispaCy output → canonical form
#
# Format: alias (lowercased, normalised) → canonical name
#
# After applying this, many entities collapse to the same canonical key.
# We then re-aggregate at the (canonical_norm, type) level.
# ---------------------------------------------------------------------------

ALIASES: dict[str, str] = {
    # ---- Drug classes / fragments → "tyrosine kinase inhibitor" ----
    "tki": "tyrosine kinase inhibitor",
    "tkis": "tyrosine kinase inhibitor",
    "egfr-tki": "egfr tyrosine kinase inhibitor",
    "egfr-tkis": "egfr tyrosine kinase inhibitor",
    "egfr tki": "egfr tyrosine kinase inhibitor",
    "egfr tkis": "egfr tyrosine kinase inhibitor",
    "egfr-tyrosine kinase inhibitor": "egfr tyrosine kinase inhibitor",
    "egfr-tyrosine": "egfr tyrosine kinase inhibitor",
    "multi-kinase": "multi-kinase inhibitor",
    "multikinase": "multi-kinase inhibitor",
    "btki": "btk inhibitor",
    "btkis": "btk inhibitor",
    "bruton tyrosine kinase inhibitor": "btk inhibitor",
    "bruton's tyrosine kinase inhibitor": "btk inhibitor",
    "cdk4/6": "cdk4/6 inhibitor",
    "cdk4/6i": "cdk4/6 inhibitor",
    "ckis": "checkpoint inhibitor",
    "checkpoint inhibitors": "checkpoint inhibitor",

    # ---- Genes/proteins: acronym ↔ expansion ----
    "epidermal growth factor receptor": "egfr",
    "epidermal growth factor receptor tyrosine kinase": "egfr",
    "egfr tyrosine kinase": "egfr",
    "bruton tyrosine kinase": "btk",
    "bruton's tyrosine kinase": "btk",
    "anaplastic lymphoma kinase": "alk",
    "human epidermal growth factor receptor 2": "her2",
    "erbb2": "her2",
    "her-2": "her2",
    "vascular endothelial growth factor receptor": "vegfr",
    "platelet-derived growth factor receptor": "pdgfr",
    "fibroblast growth factor receptor": "fgfr",
    "janus kinase": "jak",
    "mitogen-activated protein kinase": "mapk",
    "extracellular signal-regulated kinase": "erk",
    "phosphoinositide 3-kinase": "pi3k",
    "mammalian target of rapamycin": "mtor",
    "abelson tyrosine kinase": "abl",
    "abl1": "abl",
    "breakpoint cluster region": "bcr",
    "bcr-abl1": "bcr-abl",

    # ---- Diseases: acronym ↔ expansion ----
    "non-small cell lung cancer": "nsclc",
    "non small cell lung cancer": "nsclc",
    "non-small-cell lung cancer": "nsclc",
    "non-small cell lung carcinoma": "nsclc",
    "small cell lung cancer": "sclc",
    "lung adenocarcinoma": "luad",
    "hepatocellular carcinoma": "hcc",
    "renal cell carcinoma": "rcc",
    "clear cell renal cell carcinoma": "ccrcc",
    "cell renal cell carcinoma": "ccrcc",
    "chronic myeloid leukemia": "cml",
    "chronic myelogenous leukemia": "cml",
    "acute myeloid leukemia": "aml",
    "acute myelogenous leukemia": "aml",
    "acute lymphoblastic leukemia": "all",
    "chronic lymphocytic leukemia": "cll",
    "chronic lymphocytic leukemia/small lymphocytic lymphoma": "cll",
    "cll/sll": "cll",
    "small lymphocytic lymphoma": "sll",
    "diffuse large b-cell lymphoma": "dlbcl",
    "diffuse large b cell lymphoma": "dlbcl",
    "b-cell lymphoma": "dlbcl",  # imperfect but most common usage
    "mantle cell lymphoma": "mcl",
    "follicular lymphoma": "fl",
    "hodgkin lymphoma": "hl",
    "hodgkin's lymphoma": "hl",
    "non-hodgkin lymphoma": "nhl",
    "non-hodgkin's lymphoma": "nhl",
    "multiple myeloma": "mm",
    "myelodysplastic syndrome": "mds",
    "myelodysplastic syndromes": "mds",
    "pancreatic ductal adenocarcinoma": "pdac",
    "pancreatic cancer": "pdac",   # most pancreatic cancer in this lit = PDAC
    "gastrointestinal stromal tumor": "gist",
    "gastrointestinal stromal tumors": "gist",
    "gists": "gist",
    "interstitial lung disease": "ild",
    "central nervous system": "cns",
    "myeloproliferative neoplasm": "mpn",
    "myeloproliferative neoplasms": "mpn",

    # ---- Drugs: brand → generic / spelling variants ----
    "imatinib mesylate": "imatinib",
    "gleevec": "imatinib",
    "glivec": "imatinib",
    "sprycel": "dasatinib",
    "tasigna": "nilotinib",
    "tagrisso": "osimertinib",
    "iressa": "gefitinib",
    "tarceva": "erlotinib",
    "ibrance": "palbociclib",
    "kisqali": "ribociclib",
    "verzenio": "abemaciclib",
    "imbruvica": "ibrutinib",
    "calquence": "acalabrutinib",
    "brukinsa": "zanubrutinib",
    "xalkori": "crizotinib",
    "alecensa": "alectinib",
    "lorbrena": "lorlatinib",
    "nexavar": "sorafenib",
    "lenvima": "lenvatinib",
    "cabometyx": "cabozantinib",
    "cometriq": "cabozantinib",
    "stivarga": "regorafenib",
    "jakafi": "ruxolitinib",
    "jakavi": "ruxolitinib",

    # ---- Cell lines: collapse trivial duplicates ----
    "a549 cells": "a549",
    "hcc cells": "hcc cells",  # keep — different from disease "HCC"
}


# ---------------------------------------------------------------------------
# Application
# ---------------------------------------------------------------------------

def canonicalise(norm: str) -> str:
    """Apply alias map; return canonical form (or input unchanged)."""
    return ALIASES.get(norm, norm)


def is_blocked(norm: str) -> bool:
    """Return True if the normalised text should be dropped."""
    return norm in EXTENDED_STOPLIST


def apply_normalisation(entities: pd.DataFrame) -> pd.DataFrame:
    """Filter + canonicalise a DataFrame of entity occurrences.

    Returns a new DataFrame with:
    - Stop-listed entities removed.
    - Aliased entities collapsed to canonical names (norm overwritten).
    - Original raw text preserved in `text_raw` for traceability.

    Resulting columns: [pmid, model, type, raw_type, text, text_raw, norm, start, end].
    """
    df = entities.copy()

    # Drop stop-listed entries (post-NER stop-list, more aggressive than ner.py's)
    df = df.loc[~df["norm"].isin(EXTENDED_STOPLIST)].copy()

    # Canonicalise — preserve original text for inspection
    df["text_raw"] = df["text"]
    df["norm"] = df["norm"].map(canonicalise)

    # Re-apply stop-list in case a canonical form is on it (defensive)
    df = df.loc[~df["norm"].isin(EXTENDED_STOPLIST)].copy()

    return df.reset_index(drop=True)


def normalisation_summary(before: pd.DataFrame, after: pd.DataFrame) -> dict:
    """Numbers for the README: how much did normalisation prune?"""
    def n_distinct(d):
        return d.groupby(["norm", "type"]).ngroups
    return {
        "occurrences_before": len(before),
        "occurrences_after": len(after),
        "occurrences_dropped": len(before) - len(after),
        "distinct_before": n_distinct(before),
        "distinct_after": n_distinct(after),
        "distinct_collapsed": n_distinct(before) - n_distinct(after),
    }
