"""Run scispaCy NER over the corpus and emit a flat entity table.

Two specialised models are used and their outputs merged:

* `en_ner_bc5cdr_md` — DISEASE, CHEMICAL
* `en_ner_bionlp13cg_md` — GENE_OR_GENE_PRODUCT, CELL, ORGAN, TISSUE,
  CANCER, SIMPLE_CHEMICAL, AMINO_ACID, ORGANISM, ORGANISM_SUBSTANCE,
  PATHOLOGICAL_FORMATION, ANATOMICAL_SYSTEM, IMMATERIAL_ANATOMICAL_ENTITY,
  CELLULAR_COMPONENT, MULTI-TISSUE_STRUCTURE, ORGANISM_SUBDIVISION,
  DEVELOPING_ANATOMICAL_STRUCTURE

We collapse this into a small, useful set of canonical types — see
TYPE_MAP below.

Output: `data/processed/entities.jsonl`. One JSON object per detected
entity occurrence:

    {"pmid": "12345", "model": "bc5cdr", "type": "CHEMICAL",
     "text": "imatinib", "norm": "imatinib", "start": 412, "end": 420}

CLI:

    python -m src.ner                          # uses both models
    python -m src.ner --models bc5cdr          # only one
    python -m src.ner --batch-size 64
    python -m src.ner --limit 100              # for testing
"""

from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path
from typing import Iterable, Iterator

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "data" / "raw"
PROCESSED_DIR = ROOT / "data" / "processed"

# scispaCy model name → canonical entity-type mapping. We collapse the
# 18 raw scispaCy labels into 5 useful buckets for downstream analysis.
TYPE_MAP = {
    # bc5cdr labels
    "DISEASE": "DISEASE",
    "CHEMICAL": "CHEMICAL",
    # bionlp13cg labels we keep
    "GENE_OR_GENE_PRODUCT": "GENE",
    "SIMPLE_CHEMICAL": "CHEMICAL",
    "AMINO_ACID": "GENE",                # treat as gene/protein-adjacent
    "CANCER": "DISEASE",                  # bionlp13cg has its own cancer label
    "CELL": "CELL",
    "TISSUE": "ANATOMY",
    "ORGAN": "ANATOMY",
    "ANATOMICAL_SYSTEM": "ANATOMY",
    "MULTI-TISSUE_STRUCTURE": "ANATOMY",
    "PATHOLOGICAL_FORMATION": "DISEASE",
    # We deliberately drop several bionlp13cg labels that add noise
    # (ORGANISM, ORGANISM_SUBSTANCE, ORGANISM_SUBDIVISION,
    #  CELLULAR_COMPONENT, IMMATERIAL_ANATOMICAL_ENTITY,
    #  DEVELOPING_ANATOMICAL_STRUCTURE)
}

MODELS = {
    "bc5cdr": "en_ner_bc5cdr_md",
    "bionlp13cg": "en_ner_bionlp13cg_md",
}


# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------

# Phrases we never want to count as entities, even if scispaCy emits them.
# These are weak / generic terms that turn into spurious hubs in the
# co-occurrence network.
STOPLIST = {
    "patient", "patients", "treatment", "therapy", "study", "studies",
    "disease", "diseases", "drug", "drugs", "tumor", "tumors", "tumour",
    "tumours", "cancer", "cancers", "cell", "cells", "tissue", "tissues",
    "gene", "genes", "protein", "proteins", "kinase", "kinases",
    "inhibitor", "inhibitors", "compound", "compounds", "molecule",
    "molecules", "agent", "agents", "growth", "factor", "factors",
}


def normalise(text: str) -> str:
    """Light normalisation suitable for downstream string matching.

    - Lowercase
    - Strip surrounding punctuation/whitespace
    - Collapse internal whitespace
    - Drop trailing parenthetical content like "EGFR (gene)"

    Deliberately conservative: we do NOT lemmatise, stem, or do
    biomedical-specific synonym mapping. A real production system
    would link entities to UMLS/MeSH; we leave that as an honest
    limitation in the README.
    """
    s = text.strip().lower()
    s = re.sub(r"\s*\([^)]*\)\s*$", "", s)        # drop trailing "(gene)"
    s = re.sub(r"^[\s\-\u2013\u2014.,;:'\"]+", "", s)
    s = re.sub(r"[\s\-\u2013\u2014.,;:'\"]+$", "", s)
    s = re.sub(r"\s+", " ", s)
    return s


def is_useful(norm: str) -> bool:
    """Filter rules applied after normalisation."""
    if not norm:
        return False
    if len(norm) < 2:
        return False
    if norm in STOPLIST:
        return False
    if norm.isdigit():
        return False
    return True


# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------

def load_models(model_keys: list[str]) -> dict[str, "spacy.Language"]:
    import spacy
    out = {}
    for k in model_keys:
        if k not in MODELS:
            raise ValueError(f"Unknown model {k!r}. Available: {sorted(MODELS)}")
        print(f"Loading {MODELS[k]} ...")
        out[k] = spacy.load(MODELS[k])
    return out


def extract_from_records(
    records: Iterable[dict],
    nlps: dict[str, "spacy.Language"],
    batch_size: int = 64,
) -> Iterator[dict]:
    """Yield one entity dict per detected occurrence across all models."""
    records = list(records)
    n = len(records)
    print(f"Running {len(nlps)} model(s) over {n:,} abstracts...")

    for model_key, nlp in nlps.items():
        t0 = time.time()
        # We process title + abstract together, separated by a period to
        # ensure scispaCy treats them as separate sentences.
        texts = [(r["pmid"], (r.get("title", "") + ". " + r["abstract"]).strip())
                 for r in records]
        pmids = [t[0] for t in texts]
        docs = nlp.pipe([t[1] for t in texts], batch_size=batch_size)

        ent_count = 0
        for pmid, doc in zip(pmids, docs):
            for ent in doc.ents:
                raw_type = ent.label_
                ctype = TYPE_MAP.get(raw_type)
                if ctype is None:
                    continue
                norm = normalise(ent.text)
                if not is_useful(norm):
                    continue
                yield {
                    "pmid": pmid,
                    "model": model_key,
                    "type": ctype,
                    "raw_type": raw_type,
                    "text": ent.text,
                    "norm": norm,
                    "start": ent.start_char,
                    "end": ent.end_char,
                }
                ent_count += 1

        elapsed = time.time() - t0
        print(f"  {model_key}: {ent_count:,} entities in {elapsed:.1f}s "
              f"({n / max(elapsed, 0.01):.0f} abstracts/sec)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-path", default=None,
                        help="Input JSONL (default data/raw/abstracts.jsonl)")
    parser.add_argument("--out-path", default=None,
                        help="Output JSONL (default data/processed/entities.jsonl)")
    parser.add_argument("--models", nargs="+", default=["bc5cdr", "bionlp13cg"],
                        help="Which scispaCy models to run")
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--limit", type=int, default=None,
                        help="Process only the first N abstracts (for testing)")
    args = parser.parse_args()

    in_path = Path(args.in_path) if args.in_path else RAW_DIR / "abstracts.jsonl"
    out_path = Path(args.out_path) if args.out_path else PROCESSED_DIR / "entities.jsonl"

    if not in_path.exists():
        raise SystemExit(
            f"Input not found: {in_path}\n"
            f"Run: python -m src.ingest --email you@example.com && python -m src.parse"
        )

    records = [json.loads(line) for line in in_path.read_text(encoding="utf-8").splitlines()]
    if args.limit:
        records = records[:args.limit]
        print(f"Processing first {len(records):,} abstracts (--limit)")

    nlps = load_models(args.models)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_entities = 0
    with out_path.open("w", encoding="utf-8") as f:
        for ent in extract_from_records(records, nlps, batch_size=args.batch_size):
            f.write(json.dumps(ent, ensure_ascii=False) + "\n")
            n_entities += 1

    print(f"\nWrote {n_entities:,} entity occurrences → {out_path}")


if __name__ == "__main__":
    main()
