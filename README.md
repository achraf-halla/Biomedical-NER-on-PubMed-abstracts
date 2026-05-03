# Biomedical Entity Network

Extract biomedical entities (drugs, diseases, genes/proteins, cell types)
from PubMed abstracts and analyze their **co-occurrence patterns** to
surface implicit knowledge — which drugs are mentioned together with
which targets, which diseases cluster around which genes, where
research attention is concentrated.

## Why this project

A literature search returns a list of papers. A literature *map* shows
the relationships between entities across thousands of papers at once.
For a researcher scoping a new area, the second is much more useful.

The technical interest is in the unglamorous parts of the pipeline:

1. **Entity extraction at scale** — using domain-specialised NER
   (scispaCy) rather than a generic large language model. Faster,
   cheaper, and the "right tool for the job" in biomedical NLP.
2. **Co-occurrence as a signal** — pairs of entities that appear in the
   same abstract are weak evidence of association. Aggregated over
   thousands of abstracts, the signal becomes informative — and the
   noise structure (degree distribution, central nodes) is interesting
   in its own right.
3. **Honest entity reconciliation** — the same drug appears as
   "trastuzumab", "Herceptin", "anti-HER2 antibody" across abstracts.
   We deal with this explicitly rather than pretending it doesn't
   matter.

## Default corpus

`"kinase inhibitor"[Title/Abstract] AND cancer` — ~1,500 abstracts.
Broad enough to capture multiple drug classes (TKIs, CDK inhibitors,
Aurora kinase inhibitors, etc.), targets (EGFR, ALK, BRAF, ...), and
cancer types. Configurable via CLI flag — see "Custom queries" below.

## Repo layout

```
biomedical-ner-network/
├── data/
│   ├── raw/              # downloaded XML / parsed JSONL (gitignored)
│   └── processed/        # entity annotations + co-occurrence (gitignored)
├── notebooks/
│   ├── 01_corpus_eda.ipynb            # corpus structure & coverage
│   ├── 02_entities_and_network.ipynb  # NER + network analysis (raw)
│   └── 02b_normalisation.ipynb        # before/after normalisation
├── src/
│   ├── ingest.py         # Entrez fetch (esearch + efetch via Biopython)
│   ├── parse.py          # PubMed XML → structured records
│   ├── ner.py            # scispaCy NER pipeline (CLI: python -m src.ner)
│   ├── normalize.py      # alias map + extended stop-list
│   ├── network.py        # co-occurrence + graph utilities
│   └── load.py           # convenience loader
├── reports/figures/
├── requirements.txt
└── README.md
```

## Setup

```bash
python -m venv .venv
source .venv/bin/activate         # Windows: .venv\Scripts\activate
pip install -r requirements.txt

# scispaCy models (separate install — wheels hosted at AllenAI):
pip install https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_bc5cdr_md-0.5.4.tar.gz
pip install https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_bionlp13cg_md-0.5.4.tar.gz

# Pipeline:
python -m src.ingest --email you@example.com  # ~30 sec for 1,500 abstracts
python -m src.parse                            # XML → JSONL
python -m src.ner                              # ~30–60 sec on CPU
jupyter lab notebooks/02_entities_and_network.ipynb
```

NCBI requires (politely) that you supply a contact email. No API key
needed — we stay under the 3 req/s rate limit.

## Custom queries

```bash
# Smaller, more focused
python -m src.ingest --email you@example.com \
    --query '"EGFR mutation"[Title/Abstract] AND lung cancer' --max 800

# A different research area entirely
python -m src.ingest --email you@example.com \
    --query '"CRISPR" AND "genome editing"' --max 1500

# More abstracts
python -m src.ingest --email you@example.com --max 3000
```

The full PubMed query syntax is [documented here](https://pubmed.ncbi.nlm.nih.gov/help/#search-tags).

## Roadmap

- [x] Entrez ingestion (esearch + efetch with batched history)
- [x] PubMed XML parser → structured JSONL (PMID, title, abstract, year, journal, MeSH)
- [x] Corpus EDA (year distribution, journal frequency, abstract length, MeSH overview)
- [x] NER pipeline using scispaCy (Disease + Chemical + Gene/Protein + Cell)
- [x] Entity normalisation — alias map + extended stop-list (notebook 02b: before/after analysis)
- [x] Co-occurrence matrix at the abstract level
- [x] Network construction (NetworkX) + interactive HTML (pyvis)
- [x] Centrality, community detection, type-specific sub-networks
- [ ] UMLS entity linking (left as honest limitation — see notebook 02b)

## Approach

The pipeline runs in three steps, each producing a stable artifact:

1. **Ingest + Parse.** Entrez fetches PubMed XML; the parser strips it
   to one JSON record per abstract (PMID, text, year, journal, MeSH).
2. **NER.** Two scispaCy models run over each abstract:
   `en_ner_bc5cdr_md` (Disease + Chemical) and
   `en_ner_bionlp13cg_md` (Gene/Protein + Cell + Anatomy + Cancer +
   Pathology). Their 18 raw labels collapse into 5 canonical buckets.
   Light normalisation (case-fold, strip parentheticals) plus a
   stop-list of weak generic terms.
3. **Normalise.** A second pass via `src/normalize.py` applies a
   hand-curated alias map (acronym ↔ expansion: NSCLC ↔ non-small cell
   lung cancer, EGFR ↔ epidermal growth factor receptor, brand names →
   generic, etc.) and an extended stop-list to remove descriptors
   ("toxicity", "death") and fragment artifacts (standalone
   "tyrosine"). Notebook 02b shows the before/after impact explicitly.
4. **Network.** Document-frequency filter on entities, pairwise
   co-occurrence per abstract, edge-weight filter on the graph. Two
   edge weights kept side-by-side: raw count (popularity) and PMI
   (specificity above chance). NetworkX for analysis, pyvis for
   interactive HTML.

## License

MIT — see [LICENSE](LICENSE).

PubMed abstracts may include copyrighted material; we follow NCBI's
[reuse guidance](https://www.ncbi.nlm.nih.gov/home/about/policies/) and
do not redistribute the raw text in this repo. The pipeline downloads
on demand.
