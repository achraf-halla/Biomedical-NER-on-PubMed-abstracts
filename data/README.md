# Data

Raw PubMed records and parsed abstracts are **not committed** to this
repo. They are downloaded on demand from NCBI Entrez.

## Fetch

```bash
# Default: ~1,500 abstracts on "kinase inhibitor cancer"
python -m src.ingest --email you@example.com

# Then parse the XML into one JSON record per abstract:
python -m src.parse
```

This produces:

- `data/raw/pubmed_records.xml` — raw Entrez efetch output
- `data/raw/abstracts.jsonl` — one JSON record per abstract

## Source

PubMed via NCBI Entrez E-utilities. NCBI's
[usage policy](https://www.ncbi.nlm.nih.gov/books/NBK25497/) requires:

- Identifying the tool and email contact on each request (handled by
  `src/ingest.py`).
- Staying ≤ 3 requests/second without an API key (we use 2 req/s).
- Running large jobs off-peak when possible.

PubMed abstracts may incorporate material protected by copyright. This
repo's pipeline downloads abstracts for analysis; the raw text is not
redistributed here.
