"""Download PubMed abstracts via NCBI Entrez E-utilities.

Uses the recommended history-server pattern: a single esearch with
`usehistory="y"` returns a WebEnv + QueryKey, which subsequent efetch
calls reference to retrieve results in batches without re-running the
search.

CLI:

    python -m src.ingest --email you@example.com
    python -m src.ingest --email you@example.com --query 'EGFR lung cancer' --max 800
    python -m src.ingest --email you@example.com --api-key XXXX --max 5000
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

from Bio import Entrez

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "data" / "raw"
DEFAULT_QUERY = '"kinase inhibitor"[Title/Abstract] AND cancer'
DEFAULT_MAX = 1500
BATCH_SIZE = 200      # PubMed allows up to several hundred records per efetch
SLEEP_BETWEEN = 0.5   # 2 req/s — comfortably under the 3 req/s no-key limit


def fetch_pubmed(
    query: str,
    email: str,
    max_records: int = DEFAULT_MAX,
    api_key: str | None = None,
    out_path: Path | None = None,
    tool: str = "biomedical-ner-network",
) -> Path:
    """Download PubMed records matching `query`. Returns path to the XML file."""
    if not email or "@" not in email:
        raise ValueError(
            "NCBI requires a valid contact email. "
            "Pass --email you@example.com on the CLI."
        )

    Entrez.email = email
    Entrez.tool = tool
    if api_key:
        Entrez.api_key = api_key

    out_path = out_path or (RAW_DIR / "pubmed_records.xml")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Step 1: esearch with usehistory to capture the result set on NCBI's side.
    print(f"esearch: {query!r} (max {max_records:,})")
    with Entrez.esearch(
        db="pubmed",
        term=query,
        retmax=max_records,
        usehistory="y",
        sort="pub_date",
    ) as h:
        result = Entrez.read(h)

    total_found = int(result["Count"])
    pmids = result["IdList"]
    webenv = result["WebEnv"]
    query_key = result["QueryKey"]

    n_to_fetch = min(len(pmids), max_records)
    print(f"  found {total_found:,} matching records, fetching {n_to_fetch:,}")

    if n_to_fetch == 0:
        raise RuntimeError("No records match this query.")

    # Step 2: efetch in batches, concatenating the XML.
    # We keep the outer <PubmedArticleSet> wrapper from the first batch
    # and strip it from subsequent ones.
    with out_path.open("w", encoding="utf-8") as out_f:
        for start in range(0, n_to_fetch, BATCH_SIZE):
            end = min(start + BATCH_SIZE, n_to_fetch)
            print(f"  efetch: records {start:,}..{end - 1:,}")
            with Entrez.efetch(
                db="pubmed",
                rettype="xml",
                retmode="xml",
                retstart=start,
                retmax=end - start,
                webenv=webenv,
                query_key=query_key,
            ) as h:
                chunk = h.read()
                if isinstance(chunk, bytes):
                    chunk = chunk.decode("utf-8")

            if start == 0:
                out_f.write(chunk)
            else:
                # Concat: drop the <?xml?> prolog and the <PubmedArticleSet>
                # opening tag, and the closing </PubmedArticleSet> from the
                # previous content. Simplest robust approach:
                # extract just the <PubmedArticle>...</PubmedArticle> entries.
                inner = _extract_inner_articles(chunk)
                # Re-open the file in r+ mode to insert before the closing tag
                out_f.seek(0, 2)  # end of file
                # We instead rewrite the trailing tag handling by buffering:
                # easier: strip closing tag from accumulated file, append inner, re-add close.
                _append_articles(out_path, inner)

            time.sleep(SLEEP_BETWEEN)

    size_kb = out_path.stat().st_size / 1024
    print(f"Saved {out_path} ({size_kb:,.0f} KB)")
    return out_path


def _extract_inner_articles(xml_chunk: str) -> str:
    """Strip the outer <PubmedArticleSet> wrapper, keep <PubmedArticle> entries."""
    start_tag = "<PubmedArticleSet>"
    end_tag = "</PubmedArticleSet>"
    s = xml_chunk.find(start_tag)
    e = xml_chunk.rfind(end_tag)
    if s == -1 or e == -1:
        return xml_chunk
    return xml_chunk[s + len(start_tag):e]


def _append_articles(path: Path, inner: str) -> None:
    """Insert `inner` content before the closing </PubmedArticleSet> tag."""
    text = path.read_text(encoding="utf-8")
    end_tag = "</PubmedArticleSet>"
    idx = text.rfind(end_tag)
    if idx == -1:
        path.write_text(text + inner, encoding="utf-8")
    else:
        path.write_text(text[:idx] + inner + text[idx:], encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--email", required=True,
                        help="Contact email for NCBI (required by their policy)")
    parser.add_argument("--query", default=DEFAULT_QUERY,
                        help=f"PubMed query string (default: {DEFAULT_QUERY!r})")
    parser.add_argument("--max", type=int, default=DEFAULT_MAX, dest="max_records",
                        help=f"Max records to fetch (default: {DEFAULT_MAX})")
    parser.add_argument("--api-key", default=None,
                        help="Optional NCBI API key (allows 10 req/s instead of 3)")
    args = parser.parse_args()

    fetch_pubmed(
        query=args.query,
        email=args.email,
        max_records=args.max_records,
        api_key=args.api_key,
    )


if __name__ == "__main__":
    main()
