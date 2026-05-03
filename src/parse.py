"""Parse PubMed XML into a JSONL file of structured abstract records.

Output schema (one JSON object per line):

    {
        "pmid": "12345678",
        "title": "...",
        "abstract": "...",                 # full abstract text, joined if structured
        "abstract_sections": {              # only present for structured abstracts
            "BACKGROUND": "...",
            "METHODS": "..."
        },
        "year": 2020,                       # may be null
        "journal": "Nature",
        "mesh_terms": ["Antineoplastic Agents", "...", ...],
        "publication_types": ["Journal Article", "Review", ...]
    }

We deliberately skip authors, affiliations, references — they aren't
needed for the entity-level analysis and would bloat the file.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Iterator
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "data" / "raw"


def text_of(elem: ET.Element | None) -> str:
    """Return the full text content of an element, including children, stripped."""
    if elem is None:
        return ""
    # Use itertext to get text within any nested formatting (italic, sub, etc.)
    return re.sub(r"\s+", " ", " ".join(elem.itertext())).strip()


def parse_year(article: ET.Element) -> int | None:
    """Best-effort year extraction from PubDate / MedlineDate."""
    pubdate = article.find(".//Journal/JournalIssue/PubDate")
    if pubdate is None:
        return None

    y = pubdate.findtext("Year")
    if y and y.isdigit():
        return int(y)

    # MedlineDate is a free-text field like "2020 Jan-Feb" or "2019 Spring"
    medline = pubdate.findtext("MedlineDate") or ""
    match = re.search(r"\b(19|20)\d{2}\b", medline)
    if match:
        return int(match.group(0))

    return None


def parse_abstract(article: ET.Element) -> tuple[str, dict[str, str]]:
    """Extract abstract text. Returns (joined_text, sections_dict).

    Structured abstracts have multiple AbstractText elements with @Label
    (BACKGROUND, METHODS, ...). We capture both the joined version (for
    NER) and the labelled sections (for filtering / analysis).
    """
    abstract_elem = article.find("Abstract")
    if abstract_elem is None:
        return "", {}

    sections: dict[str, str] = {}
    parts: list[str] = []

    for at in abstract_elem.findall("AbstractText"):
        text = text_of(at)
        if not text:
            continue
        parts.append(text)
        label = at.attrib.get("Label")
        if label:
            sections[label.upper()] = text

    return " ".join(parts), sections


def parse_mesh(citation: ET.Element) -> list[str]:
    """Return list of MeSH descriptor names (without qualifiers)."""
    out: list[str] = []
    for mh in citation.findall(".//MeshHeadingList/MeshHeading"):
        name = mh.findtext("DescriptorName")
        if name:
            out.append(name.strip())
    return out


def parse_publication_types(article: ET.Element) -> list[str]:
    return [
        pt.text.strip()
        for pt in article.findall(".//PublicationTypeList/PublicationType")
        if pt.text and pt.text.strip()
    ]


def parse_one_article(pa: ET.Element) -> dict | None:
    citation = pa.find("MedlineCitation")
    if citation is None:
        return None

    pmid = citation.findtext("PMID")
    if not pmid:
        return None

    article = citation.find("Article")
    if article is None:
        return None

    title = text_of(article.find("ArticleTitle"))
    abstract_text, abstract_sections = parse_abstract(article)

    # Skip records with no abstract — useless for NER.
    if not abstract_text:
        return None

    record = {
        "pmid": pmid,
        "title": title,
        "abstract": abstract_text,
        "year": parse_year(article),
        "journal": text_of(article.find(".//Journal/Title")),
        "mesh_terms": parse_mesh(citation),
        "publication_types": parse_publication_types(article),
    }
    if abstract_sections:
        record["abstract_sections"] = abstract_sections
    return record


def iter_articles(xml_path: Path) -> Iterator[dict]:
    """Stream-parse the XML file, yielding one record per <PubmedArticle>."""
    # iterparse is memory-efficient for large XML files.
    context = ET.iterparse(str(xml_path), events=("end",))
    for event, elem in context:
        if elem.tag == "PubmedArticle":
            rec = parse_one_article(elem)
            if rec is not None:
                yield rec
            elem.clear()  # free memory


def parse_file(xml_path: Path, out_path: Path) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with out_path.open("w", encoding="utf-8") as f:
        for rec in iter_articles(xml_path):
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
            n += 1
    return n


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-path", default=None,
                        help="Input XML path (default data/raw/pubmed_records.xml)")
    parser.add_argument("--out-path", default=None,
                        help="Output JSONL path (default data/raw/abstracts.jsonl)")
    args = parser.parse_args()

    in_path = Path(args.in_path) if args.in_path else RAW_DIR / "pubmed_records.xml"
    out_path = Path(args.out_path) if args.out_path else RAW_DIR / "abstracts.jsonl"

    if not in_path.exists():
        raise SystemExit(
            f"Input XML not found: {in_path}\n"
            f"Run: python -m src.ingest --email you@example.com"
        )

    n = parse_file(in_path, out_path)
    print(f"Parsed {n:,} records (with abstract) → {out_path}")


if __name__ == "__main__":
    main()
