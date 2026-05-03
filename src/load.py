"""Load helpers for parsed PubMed abstracts."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "data" / "raw"


def load_abstracts(path: Path | None = None) -> pd.DataFrame:
    """Return a DataFrame of parsed abstracts."""
    path = path or (RAW_DIR / "abstracts.jsonl")
    if not path.exists():
        raise FileNotFoundError(
            f"No parsed abstracts at {path}.\n"
            f"Run: python -m src.ingest --email you@example.com && python -m src.parse"
        )
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
    return pd.DataFrame(rows)
