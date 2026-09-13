#!/usr/bin/env python3
"""Fetch a corpus from NCBI for any pathogen, and parse it.

REPLACES a genuinely dangerous script. ``parse_fmdv_corpus_metadata.py``
took its paths from module-level constants and accepted no arguments at
all — but ``workflow/modules/acquisition.nf`` invokes it with
``--accessions``, ``--genbank-out``, ``--metadata-out`` and
``--fasta-out``. Python ignores unrecognised ``sys.argv`` unless something
reads it, so the script ran, silently discarded every flag, re-parsed the
FMDV GenBank file and overwrote the FMDV corpus. Running
``-entry ACQUISITION --pathogen ndv`` would therefore have reported
success while rewriting a different pathogen's corpus and fetching
nothing.

Nothing here is pathogen-specific. Lineage lives in the pathogen's own
config (``corpus.lineage_field``), so the metadata table carries every
source qualifier and the config decides which column is the lineage.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

#: Source qualifiers copied into the metadata table. Deliberately broad:
#: which one is "lineage" differs per pathogen and is the config's call,
#: so the fetcher must not decide it here.
QUALIFIERS = (
    "collection_date", "country", "geo_loc_name", "host", "strain",
    "isolate", "serotype", "genotype", "segment", "note", "lat_lon",
)


def _get(endpoint: str, params: dict, *, retries: int = 4) -> bytes:
    """One E-utilities call, with backoff.

    NCBI rate-limits and intermittently 500s under load. A fetch that gave
    up on the first failure would silently produce a short corpus, which
    is worse than taking longer.
    """
    url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}"
    last: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=180) as response:
                return response.read()
        except Exception as exc:  # noqa: BLE001 - retried, then re-raised
            last = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"{endpoint} failed after {retries} attempts: {last}")


def search(term: str, *, retmax: int) -> list[str]:
    import re

    payload = _get("esearch.fcgi", {
        "db": "nuccore", "term": term, "retmax": retmax, "usehistory": "y",
    }).decode()
    return re.findall(r"<Id>(\d+)</Id>", payload)


def fetch_genbank(ids: list[str], out_path: Path, *, batch: int = 200) -> int:
    """Stream GenBank records to disk in batches.

    Written straight to the file rather than accumulated in memory: an
    1,800-genome corpus is tens of megabytes and there is no reason to
    hold it all before writing a byte.
    """
    written = 0
    with out_path.open("wb") as handle:
        for start in range(0, len(ids), batch):
            chunk = ids[start:start + batch]
            handle.write(_get("efetch.fcgi", {
                "db": "nuccore", "id": ",".join(chunk),
                "rettype": "gbwithparts", "retmode": "text",
            }))
            written += len(chunk)
            print(f"  fetched {written}/{len(ids)}", file=sys.stderr)
            time.sleep(0.4)
    return written


def parse(genbank_path: Path, metadata_out: Path, fasta_out: Path) -> int:
    from Bio import SeqIO

    rows: list[dict] = []
    seen: set[str] = set()
    with fasta_out.open("w") as fasta:
        for record in SeqIO.parse(str(genbank_path), "genbank"):
            accession = record.id
            # Duplicate accessions would silently collapse downstream
            # (g4watch.io.fasta raises on them now); skip and report
            # rather than write a file that cannot be read.
            if accession in seen:
                print(f"  duplicate accession skipped: {accession}", file=sys.stderr)
                continue
            seen.add(accession)

            source = next((f for f in record.features if f.type == "source"), None)
            row = {"accession": accession, "length": len(record.seq),
                   "organism": record.annotations.get("organism", "")}
            for key in QUALIFIERS:
                value = ""
                if source is not None:
                    values = source.qualifiers.get(key)
                    value = values[0] if values else ""
                row[key] = value
            # NCBI migrated /country to /geo_loc_name; carry one column so
            # a config does not have to know which era a record is from.
            #
            # And keep the COUNTRY, not the locality. The qualifier is
            # formatted "Country:region,town", so "Japan:Gifu,Gifu" and
            # "Japan:Gifu,Seki" are the same country — but a config using
            # country as its lineage field would see two lineages, and
            # every such split pushes real groups under the Appendix C
            # per-lineage floor. The full value is preserved in
            # geo_loc_name for anyone who wants the locality.
            raw_country = row.get("country") or row.get("geo_loc_name") or ""
            row["geo_loc_name"] = raw_country
            row["country"] = raw_country.split(":")[0].strip()
            rows.append(row)
            fasta.write(f">{accession}\n{str(record.seq)}\n")

    if not rows:
        raise SystemExit(f"{genbank_path}: no records parsed — refusing to write an empty corpus")

    fields = ["accession", "length", "organism", "country", *[q for q in QUALIFIERS if q != "country"]]
    with metadata_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--term", help='NCBI query, e.g. \'"Newcastle disease virus"[Organism] AND 14500:16000[SLEN]\'')
    parser.add_argument("--accessions", help="file of accessions, one per line (instead of --term)")
    parser.add_argument("--genbank-out", required=True)
    parser.add_argument("--metadata-out", required=True)
    parser.add_argument("--fasta-out", required=True)
    parser.add_argument("--retmax", type=int, default=5000)
    parser.add_argument("--reuse-genbank", action="store_true",
                        help="parse an existing GenBank file instead of re-fetching")
    args = parser.parse_args()

    genbank = Path(args.genbank_out)
    genbank.parent.mkdir(parents=True, exist_ok=True)
    Path(args.metadata_out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.fasta_out).parent.mkdir(parents=True, exist_ok=True)

    if not args.reuse_genbank:
        if args.accessions:
            ids = [line.strip() for line in Path(args.accessions).read_text().splitlines() if line.strip()]
        elif args.term:
            ids = search(args.term, retmax=args.retmax)
        else:
            parser.error("one of --term or --accessions is required")
        if not ids:
            raise SystemExit("query matched no records — refusing to write an empty corpus")
        print(f"  {len(ids)} record(s) to fetch", file=sys.stderr)
        fetch_genbank(ids, genbank)

    n = parse(genbank, Path(args.metadata_out), Path(args.fasta_out))
    print(f"  parsed {n} record(s)")
    print(f"  metadata: {args.metadata_out}")
    print(f"  fasta   : {args.fasta_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
