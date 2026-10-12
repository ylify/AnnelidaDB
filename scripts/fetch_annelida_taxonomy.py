#!/usr/bin/env python3
"""Export NCBI Annelida species-level taxonomy with GBIF and WoRMS IDs.

The exporter downloads NCBI's ``new_taxdump`` once, then caches all GBIF and
WoRMS lookups so interrupted runs can resume without repeating requests.
"""

import argparse
import csv
from dataclasses import dataclass
import hashlib
import io
import json
from pathlib import Path
import tarfile
import tempfile
import time
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode, urlsplit
from urllib.request import Request, urlopen


BASE_DIR = Path(__file__).resolve().parents[1]
DEFAULT_TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
GBIF_MATCH_URL = "https://api.gbif.org/v1/species/match"
WORMS_NAME_URL = "https://www.marinespecies.org/rest/AphiaRecordsByName/{name}"
INCLUDED_RANKS = {"species", "subspecies", "varietas", "forma"}
CSV_HEADER = [
    "Phylum", "Class", "Order", "Family", "Geneus", "Species",
    "NCBI taxonomy id", "gbif_id", "worms_id",
]
AUDIT_HEADER = ["NCBI taxonomy id", "Species", "source", "reason"]


@dataclass(frozen=True)
class TaxonRow:
    taxid: str
    name: str
    rank: str
    phylum: str
    class_name: str
    order: str
    family: str
    genus: str


@dataclass(frozen=True)
class LookupResult:
    identifier: str
    reason: str

    def __post_init__(self) -> None:
        if bool(self.identifier) == bool(self.reason):
            raise ValueError("LookupResult requires exactly one of identifier or reason")


@dataclass(frozen=True)
class ExportSummary:
    written: int
    gbif_matched: int
    worms_matched: int
    unmatched: int


class LookupCache:
    """File cache keyed by complete URL and query parameters."""

    def __init__(self, directory: Union[str, Path]):
        self.directory = Path(directory)

    @staticmethod
    def key(url: str, params: Dict[str, str]) -> str:
        encoded = json.dumps([url, sorted(params.items())], separators=(",", ":"), ensure_ascii=False)
        return hashlib.sha256(encoded.encode("utf-8")).hexdigest()

    def get(self, url: str, params: Dict[str, str]) -> Optional[Any]:
        path = self.directory / (self.key(url, params) + ".json")
        if not path.exists():
            return None
        with path.open(encoding="utf-8") as handle:
            return json.load(handle)

    def set(self, url: str, params: Dict[str, str], value: Any) -> None:
        self.directory.mkdir(parents=True, exist_ok=True)
        destination = self.directory / (self.key(url, params) + ".json")
        with tempfile.NamedTemporaryFile(
            "w", encoding="utf-8", dir=str(self.directory), delete=False
        ) as handle:
            json.dump(value, handle, ensure_ascii=False)
            temporary = Path(handle.name)
        temporary.replace(destination)


class JsonHttpClient:
    """Small retrying JSON client with a per-host request delay."""

    def __init__(self, request_delay: float = 0.25, timeout: int = 30):
        if request_delay < 0:
            raise ValueError("request_delay must not be negative")
        self.request_delay = request_delay
        self.timeout = timeout
        self._last_request_by_host: Dict[str, float] = {}

    def get_json(self, url: str, params: Dict[str, str]) -> Optional[Any]:
        host = urlsplit(url).netloc
        request_url = url + ("?" + urlencode(params) if params else "")
        for retry_delay in (1, 2, 4, None):
            elapsed = time.monotonic() - self._last_request_by_host.get(host, 0.0)
            if elapsed < self.request_delay:
                time.sleep(self.request_delay - elapsed)
            self._last_request_by_host[host] = time.monotonic()
            try:
                request = Request(request_url, headers={"User-Agent": "AnnelidaDB-taxonomy-export/1.0"})
                with urlopen(request, timeout=self.timeout) as response:
                    return json.load(response)
            except HTTPError as error:
                retryable = error.code == 429 or 500 <= error.code <= 599
                if not retryable:
                    return None
            except (URLError, TimeoutError, json.JSONDecodeError):
                pass
            if retry_delay is None:
                return None
            time.sleep(retry_delay)
        return None


def _dmp_fields(line: str) -> List[str]:
    return [field.strip() for field in line.rstrip("\n").split("\t|")]


def _parse_nodes(lines: Iterable[str]) -> Dict[str, Tuple[str, str]]:
    details: Dict[str, Tuple[str, str]] = {}
    for line in lines:
        fields = _dmp_fields(line)
        if len(fields) >= 3 and fields[0]:
            details[fields[0]] = (fields[1], fields[2])
    return details


def _parse_names(
    lines: Iterable[str], node_details: Dict[str, Tuple[str, str]]
) -> Dict[str, Tuple[str, str]]:
    names: Dict[str, Tuple[str, str]] = {}
    for line in lines:
        fields = _dmp_fields(line)
        if len(fields) >= 4 and fields[3] == "scientific name" and fields[0] in node_details:
            names[fields[0]] = (fields[1], node_details[fields[0]][1])
    return names


def parse_taxdump(nodes_text: str, names_text: str) -> Tuple[Dict[str, str], Dict[str, Tuple[str, str]]]:
    """Parse NCBI ``nodes.dmp`` and ``names.dmp`` text for tests or small inputs."""
    node_details = _parse_nodes(io.StringIO(nodes_text))
    names = _parse_names(io.StringIO(names_text), node_details)
    return {taxid: parent for taxid, (parent, _rank) in node_details.items()}, names


def load_taxdump(archive_path: Union[str, Path]) -> Tuple[Dict[str, str], Dict[str, Tuple[str, str]]]:
    """Read nodes and scientific names from an NCBI ``new_taxdump`` archive."""
    with tarfile.open(archive_path, "r:*") as archive:
        members = {Path(member.name).name: member for member in archive.getmembers() if member.isfile()}
        try:
            nodes_file = archive.extractfile(members["nodes.dmp"])
            names_file = archive.extractfile(members["names.dmp"])
        except KeyError as error:
            raise ValueError("Taxdump must contain nodes.dmp and names.dmp") from error
        if nodes_file is None or names_file is None:
            raise ValueError("Taxdump contains an unreadable nodes.dmp or names.dmp")
        with nodes_file, names_file:
            node_details = _parse_nodes(io.TextIOWrapper(nodes_file, encoding="utf-8"))
            names = _parse_names(io.TextIOWrapper(names_file, encoding="utf-8"), node_details)
    return {taxid: parent for taxid, (parent, _rank) in node_details.items()}, names


def select_annelida_taxa(
    nodes: Dict[str, str], names: Dict[str, Tuple[str, str]]
) -> List[TaxonRow]:
    annelida = [
        taxid for taxid, (name, rank) in names.items()
        if name == "Annelida" and rank == "phylum"
    ]
    if len(annelida) != 1:
        raise ValueError("Expected one scientific-name Annelida phylum node")

    children: Dict[str, List[str]] = {}
    for taxid, parent in nodes.items():
        if taxid != parent:
            children.setdefault(parent, []).append(taxid)

    descendants = set()
    pending = [annelida[0]]
    while pending:
        taxid = pending.pop()
        if taxid in descendants:
            continue
        descendants.add(taxid)
        pending.extend(children.get(taxid, []))

    rows: List[TaxonRow] = []
    for taxid in descendants:
        if taxid not in names:
            continue
        name, rank = names[taxid]
        if rank not in INCLUDED_RANKS:
            continue
        lineage = {"phylum": "", "class": "", "order": "", "family": "", "genus": ""}
        current = taxid
        while current in nodes:
            if current in names:
                lineage_rank = names[current][1]
                if lineage_rank in lineage and not lineage[lineage_rank]:
                    lineage[lineage_rank] = names[current][0]
            parent = nodes[current]
            if parent == current:
                break
            current = parent
        rows.append(TaxonRow(
            taxid=taxid,
            name=name,
            rank=rank,
            phylum=lineage["phylum"],
            class_name=lineage["class"],
            order=lineage["order"],
            family=lineage["family"],
            genus=lineage["genus"],
        ))
    return sorted(rows, key=lambda row: (row.name.casefold(), int(row.taxid)))


def to_csv_row(row: TaxonRow, gbif_id: str, worms_id: str) -> List[str]:
    value = lambda item: item or "Unassigned"
    return [
        value(row.phylum), value(row.class_name), value(row.order), value(row.family),
        value(row.genus), row.name, row.taxid, gbif_id, worms_id,
    ]


def _get_json(
    client: JsonHttpClient, url: str, params: Dict[str, str], cache: Optional[LookupCache]
) -> Optional[Any]:
    if cache is not None:
        cached = cache.get(url, params)
        if cached is not None:
            return cached
    value = client.get_json(url, params)
    if cache is not None and value is not None:
        cache.set(url, params, value)
    return value


def is_exact_name_match(candidate: str, requested: str) -> bool:
    return candidate.casefold().strip() == requested.casefold().strip()


def resolve_gbif_id(
    name: str, client: JsonHttpClient, cache: Optional[LookupCache] = None
) -> LookupResult:
    record = _get_json(client, GBIF_MATCH_URL, {"name": name}, cache)
    if not isinstance(record, dict):
        return LookupResult("", "missing_or_request_failed")
    if record.get("matchType") != "EXACT":
        return LookupResult("", "non_exact_match")
    if record.get("alternatives") or "ambiguous" in str(record.get("note", "")).casefold():
        return LookupResult("", "ambiguous_match")
    confidence = record.get("confidence")
    if isinstance(confidence, (int, float)) and confidence < 90:
        return LookupResult("", "low_confidence_match")
    identifier = record.get("usageKey")
    if not identifier:
        return LookupResult("", "missing_usage_key")
    return LookupResult(str(identifier), "")


def resolve_worms_id(
    name: str, client: JsonHttpClient, cache: Optional[LookupCache] = None
) -> LookupResult:
    url = WORMS_NAME_URL.format(name=quote(name, safe=""))
    records = _get_json(client, url, {"marine_only": "false"}, cache)
    if not isinstance(records, list):
        return LookupResult("", "missing_or_request_failed")
    exact = [
        record for record in records
        if isinstance(record, dict) and is_exact_name_match(str(record.get("scientificname", "")), name)
    ]
    if len(exact) != 1:
        return LookupResult("", "ambiguous_or_missing")
    identifier = exact[0].get("valid_AphiaID") or exact[0].get("AphiaID")
    if not identifier:
        return LookupResult("", "missing_aphia_id")
    return LookupResult(str(identifier), "")


def _write_csv_atomic(path: Path, rows: Iterable[Sequence[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8-sig", newline="", dir=str(path.parent), delete=False
    ) as handle:
        writer = csv.writer(handle)
        writer.writerows(rows)
        temporary = Path(handle.name)
    temporary.replace(path)


def export_rows(
    rows: Iterable[TaxonRow], output: Union[str, Path], client: JsonHttpClient,
    cache: Optional[LookupCache] = None,
) -> ExportSummary:
    output_path = Path(output)
    audit_path = output_path.with_name(output_path.stem + ".unmatched.csv")
    csv_rows: List[List[str]] = [CSV_HEADER]
    audit_rows: List[List[str]] = [AUDIT_HEADER]
    gbif_matched = worms_matched = 0

    for row in rows:
        gbif = resolve_gbif_id(row.name, client, cache)
        worms = resolve_worms_id(row.name, client, cache)
        if gbif.identifier:
            gbif_matched += 1
        else:
            audit_rows.append([row.taxid, row.name, "GBIF", gbif.reason])
        if worms.identifier:
            worms_matched += 1
        else:
            audit_rows.append([row.taxid, row.name, "WoRMS", worms.reason])
        csv_rows.append(to_csv_row(row, gbif.identifier, worms.identifier))

    _write_csv_atomic(output_path, csv_rows)
    _write_csv_atomic(audit_path, audit_rows)
    return ExportSummary(
        written=len(csv_rows) - 1,
        gbif_matched=gbif_matched,
        worms_matched=worms_matched,
        unmatched=len(audit_rows) - 1,
    )


def download_taxdump(url: str, cache_dir: Union[str, Path], refresh: bool = False) -> Path:
    cache_path = Path(cache_dir)
    destination = cache_path / "new_taxdump.tar.gz"
    if destination.exists() and not refresh:
        return destination
    cache_path.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers={"User-Agent": "AnnelidaDB-taxonomy-export/1.0"})
    for retry_delay in (1, 2, 4, None):
        try:
            with urlopen(request, timeout=120) as response:
                with tempfile.NamedTemporaryFile("wb", dir=str(cache_path), delete=False) as handle:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        handle.write(chunk)
                    temporary = Path(handle.name)
            temporary.replace(destination)
            return destination
        except (HTTPError, URLError, TimeoutError):
            if retry_delay is None:
                raise
            time.sleep(retry_delay)
    raise RuntimeError("Unable to download NCBI taxdump")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=BASE_DIR / "taxonomy.ncbi.csv")
    parser.add_argument("--cache-dir", type=Path, default=BASE_DIR / ".cache" / "annelida-taxonomy")
    parser.add_argument("--refresh", action="store_true", help="Download a fresh NCBI taxdump.")
    parser.add_argument("--request-delay", type=float, default=0.25, help="Minimum seconds between requests to each host.")
    parser.add_argument("--taxdump-url", default=DEFAULT_TAXDUMP_URL)
    args = parser.parse_args(argv)

    if args.request_delay < 0:
        parser.error("--request-delay must not be negative")
    if args.output.resolve() == (BASE_DIR / "taxonomy.csv").resolve():
        parser.error("Refusing to replace taxonomy.csv; choose a new --output path.")

    taxdump = download_taxdump(args.taxdump_url, args.cache_dir, args.refresh)
    nodes, names = load_taxdump(taxdump)
    rows = select_annelida_taxa(nodes, names)
    client = JsonHttpClient(request_delay=args.request_delay)
    summary = export_rows(rows, args.output, client, LookupCache(args.cache_dir / "lookups"))
    print(
        "Wrote {0.written} rows; GBIF matched: {0.gbif_matched}; "
        "WoRMS matched: {0.worms_matched}; unresolved lookups: {0.unmatched}.".format(summary)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
