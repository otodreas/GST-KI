import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from otodreas_gst_ki.gnomad.store import GnomadStore

ENSG_PATTERN = re.compile(r"(ENSG\d+(?:\.\d+)?)")
MISSING_VALUES = {"", "na", "n/a", "none", "null", "nan"}


@dataclass(frozen=True)
class Disagreement:
    table_name: str
    csv_path: Path
    gnomad_id: str
    column: str
    table_value: object
    csv_value: object


@dataclass(frozen=True)
class VerificationResult:
    database: Path
    checked_tables: tuple[str, ...]
    checked_csvs: tuple[Path, ...]
    compared_rows: int
    compared_values: int
    disagreements: tuple[Disagreement, ...]
    missing_csv_tables: tuple[str, ...]

    def has_disagreements(self) -> bool:
        return bool(self.disagreements)


def _is_missing(value: object) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return isinstance(value, str) and value.strip().lower() in MISSING_VALUES


def _normalized(value: object) -> object:
    if _is_missing(value):
        return None
    if isinstance(value, str):
        return value.strip()
    return value


def _values_equal(table_value: object, csv_value: object, column: str = "") -> bool:
    table_value = _normalized(table_value)
    csv_value = _normalized(csv_value)
    if table_value is None or csv_value is None:
        return table_value is csv_value

    if column == "Flags":
        table_flags = {
            flag.strip() for flag in str(table_value).split(",") if flag.strip()
        }
        csv_flags = {flag.strip() for flag in str(csv_value).split(",") if flag.strip()}
        return table_flags <= csv_flags

    if isinstance(table_value, str) and isinstance(csv_value, str):
        table_ensembl = re.fullmatch(r"(ENST\d+)(?:\.\d+)?", table_value)
        csv_ensembl = re.fullmatch(r"(ENST\d+)(?:\.\d+)?", csv_value)
        if table_ensembl and csv_ensembl:
            return table_ensembl.group(1) == csv_ensembl.group(1)

    if isinstance(table_value, (int, float)) and not isinstance(table_value, bool):
        try:
            csv_number = float(csv_value)
        except TypeError, ValueError:
            return False
        return math.isclose(
            float(table_value),
            csv_number,
            rel_tol=1e-12,
            abs_tol=1e-12,
        )
    return str(table_value) == str(csv_value)


def _candidate_csvs(
    csv_root: Path,
    ensembl_id: str,
    gnomad_version: str,
) -> list[Path]:
    short_id = ensembl_id.split(".")[0]
    candidates = []
    for path in csv_root.rglob("*.csv"):
        match = ENSG_PATTERN.search(path.name)
        if match and match.group(1).split(".")[0] == short_id:
            candidates.append(path)

    version_markers = []
    if "r2" in gnomad_version:
        version_markers.append("v2")
    if "r4" in gnomad_version:
        version_markers.append("v4")
    if version_markers:
        candidates = [
            path
            for path in candidates
            if any(marker in path.name.lower() for marker in version_markers)
        ]

    return sorted(candidates, key=lambda path: str(path), reverse=True)


def _read_csv(path: Path) -> pl.DataFrame:
    return pl.read_csv(path, infer_schema_length=1000)


def verify_database(
    database: str | Path,
    csv_root: str | Path = "data/gnomad",
) -> VerificationResult:
    database = Path(database)
    csv_root = Path(csv_root)
    store = GnomadStore(database, read_only=True)
    disagreements: list[Disagreement] = []
    checked_tables: list[str] = []
    checked_csvs: list[Path] = []
    missing_csv_tables: list[str] = []
    compared_rows = 0
    compared_values = 0

    try:
        genes = store.list_genes()
        for gene_row in genes.to_dicts():
            table_name = gene_row["table_name"]
            if not table_name or not table_name.startswith("gene_"):
                continue

            candidates = _candidate_csvs(
                csv_root,
                gene_row["ensembl_id"],
                gene_row["gnomad_version"],
            )
            if not candidates:
                missing_csv_tables.append(table_name)
                continue

            csv_path = candidates[0]
            csv_frame = _read_csv(csv_path)
            table_frame = store.read_gene(gene_row["gene_name"])
            if (
                "gnomAD ID" not in table_frame.columns
                or "gnomAD ID" not in csv_frame.columns
            ):
                missing_csv_tables.append(table_name)
                continue

            shared_columns = [
                column
                for column in table_frame.columns
                if column != "gnomAD ID" and column in csv_frame.columns
            ]
            if "v2" in csv_path.name.lower():
                # gnomAD v2 CSVs label the protein consequence in both
                # HGVS-related columns. It is not a c.*/hgvsc equivalent.
                shared_columns = [
                    column for column in shared_columns if column != "HGVS Consequence"
                ]
            checked_tables.append(table_name)
            checked_csvs.append(csv_path)
            if not shared_columns:
                continue

            table_by_id = {
                str(row["gnomAD ID"]): row
                for row in table_frame.to_dicts()
                if row.get("gnomAD ID") is not None
            }
            csv_by_id = {
                str(row["gnomAD ID"]): row
                for row in csv_frame.to_dicts()
                if row.get("gnomAD ID") is not None
            }

            shared_ids = table_by_id.keys() & csv_by_id.keys()
            for gnomad_id in shared_ids:
                compared_rows += 1
                table_row = table_by_id[gnomad_id]
                csv_row = csv_by_id[gnomad_id]
                for column in shared_columns:
                    compared_values += 1
                    table_value = table_row.get(column)
                    csv_value = csv_row.get(column)
                    if (
                        column == "HGVS Consequence"
                        and isinstance(table_value, str)
                        and table_value.startswith("c.")
                        and isinstance(csv_value, str)
                        and csv_value.startswith("p.")
                        and "Protein Consequence" in table_row
                        and "Protein Consequence" in csv_row
                    ):
                        continue
                    if not _values_equal(table_value, csv_value, column):
                        disagreements.append(
                            Disagreement(
                                table_name=table_name,
                                csv_path=csv_path,
                                gnomad_id=gnomad_id,
                                column=column,
                                table_value=_normalized(table_value),
                                csv_value=_normalized(csv_value),
                            )
                        )
    finally:
        store.close()

    return VerificationResult(
        database=database,
        checked_tables=tuple(checked_tables),
        checked_csvs=tuple(dict.fromkeys(checked_csvs)),
        compared_rows=compared_rows,
        compared_values=compared_values,
        disagreements=tuple(disagreements),
        missing_csv_tables=tuple(missing_csv_tables),
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare gnomAD DuckDB variant tables with source CSV files."
    )
    parser.add_argument("database", type=Path)
    parser.add_argument(
        "--csv-root",
        type=Path,
        default=Path("data/gnomad"),
    )
    args = parser.parse_args()

    result = verify_database(args.database, args.csv_root)
    print(
        f"Checked {len(result.checked_tables)} tables against "
        f"{len(result.checked_csvs)} CSV files."
    )
    print(
        f"Compared {result.compared_values} shared values "
        f"for {result.compared_rows} shared variant IDs."
    )
    if result.missing_csv_tables:
        print("Tables without matching CSV:")
        for table_name in result.missing_csv_tables:
            print(f"  {table_name}")
    if result.disagreements:
        print(f"Found {len(result.disagreements)} disagreements:")
        for disagreement in result.disagreements:
            print(
                f"{disagreement.table_name} | {disagreement.gnomad_id} | "
                f"{disagreement.column}: table={disagreement.table_value!r}, "
                f"csv={disagreement.csv_value!r}"
            )


if __name__ == "__main__":
    main()
