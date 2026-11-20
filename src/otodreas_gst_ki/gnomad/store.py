import re
from datetime import UTC, datetime
from pathlib import Path

import duckdb
import polars as pl

from otodreas_gst_ki.gnomad import structs

GENE_TABLE_PREFIX = "gene_"
_GENE_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


def gene_table_name(gene: str) -> str:
    if not _GENE_NAME.fullmatch(gene):
        raise ValueError(f"Invalid gene name: {gene!r}")
    return f"{GENE_TABLE_PREFIX}{gene}"


class GnomadStore:
    def __init__(
        self,
        database: str | Path = ":memory:",
        database_version: str | None = None,
    ):
        self.database_version = database_version
        self.connection = duckdb.connect(str(database))
        self._create_metadata_table()

    def close(self) -> None:
        self.connection.close()

    def _create_metadata_table(self) -> None:
        self.connection.execute(
            """
            CREATE TABLE IF NOT EXISTS genes (
                gene_name VARCHAR PRIMARY KEY,
                ensembl_id VARCHAR NOT NULL,
                gnomad_version VARCHAR NOT NULL,
                canonical_transcript VARCHAR NOT NULL,
                chromosomes INTEGER[],
                n_variants INTEGER NOT NULL,
                n_variants_in_canonical_transcript INTEGER NOT NULL,
                source VARCHAR,
                fetched_at TIMESTAMP NOT NULL,
                table_name VARCHAR NOT NULL UNIQUE,
                database_version VARCHAR
            )
            """
        )
        columns = {
            row[0]
            for row in self.connection.execute(
                "SELECT column_name FROM information_schema.columns "
                "WHERE table_name = 'genes'"
            ).fetchall()
        }
        if "database_version" not in columns:
            self.connection.execute(
                "ALTER TABLE genes ADD COLUMN database_version VARCHAR"
            )

    def write_dataset(self, dataset: structs.GnomadDataset) -> str:
        table_name = gene_table_name(dataset.gene)
        self.connection.begin()
        try:
            self.connection.register("_gene_frame", dataset.dataMatrix.to_arrow())
            self.connection.execute(
                f'CREATE OR REPLACE TABLE "{table_name}" AS SELECT * FROM "_gene_frame"'
            )
            self.connection.execute(
                """
                INSERT INTO genes VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT (gene_name) DO UPDATE SET
                    ensembl_id = excluded.ensembl_id,
                    gnomad_version = excluded.gnomad_version,
                    canonical_transcript = excluded.canonical_transcript,
                    chromosomes = excluded.chromosomes,
                    n_variants = excluded.n_variants,
                    n_variants_in_canonical_transcript =
                        excluded.n_variants_in_canonical_transcript,
                    source = excluded.source,
                    fetched_at = excluded.fetched_at,
                    table_name = excluded.table_name,
                    database_version = excluded.database_version
                """,
                [
                    dataset.gene,
                    dataset.ensemblId,
                    dataset.gnomadVersion,
                    dataset.canonicalTranscript,
                    dataset.chromosomes,
                    dataset.nVars,
                    dataset.nVarsInCanonicalTranscript,
                    dataset.source,
                    datetime.now(UTC).replace(tzinfo=None),
                    table_name,
                    self.database_version,
                ],
            )
            self.connection.unregister("_gene_frame")
            self.connection.commit()
        except Exception:
            self.connection.rollback()
            raise
        return table_name

    def list_genes(self) -> pl.DataFrame:
        gene_list = pl.from_arrow(
            self.connection.sql("SELECT * FROM genes").to_arrow_table()
        )
        if not type(gene_list) is pl.DataFrame:
            raise TypeError(f"Expected polars DataFrame, got {type(gene_list)}")

        return gene_list

    def read_gene(self, gene: str) -> pl.DataFrame:
        table_name = gene_table_name(gene)
        cursor = self.connection.cursor()
        exists = cursor.execute(
            "SELECT 1 FROM information_schema.tables WHERE table_name = ?",
            [table_name],
        ).fetchone()
        cursor.close()
        if not exists:
            raise KeyError(f"Gene not found: {gene}")
        gene_table = pl.from_arrow(
            self.connection.sql(f'SELECT * FROM "{table_name}"').to_arrow_table()
        )
        if type(gene_table) is not pl.DataFrame:
            raise TypeError(f"Expected polars DataFrame, got {type(gene_table)}")

        return gene_table
