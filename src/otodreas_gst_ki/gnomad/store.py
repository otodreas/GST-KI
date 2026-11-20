import re
from datetime import UTC, datetime
from pathlib import Path

import duckdb
import polars as pl

from otodreas_gst_ki.gnomad import structs

GENE_TABLE_PREFIX = "gene_"
SV_TABLE_PREFIX = "sv_"
_GENE_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


def gene_table_name(gene: str) -> str:
    if not _GENE_NAME.fullmatch(gene):
        raise ValueError(f"Invalid gene name: {gene!r}")
    return f"{GENE_TABLE_PREFIX}{gene}"


def sv_table_name(gene: str) -> str:
    if not _GENE_NAME.fullmatch(gene):
        raise ValueError(f"Invalid gene name: {gene!r}")
    return f"{SV_TABLE_PREFIX}{gene}"


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
                database_version VARCHAR,
                exp_lof DOUBLE,
                exp_mis DOUBLE,
                exp_syn DOUBLE,
                obs_lof INTEGER,
                obs_mis INTEGER,
                obs_syn INTEGER,
                oe_lof DOUBLE,
                oe_lof_lower DOUBLE,
                oe_lof_upper DOUBLE,
                oe_lof_percentile INTEGER,
                oe_mis DOUBLE,
                oe_mis_lower DOUBLE,
                oe_mis_upper DOUBLE,
                oe_syn DOUBLE,
                oe_syn_lower DOUBLE,
                oe_syn_upper DOUBLE,
                lof_z DOUBLE,
                mis_z DOUBLE,
                syn_z DOUBLE,
                pli DOUBLE,
                constraint_flags VARCHAR[]
            )
            """
        )
        existing_columns = {
            row[0]
            for row in self.connection.execute(
                "SELECT column_name FROM information_schema.columns "
                "WHERE table_name = 'genes'"
            ).fetchall()
        }
        constraint_columns = {
            "database_version": "VARCHAR",
            "exp_lof": "DOUBLE",
            "exp_mis": "DOUBLE",
            "exp_syn": "DOUBLE",
            "obs_lof": "INTEGER",
            "obs_mis": "INTEGER",
            "obs_syn": "INTEGER",
            "oe_lof": "DOUBLE",
            "oe_lof_lower": "DOUBLE",
            "oe_lof_upper": "DOUBLE",
            "oe_lof_percentile": "INTEGER",
            "oe_mis": "DOUBLE",
            "oe_mis_lower": "DOUBLE",
            "oe_mis_upper": "DOUBLE",
            "oe_syn": "DOUBLE",
            "oe_syn_lower": "DOUBLE",
            "oe_syn_upper": "DOUBLE",
            "lof_z": "DOUBLE",
            "mis_z": "DOUBLE",
            "syn_z": "DOUBLE",
            "pli": "DOUBLE",
            "constraint_flags": "VARCHAR[]",
        }
        for column_name, column_type in constraint_columns.items():
            if column_name not in existing_columns:
                self.connection.execute(
                    f"ALTER TABLE genes ADD COLUMN {column_name} {column_type}"
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
                INSERT INTO genes VALUES (
                    ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                    ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                    ?, ?, ?, ?, ?, ?, ?, ?
                )
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
                    database_version = excluded.database_version,
                    exp_lof = excluded.exp_lof,
                    exp_mis = excluded.exp_mis,
                    exp_syn = excluded.exp_syn,
                    obs_lof = excluded.obs_lof,
                    obs_mis = excluded.obs_mis,
                    obs_syn = excluded.obs_syn,
                    oe_lof = excluded.oe_lof,
                    oe_lof_lower = excluded.oe_lof_lower,
                    oe_lof_upper = excluded.oe_lof_upper,
                    oe_lof_percentile = excluded.oe_lof_percentile,
                    oe_mis = excluded.oe_mis,
                    oe_mis_lower = excluded.oe_mis_lower,
                    oe_mis_upper = excluded.oe_mis_upper,
                    oe_syn = excluded.oe_syn,
                    oe_syn_lower = excluded.oe_syn_lower,
                    oe_syn_upper = excluded.oe_syn_upper,
                    lof_z = excluded.lof_z,
                    mis_z = excluded.mis_z,
                    syn_z = excluded.syn_z,
                    pli = excluded.pli,
                    constraint_flags = excluded.constraint_flags
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
                    dataset.exp_lof,
                    dataset.exp_mis,
                    dataset.exp_syn,
                    dataset.obs_lof,
                    dataset.obs_mis,
                    dataset.obs_syn,
                    dataset.oe_lof,
                    dataset.oe_lof_lower,
                    dataset.oe_lof_upper,
                    dataset.oe_lof_percentile,
                    dataset.oe_mis,
                    dataset.oe_mis_lower,
                    dataset.oe_mis_upper,
                    dataset.oe_syn,
                    dataset.oe_syn_lower,
                    dataset.oe_syn_upper,
                    dataset.lof_z,
                    dataset.mis_z,
                    dataset.syn_z,
                    dataset.pli,
                    dataset.constraint_flags,
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
