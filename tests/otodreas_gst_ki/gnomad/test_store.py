from pathlib import Path

import polars as pl
import pytest

from otodreas_gst_ki.gnomad import database as database_module
from otodreas_gst_ki.gnomad import structs
from otodreas_gst_ki.gnomad.database import build_database, rebuild_database
from otodreas_gst_ki.gnomad.store import GnomadStore
from otodreas_gst_ki.gnomad.verify import verify_database


def _dataset(gene: str = "GSTT2B") -> structs.GnomadDataset:
    frame = pl.DataFrame(
        {
            "gnomAD ID": ["22-1-A-G"],
            "Chromosome": [22],
            "Position": [1],
            "Gene": [gene],
            "CanonicalTranscript": [True],
        }
    )
    return structs.GnomadDataset(
        gnomadVersion="v4",
        gene=gene,
        ensemblId="ENSG00000133433.11",
        canonicalTranscript="ENST00000290765.9",
        chromosomes=[22],
        nVars=1,
        nVarsInCanonicalTranscript=1,
        source="test",
        endpoint_queried="https://gnomad.example/api",
        gnomad_dataset="gnomad_r4",
        reference_genome="GRCh38",
        dataMatrix=frame,
    )


def test_store_writes_metadata_and_variants_table():
    store = GnomadStore()
    try:
        dataset = _dataset()
        dataset.obs_mis = 25
        dataset.oe_mis = 0.495
        dataset.obs_syn = 18
        dataset.oe_syn = 0.9
        assert store.write_dataset(dataset) == "variants"
        genes = store.list_genes()
        assert genes["gene_name"].to_list() == ["GSTT2B"]
        assert genes["obs_mis"].to_list() == [25]
        assert genes["oe_mis"].to_list() == [0.495]
        assert genes["obs_syn"].to_list() == [18]
        assert genes["oe_syn"].to_list() == [0.9]
        assert store.read_gene("GSTT2B").shape == (1, 9)
        variants = store.connection.sql("SELECT * FROM variants").pl()
        assert variants.columns[-4:] == [
            "fetched_at",
            "endpoint_queried",
            "gnomad_dataset",
            "reference_genome",
        ]
        assert variants["endpoint_queried"].to_list() == ["https://gnomad.example/api"]
    finally:
        store.close()


def test_store_replaces_existing_gene_data():
    store = GnomadStore()
    try:
        store.write_dataset(_dataset())
        replacement = _dataset()
        replacement.dataMatrix = replacement.dataMatrix.vstack(replacement.dataMatrix)
        replacement.nVars = 2
        replacement.nVarsInCanonicalTranscript = 2
        store.write_dataset(replacement)
        assert store.read_gene("GSTT2B").height == 2
        assert store.list_genes()["n_variants"].to_list() == [2]
    finally:
        store.close()


def test_store_keeps_multiple_genes_in_one_variants_table():
    store = GnomadStore()
    try:
        store.write_dataset(_dataset("GSTT2B"))
        store.write_dataset(_dataset("GSTM5"))
        cursor = store.connection.cursor()
        try:
            assert cursor.execute("SELECT COUNT(*) FROM variants").fetchone() == (2,)
        finally:
            cursor.close()
        assert store.read_gene("GSTT2B").height == 1
        assert store.read_gene("GSTM5").height == 1
    finally:
        store.close()


def test_verify_database_reports_only_shared_value_disagreements(tmp_path: Path):
    csv_dir = tmp_path / "gnomad"
    csv_dir.mkdir()
    csv_path = csv_dir / "gnomAD_v4.1.2_ENSG00000133433_2026_11_20.csv"
    csv_path.write_text(
        "gnomAD ID,Chromosome,Position,Reference,Alternate,Extra\n"
        "22-1-A-G,22,1,A,G,csv-only\n"
        "22-2-A-C,22,999,A,C,csv-only\n",
    )

    dataset = _dataset()
    dataset.gnomadVersion = "gnomad_r4"
    dataset.dataMatrix = pl.DataFrame(
        {
            "gnomAD ID": ["22-1-A-G", "22-2-A-C", "22-3-A-T"],
            "Chromosome": [22, 22, 22],
            "Position": [1, 2, 3],
            "Reference": ["A", "A", "A"],
            "Alternate": ["G", "C", "T"],
        }
    )
    dataset.nVars = 3

    database = tmp_path / "gnomad.duckdb"
    store = GnomadStore(database)
    try:
        store.write_dataset(dataset)
    finally:
        store.close()

    result = verify_database(database, csv_dir)

    assert result.checked_tables == ("variants",)
    assert result.checked_csvs == (csv_path,)
    assert result.compared_rows == 2
    assert result.compared_values == 8
    assert len(result.disagreements) == 1
    disagreement = result.disagreements[0]
    assert disagreement.gnomad_id == "22-2-A-C"
    assert disagreement.column == "Position"
    assert disagreement.table_value == 2
    assert disagreement.csv_value == 999


def test_rebuild_database_creates_snapshot_and_latest(tmp_path: Path):
    snapshot, latest = rebuild_database(
        [_dataset()],
        version="v4-test",
        database_dir=tmp_path,
    )

    assert snapshot == tmp_path / "snapshots" / "gnomad_v4-test.duckdb"
    assert latest == tmp_path / "gnomad_latest.duckdb"
    assert snapshot.is_file()
    assert latest.is_file()
    assert latest.is_symlink()
    assert latest.resolve() == snapshot

    store = GnomadStore(latest)
    try:
        assert store.list_genes()["database_version"].to_list() == ["v4-test"]
        table_names = {
            row[0]
            for row in store.connection.execute(
                "SELECT table_name FROM information_schema.tables "
                "WHERE table_schema = 'main'"
            ).fetchall()
        }
        assert "variants" in table_names
        assert not any(name.startswith(("gene_", "sv_")) for name in table_names)
    finally:
        store.close()


def test_failed_verification_preserves_latest_and_wal(tmp_path: Path, monkeypatch):
    latest = tmp_path / "gnomad_latest.duckdb"
    latest.write_bytes(b"existing latest")
    latest_wal = tmp_path / "gnomad_latest.duckdb.wal"
    latest_wal.write_bytes(b"existing wal")

    real_store = GnomadStore

    class UnverifiedStore:
        def __init__(
            self,
            database=":memory:",
            database_version=None,
            read_only=False,
        ):
            self.store = real_store(database, database_version, read_only)
            self.read_only = read_only

        def close(self):
            self.store.close()

        def write_dataset(self, dataset):
            return self.store.write_dataset(dataset)

        def list_genes(self):
            if self.read_only:
                return pl.DataFrame({"gene_name": ["WRONG"]})
            return self.store.list_genes()

        def read_gene(self, gene):
            return self.store.read_gene(gene)

    monkeypatch.setattr(database_module, "GnomadStore", UnverifiedStore)

    with pytest.raises(RuntimeError, match="Database verification failed"):
        rebuild_database(
            [_dataset()],
            version="v4-failed",
            database_dir=tmp_path,
        )

    assert latest.read_bytes() == b"existing latest"
    assert latest_wal.read_bytes() == b"existing wal"
    assert not (tmp_path / "snapshots" / "gnomad_v4-failed.duckdb").exists()


def test_rebuild_database_removes_stale_latest_wal(tmp_path: Path):
    snapshot, latest = rebuild_database(
        [_dataset()],
        version="v4-old",
        database_dir=tmp_path,
    )
    assert snapshot.is_file()
    assert latest.is_file()

    latest_wal = tmp_path / "gnomad_latest.duckdb.wal"
    latest_wal.write_bytes(b"stale wal")

    rebuild_database(
        [_dataset("GSTM5")],
        version="v4-new",
        database_dir=tmp_path,
    )

    assert not latest_wal.exists()
    store = GnomadStore(latest)
    try:
        assert store.list_genes()["gene_name"].to_list() == ["GSTM5"]
    finally:
        store.close()


def test_build_database_fetches_configured_genes(tmp_path: Path):
    configs = [
        structs.GnomadGeneConfig(
            gene_name="GSTT1",
            ensembl_id="ENSG00000184674.8",
            gnomad_database_version="gnomad_r2_1",
            reference_genome="GRCh37",
        ),
        structs.GnomadGeneConfig(
            gene_name="GSTM5",
            ensembl_id="ENSG00000134201.12",
            gnomad_database_version="gnomad_r4",
        ),
    ]

    class FakeClient:
        def __init__(self):
            self.calls = []

        def fetch_gene(self, ensembl_id, **kwargs):
            self.calls.append((ensembl_id, kwargs))
            return _dataset(kwargs["gene_name"])

    client = FakeClient()
    snapshot, latest = build_database(
        configs,
        version="v-test",
        database_dir=tmp_path,
        client=client,
    )

    assert snapshot.is_file()
    assert latest.is_file()
    assert client.calls == [
        (
            "ENSG00000184674.8",
            {
                "dataset": "gnomad_r2_1",
                "reference_genome": "GRCh37",
                "gene_name": "GSTT1",
                "gnomad_version": "gnomad_r2_1",
            },
        ),
        (
            "ENSG00000134201.12",
            {
                "dataset": "gnomad_r4",
                "reference_genome": "GRCh38",
                "gene_name": "GSTM5",
                "gnomad_version": "gnomad_r4",
            },
        ),
    ]
