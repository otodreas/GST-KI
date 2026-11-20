from pathlib import Path

import polars as pl
import pytest

from otodreas_gst_ki.gnomad import structs
from otodreas_gst_ki.gnomad.database import build_database, rebuild_database
from otodreas_gst_ki.gnomad.store import GnomadStore, gene_table_name, sv_table_name


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
        dataMatrix=frame,
    )


def test_gene_table_name_uses_prefix():
    assert gene_table_name("GSTT2B") == "gene_GSTT2B"


def test_gene_table_name_rejects_identifiers():
    with pytest.raises(ValueError):
        gene_table_name("GST; DROP TABLE genes")


def test_sv_table_name_uses_prefix():
    assert sv_table_name("GSTT2B") == "sv_GSTT2B"


def test_store_writes_metadata_and_gene_table():
    store = GnomadStore()
    try:
        dataset = _dataset()
        dataset.obs_mis = 25
        dataset.oe_mis = 0.495
        dataset.obs_syn = 18
        dataset.oe_syn = 0.9
        assert store.write_dataset(dataset) == "gene_GSTT2B"
        genes = store.list_genes()
        assert genes["gene_name"].to_list() == ["GSTT2B"]
        assert genes["obs_mis"].to_list() == [25]
        assert genes["oe_mis"].to_list() == [0.495]
        assert genes["obs_syn"].to_list() == [18]
        assert genes["oe_syn"].to_list() == [0.9]
        assert store.read_gene("GSTT2B").shape == (1, 5)
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

    store = GnomadStore(latest)
    try:
        assert store.list_genes()["database_version"].to_list() == ["v4-test"]
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
