from pathlib import Path

import polars as pl
import pytest

from otodreas_gst_ki.gnomad import structs
from otodreas_gst_ki.gnomad.database import rebuild_database
from otodreas_gst_ki.gnomad.store import GnomadStore, gene_table_name


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


def test_store_writes_metadata_and_gene_table():
    store = GnomadStore()
    try:
        assert store.write_dataset(_dataset()) == "gene_GSTT2B"
        assert store.list_genes()["gene_name"].to_list() == ["GSTT2B"]
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
