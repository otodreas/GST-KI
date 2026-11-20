import re
import uuid
from pathlib import Path
from time import sleep

from otodreas_gst_ki.gnomad.fetch import GnomadGraphQLClient
from otodreas_gst_ki.gnomad.store import GnomadStore
from otodreas_gst_ki.gnomad.structs import GnomadDataset, GnomadGeneConfig


def build_database(
    gene_configs: list[GnomadGeneConfig],
    version: str,
    database_dir: str | Path = "data/gnomad",
    latest_name: str = "gnomad_latest.duckdb",
    client: GnomadGraphQLClient | None = None,
    pause_seconds: float = 0,
) -> tuple[Path, Path]:
    """Fetch configured genes and rebuild the versioned DuckDB database."""
    client = client or GnomadGraphQLClient()
    datasets = []
    for config in gene_configs:
        datasets.append(
            client.fetch_gene(
                config.ensembl_id,
                dataset=config.gnomad_database_version,
                reference_genome=config.reference_genome,
                gene_name=config.gene_name,
                gnomad_version=config.gnomad_database_version,
            )
        )
        if pause_seconds:
            sleep(pause_seconds)

    return rebuild_database(
        datasets,
        version=version,
        database_dir=database_dir,
        latest_name=latest_name,
    )


def rebuild_database(
    datasets: list[GnomadDataset],
    version: str,
    database_dir: str | Path = "data/gnomad",
    latest_name: str = "gnomad_latest.duckdb",
) -> tuple[Path, Path]:
    """Build a versioned DuckDB snapshot and atomically promote it as latest."""
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", version):
        raise ValueError(f"Invalid database version: {version!r}")
    database_dir = Path(database_dir)
    snapshots_dir = database_dir / "snapshots"
    snapshots_dir.mkdir(parents=True, exist_ok=True)

    snapshot_path = snapshots_dir / f"gnomad_{version}.duckdb"
    staging_snapshot = snapshots_dir / f".{snapshot_path.name}.{uuid.uuid4().hex}.tmp"
    latest_path = database_dir / latest_name
    staging_latest = database_dir / f".{latest_name}.{uuid.uuid4().hex}.tmp"

    try:
        store = GnomadStore(str(staging_snapshot), database_version=version)
        try:
            for dataset in datasets:
                store.write_dataset(dataset)
        finally:
            store.close()

        verifier = GnomadStore(str(staging_snapshot), read_only=True)
        try:
            gene_rows = verifier.list_genes()
            actual_genes = set(gene_rows["gene_name"].to_list())
            expected_genes = {dataset.gene for dataset in datasets}
            if actual_genes != expected_genes:
                raise RuntimeError(
                    "Database verification failed: "
                    f"expected genes {sorted(expected_genes)}, "
                    f"found {sorted(actual_genes)}"
                )
            for dataset in datasets:
                gene_rows = verifier.read_gene(dataset.gene)
                if gene_rows.height != dataset.nVars:
                    raise RuntimeError(
                        "Database verification failed: "
                        f"{dataset.gene} has {gene_rows.height} variants, "
                        f"expected {dataset.nVars}"
                    )
        finally:
            verifier.close()

        staging_snapshot.replace(snapshot_path)

        Path(f"{latest_path}.wal").unlink(missing_ok=True)
        staging_latest.symlink_to(Path("snapshots") / snapshot_path.name)
        staging_latest.replace(latest_path)
    finally:
        staging_snapshot.unlink(missing_ok=True)
        staging_latest.unlink(missing_ok=True)

    return snapshot_path, latest_path
