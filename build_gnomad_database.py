# %%
"""Fetch and rebuild the gnomAD DuckDB database.

Run this file in IPython/Jupyter. Each # %% section can be executed
independently, so the latest database can be read without rebuilding it.
"""

from datetime import UTC, datetime
from pathlib import Path

from otodreas_gst_ki.gnomad import structs
from otodreas_gst_ki.gnomad.database import build_database
from otodreas_gst_ki.gnomad.store import GnomadStore

DATABASE_DIR = Path("data/gnomad")
LATEST_DATABASE = DATABASE_DIR / "gnomad_latest.duckdb"
GENE_CONFIGS = [
    structs.GnomadGeneConfig(
        gene_name=gene_name,
        ensembl_id=gene_info["geneId"],
        gnomad_database_version=(
            "gnomad_r2_1" if gene_name == "GSTT1" else "gnomad_r4"
        ),
        reference_genome="GRCh37" if gene_name == "GSTT1" else "GRCh38",
    )
    for gene_name, gene_info in structs.gene_to_ensembl.items()
]


# %% Read the latest database
store = GnomadStore(LATEST_DATABASE)

genes = store.list_genes()
print(genes)

# Run this in a later cell to read one gene:
# gst_variants = store.read_gene("GSTT2B")


# %% Build a new database version
database_version = datetime.now(UTC).strftime("v4-%Y-%m-%d")
snapshot, latest = build_database(
    GENE_CONFIGS,
    version=database_version,
    database_dir=DATABASE_DIR,
    pause_seconds=30,
)

print(f"Snapshot: {snapshot}")
print(f"Latest: {latest}")
# %%
