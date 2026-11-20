# %%
"""Fetch and rebuild the gnomAD DuckDB database.

Run this file in IPython/Jupyter. Each # %% section can be executed
independently, so the latest database can be read without rebuilding it.
"""

from datetime import UTC, datetime
from pathlib import Path

import polars as pl

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

# %% Interact with the database
# Read one gene and summarize its stored metadata.
gst_variants = store.read_gene("GSTT2B")
gst_metadata = genes.filter(pl.col("gene_name") == "GSTT2B")

print(gst_metadata)
print(gst_variants.head())

# Summarize canonical-transcript coverage for each stored gene.
coverage = genes.select(
    "gene_name",
    "n_variants",
    "n_variants_in_canonical_transcript",
).with_columns(
    (pl.col("n_variants_in_canonical_transcript") / pl.col("n_variants") * 100)
    .fill_nan(0)
    .round(1)
    .alias("canonical_coverage_percent")
)
print(coverage)

# Read all genes into a dictionary, matching the main notebook workflow.
gene_frames = {
    row["gene_name"]: store.read_gene(row["gene_name"]) for row in genes.to_dicts()
}

# Concatenate all gene variant tables for cross-gene analysis.
full_data = pl.concat(gene_frames.values(), how="diagonal")
print(full_data.head())

# Find variant IDs shared by more than one gene.
shared_variants = (
    full_data.select("gnomAD ID", "Gene")
    .unique()
    .group_by("gnomAD ID")
    .agg(
        pl.col("Gene").unique().sort().alias("Genes"),
        pl.len().alias("nGenes"),
    )
    .filter(pl.col("nGenes") > 1)
    .sort("nGenes", descending=True)
)
print(shared_variants)


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
