# %%
"""Fetch and rebuild the gnomAD DuckDB database.

Run this file in IPython/Jupyter. Each # %% section can be executed
independently, so the latest database can be read without rebuilding it.
"""

from datetime import UTC, datetime
from pathlib import Path
from time import sleep

from otodreas_gst_ki.gnomad import structs
from otodreas_gst_ki.gnomad.database import rebuild_database
from otodreas_gst_ki.gnomad.fetch import GnomadGraphQLClient
from otodreas_gst_ki.gnomad.store import GnomadStore

DATABASE_DIR = Path("data/gnomad")
LATEST_DATABASE = DATABASE_DIR / "gnomad_latest.duckdb"
GNOMAD_DATASET = "gnomad_r4"
REFERENCE_GENOME = "GRCh38"


# %% Read the latest database
store = GnomadStore(LATEST_DATABASE)

genes = store.list_genes()
print(genes)

# Run this in a later cell to read one gene:
# gst_variants = store.read_gene("GSTT2B")


# %% Build a new database version
database_version = datetime.now(UTC).strftime("v4-%Y-%m-%d")
client = GnomadGraphQLClient()
datasets = []

for gene_name, gene_info in structs.gene_to_ensembl.items():
    ensembl_id = gene_info["geneId"]
    print(f"Fetching {gene_name} ({ensembl_id})")
    datasets.append(
        client.fetch_gene(
            ensembl_id,
            dataset=GNOMAD_DATASET,
            reference_genome=REFERENCE_GENOME,
            gene_name=gene_name,
            canonical_transcript=gene_info["canonicalTranscript"],
        )
    )
    print(f"Finished fetching {gene_name}")
    t = 30
    print(f"sleeping for {t} seconds at {datetime.now().time().strftime('%H:%M:%S')}")
    sleep(t)

snapshot, latest = rebuild_database(
    datasets,
    version=database_version,
    database_dir=DATABASE_DIR,
)

print(f"Snapshot: {snapshot}")
print(f"Latest: {latest}")
# %%
