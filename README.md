# GST-KI

Analytical code for KI GST project. Master's thesis, Lund University MSc in Bioinformatics.

# Repo structure

This repo contains the Python library [otodreas-gst-ki](src/otodreas_gst_ki/), [data](data), and [tests](tests). The notebook [main.ipynb](main.ipynb) runs the main analysis using the functions surfaced in the library.

# Reproducibility

## Reproduce results

Currently, the only way to reproduce the results is to install the project and run the notebooks, which is outlined in the **Audit code** section below.

I will provide a Jupyter Notebook via Colab (not published yet). It will pip install this package and let you run my analyses quickly.

## Audit code

If you prefer to audit the code, work on this repo as a developer.

### Install required software

The repo is a Python-only uv project. Required software and download links are listed below:

- [uv](https://docs.astral.sh/uv/getting-started/installation/)

> **Note:** the appropriate Python version will be installed by uv if it's not already present when you initialize the environment.

### Build the environment

1. Fork the repo on GitHub
2. Clone and enter the repo
```sh
# Using SSH
git clone git@github.com:<YOUR_GITHUB_USERNAME>/GST-KI.git
cd GST-KI
```
3. Initialize the environment
```sh
uv sync --locked
```

### Run notebooks

To edit and run Jupyter notebooks, either open them in your IDE or run

```sh
uv run jupyter lab
```

And open your notebook of choice in the Jupyter frontend.

### Run tests

If proposing changes, run the test suite:

```sh
uv run pytest
```

Add at least one test for every unit added.

## gnomAD data storage

The library can fetch a gene directly from the gnomAD GraphQL API and store it
in DuckDB:

```python
from otodreas_gst_ki.gnomad.fetch import GnomadGraphQLClient
from otodreas_gst_ki.gnomad.store import GnomadStore

dataset = GnomadGraphQLClient().fetch_gene("ENSG00000133433.11")

store = GnomadStore("data/gnomad.duckdb")
store.write_dataset(dataset)
variants = store.read_gene("GSTT2B")
store.close()
```

To preserve a rebuildable database history, create a named snapshot and promote
it to the latest database:

```python
from otodreas_gst_ki.gnomad.database import rebuild_database

snapshot, latest = rebuild_database(
    [dataset],
    version="v4.1.2-2026-11-19",
)
```

This writes `data/gnomad/snapshots/gnomad_v4.1.2-2026-11-19.duckdb` and
`data/gnomad/gnomad_latest.duckdb`. The latest database is replaced only after
the complete snapshot is written. A later rebuild leaves earlier snapshots
available if gnomAD or the live API is unavailable.

Each gene is stored in a table named `gene_<gene symbol>`, while the `genes`
table stores gene metadata and provenance. The default API dataset is gnomAD
v4 on GRCh38; both can be overridden when fetching.