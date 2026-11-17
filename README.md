# GST-KI

Analytical code for KI GST project. Master's thesis, Lund University MSc in Bioinformatics.

# Reproducibility

## Reproduce results

Currently, the only way to reproduce the results is to install the project and run the notebooks, which is outlined in the [Audit code](audit-code) section.

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