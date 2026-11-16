# GST-KI

Analytical code for KI GST project. Master's thesis, Lund University MSc in Bioinformatics.

# Reproducibility

## Reproduce results

I will provide a Jupyter Notebook via Colab (not published yet). It will pip install this package and let you run my analyses quickly.

## Audit code

If you prefer to audit the code, work on this repo as a developer.

### Install required software

Currently, the repo is a Python-only uv project. Required software is listed below:

- [uv](https://docs.astral.sh/uv/getting-started/installation/)

> Note: the appropriate Python version will be installed by uv if it's not already present when you initialize the environment.

### Build the environment

1. Fork the repo on GitHub
2. Clone the repo
```sh
cd ~
git clone git@github.com:<YOUR_GITHUB_USERNAME>/GST-KI.git
cd GST-KI
```
3. Initialize the environment
```sh
uv sync --locked
```
### Run tests

If proposing changes, run the test suite:

```sh
uv run pytest
```

Add at least one test for every unit added.