# GST-KI

Analytical code for KI GST project. Master's thesis, Lund University MSc in Bioinformatics.

## Reproducibility

### Install required software

Currently, the repo is a Python-only uv project. Required software is listed below:

- [uv](https://docs.astral.sh/uv/getting-started/installation/)

> Note: the appropriate Python version will be installed by uv if it's not already present when you initialize the environment.

### Build the environment

1. [Install uv](https://docs.astral.sh/uv/getting-started/installation/)
2. Fork the repo on GitHub
3. Clone the repo
```sh
cd ~
git clone git@github.com:<YOUR_GITHUB_USERNAME>/GST-KI.git
cd GST-KI
```
4. Initialize the environment
```sh
uv sync --locked
```
### Run tests

If proposing changes, run the test suite:

```sh
uv run pytest
```

Add at least one test for every unit added.