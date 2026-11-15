# GST-KI

Analytical code for KI GST project. Master's thesis, Lund University MSc in Bioinformatics

## Reproducibility

Currently, the repo is a Python-only uv project. To reproduce:

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

> Note: the appropriate Python version will be installed by uv if it's not already present