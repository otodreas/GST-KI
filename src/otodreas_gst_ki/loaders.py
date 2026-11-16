import warnings
from pathlib import Path

import polars as pl
from pydantic import BaseModel, ConfigDict

from otodreas_gst_ki import structs


class GnomadDataset(BaseModel):
    gene: str
    ensemblId: str
    canonicalTranscript: str
    chromosomes: list[int]
    nVars: int
    nVarsInCanonicalTranscript: int
    filePath: Path
    dataMatrix: pl.DataFrame

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def summarize(self):
        print(
            f"{self.gene}: {self.nVarsInCanonicalTranscript} variants are in the canonical transcript {self.canonicalTranscript} ({round(self.nVarsInCanonicalTranscript / self.nVars * 100, 1)}%)"
        )


def load_gnomad(
    filename: str | Path,
    ensembl_to_gene_lookup: dict[str, str] = structs.ensembl_to_gene,
    gene_to_ensembl_lookup: dict[str, dict[str, str]] = structs.gene_to_ensembl,
    gnomad_schema_overrides: dict[str, pl.Any] = structs.gnomad_schema_overrides,
) -> GnomadDataset:
    filename = Path(filename)

    if not filename.is_file():
        raise FileNotFoundError(f"{filename} not found")

    # Read the first line to check for commas
    with open(filename, "r") as head:
        head = head.readline()

    if filename.suffix != ".csv" or "," not in head:
        raise TypeError(f"{filename} is not a CSV file")

    ensembl_id_short = filename.stem.split("_")[2]

    if ensembl_id_short == [filename.stem]:
        raise TypeError(f"Failed to extract Ensembl ID from filename {filename}")

    # Lookups
    try:
        gene = ensembl_to_gene_lookup[ensembl_id_short]
    except KeyError:
        print(f"Ensembl ID {ensembl_id_short} not found in lookup table")

    try:
        ensembl_id_long = gene_to_ensembl_lookup[gene]["geneId"]
        canonical_transcript = gene_to_ensembl_lookup[gene]["canonicalTranscript"]
    except KeyError:
        print(
            f"Gene {gene} not present in lookup table or missing geneID and canonicalTranscript keys"
        )

    df = pl.read_csv(filename, schema_overrides=gnomad_schema_overrides)
    df = df.with_columns(Gene=pl.lit(gene, dtype=pl.String))

    if df.get_column("Chromosome").n_unique() > 1:
        warnings.warn(f"{filename} contains variants from multiple chromosomes")

    n_vars = df.shape[0]
    df = df.with_columns(
        CanonicalTranscript=pl.col("Transcript") == canonical_transcript
    )
    n_vars_in_canonical_transcript = df.get_column("CanonicalTranscript").sum()

    return GnomadDataset(
        gene=gene,
        ensemblId=ensembl_id_long,
        canonicalTranscript=canonical_transcript,
        chromosomes=list(df["Chromosome"].unique()),
        nVars=n_vars,
        nVarsInCanonicalTranscript=n_vars_in_canonical_transcript,
        filePath=filename.resolve(),
        dataMatrix=df,
    )


def concat_gnomad(datasets: list[GnomadDataset]) -> pl.DataFrame:

    if len(datasets) < 2:
        raise ValueError("Pass a list of at least 2 datasets to concatenate")

    return pl.concat([d.dataMatrix for d in datasets])
