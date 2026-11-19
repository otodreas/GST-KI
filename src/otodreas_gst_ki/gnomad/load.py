import warnings
from pathlib import Path

import polars as pl
from pydantic import BaseModel, ConfigDict

from otodreas_gst_ki import structs


class GnomadDataset(BaseModel):
    gnomadVersion: str
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


def _check_csv(filename: str | Path) -> Path:
    filename = Path(filename)

    if not filename.is_file():
        raise FileNotFoundError(f"{filename} not found")

    # Read the first line to check for commas
    with open(filename, "r") as head:
        head = head.readline()

    if filename.suffix != ".csv" or "," not in head:
        raise TypeError(f"{filename} is not a CSV file")

    return filename


def load_gnomad(
    filename: str | Path,
    ensembl_to_gene_lookup: dict[str, str] = structs.ensembl_to_gene,
    gene_to_ensembl_lookup: dict[str, dict[str, str]] = structs.gene_to_ensembl,
    gnomad_schema: dict[str, pl.Any] = structs.gnomad_v4_schema,
) -> GnomadDataset:
    """Load a gnomAD CSV. Check for provenance between fields derived from hand-written lookup tables
    Returns a GnomadDataset object
    """

    filename = _check_csv(filename)

    filenameinfo = filename.stem.split("_")
    gnomad_version = filenameinfo[1]
    if not gnomad_version.startswith("v4."):
        warnings.warn(
            "gnomAD version < 4 detected. Matrix will be matched to gnomAD v4 schema"
        )

    # Extract short ensembl id
    ensembl_id_short = filenameinfo[2]

    if ensembl_id_short == [filename.stem]:
        raise TypeError(f"Failed to extract Ensembl ID from filename {filename}")

    # Lookups
    try:
        gene = ensembl_to_gene_lookup[ensembl_id_short]
    except KeyError:
        print(f"Ensembl ID {ensembl_id_short} not found in lookup table")

    try:
        ensembl_id_long = gene_to_ensembl_lookup[gene]["geneId"]
    except KeyError:
        print(f"Gene {gene} not present in lookup table or missing geneID key")
    try:
        canonical_transcript = gene_to_ensembl_lookup[gene]["canonicalTranscript"]
    except KeyError:
        print(
            f"Gene {gene} not present in lookup table or missing canonicalTranscript key"
        )

    if ensembl_id_short != ensembl_id_long.split(".")[0]:
        raise ValueError(
            f"Provenance failure for {gene}: {ensembl_id_short} not like {ensembl_id_long}"
        )

    # Currently, all data must follow v4 schema
    if gnomad_version.startswith("v4."):
        df = pl.read_csv(filename, schema_overrides=gnomad_schema)
    else:
        df = pl.read_csv(
            filename, schema_overrides={"ClinVar Variation ID": pl.String}
        )  # always force clinvar variant id
        df = df.match_to_schema(
            gnomad_schema, missing_columns="insert", extra_columns="ignore"
        )

    df = df.with_columns(Gene=pl.lit(gene, dtype=pl.String))
    df = df.with_columns(
        pl.when(pl.col("VEP Annotation") == "intron_variant")
        .then(False)
        .otherwise(True)
        .alias("Exonic")
    )

    vep_anns = df.filter(pl.col("Exonic")).get_column("VEP Annotation")
    subset = vep_anns.filter(vep_anns.str.contains("(?i)intron"))
    if not subset.is_empty():
        warnings.warn(f"Multiple VEP Annotations containing 'intron' found: {subset}")

    if df.get_column("Chromosome").n_unique() > 1:
        warnings.warn(f"{filename} contains variants from multiple chromosomes")

    n_vars = df.shape[0]
    df = df.with_columns(
        CanonicalTranscript=pl.col("Transcript") == canonical_transcript
    )
    n_vars_in_canonical_transcript = df.get_column("CanonicalTranscript").sum()

    return GnomadDataset(
        gnomadVersion=gnomad_version,
        gene=gene,
        ensemblId=ensembl_id_long,
        canonicalTranscript=canonical_transcript,
        chromosomes=list(df.get_column("Chromosome").unique()),
        nVars=n_vars,
        nVarsInCanonicalTranscript=n_vars_in_canonical_transcript,
        filePath=filename,
        dataMatrix=df,
    )


def concat_gnomad(datasets: list[GnomadDataset]) -> pl.DataFrame:

    if len(datasets) < 2:
        raise ValueError("Pass a list of at least 2 datasets to concatenate")

    return pl.concat([d.dataMatrix for d in datasets])


def find_shared_variants(datasets: list[GnomadDataset]) -> pl.DataFrame:

    combined = concat_gnomad(datasets).with_columns(
        VariantId=pl.concat_str(
            ["Chromosome", "Position", "Reference", "Alternate"], separator="-"
        )
    )

    shared = (
        combined.group_by("VariantId")
        .agg(pl.col("Gene").n_unique().alias("nGenes"))
        .filter(pl.col("nGenes") > 1)
    )

    if shared.shape[0] > 0:
        warnings.warn("Overlapping variants found")

    return shared
