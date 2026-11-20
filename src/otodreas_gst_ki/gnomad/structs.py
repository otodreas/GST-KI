from pathlib import Path

import polars as pl
from pydantic import BaseModel, ConfigDict


class GnomadDataset(BaseModel):
    gnomadVersion: str
    gene: str
    ensemblId: str
    canonicalTranscript: str
    chromosomes: list[int]
    nVars: int
    nVarsInCanonicalTranscript: int
    filePath: Path | None = None
    source: str | None = None
    dataMatrix: pl.DataFrame

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def summarize(self):
        percentage = (
            self.nVarsInCanonicalTranscript / self.nVars * 100 if self.nVars else 0
        )
        print(
            f"{self.gene}: {self.nVarsInCanonicalTranscript} variants are in the "
            f"canonical transcript {self.canonicalTranscript} ({percentage:.1f}%)"
        )


gnomad_v4_schema = {
    "gnomAD ID": pl.String,
    "Chromosome": pl.Int64,
    "Position": pl.Int64,
    "rsIDs": pl.String,
    "Reference": pl.String,
    "Alternate": pl.String,
    "Source": pl.String,
    "Filters - exomes": pl.String,
    "Filters - genomes": pl.String,
    "Transcript": pl.String,
    "HGVS Consequence": pl.String,
    "Protein Consequence": pl.String,
    "Transcript Consequence": pl.String,
    "VEP Annotation": pl.String,
    "ClinVar Germline Classification": pl.String,
    "ClinVar Variation ID": pl.String,
    "Flags": pl.String,
    "Allele Count": pl.Int64,
    "Allele Number": pl.Int64,
    "Allele Frequency": pl.Float64,
    "Homozygote Count": pl.Int64,
    "Hemizygote Count": pl.Int64,
    "Filters - joint": pl.String,
    "GroupMax FAF group": pl.String,
    "GroupMax FAF frequency": pl.Float64,
    "cadd": pl.Float64,
    "revel_max": pl.Float64,
    "spliceai_ds_max": pl.Float64,
    "pangolin_largest_ds": pl.Float64,
    "phylop": pl.Float64,
    "sift_max": pl.Float64,
    "polyphen_max": pl.Float64,
    "Allele Count African/African American": pl.Int64,
    "Allele Number African/African American": pl.Int64,
    "Homozygote Count African/African American": pl.Int64,
    "Hemizygote Count African/African American": pl.Int64,
    "Allele Count Admixed American": pl.Int64,
    "Allele Number Admixed American": pl.Int64,
    "Homozygote Count Admixed American": pl.Int64,
    "Hemizygote Count Admixed American": pl.Int64,
    "Allele Count Ashkenazi Jewish": pl.Int64,
    "Allele Number Ashkenazi Jewish": pl.Int64,
    "Homozygote Count Ashkenazi Jewish": pl.Int64,
    "Hemizygote Count Ashkenazi Jewish": pl.Int64,
    "Allele Count East Asian": pl.Int64,
    "Allele Number East Asian": pl.Int64,
    "Homozygote Count East Asian": pl.Int64,
    "Hemizygote Count East Asian": pl.Int64,
    "Allele Count European (Finnish)": pl.Int64,
    "Allele Number European (Finnish)": pl.Int64,
    "Homozygote Count European (Finnish)": pl.Int64,
    "Hemizygote Count European (Finnish)": pl.Int64,
    "Allele Count Middle Eastern": pl.Int64,
    "Allele Number Middle Eastern": pl.Int64,
    "Homozygote Count Middle Eastern": pl.Int64,
    "Hemizygote Count Middle Eastern": pl.Int64,
    "Allele Count European (non-Finnish)": pl.Int64,
    "Allele Number European (non-Finnish)": pl.Int64,
    "Homozygote Count European (non-Finnish)": pl.Int64,
    "Hemizygote Count European (non-Finnish)": pl.Int64,
    "Allele Count Amish": pl.Int64,
    "Allele Number Amish": pl.Int64,
    "Homozygote Count Amish": pl.Int64,
    "Hemizygote Count Amish": pl.Int64,
    "Allele Count South Asian": pl.Int64,
    "Allele Number South Asian": pl.Int64,
    "Homozygote Count South Asian": pl.Int64,
    "Hemizygote Count South Asian": pl.Int64,
    "Allele Count Remaining": pl.Int64,
    "Allele Number Remaining": pl.Int64,
    "Homozygote Count Remaining": pl.Int64,
    "Hemizygote Count Remaining": pl.Int64,
}

# Full Ensembl gene IDs and Ensembl canonical transcripts
gene_to_ensembl = {
    "GSTM1": {
        "geneId": "ENSG00000134184.13",
        "canonicalTranscript": "ENST00000309851.10",
    },
    "GSTM3": {
        "geneId": "ENSG00000134202.12",
        "canonicalTranscript": "ENST00000361066.7",
    },
    "GSTM5": {
        "geneId": "ENSG00000134201.12",
        "canonicalTranscript": "ENST00000256593.8",
    },
    "GSTO1": {
        "geneId": "ENSG00000148834.13",
        "canonicalTranscript": "ENST00000369713.10",
    },
    "GSTT2B": {
        "geneId": "ENSG00000133433.11",
        "canonicalTranscript": "ENST00000290765.9",
    },
    "GSTA1": {
        "geneId": "ENSG00000243955.6",
        "canonicalTranscript": "ENST00000334575.6",
    },
    "GSTA3": {
        "geneId": "ENSG00000174156.15",
        "canonicalTranscript": "ENST00000211122.4",
    },
    "GSTA4": {
        "geneId": "ENSG00000170899.11",
        "canonicalTranscript": "ENST00000370963.9",
    },
    "GSTA5": {
        "geneId": "ENSG00000182793.12",
        "canonicalTranscript": "ENST00000370989.6",
    },
    # "GSTK1": {
    #     "geneId": "ENSG00000197448.14",
    #     "canonicalTranscript": "ENST00000358406.10",
    # },
    "GSTM2": {
        "geneId": "ENSG00000213366.13",
        "canonicalTranscript": "ENST00000241337.9",
    },
    "GSTM4": {
        "geneId": "ENSG00000168765.18",
        "canonicalTranscript": "ENST00000369836.9",
    },
    "GSTO2": {
        "geneId": "ENSG00000065621.15",
        "canonicalTranscript": "ENST00000338595.7",
    },
    "GSTP1": {
        "geneId": "ENSG00000084207.18",
        "canonicalTranscript": "ENST00000398606.10",
    },
    "GSTT2": {
        "geneId": "ENSG00000099984.11",
        "canonicalTranscript": "ENST00000621118.4",
    },
    "GSTZ1": {
        "geneId": "ENSG00000100577.19",
        "canonicalTranscript": "ENST00000216465.10",
    },
    "GSTT1": {
        "geneId": "ENSG00000184674.8",
        "canonicalTranscript": "ENST00000248935.5",
    },
}

# Reverse lookup of gene names from short Ensembl gene IDs
ensembl_to_gene = {
    "ENSG00000134184": "GSTM1",
    "ENSG00000134202": "GSTM3",
    "ENSG00000134201": "GSTM5",
    "ENSG00000148834": "GSTO1",
    "ENSG00000133433": "GSTT2B",
    "ENSG00000243955": "GSTA1",
    "ENSG00000174156": "GSTA3",
    "ENSG00000170899": "GSTA4",
    "ENSG00000182793": "GSTA5",
    # "ENSG00000197448": "GSTK1",
    "ENSG00000213366": "GSTM2",
    "ENSG00000168765": "GSTM4",
    "ENSG00000065621": "GSTO2",
    "ENSG00000084207": "GSTP1",
    "ENSG00000099984": "GSTT2",
    "ENSG00000100577": "GSTZ1",
    "ENSG00000184674": "GSTT1",
}
