import polars as pl

# These columns can have unclear types a long way down
gnomad_schema_overrides = {"ClinVar Variation ID": pl.String}

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
}

# Reverse lookup of gene names from short Ensembl gene IDs
ensembl_to_gene = {
    "ENSG00000134184": "GSTM1",
    "ENSG00000134202": "GSTM3",
    "ENSG00000134201": "GSTM5",
    "ENSG00000148834": "GSTO1",
    "ENSG00000133433": "GSTT2B",
}
