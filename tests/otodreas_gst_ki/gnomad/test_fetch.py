import polars as pl
from gql.transport.exceptions import TransportQueryError

from otodreas_gst_ki.gnomad.fetch import GnomadGraphQLClient


def test_fetch_gene_normalizes_graphql_results():
    responses = [
        {
            "gene": {
                "gene_id": "ENSG00000133433",
                "gene_version": "11",
                "symbol": "GSTT2B",
                "chrom": "22",
                "start": 1,
                "stop": 2,
                "canonical_transcript_id": "ENST00000290765.9",
            }
        },
        {
            "region": {
                "variants": [
                    {
                        "variant_id": "22-1-A-G",
                        "chrom": "22",
                        "pos": 1,
                        "ref": "A",
                        "alt": "G",
                        "rsid": "rs1",
                        "gene_id": "ENSG00000133433",
                        "gene_symbol": "GSTT2B",
                        "transcript_id": "ENST00000290765.9",
                        "hgvsc": "c.1A>G",
                        "hgvsp": None,
                        "consequence": "missense_variant",
                        "consequence_in_canonical_transcript": True,
                        "flags": ["lc lof"],
                    }
                ]
            }
        },
    ]

    def execute(_document, *, variable_values):
        assert variable_values["referenceGenome"] == "GRCh38"
        return responses.pop(0)

    dataset = GnomadGraphQLClient(execute=execute).fetch_gene("ENSG00000133433.11")

    assert dataset.gene == "GSTT2B"
    assert dataset.ensemblId == "ENSG00000133433.11"
    assert dataset.dataMatrix.select("Chromosome").to_series().to_list() == [22]
    assert dataset.dataMatrix.select("CanonicalTranscript").to_series().to_list() == [
        True
    ]
    assert isinstance(dataset.dataMatrix, pl.DataFrame)


def test_fetch_gene_splits_region_when_api_reports_too_many_variants():
    responses = iter(
        [
            {
                "gene": {
                    "gene_id": "ENSG00000134201",
                    "gene_version": "12",
                    "symbol": "GSTM5",
                    "chrom": "1",
                    "start": 100,
                    "stop": 104,
                    "canonical_transcript_id": "ENST00000256593",
                }
            }
        ]
    )
    calls = []

    def execute(document, *, variable_values):
        if "geneId" in variable_values:
            return next(responses)
        calls.append(variable_values)
        if variable_values["start"] == 100 and variable_values["stop"] == 104:
            raise TransportQueryError("This region has too many variants to display.")
        midpoint = (variable_values["start"], variable_values["stop"])
        variant_id = f"{midpoint[0]}-{midpoint[1]}"
        return {
            "region": {
                "variants": [
                    {
                        "variant_id": variant_id,
                        "chrom": "1",
                        "pos": midpoint[0],
                        "ref": "A",
                        "alt": "G",
                        "rsid": None,
                        "gene_id": "ENSG00000134201",
                        "gene_symbol": "GSTM5",
                        "transcript_id": "ENST00000256593",
                        "hgvsc": None,
                        "hgvsp": None,
                        "consequence": "missense_variant",
                        "consequence_in_canonical_transcript": True,
                        "flags": [],
                    }
                ]
            }
        }

    dataset = GnomadGraphQLClient(execute=execute).fetch_gene("ENSG00000134201.12")

    assert dataset.gene == "GSTM5"
    assert dataset.nVars == 2
    assert len(calls) == 3
