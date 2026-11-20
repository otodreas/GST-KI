from collections.abc import Callable
from typing import Any

import polars as pl
from gql import Client, gql
from gql.transport.exceptions import TransportQueryError
from gql.transport.requests import RequestsHTTPTransport

from otodreas_gst_ki.gnomad import structs

GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"

GENE_REGION_QUERY = gql(
    """
    query GeneRegion(
        $geneId: String!
        $referenceGenome: ReferenceGenomeId!
    ) {
        gene(gene_id: $geneId, reference_genome: $referenceGenome) {
            gene_id
            gene_version
            symbol
            chrom
            start
            stop
            canonical_transcript_id
            gnomad_constraint {
                exp_lof
                exp_mis
                exp_syn
                obs_lof
                obs_mis
                obs_syn
                oe_lof
                oe_lof_lower
                oe_lof_upper
                oe_lof_percentile
                oe_mis
                oe_mis_lower
                oe_mis_upper
                oe_syn
                oe_syn_lower
                oe_syn_upper
                lof_z
                mis_z
                syn_z
                pli
                flags
            }
        }
    }
    """
)

REGION_VARIANTS_QUERY = gql(
    """
    query RegionVariants(
        $chrom: String!
        $start: Int!
        $stop: Int!
        $referenceGenome: ReferenceGenomeId!
        $dataset: DatasetId!
    ) {
        region(
            chrom: $chrom
            start: $start
            stop: $stop
            reference_genome: $referenceGenome
        ) {
            variants(dataset: $dataset) {
                variant_id
                chrom
                pos
                ref
                alt
                rsid
                gene_id
                gene_symbol
                transcript_id
                hgvsc
                hgvsp
                consequence
                consequence_in_canonical_transcript
                flags
            }
        }
    }
    """
)


class GnomadGraphQLClient:
    def __init__(
        self,
        endpoint: str = GNOMAD_API_URL,
        execute: Callable[..., dict[str, Any]] | None = None,
    ):
        self._execute = execute
        if execute is None:
            transport = RequestsHTTPTransport(
                url=endpoint,
                use_json=True,
            )
            self._client = Client(
                transport=transport, fetch_schema_from_transport=False
            )
        else:
            self._client = None

    def _execute_query(
        self,
        document: Any,
        variables: dict[str, Any],
        *,
        database: str | None = None,
    ) -> dict[str, Any]:
        try:
            if self._execute is not None:
                return self._execute(document, variable_values=variables)
            assert self._client is not None
            with self._client as client:
                return client.execute(document, variable_values=variables)
        except TransportQueryError as error:
            if "gene not found" in str(error).lower() and "geneId" in variables:
                raise ValueError(
                    f"Gene {variables['geneId']} not found in database "
                    f"{database or 'unknown'} on reference genome "
                    f"{variables.get('referenceGenome', 'unknown')}"
                ) from error
            raise

    def _fetch_region_variants(
        self,
        *,
        chrom: str,
        start: int,
        stop: int,
        reference_genome: str,
        dataset: str,
    ) -> list[dict[str, Any]]:
        variables = {
            "chrom": chrom,
            "start": start,
            "stop": stop,
            "referenceGenome": reference_genome,
            "dataset": dataset,
        }
        try:
            result = self._execute_query(REGION_VARIANTS_QUERY, variables)
        except TransportQueryError as error:
            message = str(error)
            if "too many variants" not in message.lower():
                raise
            if stop - start <= 1:
                raise

            midpoint = start + (stop - start) // 2
            left = self._fetch_region_variants(
                chrom=chrom,
                start=start,
                stop=midpoint,
                reference_genome=reference_genome,
                dataset=dataset,
            )
            right = self._fetch_region_variants(
                chrom=chrom,
                start=midpoint,
                stop=stop,
                reference_genome=reference_genome,
                dataset=dataset,
            )

            variants_by_id = {
                variant["variant_id"]: variant for variant in (*left, *right)
            }
            return list(variants_by_id.values())

        return result.get("region", {}).get("variants", [])

    def fetch_gene(
        self,
        gene_id: str,
        *,
        dataset: str = "gnomad_r4",
        reference_genome: str = "GRCh38",
        gene_name: str | None = None,
        canonical_transcript: str | None = None,
        gnomad_version: str = "v4",
    ) -> structs.GnomadDataset:
        gene_result = self._execute_query(
            GENE_REGION_QUERY,
            {
                "geneId": gene_id.split(".")[0],
                "referenceGenome": reference_genome,
            },
            database=dataset,
        )
        gene = gene_result.get("gene")
        if gene is None:
            raise ValueError(f"gnomAD returned no gene for {gene_id}")

        actual_gene_id = gene["gene_id"]
        if actual_gene_id != gene_id.split(".")[0]:
            raise ValueError(
                f"gnomAD returned {actual_gene_id} for requested gene {gene_id}"
            )
        gene_name = gene_name or gene["symbol"]
        canonical_transcript = canonical_transcript or gene["canonical_transcript_id"]
        variants = self._fetch_region_variants(
            chrom=gene["chrom"],
            start=gene["start"],
            stop=gene["stop"],
            reference_genome=reference_genome,
            dataset=dataset,
        )
        data = _variants_to_frame(variants, gene_name, canonical_transcript)
        constraint = gene.get("gnomad_constraint") or {}
        return structs.GnomadDataset(
            gnomadVersion=gnomad_version,
            gene=gene_name,
            ensemblId=f"{actual_gene_id}.{gene['gene_version']}",
            canonicalTranscript=canonical_transcript,
            chromosomes=sorted(data.get_column("Chromosome").unique().to_list())
            if data.height
            else [],
            nVars=data.height,
            nVarsInCanonicalTranscript=round(
                data.get_column("CanonicalTranscript").sum()
            )
            if data.height
            else 0,
            source=f"{GNOMAD_API_URL} ({dataset}, {reference_genome})",
            endpoint_queried=GNOMAD_API_URL,
            gnomad_dataset=dataset,
            reference_genome=reference_genome,
            dataMatrix=data,
            exp_lof=constraint.get("exp_lof"),
            exp_mis=constraint.get("exp_mis"),
            exp_syn=constraint.get("exp_syn"),
            obs_lof=constraint.get("obs_lof"),
            obs_mis=constraint.get("obs_mis"),
            obs_syn=constraint.get("obs_syn"),
            oe_lof=constraint.get("oe_lof"),
            oe_lof_lower=constraint.get("oe_lof_lower"),
            oe_lof_upper=constraint.get("oe_lof_upper"),
            oe_lof_percentile=constraint.get("oe_lof_percentile"),
            oe_mis=constraint.get("oe_mis"),
            oe_mis_lower=constraint.get("oe_mis_lower"),
            oe_mis_upper=constraint.get("oe_mis_upper"),
            oe_syn=constraint.get("oe_syn"),
            oe_syn_lower=constraint.get("oe_syn_lower"),
            oe_syn_upper=constraint.get("oe_syn_upper"),
            lof_z=constraint.get("lof_z"),
            mis_z=constraint.get("mis_z"),
            syn_z=constraint.get("syn_z"),
            pli=constraint.get("pli"),
            constraint_flags=constraint.get("flags"),
        )


def _variants_to_frame(
    variants: list[dict[str, Any]], gene: str, canonical_transcript: str
) -> pl.DataFrame:
    rows = []
    for variant in variants:
        rows.append(
            {
                "gnomAD ID": variant["variant_id"],
                "Chromosome": int(variant["chrom"]),
                "Position": variant["pos"],
                "rsIDs": variant.get("rsid"),
                "Reference": variant["ref"],
                "Alternate": variant["alt"],
                "Transcript": variant.get("transcript_id"),
                "HGVS Consequence": variant.get("hgvsc"),
                "Protein Consequence": variant.get("hgvsp"),
                "VEP Annotation": variant.get("consequence"),
                "Flags": ",".join(variant.get("flags") or []),
                "Gene": gene,
                "Exonic": variant.get("consequence") != "intron_variant",
                "CanonicalTranscript": variant.get("transcript_id")
                == canonical_transcript,
            }
        )
    if not rows:
        return pl.DataFrame(
            schema={
                "gnomAD ID": pl.String,
                "Chromosome": pl.Int64,
                "Position": pl.Int64,
                "rsIDs": pl.String,
                "Reference": pl.String,
                "Alternate": pl.String,
                "Transcript": pl.String,
                "HGVS Consequence": pl.String,
                "Protein Consequence": pl.String,
                "VEP Annotation": pl.String,
                "Flags": pl.String,
                "Gene": pl.String,
                "Exonic": pl.Boolean,
                "CanonicalTranscript": pl.Boolean,
            }
        )
    return pl.DataFrame(rows, infer_schema_length=None)
