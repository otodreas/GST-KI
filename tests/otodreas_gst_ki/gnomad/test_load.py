from pathlib import Path

import pytest

from otodreas_gst_ki.gnomad import load


class TestGnomad:
    def test_gnomad_input_not_file(self):
        with pytest.raises(FileNotFoundError):
            load.load_gnomad("tests/data/gnomad/non_existent_file.csv")

    def test_gnomad_no_csv(self):
        with pytest.raises(TypeError):
            load.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short_tab.csv"
            )

    def test_gnomad_multichromo(self):
        with pytest.warns(UserWarning):
            load.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short_multichromo.csv"
            )

    def test_gnomad_load_success(self):
        dataset = load.load_gnomad(
            "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
        )
        assert dataset.gene == "GSTT2B"
        assert dataset.dataMatrix.shape == (9, 75)
        assert dataset.ensemblId == "ENSG00000133433.11"
        assert dataset.canonicalTranscript == "ENST00000290765.9"
        assert dataset.chromosomes == [22]
        assert dataset.filePath == Path(
            "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
        )

    def test_gnomad_concat_too_few_datasets(self):
        dataset = load.load_gnomad(
            "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
        )
        with pytest.raises(ValueError):
            load.concat_gnomad([dataset])

    def test_gnomad_concat_success(self):
        datasets = [
            load.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
            ),
            load.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
            ),
        ]
        assert load.concat_gnomad(datasets).shape == (18, 75)

    def test_shared_variants(self):
        # Since gene names are parsed from the file name, a test file for ENSG00000133433 gets renamed to ENSG00000148834.
        # This makes the loader "think" that the two identical datasets belong to different genes, thus having 100% overlap
        datasets_shared = [
            load.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
            ),
            load.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000148834_2026_11_12_16_00_24_short_difgene.csv"
            ),
        ]

        # These two files are headers of two completely non-overlapping files
        datasets_not_shared = [
            load.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
            ),
            load.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000148834_2026_11_12_16_01_30_short.csv"
            ),
        ]

        with pytest.warns(UserWarning):
            shared_vars_frame = load.find_shared_variants(datasets_shared)
        assert shared_vars_frame.shape == (9, 2)
        assert load.find_shared_variants(datasets_not_shared).shape == (0, 2)

    def test_provenance_check(self):
        bad_table = {
            "ENSG00000134184": "GSTM1",
            "ENSG00000134202": "GSTM3",
            "ENSG00000134201": "GSTM5",
            "ENSG00000148834": "GSTO1",
            "ENSG00000133432": "GSTT2B",  # the ensembl id is incorrect here
            "ENSG00000243955": "GSTA1",
            "ENSG00000174156": "GSTA3",
            "ENSG00000170899": "GSTA4",
            "ENSG00000182793": "GSTA5",
            "ENSG00000197448": "GSTK1",
            "ENSG00000213366": "GSTM2",
            "ENSG00000168765": "GSTM4",
            "ENSG00000065621": "GSTO2",
            "ENSG00000084207": "GSTP1",
            "ENSG00000099984": "GSTT2",
            "ENSG00000100577": "GSTZ1",
        }
        with pytest.raises(ValueError):
            load.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133432_2026_11_12_16_00_24_short_provenance_mismatch.csv",
                ensembl_to_gene_lookup=bad_table,
            )

    def test_gnomad_version_extraction(self):
        assert (
            load.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
            ).gnomadVersion
            == "v4.1.2"
        )

    def test_multiple_intron_veps(self):
        with pytest.warns(UserWarning):  # case insensitive
            load.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short_multiple_vep_intron.csv"
            )
