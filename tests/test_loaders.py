from pathlib import Path

import pytest

from otodreas_gst_ki import loaders


class TestGnomad:
    def test_gnomad_input_not_file(self):
        with pytest.raises(FileNotFoundError):
            loaders.load_gnomad("tests/data/gnomad/non_existent_file.csv")

    def test_gnomad_no_csv(self):
        with pytest.raises(TypeError):
            loaders.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short_tab.csv"
            )

    def test_gnomad_multichromo(self):
        with pytest.warns(UserWarning):
            loaders.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short_multichromo.csv"
            )

    def test_gnomad_load_success(self):
        dataset = loaders.load_gnomad(
            "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
        )
        assert dataset.gene == "GSTT2B"
        assert dataset.dataMatrix.shape == (9, 74)
        assert dataset.ensemblId == "ENSG00000133433.11"
        assert dataset.canonicalTranscript == "ENST00000290765.9"
        assert dataset.chromosomes == [22]
        assert (
            dataset.filePath
            == Path(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
            ).resolve()
        )

    def test_gnomad_concat_too_few_datasets(self):
        dataset = loaders.load_gnomad(
            "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
        )
        with pytest.raises(ValueError):
            loaders.concat_gnomad([dataset])

    def test_gnomad_concat_success(self):
        datasets = [
            loaders.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
            ),
            loaders.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
            ),
        ]
        assert loaders.concat_gnomad(datasets).shape == (18, 74)

    def test_shared_variants(self):
        # Since gene names are parsed from the file name, a test file for ENSG00000133433 gets renamed to ENSG00000148834.
        # This makes the loader "think" that the two identical datasets belong to different genes, thus having 100% overlap
        datasets_shared = [
            loaders.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
            ),
            loaders.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000148834_2026_11_12_16_00_24_short_difgene.csv"
            ),
        ]

        # These two files are headers of two completely non-overlapping files
        datasets_not_shared = [
            loaders.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
            ),
            loaders.load_gnomad(
                "tests/data/gnomad/gnomAD_v4.1.2_ENSG00000148834_2026_11_12_16_01_30_short.csv"
            ),
        ]

        with pytest.warns(UserWarning):
            shared_vars_frame = loaders.find_shared_variants(datasets_shared)
        assert shared_vars_frame.shape == (9, 2)
        assert loaders.find_shared_variants(datasets_not_shared).shape == (0, 2)
