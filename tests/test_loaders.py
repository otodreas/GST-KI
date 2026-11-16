from pathlib import Path

import pytest

from otodreas_gst_ki import loaders


class TestGnomad:
    def test_gnomad_input_not_file(self):
        with pytest.raises(FileNotFoundError):
            loaders.load_gnomad("tests/non_existent_file.csv")

    def test_gnomad_no_csv(self):
        with pytest.raises(TypeError):
            loaders.load_gnomad(
                "tests/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short_tab.csv"
            )

    def test_gnomad_multichromo(self):
        with pytest.warns(UserWarning):
            loaders.load_gnomad(
                "tests/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short_multichromo.csv"
            )

    def test_gnomad_load_success(self):
        dataset = loaders.load_gnomad(
            "tests/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
        )
        assert dataset.gene == "GSTT2B"
        assert dataset.dataMatrix.shape == (9, 74)
        assert dataset.ensemblId == "ENSG00000133433.11"
        assert dataset.canonicalTranscript == "ENST00000290765.9"
        assert dataset.chromosomes == [22]
        assert (
            dataset.filePath
            == Path(
                "tests/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
            ).resolve()
        )

    def test_gnomad_concat_too_few_datasets(self):
        dataset = loaders.load_gnomad(
            "tests/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
        )
        with pytest.raises(ValueError):
            loaders.concat_gnomad([dataset])

    def test_gnomad_concat_success(self):
        datasets = [
            loaders.load_gnomad(
                "tests/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
            ),
            loaders.load_gnomad(
                "tests/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
            ),
        ]
        assert loaders.concat_gnomad(datasets).shape == (18, 74)
