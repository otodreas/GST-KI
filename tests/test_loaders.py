from pathlib import Path

import pytest

from gst_ki import loaders


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
        assert next(iter(dataset.keys())) == "GSTT2B"
        assert dataset["GSTT2B"].dataMatrix.shape == (9, 73)
        assert dataset["GSTT2B"].ensemblId == "ENSG00000133433.11"
        assert dataset["GSTT2B"].canonicalTranscript == "ENST00000290765.9"
        assert dataset["GSTT2B"].chromosomes == [22]
        assert (
            dataset["GSTT2B"].filePath
            == Path(
                "tests/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
            ).resolve()
        )
