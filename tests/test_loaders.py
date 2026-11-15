import pytest

from gst_ki import loaders


def test_gnomad_input_not_file():
    with pytest.raises(FileNotFoundError):
        loaders.load_gnomad("tests/non_existent_file.csv")


def test_gnomad_load_success():
    dataset = loaders.load_gnomad(
        "tests/gnomAD_v4.1.2_ENSG00000133433_2026_11_12_16_00_24_short.csv"
    )
    assert next(iter(dataset.keys())) == "GSTT2B"
    assert dataset["GSTT2B"].dataMatrix.shape == (9, 73)
