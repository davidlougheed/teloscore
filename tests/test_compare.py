import io
from pathlib import Path
from teloscore import compare, scoring

TEST_DATA_DIR = Path(__file__).parent.parent / "test_data"

HG002 = TEST_DATA_DIR / "tlens_by_allele_HG002.tsv"
HG003 = TEST_DATA_DIR / "tlens_by_allele_HG003.tsv"
HG004 = TEST_DATA_DIR / "tlens_by_allele_HG004.tsv"

OUT_0_1_0_HG002_HG003 = TEST_DATA_DIR / "out_0_1_0_HG002_HG003.tsv"
OUT_0_1_0_HG002_HG004 = TEST_DATA_DIR / "out_0_1_0_HG002_HG004.tsv"

with open(OUT_0_1_0_HG002_HG003, "r") as ofh:
    OUT_0_1_0_HG002_HG003_DATA = ofh.read()

with open(OUT_0_1_0_HG002_HG004, "r") as ofh:
    OUT_0_1_0_HG002_HG004_DATA = ofh.read()


def test_loading():
    hg002_alleles, hg003_alleles = compare._read_sample_files(scoring.Scoring1(), HG002, HG003)

    assert len(hg002_alleles) == 65
    assert len(hg003_alleles) == 43


def test_compare_s1_hg002_hg003():
    with io.StringIO() as fh:
        compare.compare_samples(scoring.Scoring1(), HG002, HG003, fh)
        fh.seek(0)
        matrix_out = fh.read()
    assert matrix_out == OUT_0_1_0_HG002_HG003_DATA


def test_compare_s1_hg002_hg004():
    with io.StringIO() as fh:
        compare.compare_samples(scoring.Scoring1(), HG002, HG004, fh)
        fh.seek(0)
        matrix_out = fh.read()
    assert matrix_out == OUT_0_1_0_HG002_HG004_DATA


def test_compare_s2():
    with io.StringIO() as fh:
        compare.compare_samples(scoring.Scoring2(), HG002, HG004, fh)
        fh.seek(0)
        # TODO
