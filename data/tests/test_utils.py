import pytest
from utils import dataclasses


@pytest.fixture
def sequence():
    return """TAAATTATATCAAATATGCTGTTTTGAATTTTATTTTTCCCTTTAAACAT
AAAGACACACATTCAGTTCATTGTGCTAGATAAATTACCAGTGCGATCAC
AAATTAAGAAATGCAATTCAAAGAATTTTGCATACAAGGAGTCCTGAAAG
TGTTAATAACTTTTGATGCAAAGATAATTTTATGAAAGTAATAGAAGACT
AAAAAAGGTACAAAATAACTATTATGTAAGTATTTTCCTTTTTCTGAATC
ACCCATGATTACTTTTTCCACCAAGCAAAAACTAACTGCATACTTCAGAC
CTGTCTCAAATCTCCCCAGCCTCTTTTCCTAAACCTCCCCAGCCTCTCAG
GACAGACAGGCTGCTCCTGTATTTTGTGCATTCTGCTATTTTTAGCAAGA
GGCCTATTTTGTCAGTGTTGTCTGAATAGTATTTGCCAACTCTCAGACTT
TCAGTCACTTATTTGTTTATTTATTTATTTATTTGTCTCCTTTTCTTGTA
TTTCTCTTTTCCTTTTCTTTCCTTTCTTTTTCCCTTTCCTCCTCCCTTCC
TTTGCTTACTTATTTTTTTTCCCTTTAATTCCCATTCACTATTTCCATGA"""


def test_seqrecord(sequence):
    sr = dataclasses.SeqRecord(id="seq1", name="seq1", seq=sequence)
    print(sr)
