import pytest
from fastqandfurious import fastqandfurious, _fastqandfurious
from Bio import SeqIO

filename = "data/test.fq"

def _test_readfastq_iter(bufsize, entrypos):
    with open(filename, 'rb') as fh, open(filename, 'rt') as fh_bp:
        for entry, entry_bp in zip(fastqandfurious.readfastq_iter(fh, bufsize, _entrypos=entrypos),
                                   SeqIO.parse(fh_bp, "fastq")):
            assert entry.header == b'@'+entry_bp.description.encode('ascii')
            assert entry.sequence == str(entry_bp.seq).encode('ascii')

def test_readfastq_iter():
    bufsize = 600
    _test_readfastq_iter(bufsize, fastqandfurious._entrypos)
    # buffer too low
    with pytest.raises(RuntimeError):
        bufsize = 100
        _test_readfastq_iter(bufsize, fastqandfurious._entrypos)
        
def test_readfastq_c_iter():
    bufsize = 600
    _test_readfastq_iter(bufsize, _fastqandfurious.entrypos)
    # buffer too low
    with pytest.raises(RuntimeError):
        bufsize = 100
        _test_readfastq_iter(bufsize, _fastqandfurious.entrypos)
