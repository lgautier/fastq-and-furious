import pytest
from fastqandfurious import fastqandfurious, _fastqandfurious
from Bio import SeqIO


def _test_readfastq_iter(filename, bufsize, entrypos):
    with open(filename, 'rb') as fh, open(filename, 'rt') as fh_bp:
        for entry, entry_bp in zip(fastqandfurious.readfastq_iter(fh, bufsize, _entrypos=entrypos),
                                   SeqIO.parse(fh_bp, "fastq")):
            header, sequence, quality = entry
            assert header == entry_bp.description.encode('ascii')
            assert sequence == str(entry_bp.seq).encode('ascii')

def test_readfastq_iter():
    for filename in ("data/test.fq", "data/test_longqualityheader.fq"):
        bufsize = 600
        _test_readfastq_iter(filename, bufsize, fastqandfurious._entrypos)
        # buffer too low
        with pytest.raises(RuntimeError):
            bufsize = 100
            _test_readfastq_iter(filename, bufsize, fastqandfurious._entrypos)
    bufsize = 600
    with pytest.raises(ValueError):
        _test_readfastq_iter("data/test_tricky.fq", bufsize, fastqandfurious._entrypos)
        
def test_readfastq_c_iter():
    for filename in ("data/test.fq", ):
        bufsize = 600
        _test_readfastq_iter(filename, bufsize, _fastqandfurious.entrypos)
        # buffer too low
        with pytest.raises(RuntimeError):
            bufsize = 100
            _test_readfastq_iter(filename, bufsize, _fastqandfurious.entrypos)
    bufsize = 600
    with pytest.raises(ValueError):
        _test_readfastq_iter("data/test_tricky.fq", bufsize, fastqandfurious._entrypos)


def _test_readfastq_abspos(filename, bufsize, entrypos):
    with open(filename, 'rb') as fh, \
         open(filename, 'rt') as fh_bp, \
         open(filename, 'rb') as fh2:
        data = fh2.read()
        for i, (posarray, entry_bp) in enumerate(zip(fastqandfurious.readfastq_iter(fh, bufsize,
                                                                                    entryfunc=fastqandfurious.entryfunc_abspos,
                                                                                    _entrypos=entrypos),
                                                     SeqIO.parse(fh_bp, "fastq"))):
            header = data[(posarray[0]+1):posarray[1]]
            sequence = data[posarray[2]:posarray[3]]
            assert header == entry_bp.description.encode('ascii')
            assert sequence == str(entry_bp.seq).encode('ascii')

def test_readfastq_abspos():
    for filename in ("data/test.fq", "data/test_longqualityheader.fq"):
        bufsize = 600
        _test_readfastq_abspos(filename, bufsize, fastqandfurious._entrypos)
        # buffer too low
        with pytest.raises(RuntimeError):
            bufsize = 100
            _test_readfastq_abspos(filename, bufsize, fastqandfurious._entrypos)
    bufsize = 600
    with pytest.raises(ValueError):
        _test_readfastq_iter("data/test_tricky.fq", bufsize, fastqandfurious._entrypos)
