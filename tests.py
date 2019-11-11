import pytest
from fastqandfurious import fastqandfurious, _fastqandfurious
from Bio import SeqIO


def _test_readfastq_iter(filename, bufsize, entrypos):
    with open(filename, 'rb') as fh, open(filename, 'rt') as fh_bp:
        entry_iter = zip(fastqandfurious.readfastq_iter(fh, bufsize,
                                                        _entrypos=entrypos),
                         SeqIO.parse(fh_bp, "fastq"))
        for entry, entry_bp in entry_iter:
            header, sequence, quality = entry
            assert header == entry_bp.description.encode('ascii')
            assert sequence == str(entry_bp.seq).encode('ascii')


@pytest.mark.parametrize('filename',
                         ('data/test.fq', 'data/test_longqualityheader.fq'))
@pytest.mark.parametrize('func',
                         (fastqandfurious._entrypos,
                          _fastqandfurious.entrypos))
def test_readfastq_iter(filename, func):
    bufsize = 600
    _test_readfastq_iter(filename, bufsize, func)


@pytest.mark.parametrize('filename',
                         ('data/test.fq', 'data/test_longqualityheader.fq'))
def test_readfastq_iter_buftoosmall(filename):
    bufsize = 100
    with pytest.raises(RuntimeError):
        _test_readfastq_iter(filename, bufsize, fastqandfurious._entrypos)


@pytest.mark.parametrize('func',
                         (fastqandfurious._entrypos,
                          _fastqandfurious.entrypos))
def test_readfastq_iter_error(func):
    bufsize = 600
    with pytest.raises(ValueError):
        _test_readfastq_iter("data/test_tricky.fq", bufsize, func)


@pytest.mark.parametrize('filename',
                         ('data/test.fq', 'data/test_longqualityheader.fq'))
def test_readfastq_c_iter_buftoosmall(filename):
    bufsize = 100
    with pytest.raises(RuntimeError):
        _test_readfastq_iter(filename, bufsize, _fastqandfurious.entrypos)


def _test_readfastq_abspos(filename, bufsize, entrypos):
    with open(filename, 'rb') as fh, \
         open(filename, 'rt') as fh_bp, \
         open(filename, 'rb') as fh2:
        data = fh2.read()
        entry_iter = zip(fastqandfurious
                         .readfastq_iter(
                             fh, bufsize,
                             entryfunc=fastqandfurious.entryfunc_abspos,
                             _entrypos=entrypos
                         ),
                         SeqIO.parse(fh_bp, "fastq"))
        for i, (posarray, entry_bp) in enumerate(entry_iter):
            header = data[(posarray[0]+1):posarray[1]]
            sequence = data[posarray[2]:posarray[3]]
            assert header == entry_bp.description.encode('ascii')
            assert sequence == str(entry_bp.seq).encode('ascii')


@pytest.mark.parametrize('filename',
                         ('data/test.fq', 'data/test_longqualityheader.fq'))
def test_readfastq_abspos(filename):
    bufsize = 600
    _test_readfastq_abspos(filename, bufsize, fastqandfurious._entrypos)


@pytest.mark.parametrize('filename',
                         ('data/test.fq', 'data/test_longqualityheader.fq'))
def test_readfastq_abspos_buftoosmall(filename):
    bufsize = 100
    with pytest.raises(RuntimeError):
        _test_readfastq_abspos(filename, bufsize, fastqandfurious._entrypos)


def test_readfastq_abspos_error():
    bufsize = 600
    with pytest.raises(ValueError):
        _test_readfastq_abspos("data/test_tricky.fq",
                               bufsize, fastqandfurious._entrypos)
