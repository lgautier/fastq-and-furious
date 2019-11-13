import pytest
from array import array
from fastqandfurious import fastqandfurious, _fastqandfurious
from Bio import SeqIO

HEADER = 'foo#2'
SEQUENCE = 'AATTGCCG'
QUALITY = '3425@!#!'

ENTRIES = """
@{header}
{sequence}
+
{quality}
@bar{header}
""".format(header=HEADER,
           sequence=SEQUENCE,
           quality=QUALITY).encode('ascii')

ENTRIES_FINAL = """
@{header}
{sequence}
+
{quality}
""".format(header=HEADER,
           sequence=SEQUENCE,
           quality=QUALITY).encode('ascii')

ENTRIES_QUALHEAD = """
@{header}
{sequence}
+
{quality}
@bar{header}
""".format(header=HEADER,
           sequence=SEQUENCE,
           quality=QUALITY).encode('ascii')


@pytest.mark.parametrize('func',
                         (fastqandfurious._entrypos,
                          _fastqandfurious.entrypos))
@pytest.mark.parametrize('entries,status',
                         ((ENTRIES, fastqandfurious.COMPLETE),
                          (ENTRIES_FINAL, fastqandfurious.MISSING_QUAL_END),
                          (ENTRIES_QUALHEAD, fastqandfurious.COMPLETE)))
def test_parseentry(func, entries, status):
    posbuffer = array('q', [-1, ] * 6)
    offset = 0
    globaloffset = 0
    res_status = func(entries, offset, posbuffer)
    assert res_status == status
    (header,
     sequence,
     quality) = fastqandfurious.entryfunc(entries,
                                          posbuffer,
                                          globaloffset)
    assert header == HEADER.encode('ascii')
    assert sequence == SEQUENCE.encode('ascii')
    assert quality == QUALITY.encode('ascii')
    

@pytest.mark.parametrize('func',
                         (fastqandfurious._entrypos,
                          _fastqandfurious.entrypos))
@pytest.mark.parametrize('entries',
                         (ENTRIES, ENTRIES_FINAL, ENTRIES_QUALHEAD))
@pytest.mark.parametrize('bufsize,status',
                         ((3, fastqandfurious.MISSING_SEQHEADER_END),
                          (9, fastqandfurious.MISSING_SEQ_END),
                          (20, fastqandfurious.MISSING_QUAL_END)))
def test_parseentry_incomplete(func, entries, bufsize, status):
    posbuffer = array('q', [-1, ] * 6)
    offset = 0
    globaloffset = 0
    res_status = func(entries[:bufsize], offset, posbuffer)
    assert res_status == status

    
@pytest.mark.parametrize('filename,fixmultiline',
                         (('data/test.fq', False),
                          ('data/test_longqualityheader.fq', False),
                          ('data/test_multiline.fq', True)))
@pytest.mark.parametrize('func',
                         (fastqandfurious._entrypos,
                          _fastqandfurious.entrypos))
@pytest.mark.parametrize('bufsize', (600, 700))
def _test_readfastq_iter(filename, fixmultiline, func, bufsize):
    with open(filename, 'rb') as fh, open(filename, 'rt') as fh_bp:
        entry_iter = zip(fastqandfurious.readfastq_iter(fh, bufsize,
                                                        _entrypos=entrypos),
                         SeqIO.parse(fh_bp, "fastq"))
        for entry, entry_bp in entry_iter:
            header, sequence, quality = entry
            assert header == entry_bp.description.encode('ascii')
            
            if fixmultiline:
                assert (sequence.replace(b'\n', b'')
                        ==
                        str(entry_bp.seq).encode('ascii'))
            else:
                assert sequence == str(entry_bp.seq).encode('ascii')


def _test_readfastq_abspos(filename, bufsize, entrypos,
                           fixmultiline=False):
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
            if fixmultiline:
                assert (sequence.replace(b'\n', b'')
                        ==
                        str(entry_bp.seq).encode('ascii'))
            else:
                assert sequence == str(entry_bp.seq).encode('ascii')


@pytest.mark.parametrize('filename,fixmultiline',
                         (('data/test.fq', False),
                          ('data/test_longqualityheader.fq', False),
                          ('data/test_multiline.fq', True)))
@pytest.mark.parametrize('bufsize', (100, 200, 600, 700))
def test_readfastq_abspos(filename, fixmultiline, bufsize):
    _test_readfastq_abspos(filename, bufsize, fastqandfurious._entrypos,
                           fixmultiline=fixmultiline)
