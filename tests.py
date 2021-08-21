import pytest
from array import array
from fastqandfurious import fastqandfurious, _fastqandfurious
from Bio import SeqIO

HEADER = 'foo#2'
SEQUENCE = 'AATTGCCG'
QUALITY = '3425@!#!'
MLINE_SEQUENCE = 'AATTGCCG\nGCCGTA'
MLINE_QUALITY = '3425@!#!\n255212'

ENTRIES_TPL = """
@{header}
{sequence}
+
{quality}
@bar{header}
"""

ENTRIES_FINAL_TPL = """
@{header}
{sequence}
+
{quality}
"""

ENTRIES_QUALHEAD_TPL = """
@{header}
{sequence}
+
{quality}
@bar{header}
"""


@pytest.mark.parametrize('func',
                         (fastqandfurious.entrypos,
                          _fastqandfurious.entrypos))
@pytest.mark.parametrize('entries_tpl,status',
                         ((ENTRIES_TPL, fastqandfurious.COMPLETE),
                          (ENTRIES_FINAL_TPL, fastqandfurious.MISSING_QUAL_END),
                          (ENTRIES_QUALHEAD_TPL, fastqandfurious.COMPLETE)))
@pytest.mark.parametrize('test_sequence,test_quality',
                         ((SEQUENCE, QUALITY), (MLINE_SEQUENCE, MLINE_QUALITY)))
def test_parseentry(func, entries_tpl, status, test_sequence, test_quality):
    entries = entries_tpl.format(header=HEADER,
                                 sequence=test_sequence,
                                 quality=test_quality).encode('ascii')
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
    assert sequence == test_sequence.encode('ascii')
    assert quality == test_quality.encode('ascii')
    

@pytest.mark.parametrize('func',
                         (fastqandfurious.entrypos,
                          _fastqandfurious.entrypos))
@pytest.mark.parametrize('entries_tpl',
                         (ENTRIES_TPL, ENTRIES_FINAL_TPL, ENTRIES_QUALHEAD_TPL))
@pytest.mark.parametrize(
    'func_bufsize,status',
    ((lambda l_header, l_sequence, l_quality: l_header - 2,
      fastqandfurious.MISSING_SEQHEADER_END),
     (lambda l_header, l_sequence, l_quality: l_header + l_sequence - 1,
      fastqandfurious.MISSING_SEQ_END),
     (lambda l_header, l_sequence, l_quality: l_header + l_sequence + l_quality,
      fastqandfurious.MISSING_QUAL_END))
)
@pytest.mark.parametrize('test_sequence,test_quality',
                         ((SEQUENCE, QUALITY), (MLINE_SEQUENCE, MLINE_QUALITY)))
def test_parseentry_incomplete(func, entries_tpl, func_bufsize, status, test_sequence, test_quality):
    entries = entries_tpl.format(header=HEADER,
                                 sequence=test_sequence,
                                 quality=test_quality).encode('ascii')
    bufsize = func_bufsize(len(HEADER), len(test_sequence), len(test_quality))
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
                         (fastqandfurious.entrypos,
                          _fastqandfurious.entrypos))
@pytest.mark.parametrize('bufsize', (600, 700))
def _test_readfastq_iter(filename, fixmultiline, func, bufsize):
    with open(filename, 'rb') as fh, open(filename, 'rt') as fh_bp:
        entry_iter = zip(fastqandfurious.readfastq_iter(fh, bufsize,
                                                        entrypos=entrypos),
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
                             entrypos=entrypos
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
    _test_readfastq_abspos(filename, bufsize, fastqandfurious.entrypos,
                           fixmultiline=fixmultiline)
