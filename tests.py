import pytest
from array import array
import enum
import textwrap
from fastqandfurious import fastqandfurious, _fastqandfurious
from Bio import SeqIO

HEADER = 'foo#2'
SEQUENCE = 'AATTGCCG'
QUALITY = '3425@!#!'
MLINE_SEQUENCE = 'AATTGCCG\nGCCGTA'
MLINE_QUALITY = '3425@!#!\n255212'

class ENTRIES_FQ_TPL(enum.Enum):

    FINAL = """
@{header}
{sequence}
+
{quality}
"""

    QUALHEAD = """
@{header}
{sequence}
+
{quality}
@bar{header}
"""

    NOQUAL = """
@{header}
{sequence}
+
"""

class ENTRIES_FA_TPL(enum.Enum):

    NOTFINAL = textwrap.dedent("""
    >{header}
    {sequence}
    >{header}_2
    {sequence}
    """)

    FINAL = textwrap.dedent("""
    >{header}
    {sequence}
    """)

    NOSEQ = textwrap.dedent("""
    >{header}
    """)


@pytest.mark.parametrize('func',
                         (fastqandfurious.entrypos,
                          _fastqandfurious.entrypos))
@pytest.mark.parametrize('entries_tpl,status',
                         ((ENTRIES_FQ_TPL.FINAL.value, fastqandfurious.MISSING_QUAL_END),
                          (ENTRIES_FQ_TPL.QUALHEAD.value, fastqandfurious.COMPLETE)))
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
                         (fastqandfurious.entrypos_fasta,))
@pytest.mark.parametrize(
    'entries_tpl,status,test_sequence',
    ((ENTRIES_FA_TPL.NOTFINAL.value, fastqandfurious.COMPLETE, SEQUENCE),
     (ENTRIES_FA_TPL.NOTFINAL.value, fastqandfurious.COMPLETE, MLINE_SEQUENCE),
     (ENTRIES_FA_TPL.FINAL.value, fastqandfurious.MISSING_SEQ_END, SEQUENCE),
     (ENTRIES_FA_TPL.FINAL.value, fastqandfurious.MISSING_SEQ_END, MLINE_SEQUENCE),
     (ENTRIES_FA_TPL.NOSEQ.value, fastqandfurious.MISSING_SEQ_BEG, ''),
    )
)
def test_parseentry(func, entries_tpl, status, test_sequence):
    entries = entries_tpl.format(header=HEADER,
                                 sequence=test_sequence).encode('ascii')
    posbuffer = array('q', [-1, ] * 6)
    offset = 0
    globaloffset = 0
    res_status = func(entries, offset, posbuffer)
    assert res_status == status
    (header,
     sequence) = fastqandfurious.entryfunc_fasta(entries,
                                                 posbuffer,
                                                 globaloffset)
    assert header == HEADER.encode('ascii')
    assert sequence == test_sequence.encode('ascii')


@pytest.mark.parametrize(
    'func',
    (fastqandfurious.entrypos, 
     pytest.param(
         _fastqandfurious.entrypos,
         marks=pytest.mark.xfail(reason='C parser not consistent with Python parser')))
)
@pytest.mark.parametrize('status,test_sequence',
                         ((fastqandfurious.MISSING_QUAL_BEGIN, SEQUENCE),
                          (fastqandfurious.MISSING_QUAL_BEGIN, MLINE_SEQUENCE)))
def test_parseentry_missqual(func, status, test_sequence):
    entries_tpl = ENTRIES_FQ_TPL.NOQUAL.value
    test_quality = ''
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
                         (ENTRIES_FQ_TPL.FINAL.value, ENTRIES_FQ_TPL.QUALHEAD.value))
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
