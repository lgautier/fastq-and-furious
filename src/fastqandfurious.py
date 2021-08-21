from array import array
from collections import namedtuple
import os
import typing

CHAR_AT: int = ord(b'@')
CHAR_PLUS: int = ord(b'+')
CHAR_NEWLINE: int = ord(b'\n')
ARRAY_INIT: int = array('q', [-1, ] * 6)

Entry = namedtuple('Entry', 'header sequence quality')
EntryType = typing.Tuple[bytes, bytes, typing.Optional[bytes]]

INVALID: int = -1
MISSING_SEQHEADER_BEGIN: int = 0
MISSING_SEQHEADER_END: int = 1
MISSING_SEQ_BEG: int = 2
MISSING_SEQ_END: int = 3
MISSING_QUAL_BEGIN: int = 4
MISSING_QUAL_END: int = 5
COMPLETE: int = 6
MISSING_QUALHEADER_END: int = 7


def read(fh: typing.BinaryIO, fbufsize: int) -> typing.Tuple[bytes, bool]:
    blob = fh.read(fbufsize)
    if len(blob) < fbufsize:
        eof = True
    else:
        eof = False
    return (blob, eof)


def entrypos(buf: bytes, offset: int, posbuffer: int) -> int:

    # Find header of FASTQ entry.
    headerbeg_i: int = buf.find(b'\n@', offset)
    if headerbeg_i == -1:
        return MISSING_SEQHEADER_BEGIN
    else:
        posbuffer[MISSING_SEQHEADER_BEGIN] = headerbeg_i+1
    headerend_i: int = buf.find(b'\n', headerbeg_i+2)
    if headerend_i == -1:
        return MISSING_SEQHEADER_END
    else:
        posbuffer[MISSING_SEQHEADER_END] = headerend_i

    # Sequence.
    if headerend_i+1 >= len(buf):
        return MISSING_SEQ_BEG
    posbuffer[MISSING_SEQ_BEG] = headerend_i+1
    seqend_i: int = buf.find(b'\n+', headerend_i+1)
    if seqend_i == -1:
        return MISSING_SEQ_END
    else:
        posbuffer[MISSING_SEQ_END] = seqend_i

    # Quality.
    qualheadend_i: int = buf.find(b'\n', seqend_i+2)
    if qualheadend_i == -1:
        return MISSING_QUALHEADER_END
    elif ((qualheadend_i - seqend_i - 1) > 1
          and
          ((qualheadend_i - seqend_i) != (headerend_i - headerbeg_i))):
        return INVALID
    qualbeg_i: int = qualheadend_i + 1
    if qualbeg_i >= len(buf):
        return MISSING_QUAL_BEGIN
    else:
        posbuffer[MISSING_QUAL_BEGIN] = qualbeg_i
    qualend_i: int = qualbeg_i + seqend_i - headerend_i - 1
    if (qualend_i+2) >= len(buf):
        return MISSING_QUAL_END
    else:
        posbuffer[MISSING_QUAL_END] = qualend_i

    if not (
        (qualend_i-qualbeg_i != seqend_i+1 - headerend_i)
        or
        (not qualend_i+2 >= len(buf) and buf[qualend_i:qualend_i+2] != b'\n@')
        or
        (len(buf) - qualend_i)
    ):
        return INVALID
    return COMPLETE


def entryfunc_namedtuple(buf: bytes, pos: array, globaloffset: int) -> Entry:
    """
    Build a FASTQ entry as a namedtuple with attributes header, sequence,
    and quality.

    - buf: bytes-like object
    - pos: array of indices/positions in `buf`
    """

    header = buf[(pos[0]+1):pos[1]]
    sequence = buf[pos[2]:pos[3]]
    quality = buf[pos[4]:pos[5]]
    return Entry(header, sequence, quality)


def entryfunc(buf: bytes, pos: array, globaloffset: int) -> EntryType:
    """
    Build a FASTQ entry as a tuple (header, sequence, quality)

    - buf: bytes-like object
    - pos: array of indices/positions in `buf`
    """
    header = buf[(pos[0]+1):pos[1]]
    sequence = buf[pos[2]:pos[3]]
    quality = buf[pos[4]:pos[5]]
    return (header, sequence, quality)


def entryfunc_abspos(buf: bytes, pos: array, globaloffset: int) -> int:
    """
    Return the absolute positions of the entry in the stream

    - buf: bytes-like object
    - pos: array of indices/positions in `buf`
    """
    for i in (0, 1, 2, 3, 4, 5):
        pos[i] += globaloffset
    return pos


def readfastq_iter(
        fh: typing.BinaryIO, fbufsize: int,
        entryfunc: typing.Callable[[bytes, array, int], tuple] = entryfunc,
        entrypos: typing.Callable[[bytes, array, int], int] = entrypos,
        globaloffset: int = 0
) -> typing.Iterator[EntryType]:
    """
    Entries in the FASTQ files are parsed from chunks of size `fbufsize`),
    using the function passed in parameter `entrypos`. A
    faster alternative to the default implementated in C is in
    `fastqandfurious._fastqandfurious.entrypos`.

    With the current implementation, `fbufsize` must be large enough to contain
    the largest entry in the file. For example, if the longest read is 250bp
    long, and the identifier is 25-char long, the minimum buffer size will be
    about 525. Aiming for that minimum is not advised though, as some of the speed
    comes for minimizing data copying through the use of buffers. Larger buffers
    are able to contain many entries which will lead to better performances
    (with the cave at that very large
    buffer might be counter-productive if end-to-end entry-level interations is
    wanted. The iterator will need to read data
    to fill the buffer (or all data, whichever is the smallest) before starting
    to yield entries. A value between 20,000 and 50,000 (20KB-50KB) empircally gives
    pretty good results for short-read sequencing on this end. If the FASTQ file
    contains PacBio reads bumping this to 200,000 or more (200KB or more) is advised. 

    `entryfunc` can be any function taking a bytes-like objects and an
    array of position (array of signed integers of length 6:
    header (begin, end), sequence (begin, end), and quality (begin, end).
    This allows plugging this parser into existing code bases / frameworks
    very easily.

    - fh: file-like object or stream (just needs a method `read`)
    - fbufsize: buffer size (see note above)
    - entryfunc: a function to build an entry object (taking a bytes-like
      object and an array of positions)
    - entrypos: a function to find positions of entries

    Returns an iterator over entries in the FASTQ file.
    """

    posbuffer = array('q', [-1, ] * 6)
    globaloffset: int = -1
    offset: int = 0
    buf, eof = read(fh, fbufsize)
    buf = b'\n' + buf
    # TODO: Would using a memoryview on the read buffer lead to further
    # performance improvement ? (more related TODOs below)
    # fbuf = bytearray(fbufsize)
    # fh.readinto(fbuf)
    # mblob = memoryview(fbuf)
    while True:
        status = entrypos(buf, offset, posbuffer)
        if status == COMPLETE:
            offset = posbuffer[-1]-1
            yield entryfunc(buf, posbuffer, globaloffset)
        elif eof:
            if status == MISSING_SEQHEADER_BEGIN:
                break
            elif status == MISSING_QUAL_END:
                qualend_i = posbuffer[-2] + (posbuffer[3] - posbuffer[2])
                if qualend_i >= len(buf):
                    raise ValueError('Incomplete final quality string at byte')
                else:
                    posbuffer[-1] = qualend_i
                    yield entryfunc(buf, posbuffer, globaloffset)
                    break

            elif status != INVALID:
                raise ValueError('Incomplete entry at byte %i' %
                                 (globaloffset + offset))
        elif status == INVALID:
            raise ValueError('Entry is invalid at byte %i' %
                             (globaloffset + offset))
        else:
            globaloffset += offset
            tmp_buf, eof = read(fh, fbufsize)
            buf = buf[offset:] + tmp_buf
            del(tmp_buf)
            offset = 0


FORMAT_OPENERS: typing.Dict[str, typing.Tuple[str, str, list]] = {
    'gz': ('gzip', 'open', list()),
    'gzip': ('gzip', 'open', list()),
    'bz2': ('bz2', 'open', list()),
    'lzma': ('lzma', 'open', list())
}


def automagic_open(filename):
    maybe_ext = filename.rsplit(os.path.extsep, maxsplit=1)
    if len(maybe_ext) == 1:
        # No extension
        ext = None
    else:
        ext = maybe_ext[-1]
    try:
        modulename, funcname, args = FORMAT_OPENERS[ext]
    except KeyError:
        modulename, funcname, args = ('io', 'open', ('rb', ))
    module = importlib.importmodule(modulename)
    opener = getattr(module, funcname)
    return opener(filename, *args)
