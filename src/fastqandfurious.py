from typing import Callable
from array import array
CHAR_PLUS = ord(b'+')
CHAR_NEWLINE = ord(b'\n')
ARRAY_INIT = array('q', [-1,] * 6)

from collections import namedtuple
Entry = namedtuple('Entry', 'header sequence quality')

def nextentrypos(blob, backlog) -> int:
    """
    Return the position (index) for the begining of next FASTQ entry in a "blob" of text,
    or -1 if it cannot find any.

    :param blob: a "blob" of text as a bytes-like object
    :param backlog: bytes-like object with an unfinished entry from a "blob"
                    in the input right before this one. It can be an empty
                    sequence (e.g., b'')
    
    """

    
    # look for next "@" starting a line
    headerbeg_i = blob.find(b'\n@', 0)

    lbacklog = len(backlog)

    if lbacklog == 0:
        # this is all we have to make a decision then
        return headerbeg_i
    
    # it may still be begining of a "quality" line
    if (headerbeg_i == 0):
        if lbacklog > 1:
            if (backlog[-1] == CHAR_PLUS) and (backlog[-2] == CHAR_NEWLINE):
                # if so, look for the next line starting with "@"
                headerbeg_i = blob.find(b'\n@', headerbeg_i+1)
    elif (headerbeg_i == 1):
        if (blob[headerbeg_i-1] == CHAR_PLUS) and (backlog[-1] == CHAR_NEWLINE):
            # if so, look for the next line starting with "@"
            headerbeg_i = blob.find(b'\n@', headerbeg_i+1)
    else:
        # headerbeg_i >= 2
        if (blob[headerbeg_i-1] == CHAR_PLUS) and (blob[headerbeg_i-2] == CHAR_NEWLINE):
            # if so, look for the next line starting with "@"
            headerbeg_i = blob.find(b'\n@', headerbeg_i+1)
    return headerbeg_i

def _entrypos(blob, offset, posbuffer):
    posbuffer[:] = ARRAY_INIT
    lblob = len(blob)
    # header
    headerbeg_i = blob.find(b'@', offset)
    posbuffer[0] = headerbeg_i
    if (headerbeg_i > 0) and (blob[headerbeg_i-1] == CHAR_PLUS):
        raise Exception("Lost sync' in input data.")
    if headerbeg_i == -1:
        return 0
    headerend_i = blob.find(b'\n', headerbeg_i+1)
    posbuffer[1] = headerend_i
    if headerend_i == -1:
        return 1
    
    # sequence
    seqbeg_i = headerend_i+1
    posbuffer[2] = seqbeg_i
    if seqbeg_i == lblob:
        return 2
    seqend_i = blob.find(b'\n', seqbeg_i)
    if seqend_i == -1 or (seqend_i + 1) == lblob:
        return 3
    else:
        posbuffer[3] = seqend_i
        if blob[seqend_i + 1] != ord(b'+'):
            # multi-line FASTQ :/
            raise ValueError("Multi-line FASTQ. Bye. (expected '+' but got '%s')." % \
                             blob[max(0, seqend_i-3):min(len(blob), seqend_i+4)])
        if (seqend_i + 2) >= lblob:
            return 4
        if blob[seqend_i + 2] == CHAR_NEWLINE:
            # if most-common situation where separator for quality sequence only a "+"
            qualbeg_i = seqend_i+3
        else:
            # name in header can optionally be repeated
            lheader = posbuffer[1] - posbuffer[0] + 1
            if (blob[seqend_i + lheader] == CHAR_NEWLINE) and \
               (blob[(posbuffer[0]+1):posbuffer[1]] == blob[(seqend_i + 2):(seqend_i+lheader)]):
                qualbeg_i = seqend_i + lheader + 1
            else:
                raise ValueError("Invalid quality header (sequence header is '%s' and quality header is '%s'." % \
                                 (blob[(posbuffer[0]):posbuffer[1]],
                                  blob[(seqend_i+1):(seqend_i+lheader)],))
    # quality
    if qualbeg_i >= lblob or seqend_i == -1:
        return 4
    else:
        posbuffer[4] = qualbeg_i
        qualend_i = qualbeg_i + (seqend_i - seqbeg_i)
        if qualend_i >= lblob:
            qualend_i = -1
        else:
            assert blob[qualend_i] == CHAR_NEWLINE
        posbuffer[5] = qualend_i
    if qualend_i == -1:
        return 5
    else:
        return 6

def entryfunc_namedtuple(buf, pos: array) -> Entry:
    """ 
    Build a FASTQ entry as a namedtuple with attributes header, sequence, and quality.

    - buf: bytes-like object
    - pos: array of indices/positions in `buf`
    """
    header = buf[(pos[0]+1):pos[1]]
    sequence = buf[pos[2]:pos[3]]
    quality = buf[pos[4]:pos[5]]
    return Entry(header, sequence, quality)

def entryfunc(buf, pos: array) -> tuple:
    """ 
    Build a FASTQ entry as a tuple (header, sequence, quality)

    - buf: bytes-like object
    - pos: array of indices/positions in `buf`
    """
    header = buf[(pos[0]+1):pos[1]]
    sequence = buf[pos[2]:pos[3]]
    quality = buf[pos[4]:pos[5]]
    return (header, sequence, quality)
    
def readfastq_iter(fh, fbufsize: int, entryfunc = entryfunc, _entrypos = _entrypos):
    """
    The entries in the FASTQ files are parsed from chunks of size `fbufsize`),
    using the function `_entrypos` (that be changed as a parameter - an
    faster implementation in C iblob[(posbuffer[0]+1):posbuffer[1]]s in `fastqandfurious._fastqandfurious.entrypos`).

    With the current implementation, `fbufsize` must be large enough to contain
    the largest entry in the file. For example, is the longest read is 250bp long,
    and the identifier is 25-char long, the minimum buffer size will be about
    525. Aiming for that minimum is not advised, as some of the speed comes
    for working on buffers, and larger buffers able to contain many entries will
    lead to better performances (with the cave at that very large buffer might be
    counter-productive as the iterator will need to read data to fill the buffer
    (or all data, whichever is the smallest) before starting to yield entries. A
    value between 20,000 and 50,000 (20KB-50KB) gives pretty good results on this end.

    `entryfunc` can be any function taking a bytes-like objects and an 
    array of position (array of signed integers of length 6: header (begin, end),
    sequence (begin, end), and quality (begin, end). This allows plugging this
    parser into existing code bases / frameworks very easily.
 
    - fh: file-like object or stream (just needs a method `read`)
    - fbufsize: buffer size (see note above)
    - entryfunc: a function to build an entry object (taking a bytes-like object and an array of positions)
    - _entrypos: a function to find positions of entries 

    Returns an iterator over entries in the FASTQ file.

    """
    posbuffer = array('q', [-1, ] * 6)
    fbuf = bytearray(fbufsize)
    blob = fh.read(fbufsize)
    #mblob = memoryview(blob)
    carryover = False
    offset = 0
    backlog = b''
    # FIXME
    header = b''
    lblob = len(blob)
    while blob != b'':
        posbuffer[0] = 0
        posbuffer[2] = 0
        posbuffer[3] = 0
        qualend_i = 0
        while qualend_i  < (lblob-1) and posbuffer[3] >= 0 and posbuffer[2] >= 0:
            npos = _entrypos(blob, offset, posbuffer)
            qualend_i = posbuffer[5]
            if npos < 6:
                #if headerbeg_i == -1:
                #    import pdb; pdb.set_trace()
                #    raise RuntimeError("Missing begining of entry ?")
                carryover = True
                backlog = blob[offset:] # bytes(mblob[offset:])
                break
            else:
                #(headerbeg_i, headerend_i, seqbeg_i, seqend_i, qualbeg_i, qualend_i) = posbuffer
                yield entryfunc(blob, posbuffer)
                offset = qualend_i+1
        blob = fh.read(fbufsize)
        #mblob = memoryview(blob)
        lblob = len(blob)
        offset = 0
        if carryover:
            if lblob == 0:
                if backlog == b'\n':
                    continue
                else:
                    nextentry_i = -1
            else:
                nextentry_i = nextentrypos(blob, backlog)
            if nextentry_i == -1:
                raise RuntimeError("Incomplete last entry, or buffer too small.")
            # FIXME:
            tmp = backlog
            backlog = backlog + blob[:(nextentry_i+1)]
            offset = nextentry_i+1
            npos = _entrypos(backlog, 0, posbuffer)
            if npos < 6:
                # FIXME: handle this case
                raise RuntimeError("The buffer is too small !")
            yield entryfunc(backlog, posbuffer)
            backlog = b''
            carryover = False


