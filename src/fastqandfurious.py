from array import array
CHAR_PLUS = ord(b'+')
CHAR_NEWLINE = ord(b'\n')
ARRAY_INIT = array('q', [-1,] * 6)

from collections import namedtuple
Entry = namedtuple('Entry', 'header sequence')

def _nextentrypos(blob, backlog):
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
        return 2
    else:
        posbuffer[3] = seqend_i
        if blob[seqend_i + 1] != ord(b'+'):
            # multi-line FASTQ :/
            raise ValueError("Multi-line FASTQ. Bye.")
        if (seqend_i + 2) >= lblob:
            return 4
        assert blob[seqend_i + 2] == ord(b'\n')
    # quality
    qualbeg_i = seqend_i+3
    if qualbeg_i >= lblob or seqend_i == -1:
        return 4
    else:
        posbuffer[4] = qualbeg_i
        qualend_i = qualbeg_i + (seqend_i - seqbeg_i)
        if qualend_i >= lblob:
            qualend_i = -1
        else:
            assert blob[qualend_i] == ord(b'\n')
        posbuffer[5] = qualend_i
    if qualend_i == -1:
        return 5
    else:
        return 6

def readfastq_iter(fh, fbufsize: int, _entrypos = _entrypos):
    """
    The entries in the FASTQ files are parsed from chunks of size `fbufsize`),
    using the function `_entrypos` (that be changed as a parameter - an
    faster implementation in C is in `fastqandfurious._fastqandfurious.entrypos`).

    With the current implementation, `fbufsize` must be large enough to contain
    the largest entry in the file. For example, is the longest read is 250bp long,
    and the identifier is 25-char long, the minimum buffer size will be about
    525. Aiming for that minimum is not advised, as some of the speed comes
    for working on buffers, and larger buffers able to contain many entries will
    lead to better performances (with the cave at that very large buffer might be
    counter-productive as the iterator will need to read data to fill the buffer
    (or all data, whichever is the smallest) before starting to yield entries. A
    value of 20,000 (20KB) lead pretty good results on this end.
 
    - fh: file-like object or stream (just needs a method `read`)
    - fbufsize: buffer size (see note above)
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
        headerbeg_i = 0
        seqbeg_i = 0
        seqend_i = 0
        qualend_i = 0
        while qualend_i < (lblob-1) and seqbeg_i >= 0 and seqend_i >= 0:
            npos = _entrypos(blob, offset, posbuffer)
            (headerbeg_i, headerend_i, seqbeg_i, seqend_i, qualbeg_i, qualend_i) = posbuffer
            if qualend_i == -1 or qualbeg_i == -1 or seqend_i == -1 or seqbeg_i == -1 or headerend_i == -1:
                #if headerbeg_i == -1:
                #    import pdb; pdb.set_trace()
                #    raise RuntimeError("Missing begining of entry ?")
                carryover = True
                backlog = blob[offset:] # bytes(mblob[offset:])
                break
            else:
                header = blob[headerbeg_i:headerend_i] # bytes(mblob[headerbeg_i:headerend_i])
                sequence = blob[seqbeg_i:seqend_i] # bytes(mblob[seqbeg_i:seqend_i])
                yield Entry(header, sequence)
                offset = qualend_i+1
        blob = fh.read(fbufsize)
        #mblob = memoryview(blob)
        lblob = len(blob)
        offset = 0
        if carryover:
            headerbeg_i = _nextentrypos(blob, backlog)
            if headerbeg_i == -1:
                raise RuntimeError("Incomplete last entry, or buffer too small.")
            # FIXME:
            tmp = backlog
            backlog = backlog + blob[:(headerbeg_i+1)]
            offset = headerbeg_i+1
            npos = _entrypos(backlog, 0, posbuffer)
            (headerbeg_i, headerend_i, seqbeg_i, seqend_i, qualbeg_i, qualend_i) = posbuffer
            if qualend_i == -1:
                # FIXME: handle this case
                raise RuntimeError("The buffer is too small !")
            sequence = backlog[seqbeg_i:seqend_i]
            header = backlog[headerbeg_i:headerend_i]
            yield Entry(header, sequence)
            backlog = b''
            carryover = False


