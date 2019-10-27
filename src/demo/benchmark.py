import time
import io

REFRESH_RATE=int(1E5)


def benchmark_faf(fh, bufsize: int = int(2**16)):
    from fastqandfurious import fastqandfurious
    total_seq = int(0)
    t0 = time.time()
    it = fastqandfurious.readfastq_iter(fh, bufsize)
    for i, e in enumerate(it):
        total_seq += len(e[1])
        if i % REFRESH_RATE == 0:
            t1 = time.time()
            print('\r%.2fMB/s' % (total_seq/(1E6)/(t1-t0)), end='', flush=True)
    print()
    print('%i entries in %.3f seconds.' % (i+1, time.time()-t0))


def benchmark_faf_c(fh, bufsize: int = int(2**16)):
    from fastqandfurious import fastqandfurious, _fastqandfurious
    total_seq = int(0)
    t0 = time.time()
    it = fastqandfurious.readfastq_iter(fh, bufsize, _entrypos=_fastqandfurious.entrypos)
    try:
        for i, e in enumerate(it):
            total_seq += len(e[1])
            if i % REFRESH_RATE == 0:
                t1 = time.time()
                print('\r%.2fMB/s' % (total_seq/(1E6)/(t1-t0)), end='', flush=True)
    finally:
        print()
        print('%i entries in %.3f seconds.' % (i+1, time.time()-t0))


def benchmark_faf_c_index(fh, fh_index, bufsize: int = int(2**16)):
    from fastqandfurious import fastqandfurious, _fastqandfurious
    from array import array
    total_seq = int(0)
    t0 = time.time()
    fh = io.BufferedReader(fh, buffer_size = bufsize*8)
    try:
        offset = 0
        entry_i = 0
        buf = None
        while True:
            posarray = array('q')
            try:
                posarray.fromfile(fh_index, 6)
            except EOFError as eofe:
                break
            buf = fh.read(posarray[-1]-offset+1)
            tmp = posarray[-1]
            _fastqandfurious.arrayadd_q(posarray, -offset)
            e = (buf[posarray[0]:posarray[1]],
                 buf[posarray[2]:posarray[3]],
                 buf[posarray[4]:posarray[5]])
            offset = tmp+1
            entry_i += 1
            total_seq += len(e[1])
            if entry_i % REFRESH_RATE == 0:
                t1 = time.time()
                print('\r%.2fMB/s' % (total_seq/(1E6)/(t1-t0)), end='', flush=True)
    finally:
        print()
        print('%i entries in %.3f seconds.' % (entry_i+1, time.time()-t0))


def benchmark_ngsplumbing(fh):
    import ngs_plumbing.fastq
    total_seq = int(0)
    t0 = time.time()
    it = ngs_plumbing.fastq.read_fastq(fh)
    for i, e in enumerate(it):
        total_seq += len(e.sequence)
        if i % REFRESH_RATE == 0:
            t1 = time.time()
            print('\r%.2fMB/s' % (total_seq/(1E6)/(t1-t0)), end='', flush=True)
    print()
    print('%i entries in %.3f seconds.' % (i+1, time.time()-t0))


def benchmark_screed(fn):
    import screed
    total_seq = int(0)
    t0 = time.time()
    it = screed.open(fn)
    for i, e in enumerate(it):
        total_seq += len(e.sequence)
        if i % REFRESH_RATE == 0:
            t1 = time.time()
            print('\r%.2fMB/s' % (total_seq/(1E6)/(t1-t0)), end='', flush=True)
    print()
    print('%i entries' % (i+1))


def benchmark_biopython(fh):
    from Bio import SeqIO
    total_seq = int(0)
    t0 = time.time()
    it = SeqIO.parse(fh, "fastq")
    for i, e in enumerate(it):
        total_seq += len(e.seq)
        if i % REFRESH_RATE == 0:
            t1 = time.time()
            print('\r%.2fMB/s' % (total_seq/(1E6)/(t1-t0)), end='', flush=True)
    print()
    print('%i entries in %.3f seconds.' % (i+1, time.time()-t0))


def benchmark_biopython_faster(fh):
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    total_seq = int(0)
    t0 = time.time()
    it = FastqGeneralIterator(fh)
    for i, (title, seq, qual) in enumerate(it):
        total_seq += len(seq)
        if i % REFRESH_RATE == 0:
            t1 = time.time()
            print('\r%.2fMB/s' % (total_seq/(1E6)/(t1-t0)), end='', flush=True)
    print()
    print('%i entries in %.3f seconds.' % (i+1, time.time()-t0))


def benchmark_biopython_adapter(fh):
    total_seq = int(0)
    t0 = time.time()

    from fastqandfurious import fastqandfurious
    from fastqandfurious._fastqandfurious import arrayadd_b
    from Bio.SeqRecord import SeqRecord
    from array import array

    def biopython_entryfunc(buf, posarray, globaloffset):
        name = buf[posarray[0]:posarray[1]].decode('ascii')
        quality = array('b')
        quality.frombytes(buf[posarray[4]:posarray[5]])
        arrayadd_b(quality, -33)
        entry = SeqRecord(seq=buf[posarray[2]:posarray[3]].decode('ascii'),
                          id=name,
                          name=name,
                          letter_annotations={'phred_quality': quality})
        return entry

    bufsize = 20000
    it = fastqandfurious.readfastq_iter(fh, bufsize, biopython_entryfunc)
    for i, e in enumerate(it):
        total_seq += len(e.seq)
        if i % REFRESH_RATE == 0:
            t1 = time.time()
            print('\r%.2fMB/s' % (total_seq/(1E6)/(t1-t0)), end='', flush=True)
    print()
    print('%i entries in %.3f seconds.' % (i+1, time.time()-t0))


def _opener(filename):
    if filename.endswith('.gz'):
        import gzip
        return gzip.open
    elif filename.endswith('.bz2'):
        import bz2
        return bz2.open
    elif filename.endswith('.lzma'):
        import lzma
        return lzma.open
    else:
        return open


def run_speed(args):

    print('Running benchmark on file %s' % args.filename)

    if not args.no_screed:
        print('---')
        print('screed:')
        try:
            benchmark_screed(args.filename)
        except Exception as e:
            print('Error: %s' % str(e))

    lst = list()
    if not args.no_biopython:
        lst.append(('biopython', benchmark_biopython, 'rt'))
        lst.append(('biopython_fastqiterator', benchmark_biopython_faster, 'rt'))
        if args.with_biopython_adapter:
            lst.append(('biopython_adapter', benchmark_biopython_adapter, 'rb'))
    if not args.no_ngs_plumbing:
        lst.append(('ngs_plumbing', benchmark_ngsplumbing, 'rb'))
    if not args.no_fastqandfurious_python:
        lst.append(('fastqandfurious', benchmark_faf, 'rb'))
    lst.append(('fastqandfurious (w/ C-ext)', benchmark_faf_c, 'rb'))
    lst.append(('fastqandfurious (w/ C-ext and indexing)', benchmark_faf_c_index, 'rb'))

    for name, func, mode in lst:
        print('---')
        print(name)
        openfunc = _opener(args.filename)
        if name in ('biopython', 'biopython_fastqiterator'):
            with openfunc(args.filename, mode=mode) as fh:
                try:
                    func(fh)
                except Exception as e:
                    print('Error: %s' % str(e))
        elif name == 'fastqandfurious (w/ C-ext and indexing)':
            import tempfile
            from fastqandfurious import fastqandfurious, _fastqandfurious
            bufsize = int(5E4)
            with tempfile.NamedTemporaryFile(mode='r+b') as fh_index:
                with openfunc(args.filename, mode=mode) as fh:
                    print('  building index...', end='', flush=True)
                    it = fastqandfurious.readfastq_iter(fh, bufsize,
                                                        entryfunc=fastqandfurious.entryfunc_abspos,
                                                        _entrypos=_fastqandfurious.entrypos)
                    for i, pos in enumerate(it):
                        pos.tofile(fh_index)
                    fh_index.flush()
                    fh_index.seek(0)
                    print('done.')
                with openfunc(args.filename, mode=mode) as fh:
                    #try:
                    func(fh, fh_index)
                    #except Exception as e:
                    #    print('Error: %s' % str(e))
        else:
            with open(args.filename, mode, buffering = args.io_buffersize) as f:
                with openfunc(f) as fh:
                    #try:
                    if True:
                        func(fh)
                    #except Exception as e:
                    #    print('Error: %s' % str(e))


def _screed_iter(fn):
    import screed
    it = screed.open(fn)
    for i, e in enumerate(it):
        yield (i, e.name.encode('ascii'), str(e.sequence).encode('ascii'))


def _biopython_iter(fn, mode, buffering):
    from Bio import SeqIO
    openfunc = _opener(fn)
    with openfunc(fn, mode=mode) as fh: 
        for i, e in enumerate(SeqIO.parse(fh, "fastq")):
            yield (i, e.description.encode('ascii'), str(e.seq).encode('ascii'))


def _ngs_plumbing_iter(fn, mode, buffering):
    import ngs_plumbing.fastq
    openfunc = _opener(fn)
    with open(fn, mode, buffering = buffering) as f:
        with openfunc(f) as fh: 
            it = ngs_plumbing.fastq.read_fastq(fh)
            for i, e in enumerate(it):
                yield (i, e.header[1:], e.sequence)


def _fastqandfurious_iter(fn, mode, buffering):
    from fastqandfurious import fastqandfurious
    bufsize = int(5E4)
    openfunc = _opener(fn)
    with open(fn, mode, buffering = buffering) as f:
        with openfunc(f) as fh: 
            it = fastqandfurious.readfastq_iter(fh, bufsize)
            for i, e in enumerate(it):
                yield (i, e.header, e.sequence)


def _fastqandfurious_c_iter(fn, mode, buffering):
    from fastqandfurious import fastqandfurious, _fastqandfurious
    bufsize = int(5E4)
    openfunc = _opener(fn)
    with open(fn, mode, buffering = buffering) as f:
        with openfunc(f) as fh: 
            it = fastqandfurious.readfastq_iter(fh, bufsize, _entrypos=_fastqandfurious.entrypos)
            for i, e in enumerate(it):
                yield (i, e.header, e.sequence)


def run_compare(args):
    res = list()

    for name in args.parser:
        if name == 'screed':
            res.append(_screed_iter(args.filename))
        elif name == 'biopython':
            res.append(_biopython_iter(args.filename, 'rt', 2**16))
        elif name == 'ngs_plumbing':
            res.append(_ngs_plumbing_iter(args.filename, 'rb', 2**16))
        elif name == 'fastqandfurious':
            res.append(_fastqandfurious_iter(args.filename, 'rb', 2**16))
        elif name == 'fastqandfurious_c':
            res.append(_fastqandfurious_c_iter(args.filename, 'rb', 2**16))
        else:
            raise ValueError('Unknown parser name.')
    for i, (e1, e2) in enumerate(zip(*res), 0):
        if i % 5000 == 0:
            print('\rcompared %i entries' % (i), end='', flush=True)
        assert e1[0] == e2[0] # entry number
        assert e1[1] == e2[1] # header
        assert e1[2] == e2[2] # sequence
    print('\rcompared %i entries' % (i+1))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(help='sub-command help')
    
    parser_speed = subparsers.add_parser('speed', help = "benchmark for speed")
    
    parser_speed.add_argument('--no-screed',
                              action='store_true',
                              help='Do not test "screed"')        
    parser_speed.add_argument('--no-ngs_plumbing',
                              action='store_true',
                              help='Do not test "ngs_plumbing"')        
    parser_speed.add_argument('--no-biopython',
                              action='store_true',
                              help='Do not test "biopython"')
    parser_speed.add_argument('--no-fastqandfurious_python',
                              action='store_true',
                              help='Do not test "fastqandfurious" (Python-only)')
    parser_speed.add_argument('--with-biopython-adapter',
                              action='store_true',
                              help='Test with adapter for "biopython" (unless --no-biopython specified)')
    parser_speed.add_argument('--io-buffersize',
                              type = int,
                              default = int(50E3))
    parser_speed.add_argument('filename',
                              help='name of a FASTQ file (optionally gzip, bzip2, or lzma-compressed)')
    parser_speed.set_defaults(func=run_speed)


    parser_compare = subparsers.add_parser('compare', help = "compare output")     
    parser_compare.add_argument('parser',
                                nargs = 2,
                                choices = ('screed', 'biopython', 'ngs_plumbing', 'fastqandfurious', 'fastqandfurious_c'),
                                help='Parser to use.')        
    
    parser_compare.add_argument('filename',
                              help='name of a FASTQ file (optionally gzip, bzip2, or lzma-compressed)')

    parser_compare.set_defaults(func=run_compare)
    
    
    args = parser.parse_args()
    args.func(args)
