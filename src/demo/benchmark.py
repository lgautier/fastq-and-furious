import time

def benchmark_faf(fh, bufsize: int = 50000):
    from fastqandfurious import fastqandfurious
    total = int(0)
    t0 = time.time()
    it = fastqandfurious.readfastq_iter(fh, bufsize)
    for i, e in enumerate(it):
        total += len(e.sequence)
        if i % 20000 == 0:
            t1 = time.time()
            print('\r%.2fMB/s' % (total/(1E6)/(t1-t0)), end='', flush=True)
    print()
    print('%i entries' % (i+1))

def benchmark_faf_c(fh, bufsize: int = 50000):
    from fastqandfurious import fastqandfurious, _fastqandfurious
    total = int(0)
    t0 = time.time()
    it = fastqandfurious.readfastq_iter(fh, bufsize, _entrypos=_fastqandfurious.entrypos)
    try:
        for i, e in enumerate(it):
            total += len(e.sequence)
            if i % 20000 == 0:
                t1 = time.time()
                print('\r%.2fMB/s' % (total/(1E6)/(t1-t0)), end='', flush=True)
    finally:
        print()
        print('%i entries' % (i+1))

def benchmark_ngsplumbing(fh):
    import ngs_plumbing.fastq
    total = int(0)
    t0 = time.time()
    it = ngs_plumbing.fastq.read_fastq(fh)
    for i, e in enumerate(it):
        total += len(e.sequence)
        if i % 20000 == 0:
            t1 = time.time()
            print('\r%.2fMB/s' % (total/(1E6)/(t1-t0)), end='', flush=True)
    print()
    print('%i entries' % (i+1))

def benchmark_screed(fn):
    import screed
    total = int(0)
    t0 = time.time()
    it = screed.open(fn)
    for i, e in enumerate(it):
        total += len(e.sequence)
        if i % 20000 == 0:
            t1 = time.time()
            print('\r%.2fMB/s' % (total/(1E6)/(t1-t0)), end='', flush=True)
    print()
    print('%i entries' % (i+1))

def benchmark_biopython(fh):
    from Bio import SeqIO
    total = int(0)
    t0 = time.time()
    it = SeqIO.parse(fh, "fastq")
    for i, e in enumerate(it):
        total += len(e.seq)
        if i % 20000 == 0:
            t1 = time.time()
            print('\r%.2fMB/s' % (total/(1E6)/(t1-t0)), end='', flush=True)
    print()
    print('%i entries' % (i+1))


def benchmark_biopython_adapter(fh):
    total = int(0)
    t0 = time.time()

    from fastqandfurious import fastqandfurious
    from Bio.SeqRecord import SeqRecord

    def biopython_entryfunc(buf, posarray):
        name = buf[posarray[0]:posarray[1]].decode('ascii')
        entry = SeqRecord(seq=buf[posarray[2]:posarray[3]].decode('ascii'),
                          id=name,
                          name=name)
        return entry

    bufsize = 20000
    it = fastqandfurious.readfastq_iter(fh, bufsize, biopython_entryfunc)
    for i, e in enumerate(it):
        total += len(e.seq)
        if i % 20000 == 0:
            t1 = time.time()
            print('\r%.2fMB/s' % (total/(1E6)/(t1-t0)), end='', flush=True)
    print()
    print('%i entries' % (i+1))

    
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
        lst.append(('biopython_adapter', benchmark_biopython_adapter, 'rb'))
    if not args.no_ngs_plumbing:
        lst.append(('ngs_plumbing', benchmark_ngsplumbing, 'rb'))
    lst.append(('fastqandfurious', benchmark_faf, 'rb'))
    lst.append(('fastqandfurious (w/ C-ext)', benchmark_faf_c, 'rb'))
    
    for name, func, mode in lst:
        print('---')
        print(name)
        if args.filename.endswith('.gz'):
            import gzip
            with gzip.open(args.filename, mode) as fh:
                try:
                    func(fh)
                except Exception as e:
                    print('Error: %s' % str(e))
        elif args.filename.endswith('.bz2'):
            import bz2
            with bz2.open(args.filename, mode) as fh:
                try:
                    func(fh)
                except Exception as e:
                    print('Error: %s' % str(e))        
        elif args.filename.endswith('.lzma'):
            import lzma
            with lzma.open(args.filename, mode) as fh:
                try:
                    func(fh)
                except Exception as e:
                    print('Error: %s' % str(e))        



def _opener(filename, mode):
    if filename.endswith('.gz'):
        import gzip
        return gzip.open(args.filename, mode)
    elif filename.endswith('.bz2'):
        import bz2
        return bz2.open(args.filename, mode)
    elif filename.endswith('.lzma'):
        import lzma
        return lzma.open(args.filename, mode)

def _screed_iter(fn):
    import screed
    it = screed.open(fn)
    for i, e in enumerate(it):
        yield (i, b'@'+e.name.encode('ascii'), str(e.sequence).encode('ascii'))

def _biopython_iter(fn, mode):
    from Bio import SeqIO
    with _opener(fn, mode) as fh:
        for i, e in enumerate(SeqIO.parse(fh, "fastq")):
            yield (i, b'@'+e.description.encode('ascii'), str(e.seq).encode('ascii'))

def _ngs_plumbing_iter(fn, mode):
    import ngs_plumbing.fastq
    with _opener(fn, mode) as fh:
        it = ngs_plumbing.fastq.read_fastq(fh)
        for i, e in enumerate(it):
            yield (i, e.header, e.sequence)

def _fastqandfurious_iter(fn, mode):
    from fastqandfurious import fastqandfurious
    bufsize = int(5E4)
    with _opener(fn, mode) as fh:
        it = fastqandfurious.readfastq_iter(fh, bufsize)
        for i, e in enumerate(it):
            yield (i, e.header, e.sequence)

def _fastqandfurious_c_iter(fn, mode):
    from fastqandfurious import fastqandfurious, _fastqandfurious
    bufsize = int(5E4)
    with _opener(fn, mode) as fh:
        it = fastqandfurious.readfastq_iter(fh, bufsize, _entrypos=_fastqandfurious.entrypos)
        for i, e in enumerate(it):
            yield (i, e.header, e.sequence)

def run_compare(args):
    res = list()
    for name in args.parser:
        if name == 'screed':
            res.append(_screed_iter(args.filename))
        elif name == 'biopython':
            res.append(_biopython_iter(args.filename, 'rt'))
        elif name == 'ngs_plumbing':
            res.append(_ngs_plumbing_iter(args.filename, 'rb'))
        elif name == 'fastqandfurious':
            res.append(_fastqandfurious_iter(args.filename, 'rb'))
        elif name == 'fastqandfurious_c':
            res.append(_fastqandfurious_c_iter(args.filename, 'rb'))
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
