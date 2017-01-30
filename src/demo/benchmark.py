import time

def benchmark_faf(fh, bufsize: int = 10000):
    from fastqandfurious import fastqandfurious
    total = int(0)
    t0 = time.time()
    it = fastqandfurious.readfastq_iter(fh, bufsize)
    for i, (header, sequence) in enumerate(it):
        total += len(sequence)
        if i % 10 == 0:
            t1 = time.time()
            print('\r%.2fMB/s' % (total/(1E6)/(t1-t0)), end='', flush=True)
    print()
    print('%i entries' % (i+1))

def benchmark_faf_c(fh, bufsize: int = 10000):
    from fastqandfurious import fastqandfurious, _fastqandfurious
    total = int(0)
    t0 = time.time()
    it = fastqandfurious.readfastq_iter(fh, bufsize, _entrypos=_fastqandfurious.entrypos)
    try:
        for i, (header, sequence) in enumerate(it):
            total += len(sequence)
            if i % 10 == 0:
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
        if i % 5000 == 0:
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
        if i % 5000 == 0:
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
        if i % 5000 == 0:
            t1 = time.time()
            print('\r%.2fMB/s' % (total/(1E6)/(t1-t0)), end='', flush=True)
    print()
    print('%i entries' % (i+1))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-screed',
                        action='store_true',
                        help='Do not test "screed"')        
    parser.add_argument('--no-ngs_plumbing',
                        action='store_true',
                        help='Do not test "ngs_plumbing"')        
    parser.add_argument('--no-biopython',
                        action='store_true',
                        help='Do not test "biopython"')        
    parser.add_argument('filename',
                        help='name of a FASTQ file (optionally gzip or bzip2-compressed)')

    args = parser.parse_args()

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
    if not args.no_ngs_plumbing:
        lst.append(('ngs_plumbing', benchmark_ngsplumbing, 'rb'))
    lst.append(('fastqandfurious', benchmark_faf, 'rb'))
    lst.append(('fastqandfurious (C parts)', benchmark_faf_c, 'rb'))
    
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

