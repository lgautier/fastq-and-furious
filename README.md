# fastq-and-furious

(... because it turned out that the performance bottleneck for an
[algorithm-focused implementation of bottom-sketches (MinHash sketches)](https://github.com/lgautier/mashing-pumpkins)
was the parsing of FASTQ files).

[![Build Status](https://travis-ci.org/lgautier/fastq-and-furious.svg?branch=master)](https://travis-ci.org/lgautier/fastq-and-furious)

Efficient handling of FASTQ files(*) from Python ( *: no multi-line FASTQ though, in which case an exception will be raised when parsing).

## Installation

Python >= 3.5 (with headers) and a C99-compliant compiler are the only requirements.

```bash

pip install git+https://github.com/lgautier/fastq-and-furious.git

```

Should C extensions be unavaible (no C compiling possible, Pypy, other reasons), the Python-only version can still be used
(see section Performance below).


To run the tests after installation, one will need a clone of the repository. From the root of the repository, run:

```bash

python -m pytest --cov=fastqandfurious --cov-report xml --cov-report term tests.py

```

If `coverage` is not installed / is not working for you, the tests can be run without coverage analysis with:

```bash

python -m pytest tests.py

```


## Documentation

The documentation is currently a little sparse. The docstring
for `fastqandfurious.fastqandfurious.readfastq_iter()` is a good starting point.
The code for the benchmark (see below) is also a good source of information as
it can show how to use when compared to the other parser benchmarked.

## Performance

There is a little utility to try it out on your own systems and files (there are options,
available with the flag `--help`).

The two mode are "speed" and "compare", the former benchmarking the speed of different
parsers and the second comparing the output of different parsers (not so good to be
fast if not correct).

### Speed

```bash

python -m fastqandfurious.demo.benchmark speed <FASTQ or FASTQ.gz or FASTQ.bz2 or FASTQ.lzma file>

```

Note that third-party library parsing FASTQ files are required in order to be able to run the full
benchmark.

With a gzip-compressed FASTQ file of 60MB (size compressed) with 273,639 entries,
the benchmark is
(the throughput is for the DNA sequences in the file - headers and quality strings
are not counted):


| parser | throughput |
|---|---|
| screed | 21.96MB/s |
| biopython | 9.83MB/s |
| ngs_plumbing | 31.54MB/s |
| fastqandfurious (python-only) | 47.95MB/s |
| fastqandfurious (using C extension) | 62.81MB/s |


With a gzip-compressed FASTQ file of 700MB (size compressed) with 20,853,696 entries,
the benchmark is
(the throughput is for the DNA sequences in the file - headers and quality strings
are not counted):


| parser | throughput |
|---|---|
| screed | 3.70MB/s |
| biopython | 3.70MB/s |
| ngs_plumbing | 5.68MB/s |
| fastqandfurious (python-only) | 10.58MB/s |
| fastqandfurious (with C extension) | 19.64MB/s |


### Compare

To compare the output two parsers, for example `biopython` and our parser:

```bash

python -m fastqandfurious.demo.benchmark compare biopython fastqandfurious <FASTQ | FASTQ.gz | FASTQ.bz2 | FASTQ.lzma>

```