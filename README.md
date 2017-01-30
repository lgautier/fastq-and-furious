# fastq-and-furious

(... because it turned out that the performance bottleneck for an algorithm-focused implementation
of bottom-sketches (MinHash sketches) was the parsing of FASTQ files).

[![Build Status](https://travis-ci.org/lgautier/fastq-and-furious.svg?branch=master)](https://travis-ci.org/lgautier/fastq-and-furious)

Efficient handling of FASTQ files(*) from Python ( *: no multi-line FASTQ though...)

## Installation

Python >= 3.5 (with headers) and a C compiler are the only requirements.

```

pip install git+https://github.com/lgautier/fastq-and-furious.git

```


## Documentation

The documentation is currently a little sparse. The docstring
for `fastqandfurious.fastqandfurious.readfastq_iter()` is best starting point.
The code for the benchmark (see below) is also a good source of information as
it can show how to use when compared to the other parser benchmarked.

## Performance

There is a little utility to try it out on your own files (there are options,
available with the flag `--help`):

```bash

python -m fastqandfurious.demo.benchmark <FASTQ or FASTQ.gz or FASTQ.bz2 file>

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
| fastqandfurious | 36.73MB/s |
| fastqandfurious (using C parts) | 49.73MB/s |


With a gzip-compressed FASTQ file of 700MB (size compressed) with 20,853,696 entries,
the benchmark is
(the throughput is for the DNA sequences in the file - headers and quality strings
are not counted):


| parser | throughput |
|---|---|
| screed | 3.70MB/s |
| biopython | 3.70MB/s |
| ngs_plumbing | 5.68MB/s |
| fastqandfurious | 8.86MB/s |
| fastqandfurious (using C parts) | 15.12MB/s |
