# fastq-and-furious

Efficient handling of FASTQ files(*) from Python ( *: no multi-line FASTQ though...)

The documentation is currently a little sparse. The docstring
for `fastqandfurious.fastqandfurious.readfastq_iter()` is best starting point.

There is a little utility to try it out on your own files.

```bash

python -m fastqandfurious.demo.benchmark <FASTQ or FASTQ.gz or FASTQ.bz2 file>

```

With a gzip-compressed FASTQ file of 60MB (size compressed) with 273639 entries,
the benchmark is
(the throughput is the for DNA sequences in the file - headers and quality strings
are not counted):


| parser | throughput |
|---|---|
| screed | 21.96MB/s |
| biopython | 9.83MB/ss |
| ngs_plumbing | 31.54MB/s |
| fastqandfurious | 36.73MB/s |
| fastqandfurious (C parts) | 49.73MB/s |
