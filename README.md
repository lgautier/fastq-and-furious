# fastq-and-furious

Efficient handling of FASTQ files(*) from Python ( *: no multi-line FASTQ though...)


There is a little utility to try it out on your own files.

```bash

python -m fastqandfurious.demo.benchmark <FASTQ or FASTQ.gz or FASTQ.bz2 file>

```

With a gzip-compressed FASTQ file of 60MB (size compressed), the benchmark is
(the throughput is the for DNA sequences in the file - headers and quality strings
are not counted):


```
---
screed:
21.96MB/s
273639 entries
---
biopython
9.83MB/ss
273639 entries
---
ngs_plumbing
31.54MB/s
273639 entries
---
fastqandfurious
36.73MB/s
273639 entries
---
fastqandfurious (C parts)
49.73MB/s
273639 entries
```