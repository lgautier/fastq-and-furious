# fastq-and-furious

(... because it turned out that the performance bottleneck for an
[algorithm-focused implementation of bottom-sketches (MinHash sketches)](https://github.com/lgautier/mashing-pumpkins)
was the parsing of FASTQ files).

[![Python package](https://github.com/lgautier/fastq-and-furious/actions/workflows/pythonpackage.yml/badge.svg)](https://github.com/lgautier/fastq-and-furious/actions/workflows/pythonpackage.yml)
[![Pypi release](https://img.shields.io/pypi/v/fastq-and-furious.svg)](https://img.shields.io/pypi/v/fastq-and-furious.svg)

Efficient handling of FASTQ files(*) from Python ( *: no multi-line FASTQ though, in which case an exception will be raised when parsing).

<img src="throughput.svg">


[Project page](https://lgautier.github.io/fastq-and-furious/)