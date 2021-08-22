Performance
-----------

There is a little utility to try it out on your own systems and files (there are options,
available with the flag `--help`).

The two mode are "speed" and "compare", the former benchmarking the speed of different
parsers and the second comparing the output of different parsers (not so good to be
fast if not correct).

Speed
^^^^^

.. code-block:: bash

   python -m fastqandfurious.demo.benchmark speed <FASTQ or FASTQ.gz or FASTQ.bz2 or FASTQ.lzma file>

Note that third-party library parsing FASTQ files are required in order to be able to run the full
benchmark.

With a gzip-compressed FASTQ file of 146MB (size compressed) with
1,562,120 entries,
the benchmark is
(the throughput is for the DNA sequences in the file - headers and quality strings
are not counted):

+----------------------------------+-------------------+--------------------------------------------+
|                           parser | throughput (MB/s) | notes                                      |
+==================================+===================+============================================+
|                           screed | 11.0              |                                            |
+----------------------------------+-------------------+--------------------------------------------+
|                        biopython |  6.8              |                                            |
+----------------------------------+-------------------+--------------------------------------------+
|   biopython FastqGeneralIterator | 34.5              | `Bio.SeqIO.QualityIO.FastqGeneralIterator` |
+----------------------------------+-------------------+--------------------------------------------+
|   pyfastx                        | 51.7              |                                            |
+----------------------------------+-------------------+--------------------------------------------+
|                  fastqandfurious | 32.0              | pure python                                |
+----------------------------------+-------------------+--------------------------------------------+
|         fastqandfurious w/ c-ext | 48.7              | using C extension in the package           |
+----------------------------------+-------------------+--------------------------------------------+
| fastqandfurious w/ c-ext + index | 37.7              | Like above and w/ index of entry positions |
+----------------------------------+-------------------+--------------------------------------------+


`fastqandfurious` with c-extension is 43% faster than Biopython's FastqGeneralIterator. The relatively recent `pyfastx`
is only 6% faster than `fastqandfurious` while at the cost of pretty much all Python-level flexibility in `fastqandfurious`.
For example, `fastqandfurious` can handle input from any Python `io` stream. This allows the use other compression algorithms other than gzip
(e.g., LZO-compression), or to data not in files (e.g., network streams, pipes, or interprocess communications).

Compare
^^^^^^^

To compare the output of two parsers, for example `biopython` and our parser:

.. code-block:: bash

   python -m fastqandfurious.demo.benchmark compare biopython fastqandfurious \
     <FASTQ | FASTQ.gz | FASTQ.bz2 | FASTQ.lzma>

