.. fastq-and-furious documentation master file, created by
   sphinx-quickstart on Sat Aug 21 11:38:28 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

fastq-and-furious
=================

Probably about as fast as FASTQ parsing in Python can be, but in a toolkit to build solutions.

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Installation
------------

Python >= 3.5 (with headers) and a C99-compliant compiler are the only requirements.

.. code-block:: bash

   pip install fastq-and-furious


To install the development version do

.. code-block:: bash

   pip install git+https://github.com/lgautier/fastq-and-furious.git

The package contains an optional C-extension module. That C code is how the best performance can be obtained with CPython,
but the Python-only version has the same API to allow the use of Pypy.


Running the tests
^^^^^^^^^^^^^^^^^

To run the tests after installation, one will need a clone of the repository. From the root of the repository, run:

.. code-block:: bash

   python -m pytest --cov=fastqandfurious --cov-report xml \
          --cov-report term tests.py

If :mod:`coverage` is not installed / is not working for you, the tests can be run without coverage analysis with:

.. code-block:: bash

   python -m pytest tests.py



Documentation
-------------

This documentation is work-in-progress. The docstring
for :func:`fastqandfurious.fastqandfurious.readfastq_iter` contains complementary information.
The code for the benchmark (see below) is also a good source of information as
it can show how to use when compared to other parsers benchmarked.


General idea
^^^^^^^^^^^^

In a nutshell, the package provide function that takes as input any Python IO stream,
such as one obtained by opening a file or any compressed file or network stream Python
has an IO stream for, a buffersize, a callback function
to build a entry from the elements in the stream parsed, and an optional callback
function to find the positions of FASTQ entry elements. That function will return
an iterator that will yield entries as `entryfunc` builds them.

For example, with a simple uncompressed file and an `entryfunc` defined in the package:

.. code-block:: python

   import fastqandfurious.fastqandfurious as fqf
   from fastqandfurious.fastqandfurious import entryfunc

   bufsize = 20000
   with open("a/fastq/file.fq", "rb") as fh:
       it = fqf.readfastq_iter(fh, bufsize, entryfunc)
       for sequence in it:
           # do something
	   pass

The file's compression is decoupled from parsing, which allows us to have generic code.
For example, parsing FASTQ data in a gzip-compressed file works the same way:

.. code-block:: python

   import gzip
   import fastqandfurious.fastqandfurious as fqf
   from fastqandfurious.fastqandfurious import entryfunc

   with gzip.open("a/fastq/file.fq", "rb") as fh:
       it = fqf.readfastq_iter(fh, bufsize, entryfunc)
       for entry in it:
           # do something
	   pass


Writing an extension automagic reader name-based reader: 

.. code-block:: python

   import gzip
   import fastqandfurious.fastqandfurious as fqf
   from fastqandfurious.fastqandfurious import entryfunc

   for filename in ('file_a.fq', 'file_b.fq.gz',
                    'file_c.fq.bz2', 'file_d.fq.lzma'):
       with fqf.automagic_open(filename) as fh:
           it = fqf.readfastq_iter(fh, bufsize, entryfunc)
           for entry in it:
               # do something
	       pass


.. autofunction:: fastqandfurious.fastqandfurious.readfastq_iter


Faster with C
^^^^^^^^^^^^^

The optional C-extension has the same API, but is faster. It is made optional
to allow the use of :mod:`fastqandfurious` with Pypy or in other situations
where building the C-extension in not possible.

For example, reading the same FASTQ file, first without the C-extension speeding
the code and second with the C-extension:

.. code-block:: python
	
   import fastqandfurious.fastqandfurious as fqf
   import fastqandfurious._fastqandfurious as _fqf
   from fastqandfurious.fastqandfurious import entryfunc

   bufsize = 20000

   for entrypos in (fqf.entrypos, _fqf.entrypos):
       with open("a/fastq/file.fq", "rb") as fh:
           it = fqf.readfastq_iter(fh, bufsize, entryfunc, entrypos=fqf.entrypos))
           for sequence in it:
               # do something
	       pass


Instant faster code (just add <strike>water</strike> adapter code)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

That design also lets us just drop the parser into an existing code base, or keep working
with a library you are most familiar with, writing a short adapter (and observe
immediate performance gains - see benchmark below).

For example, to existing codebase using biopython's :class:`SeqRecord` one only needs
to provide an `entryfunc` that builds entries accordingly:

.. code-block:: python

   import fastqandfurious.fastqandfurious as fqf
   from fastqandfurious._fastqandfurious import arrayadd_b
   from Bio.SeqRecord import SeqRecord
   from array import array

   def biopython_entryfunc(buf, posarray, globaloffset):
       name = buf[posarray[0]:posarray[1]].decode('ascii')
       quality = array('b')
       quality.frombytes(buf[posarray[4]:posarray[5]])
       arrayadd_b(quality, -33)
       entry = SeqRecord(
                   seq=buf[posarray[2]:posarray[3]].decode('ascii'),
                   id=name,
                   name=name,
                   letter_annotations={'phred_quality': quality}
	       )
       return entry


   bufsize = 20000
   with open("a/fastq/file.fq", "rb") as fh:
       it = fqf.readfastq_iter(fh, bufsize,
                               biopython_entryfunc)
       for entry in it:
           # do something
	   pass


Performance by design
^^^^^^^^^^^^^^^^^^^^^

The design is obviously also offering various performance gains by allowing to only build entry components
as needed. Whenever entries in a FASTQ file are filtered out is significant number this can reduce
the overhead of creating instances for each entry to delete them soon after.

For example, writing a filter on read length could be done with:

.. code-block:: python

   LENGTH_THRESHOLD = 25

   def lengthfilter_entryfunc(buf, posarray):
       if posarray[3] - posarray[2] < LENGTH_THRESHOLD:
           return buf[posarray[2]:posarray[3]]
       else:
           return None

   with open("a/fastq/file.fq", "rb") as fh:
       it = fqf.readfastq_iter(fh, bufsize,
                               lengthfilter_entryfunc)
       for sequence in it:
           if sequence is None:
	       # do nothing
	   else:
               # do something
	       pass

Fetching the positions for the elements of an entry (name/ID, sequence, quality) is also allowing us
to store the positions associated with a FASTQ for future use (see the code for `fastqandfurious_c_index`
in `fastqandfurious.demo.benchmark`).

This is essentially like storing a table of positions:

+----------+----------+---------+---------+-------------+-------------+
| name_beg | name_end | seq_beg | seq_end | quality_beg | quality_end |
+==========+==========+=========+=========+=============+=============+
| 0        | 20       | 22      | 172     | 24          | 274         |
+----------+----------+---------+---------+-------------+-------------+
| 275      | 295      | 296     | 446     | 448         | 598         |
+----------+----------+---------+---------+-------------+-------------+

Whenever the FASTQ must be used again, that table can be used to quickly extract
data elements without having to parse them. Where :mod:`fastqandfurious` shines
is that there is complete flexibility about how to store such indexes.

That approach opens the door for implementing masking strategies to avoid
saving a FASTQ file after each filtering or read-trimming step. Excluding reads or trimming either end of the read
could be done by only deleting rows or modifying the values in a table of indices. Assuming that int64 is required for the indices, each entry (read) takes 6x8=48 bytes uncompressed.
In comparison each 120 base read (sequence + quality scores) takes
over 250 bytes uncompressed.

If having the quality as a sequence of integer ajusted for the eventual offset of 33, there is a also a C utility:

.. code-block:: python

   from fastqandfurious._fastqandfurious import arrayadd_b
   from array import array
   quality = array('b')
   quality.frombytes(buf[posarray[4]:posarray[5]])
   arrayadd_b(quality, -33)

The design is obviously also offering various performance gains by allowing to only build entry components
as needed. For example, writing a filter on read length could be done with:

.. code-block:: python

   def lengthfilter_entryfunc(buf, posarray):
       LENGTH_THRESHOLD = 25
       if posarray[3] - posarray[2] < LENGTH_THRESHOLD:
           return buf[posarray[2]:posarray[3]]
       else:
           return None

   with open("a/fastq/file.fq", "rb") as fh:
       it = fastqandfurious.readfastq_iter(fh, bufsize,
                                           lengthfilter_entryfunc)
       for sequence in it:
           if sequence is None:
	       # do nothing
	       pass
	   else:
               # do something
	       pass

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


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
