
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


