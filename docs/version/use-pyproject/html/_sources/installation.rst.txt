
Installation
------------

Python >= 3.5 (with headers) and a C99-compliant compiler are the only requirements.

.. code-block:: bash

   pip install fastq-and-furious


To install the development version do

.. code-block:: bash

   pip install git+https://github.com/lgautier/fastq-and-furious.git


Running the tests
^^^^^^^^^^^^^^^^^

To run the tests after installation, one will need a clone of the repository. From the root of the repository, run:

.. code-block:: bash

   python -m pytest --cov=fastqandfurious --cov-report xml \
          --cov-report term tests.py

If :mod:`coverage` is not installed or is not working for you, the tests can be run without coverage analysis with:

.. code-block:: bash

   python -m pytest tests.py

