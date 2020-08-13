============
Installation
============

Requirements
------------

* python>=3.7


pip
---

.. code-block:: bash

    git clone git@github.com:lan496/dsenum.git
    cd dsenum
    pip install -r requirements.txt
    pip install -e .


conda
-----

.. code-block:: bash

    git clone git@github.com:lan496/dsenum.git
    cd dsenum
    conda env create -f environment.yml
    conda activate dsenum
    pip install -e .

for development
---------------

.. code-block:: bash

    git clone git@github.com:lan496/dsenum.git
    cd dsenum
    conda env create -f environment.yml
    pip install -r requirements-dev.txt
    pip install -e .
    pre-commit install
