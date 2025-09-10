Installation
========================

Install the Mamba environment
-----------------------------

You may skip this step if your Mamba environment has been installed already.

Step 1: Download the installation script for Miniconda3
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

macOS (Intel)
'''''''''''''

.. code-block:: bash

  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

macOS (Apple Silicon)
'''''''''''''''''''''''

.. code-block:: bash

  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh

Linux
''''''''

.. code-block:: bash

  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

Step 2: Install Miniconda3
"""""""""""""""""""""""""""

.. code-block:: bash

  chmod +x Miniconda3-latest-*.sh && ./Miniconda3-latest-*.sh

During the installation, a path :code:`<base-path>` needs to be specified as the base location of the Python environment.
After the installation is done, add the following lines into your shell environment (e.g., :code:`~/.bashrc` or :code:`~/.zshrc`) to enable the :code:`mamba` package manager (remember to replace :code:`<base-path>` with your actual installation path):

.. code-block:: bash

  export PATH="<base-path>/bin:$PATH"
  . <base-path>/etc/profile.d/conda.sh

Step 3: Test your Installation
"""""""""""""""""""""""""""""""

.. code-block:: bash

  source ~/.bashrc  # assume you are using Bash shell
  which python      # should return a path under <base-path>
  which mamba       # should return a path under <base-path>


Install Python version of `CRO`
--------------------------------

For a clean installation, first create a new environment named :code:`env_pyCRO` via :code:`mamba`:

.. code-block:: bash

    mamba create -n env_pyCRO python=3.12   # supports Python 3.12 and onwards
    mamba activate env_pyCRO

Then install required dependencies via :code:`mamba`:

.. code-block:: bash

    mamba install -c conda-forge numpy xarray statsmodels netcdf4

Once the above dependencies are installed, simply install pyCRO via :code:`pip`:

.. code-block:: bash

    pip install pyCRO

You are now ready to import pyCRO in Python:

.. code-block:: python

    from pyCRO import *

