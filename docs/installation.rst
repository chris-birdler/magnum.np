############
Installation
############

********************************
from Python Package Index (PyPi)
********************************
For a clean and independent system, we start with a clean virtual python environment (this step could be omitted, if you would like to install magnum.np the available python environment)

.. code-block:: bash

  mkdir venv
  python -m venv venv
  source venv/bin/activate

Finally install a release versions of magnum.np by means of pip:

.. code-block:: bash

  pip install magnumnp


*****************************
from source code (gitlab.com)
*****************************
More advanced users can also install magnum.np from source code.
It can be downloaded from https://gitlab.com/magnum.np/magnum.np .

After activating the virtual environment magnum.np can be simply installed using pip. For example
installing with the -e option also allows to modify the source code:

.. code-block:: bash

  git clone https://gitlab.com/magnum.np/magnum.np
  cd magnum.np
  pip install -e .

Note that a default version of `pytorch <http://www.pytorch.org/>`__ is included in magnum.np's dependecy list. If you would like to uses a specific pytorch version (fitting your installed CUDA library) it needs to be installed in advance.

