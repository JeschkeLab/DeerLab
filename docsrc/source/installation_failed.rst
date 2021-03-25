:orphan:

.. _installation_failed:

====================
Installation failed
====================

Known Issues #1: DLL load failed
--------------------------------


On a **Windows** computer, if you are trying to run a DeerLab function, you might get the following message:

.. code-block:: text

    ImportError: DLL load failed: The specified module could not be found.

This happens when the MKL libraries have not been properly linked in ``numpy``, ``scipy`` or ``cvxopt`` installations (typically ``numpy``). This can happen due to firewall restrictions, user rights, or internet connection issues during the DeerLab installation. To solve this, the best solution is to manually install as follows. 

1) Go to https://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy

2) Download the appropiate ``numpy`` wheels file according to your installed Python version and Windows system:

.. code-block:: text

                Python version (3.x)
    Package name          |         Windows architecture (32-64 bit)
       |                  |          |
       v                  v          v
    numpy-1.19.1+mkl-cp36-cp36m-win_amd64.whl


3) Once downloaded, open a terminal at the location of the ``.whl`` file and run the command:

.. code-block:: text

    python -m pip install "numpy-1.19.1+mkl-cp36-cp36m-win_amd64.whl"


making sure that the name of the ``.whl`` file matches the one that you downloaded.

This will install ``numpy`` and properly link all MKL DLL files. DeerLab should work now. Should the error persists, repeat this process for the ``scipy`` and ``cvxopt`` packages (in that order).
