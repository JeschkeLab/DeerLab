Contributing
=========================================

First of all, thanks for considering contributing to DeerLab! New features that satisfy the needs of dipolar EPR 
spectroscopy are most welcome. You can either `open an issue <https://github.com/JeschkeLab/DeerLab/issues/new>`_
on GitHub to discuss your idea first, or if you have working code, submit a pull-request (PR).

Guidelines
------------------------------------------
When developing new DeerLab code, please make sure to follow these guidelines to ensure a smooth development of DeerLab. While these guidelines are not strictly controlled in the automated workflows used to check PRs, if not fulfilled the PR can be denied by a reviewer.

General 
    - Changes to the numerical or mathematical structure must be scientifically justified.
    - Code must be written in Python. Avoid interfacing different programming languages.
    - All mathematical operations must be performed with the Numpy/Scipy package(s). Do not use the Pandas framework.
    - Avoid the use of object-oriented programming (i.e. classes) unless strictly necessary. DeerLab is mainly based on functional programming.
    - Avoid introducing new dependencies to DeerLab (e.g. in the form of new Python packages). 
    - When using a third-party open-source software make sure it is compatible with DeerLab's MIT license.  

Testing
    - New features must be accompanied by unit tests. Changes to existing functionality must pass through the old tests.  
    - Use the ``assert`` function as logical unit of the tests.
    - Add comments to tests to describe functionality being tested. 
    - When testing against reference values, justify or provide their origin.

Documentation
    - Make sure the code is fully documented and commented, especially parts of the code that might be difficult to understand for beginner users.
    - Comment smartly and concisely throughout the file.
    - For new functions add a docstring as header (see other DeerLab functions for examples), describing functionality, inputs and outputs.


Development workflow
------------------------------------------

Here is how to set up DeerLab in a development environent and the basic workflow for developers.

1. Fork this repo, then clone your fork locally.

2. (Optional) Create a new Python virtual environment and activate it.

3. Change to the directory of your clone. All commands are run in the top level directory where the ``setup.py`` file is located.

3. Install DeerLab in development mode with all development dependencies::

    python -m setup.py develop

4. Run the test suite with ``pytest`` to check that the fresh installation is working correctly::

    pytest

5. Start developing.

6. (Optional) Before proposing new changes run ``pytest`` again to make sure your functionality passes the original tests (and the new ones you have added).

7. Push you changes to your forked repository and open a new PR on GitHub. Your changes will be again automatically
   tested, and if all tests pass then it will require the approval of DeerLab's core development team to finally merge
   your contribution into the source code.

Contributors
-------------

This is a list of GitHub users that have contributed to DeerLab, ordered by number of commits. 

..  ghcontributors:: JeschkeLab/DeerLab