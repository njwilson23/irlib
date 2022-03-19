Documentation
=============

In addition to the basic information here, documentation can be found in `doc`.
In order to build the documentation, [Sphinx](http://sphinx-doc.org/) must be
installed, with the ``numpydoc`` extension. The extensions may be installed by

    conda install sphinx

or

    conda install numpydoc

Then, from the ``doc/`` directory, type

    make html

If LaTeX is available, the documentation can be compiled into a PDF. Type

    make latexpdf

A pdf copy of the current documentation can be found on the project github site



