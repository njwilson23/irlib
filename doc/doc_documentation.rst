Documentation:
==============

TODO - check if we can run sphinx with sphinx dependency or need numpydoc

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

A copy can be found [here](https://dl.dropboxusercontent.com/u/375008/irlib_manual.pdf).  
TODO - replace this file link - or put pdf in git repo. 



