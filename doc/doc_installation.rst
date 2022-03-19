Installation
============

Installation with Conda
-----------------------

Installing and setting up irlib is best managed in a conda environment.
Follow the conda instructions followed by Linux or Windows instructions
depending on your system. Mac should be similar to linux. The source
code is on github at https://github.com/njwilson23/irlib

There are 5 steps:

Setup a conda environment
~~~~~~~~~~~~~~~~~~~~~~~~~

These steps will set up and manage your Python environment and
dependencies for irlib.  For more on the dependencies themselves, see below.

-  Install Anaconda or Miniconda python 3 64-bit
-  Open conda prompt **as an administrator** or with **write
   permissions** to the conda directory and create an environment
   specifically to use irlib (Run one of these options):

-  This will be a bare-bones installation to run irlib

``conda create -n irlib -c conda-forge python h5py scipy matplotlib cython geopandas``

-  As above but also installs vitables, an hdf viewer and sphinx:

``conda create -n irlib -c conda-forge python h5py scipy matplotlib cython geopandas sphinx vitables numpydoc``

-  Or use the environment file in the repository:

``conda create -n environment.yml``

-  To run irlib you need to work out of a conda-aware console and type:

``conda activate irlib``


Make an irlib directory
~~~~~~~~~~~~~~~~~~~~~~~

This example shows how to create a folder called 'py' in your own home
directory with irlib files in it, however, you can install irlib wherever you
wish (including within the conda environment folders). Make a directory
called 'py' in your home folder. This will hold the irlib directory. Change 
directory into that folder.


Download irlib
~~~~~~~~~~~~~~

Go to https://github.com/njwilson23/irlib and download a
zip of latest irlib and unzip it in that folder. Alternative: if you
have git installed you can type the following in the terminal:

::

    git clone git@github.com:njwilson23/irlib.git

This makes a directory in your home folder called py/irlib Note: 
make a copy of the irlib zip file or directory for safekeeping. If you
start messing around, it's good to get the original back without any
fuss (mostly important in the field if you are not near internet).



Set the operating system path
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This step is not strictly necessary but needs to be done if you want to
type the irlib command line executables from any folder (including where
your data files are located).

**LINUX**

Instructions assume you are using Bash and you installed to
the directory from the example above.

Find the hidden .bashrc file and open it in an editor. At the
      bottom of this file type and then save the file:

| ``# Set path for irlib python scripts HERE``
| ``export PATH=$PATH:~/py/irlib``

Then, in a terminal, type the following to make the change
      permanent:

``source .bashrc``

**WINDOWS**

To add path for the current session:

``set PATH=%PATH%;C:\your\path\to\irlib\code``

To permanently add path, but not for current session:

``setx PATH=%PATH%;C:\your\path\to\irlib\code``

To view you current operating system path:

``echo %PATH%``


Set the conda environment path
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This must be done so Python can find the irlib libraries when in the
irlib conda environment.

-  Activate the conda environment and ensure that the irlib files are
   always available. In a terminal type:

 ``conda activate irlib``
 ``conda develop`` 

Where is where the irlib code is located

-  You may need to restart conda or reactivate irlib environment for
   this to take effect.
-  **NOTE** that if you are sharing your conda environment with other
   users they will be using that same version of irlib that you have
   specified!

Testing
~~~~~~~

Open a terminal, activate your irlib conda environment and type:

``h5_dumpmeta.py -h``

You should see the useage message starting like so: 

``usage: h5_dumpmeta.py [-h] [-o OUTFILE] [-c] [-w] [-l] [--clobber] 
[--swap_lon] [--swap_lat] infile``

Then see if it works with an h5 file (in this example it is called
'survey.h5'):

``h5_dumpmeta.py survey.h5``

It will output some metadata to the screen.

If that doesn't work, check your environment is activated, your paths
are set and make sure that the python files are executable.


Dependencies
------------
In this section the main irlib dependencies are listed and discussed. *If you installed 
with conda as above you should have these dependencies already and you don't need to 
read this section.* 

*radar_tools* is built upon a number of standard tools from the scientific
Python ecosystem. The following are *required*:

.. _Python: http://python.org/
.. _Numpy: http://www.numpy.org/
.. _Scipy: http://scipy.org/SciPy
.. _h5py: https://www.h5py.org/
.. _matplotlib: http://matplotlib.org/
.. _pandas: https://pandas.pydata.org/
.. _geopandas: https://geopandas.org/ 
.. _Cython: http://cython.org/
.. _Sphinx: https://www.sphinx-doc.org/en/master/
.. _Vitables: https://vitables.org/
.. _Github: http://njwilson23.github.com/radar_tools
.. _gstat: http://www.gstat.org/
.. _numpydoc: https://numpydoc.readthedocs.io/en/latest/index.html


- Python_ : Already installed for Linux/Mac OS X users

- Numpy_ : Basic array type, analogous to a matrix in MATLAB, except better

- Scipy_ : Wrappers for scientific libraries used for efficient filtering

- h5py_ : interface for HDF datasets

- matplotlib_ : Plotting library required for GUI tools

- pandas_ : Powerful Python data analysis toolkit

- geopandas_ : Python library that enables geopspatial data interchange. 

- Cython_ : Python compiler for improving performance

Finally, these are *nice to have*:

- Sphinx_ : Documentation generator library. 

- numpydoc_ : A sphinx extension containing styling.

- Vitables_ : An hdf viewer to look at the structure of h5 files more visually. (you can also use hdfview or another alternative)



Alternative installations
-------------------------

Using a package manager (e.g. APT, rpm, pacman, or Homebrew) download all the dependencies above


The latest version is on Github_. After downloading either directly or using the
command

::

    git clone git@github.com:njwilson23/irlib.git

Installation is best done with ``pip``, the Python package manager.

::

    >> cd irlib/    # or wherever it's downloaded to
    >> pip install .

Assuming that dependencies are available (see above), this will take care of
installing ``radar_tools`` properly. 

To use the *pywavelet* wavelet transform algorithms, navigate to
``irlib/external`` and follow the directions in the ``README`` file, being sure
to move the created file ``pywavelet.so`` to some place from which it can be
imported.

Installing manually
~~~~~~~~~~~~~~~~~~~

Alternatively, *irlib* can be build in place without ``pip`` by doing

::

    >> python setup.py build_ext --inplace

Path
~~~~

For convenience, programs that make up *radar\_tools* should be on the execution
``PATH``. If ``pip`` was used, this should be taken care of. Otherwise, on Linux
and Mac OS X, one can add the following line to the ``.bashrc``:

::

    export PATH=$PATH:~/python/irlib

On Windows, one should be able to modify the *Path* variable by right clicking
on **My Computer** and going to *Properties -> Advanced System Settings ->
Environment Variables*.


