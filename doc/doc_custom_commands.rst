Adding customized filters to the irlib GUI apps
===============================================

*icepick* dynamically loads commands and filters when it starts. This
happens in the lines in the ``icepick2.py`` file that look like:

.. code:: python

    console.register(irlib.app.filters)
    console.register(irlib.app.pickcommands)
    console.register(irlib.app.mapcommands)

Adding custom filters has two parts:

1. writing new ``Command`` classes and placing them in a Python module
2. registering the module

Writing a custom Command module
-------------------------------

    Terminology:

    In traditional object-oriented programming, we talk about
    **classes** and **instances**. An instance is a bundle of data,
    while a class is a category of instances. Classes can inherit from
    other classes, receiving some of their attributes. As an example, in
    real life an ``Apple`` class might be a subclass of a ``Fruit``
    class, and the apple sitting on my desk is an object, or a specific
    instance of an ``Apple``.

    In irlib, commands are represented by classes, and you can think of
    an instance being created and used whenever you run one.

Commands are defined by creating a Python class definition that inherits
from ``irlib.app.commands.Command``. If our command is a filter, it's
better to inherit from ``irlib.app.filters.FilterCommandBase``, which is
a subclass of ``Command``.

The command class has two important attributes and one important method:
- The ``cmd`` attribute is a string that defines the command that runs
the filter - The ``helpstr`` attribute is a string containing a
description of how the command is used and what it does. - The
``apply(G, args)`` method is a function that takes a radar
``LineGather`` instance and a list of ``args`` and does something to it.
This is the part that makes the command do something.

For example, the linear gain control filter (invoked with ``gc``) is
defined as (with annotations added):

.. code:: python

    class LinearGainControl(FilterCommandBase):
        # these are the two methods that make the command usable and give it
        # help documentation
        cmd = "gc"
        helpstr = """Linear gain control

        gc [n]

        Apply a time-dependent linear gain control (t^n) to each trace, with the
        exponent `n` taking the default value of 1.0. """

        # This is the method that performs an action
        #
        # The first argument (called `self` out of tradition) is a reference to the
        # calling instance, and can be ignored.
        # G is a LineGather instance, which is the irlib object that holds the radar
        # data from a single line
        # `args` is a list. If we typed "gc 1 2 3" in icepick, args would be [1, 2, 3]
        def apply(self, G, args):
            # args may be an empty list, or it may contain a number used to set the
            # exponent on the gain filter.
            if len(args) > 0:
                npow = float(args[0])
            else:
                npow = 1.0
            # This filter simply calls the LineGather method that implements gain
            # control. We could do anything here, however.
            G.DoTimeGainControl(npow=npow)
            return

Here is an example of a filter that, for illustrative purposes, reverses
the radar wave polarity and sets a maximum voltage.

.. code:: python

    # saved as myfilters.py

    import numpy as np
    from irlib.app.filters import FilterCommandBase

    class ReversePolarityAndCap(FilterCommandBase):
        cmd = "toy_filter"
        helpstr = """Toy filter illustrating how to build commands

        toy_filter arg1 arg2
        """

        def apply(self, G, args):
            for arg in args:
                # Print each argument
                print("argument: " + str(arg))

            # Reverse wave polarity
            G.data = -G.data

            # Clip data to a maximum and a minimum value
            G.data = np.clip(G.data, -1.0, 1.0)
            return

You can place as many custom commands as you like in ``myfilters.py``.

Registering the commands
------------------------

At the moment, ``icepick2.py`` needs to be modified directly to add new
filters. So to use the command in ``myfilters.py``, we would first
import it at the top of ``icepick2.py``:

.. code:: python

    # [other import statements]
    import myfilters
    ...

and then we would register the module toward the bottom, but before the
*icepick* main loop starts:

.. code:: python

    ...
    console.register(myfilters)
    console.start()

Testing it out
--------------

Opening ``icepick``, we can see our new command:

::

    >> help

    ...

        Available Filter commands

        dewow
        lowpass_td
        ringing
        lowpass
        agc
        migfk
        reverse
        power
        highpass
        toy_filter          <-- we did this!
        highpass_td
        gc

    ...

::

    >> help toy_filter

    Toy filter illustrating how to build commands

        toy_filter arg1 arg2

::

    >> f toy_filter hello world!
    argument: hello
    argument: world!

...and the radargram data gets flipped! (Totally useful.)
