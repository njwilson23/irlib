
from __future__ import print_function
import traceback

def apply_command(registry, inputs, stateobj, cmdtype):
    """ Attempt to apply a command, raising a KeyError if a suitable command cannot be found.

    Parameters
    ----------
    registry:   dictionary of commands and their associated Command classes
    inputs:     list containing first the command string, and then any arguments
    stateobj:   the stateful object that may be modified by the command
                e.g. a Console or Gather instance
    """
    if len(inputs) > 1:
        cmd = inputs[0]
        args = inputs[1:]
    else:
        cmd = inputs[0]
        args = []

    if (cmd in registry) and (cmdtype in registry[cmd]):

        cmdclass = registry[cmd][cmdtype]()

        try:
            cmdclass.apply(stateobj, args)
        except Exception as e:
            traceback.print_exc()
            print(CommandApplicationError(e))

    else:
         raise KeyError("No command definition '{0}' found".format(cmd))


def help_command(registry, cmd):
    """ Print the help documentation for *cmd*, or raise Key if it cannot be
    found. """
    if cmd in registry:
        return registry[cmd].helpstr
    else:
         raise KeyError("No command definition '{0}' found".format(cmd))

def list_filters():
    return sorted(get_commands(filters).keys())

class CommandApplicationError(Exception):
    def __init__(self, exception):
        self.e = exception
    def __str__(self):
        return self.e.__str__()

