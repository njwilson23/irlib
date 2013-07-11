
import filters

def get_commands(module):
    """ Return the commands in module """
    commands = [a.cmd for a in module.__dict__ if super(a) is filters.Command]
    return commands

def apply_filter(input_string, G):
    """ Attempt to find a Command matching an input_string, and apply. Raises
    CommandSearchError otherwise. """
    if " " in input_string:
        cmd, args = input_string.split(None, 1)
    else:
        cmd = input_string
        args = None

    if cmd in get_commands(filters):
        try:
            filters.__dict__[cmd].apply(G, args)
        except Exception as e:
            raise CommandApplicationError(e.message)
    else:
         raise CommandSearchError("No command definition '{0}' found".format(cmd))

def help_filter(cmd):
    """ Print the help documentation for *cmd*, or raise CommandSearchError if
    it cannot be found. """
    if cmd in get_commands(filters):
        print filters.__dict__[cmd].helpstr
    else:
         raise CommandSearchError("No command definition '{0}' found".format(cmd))
    return

def list_filters():
    for cmd in sorted(get_commands(filters)):
        print cmd
    return

class CommandSearchError(Exception):
    def __init__(self, message="No message"):
        self.message = message
    def __str__(self):
        return self.message

class CommandApplicationError(Exception):
    def __init__(self, message="No message"):
        self.message = message
    def __str__(self):
        return self.message

