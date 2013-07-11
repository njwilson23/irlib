
import filters

def get_commands(module):
    """ Return the commands in module """
    commands = [a for a in module.__dict__.values()
                                    if (type(a) is type)
                                        and a.__base__ is filters.Command]
    commands = {}
    for key in module.__dict__:
        val = module.__dict__[key]
        if (type(val) is type) and val.__base__ is filters.Command:
            commands[val.cmd] = val
    return commands

def apply_filter(inputs, G):
    """ Attempt to find a Command matching an input_string, and apply. Raises
    CommandSearchError otherwise. """
    if len(inputs) > 1:
        cmd = inputs[0]
        args = inputs[1:]
    else:
        cmd = inputs[0]
        args = []

    command_dict = get_commands(filters)
    if cmd in command_dict:
        try:
            cmdobj = command_dict[cmd]()
            cmdobj.apply(G, args)
        except Exception as e:
            raise CommandApplicationError(e)
    else:
         raise CommandSearchError("No command definition '{0}' found".format(cmd))

def help_filter(cmd):
    """ Print the help documentation for *cmd*, or raise CommandSearchError if
    it cannot be found. """
    command_dict = get_commands(filters)
    if cmd in command_dict:
        return command_dict[cmd].helpstr
    else:
         raise CommandSearchError("No command definition '{0}' found".format(cmd))

def list_filters():
    for cmd in sorted(get_commands(filters).keys()):
        print cmd
    return

class CommandSearchError(Exception):
    def __init__(self, message="No message"):
        self.message = message
    def __str__(self):
        return self.message

class CommandApplicationError(Exception):
    def __init__(self, exception):
        self.e = exception
    def __str__(self):
        return e.__str__

