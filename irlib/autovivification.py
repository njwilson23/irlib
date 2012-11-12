class AutoVivification(dict):
    """ Implementation of perl's autovivification feature. Copied from
    http://stackoverflow.com/questions/635483/...
    what-is-the-best-way-to-implement-nested-dictionaries-in-python
    """
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
