""" Functions for reading pulseEKKO data files. """

import os
import numpy as np
from struct import unpack

def parse_header(lines):
    """ Read a header string and return a dictionary.

    str -> dict
    """
    meta = {}
    for line in lines:
        if "=" in line:
            k,v = line.split("=", 1)
            meta[k.strip()] = v.strip()
        elif (("-" in line) or ("/" in line)) and len(line.strip()) == 8:
            meta["date"] = line.strip()

    return meta

def parse_data(s):
    """ Read a data string and return a dictionary and a data array.

    str -> (dict, array)
    """
    i = 0
    dlist = []
    meta = {}

    while True:
        if len(s) < i+128:
            break
        hdr = unpack("32f", s[i:i+128])
        nsmp = int(hdr[2])
        d = unpack("{0}h".format(nsmp), s[i+128:i+128+2*nsmp])
        meta[i] = hdr
        dlist.append(d)
        i += (128 + 2*nsmp)

    # Pad short traces with zeros
    maxlen = max([len(a) for a in dlist])
    dlist_even = map(lambda a,n: a if len(a) == n else a+(n-len(a))*[0],
                     dlist, (maxlen for _ in dlist))
    darray = np.vstack(dlist_even).T
    return meta, darray

def read_pulseEKKO(path):
    """ Search for header and data files matching path, open them, and return a
    dictionary of line metadata, a dictionary of trace metadata, and an array
    of radar data.

    str -> (dict, dict, array)
    """
    directory, nm = os.path.split(path)
    if nm + ".HD" not in os.listdir(directory):
        raise IOError("{0}.HD not found".format(nm))
    elif nm + ".DT1" not in os.listdir(directory):
        raise IOError("{0}.DT1 not found".format(nm))

    with open(path + ".HD", "r") as f:
        lnmeta = parse_header(f.readlines())

    with open(path + ".DT1", "rb") as f:
        trmeta, darray = parse_data(f.read())

    return lnmeta, trmeta, darray

