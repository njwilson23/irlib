#! /usr/bin/env python
""" Benchmark program for irlib app framework """

import sys

import irlib.app.filters
import irlib.app.pickcommands
import irlib.app.mapcommands
from irlib.app.console import Console
from irlib.app.components import Radargram

bannertext = """IcePick v0.5
Type 'help' for assistance with available commands."""

try:
    console = Console("icepick2", bannertext=bannertext)
except IOError:
    sys.exit(1)
console.register(irlib.app.filters)
console.register(irlib.app.pickcommands)
console.register(irlib.app.mapcommands)
console.start()

