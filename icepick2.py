#! /usr/bin/env python
""" Benchmark program for irlib app framework """

from irlib.app.console import Console
from irlib.app.components import Radargram

bannertext = """IcePick v0.4-dev
Type 'help' for assistance with available commands."""

console = Console("icepick2", bannertext=bannertext)

#from irlib.app.components import PickWindow
#pw = PickWindow(console.line)
#pw.connect_radargram(console.get_appwindows(Radargram)[0])
#console.appwindows.append(pw)

