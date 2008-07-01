#! /usr/bin/env python

import sys
sys.path.append("../../bin")
import Gnuplot
import cPickle as pickle
from optparse import OptionParser
from error_tags import *
import os.path
import table




errors = pickle.load(file(sys.argv[1]))
errors = errors.filter(ETAG_COMPONENT, "sum")
errors = errors.filter(ETAG_TYPE, "L2")
resolutions = errors.tag_values(ETAG_RESOLUTION)
resolutions.sort()
degrees = errors.tag_values(ETAG_PDEGREE)
degrees.sort()
timesteps = errors.tag_values(ETAG_DT)
timesteps.sort()
timesteps.reverse()

def unique_value(table):
    if len(table) == 0:
        return "-"
    steps = table.tag_values(ETAG_STEP)
    last_step = max(steps)
    list = table.filter(ETAG_STEP, last_step).values()
    if len(list) > 1:
        print table.DataTable
        assert False
    if list[0] < 100:
        return "%6.3e" % list[0]
    else:
        return "-"

errors.add_dependent_tag("ppw", 
        lambda tags: float(tags[ETAG_RESOLUTION][0])*tags[ETAG_PDEGREE]/tags["usrres1"])
ppws = errors.tag_values("ppw")
ppws.sort()
lwaves = errors.tag_values("usrres1")
lwaves.sort()

t = table.Table()
t.append_row(["PPW", "dt"] + degrees)

for ppw in ppws:
    for dt in timesteps:
        t.append_row([ppw,dt]+[
                unique_value(errors
                    .filter("ppw", ppw)
                    .filter(ETAG_DT, dt)
                    .filter(ETAG_PDEGREE, deg)
                    )
                for deg in degrees])
print t
