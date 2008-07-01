#! /usr/bin/env python

import sys
sys.path.append("../../../bin")
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

if False:
    t = table.Table()
    t.append_row(["#El", "dt"] + degrees)
    for res in resolutions:
        for dt in timesteps:
            t.append_row([res[0],dt]+[
                    unique_value(errors
                    .filter(ETAG_RESOLUTION, res)
                    .filter(ETAG_DT, dt)
                    .filter(ETAG_PDEGREE, deg)
                    )
                    for deg in degrees])
    print t
else:
    errors.add_dependent_tag("pointcount", 
            lambda tags: tags[ETAG_RESOLUTION][0]*tags[ETAG_PDEGREE])
    pointcounts = errors.tag_values("pointcount")
    pointcounts.sort()
    lwaves = errors.tag_values("usrres1")
    lwaves.sort()

    for lwave in lwaves:
        t = table.Table()
        t.append_row(["#Pts", "dt"] + degrees)

        print
        print "Wavelengths in Domain: %d" % lwave
        print

        for pointcount in pointcounts:
            for dt in timesteps:
                t.append_row([int(pointcount),dt]+[
                        unique_value(errors
                            .filter("usrres1", lwave)
                            .filter("pointcount", pointcount)
                            .filter(ETAG_DT, dt)
                            .filter(ETAG_PDEGREE, deg)
                            )
                        for deg in degrees])
        print t
        print "#Pts = #Elements * Polynomial degree"
