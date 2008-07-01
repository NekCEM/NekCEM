#! /usr/bin/env python

from __future__ import division

import sys
sys.path.append("../../bin")
import Gnuplot
import cPickle as pickle
import os
from data_tags import *

def product(iterable):
    import operator
    return reduce(operator.mul, iterable)

gp = Gnuplot.Gnuplot()
gp("set terminal postscript")

rundir = sys.argv[1]
outdir = sys.argv[2]

times = pickle.load(file(os.path.join(rundir, "times.pickle")))
steps = max(times.tag_values(TAG_STEP))
times = times.filter(TAG_STEP, steps)
errors = pickle.load(file(os.path.join(rundir, "errors.pickle")))\
        .filter(TAG_COMPONENT, "sum")\
        .filter(TAG_TYPE, "L2")\
        .filter(TAG_STEP, steps)

cases = errors.tag_values(TAG_REASTEM)
cases.sort()
degrees = errors.tag_values(TAG_PDEGREE)
degrees.sort()

def unique_value(table):
    if len(table) == 0:
        raise RuntimeError, "no value in table"
    if len(table) > 1:
        raise RuntimeError, "more than one value in table"
    return table.values()[0]




for case in cases:
    case_errors = errors.filter(TAG_REASTEM, case)

    resolutions = case_errors.tag_values(TAG_RESOLUTION)
    resolutions.sort()

    # grid-cpu graph
    graphname = "grid-cpu-%s.ps" % case
    gp("set output \"%s\"" % (os.path.join(outdir, graphname)))
    gp("set logscale xy")
    gp.xlabel("Number of grid points")
    gp.ylabel("CPU time [s]")
    gp.title("Test case %s, %d time steps" % (case, steps))
    datasets = []
    for deg in degrees:
        x_values = []
        y_values = []
        for res in resolutions:
            dim = len(res)
            gridp = (deg+1)**dim * product(res)
            time = unique_value(times\
                    .filter(TAG_REASTEM, case)
                    .filter(TAG_RESOLUTION, res)
                    .filter(TAG_PDEGREE, deg))
            x_values.append(gridp)
            y_values.append(time)
        datasets.append(Gnuplot.Data(x_values, y_values,
            title=("N=%d" % deg), with="linespoints"))
    gp.plot(*datasets)

    # grid-error graph
    graphname = "grid-error-%s.ps" % case
    gp("set output \"%s\"" % (os.path.join(outdir, graphname)))
    gp("unset logscale xy")
    gp("set logscale y")
    gp.xlabel("Number of grid points")
    gp.ylabel("L2 error")
    gp.title("Test case %s, %d time steps" % (case, steps))
    datasets = []
    for deg in degrees:
        x_values = []
        y_values = []
        for res in resolutions:
            dim = len(res)
            gridp = (deg+1)**dim * product(res)
            error = unique_value(errors\
                    .filter(TAG_REASTEM, case)
                    .filter(TAG_RESOLUTION, res)
                    .filter(TAG_PDEGREE, deg))
            x_values.append(gridp)
            y_values.append(error)
        datasets.append(Gnuplot.Data(x_values, y_values,
            title=("N=%d" % deg), with="linespoints"))
    gp.plot(*datasets)

    # error-cpu graph
    graphname = "error-cpu-%s.ps" % case
    gp("set output \"%s\"" % (os.path.join(outdir, graphname)))
    gp("unset logscale xy")
    gp("set logscale x")
    gp.xlabel("L2 error")
    gp.ylabel("CPU time [s]")
    gp.title("Test case %s, %d time steps" % (case, steps))
    datasets = []
    for deg in degrees:
        x_values = []
        y_values = []
        for res in resolutions:
            dim = len(res)
            error = unique_value(errors\
                    .filter(TAG_REASTEM, case)
                    .filter(TAG_RESOLUTION, res)
                    .filter(TAG_PDEGREE, deg))
            time = unique_value(times\
                    .filter(TAG_REASTEM, case)
                    .filter(TAG_RESOLUTION, res)
                    .filter(TAG_PDEGREE, deg))
            x_values.append(error)
            y_values.append(time)
        datasets.append(Gnuplot.Data(x_values, y_values,
            title=("N=%d" % deg), with="linespoints"))
    gp.plot(*datasets)

    # degree-cpupgpt graph
    graphname = "degree-cpupgpt-%s.ps" % case
    gp("set output \"%s\"" % (os.path.join(outdir, graphname)))
    gp("unset logscale xy")
    gp("set logscale y")
    gp.xlabel("L2 error")
    gp.ylabel("CPU time per gridpoint per timestep [s]")
    gp.title("Test case %s, %d time steps" % (case, steps))
    datasets = []
    for res in resolutions:
        x_values = []
        y_values = []
        for deg in degrees:
            dim = len(res)
            gridp = (deg+1)**dim * product(res)
            timepgpt = unique_value(times\
                    .filter(TAG_REASTEM, case)
                    .filter(TAG_RESOLUTION, res)
                    .filter(TAG_PDEGREE, deg))/gridp/steps
            x_values.append(deg)
            y_values.append(timepgpt)
        res_str = "x".join([str(x) for x in res])
        datasets.append(Gnuplot.Data(x_values, y_values,
            title=res_str, with="linespoints"))
    gp.plot(*datasets)
