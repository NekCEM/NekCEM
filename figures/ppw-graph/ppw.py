#! /usr/bin/env python
import Gnuplot
import cPickle as pickle
from optparse import OptionParser
from data_tags import *
import sys
import os.path




usage = "usage: %prog PICKLE-FILE ELEMENT-COUNT"
parser = OptionParser(usage=usage)
parser.add_option(
        "-p", "--ps", dest="postscript", action="store_true",
        help="Output to .ps instead of to screen.")
parser.add_option(
        "-z", "--zwaves", dest="zwaves", action="store_true",
        help="Draw # of waves on X axis instead of PPW.")
parser.add_option(
        "-l", "--loglog", dest="loglog", action="store_true",
        help="Use a log-log-scale for the graph")
options, args = parser.parse_args()

if len(args) != 2:
    print "wrong number of arguments"
    print "use -h for help"
    sys.exit(1)

in_filename = args[0]
out_filename = os.path.splitext(in_filename)[0] + ".ps"

gp = Gnuplot.Gnuplot()
if options.zwaves:
    gp.xlabel("Number of wavelengths in domain")
else:
    gp.xlabel("Points per Wavelength")
gp.ylabel("L2 Error")
if options.loglog:
    gp("set logscale xy")
else:
    gp("set logscale y")
if options.postscript:
    gp("set terminal postscript")
    gp("set output \"%s\"" % out_filename)

errors = pickle.load(file(args[0]))
errors = errors.filter(TAG_COMPONENT, "sum")
errors = errors.filter(TAG_TYPE, "L2")
errors = errors.filter(TAG_STEP, max(errors.tag_values(TAG_STEP)))

err = []
ppw = []

elcount = int(args[1])

datasets = []

degrees = errors.tag_values(TAG_PDEGREE)
degrees.sort()

for degree in degrees:
    deg_errors = errors.filter(TAG_PDEGREE, degree)
    zwaves_values = deg_errors.tag_values("usrres1")

    error_values = deg_errors.values_sorted_by("usrres1")
    #error_values = [err / zw for err, zw in zip(error_values, zwaves_values)]
    if options.zwaves:
        x_values = zwaves_values
    else:
        error_values.reverse() # increasing ppw is decreasing zwaves
        points = float(elcount * (degree + 1))
        x_values = [points / zwaves for zwaves in zwaves_values]

    x_values.sort()

    datasets.append(Gnuplot.Data(x_values, error_values,
        title="Degree %d" % degree,
        with="linespoints"))

gp.plot(*datasets)    
if not options.postscript:
    raw_input()
