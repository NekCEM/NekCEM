# Documentation can be found in doc/sisyphus.html. See also Section
# IV on how to generate that HTML from this code.

"""
Sisyphus is an automatic test runner. It takes a single argument,
which is a Python file (see http://docs.python.org/tut for a tutorial)
which is run in a way that a few special functions are available.

After running the tests specified by the script, sisyphus writes
a report in HTML format to a given output directory.

I. How to write a test script
-----------------------------

As a matter of convention, the extension of test scripts is `.sis'.

Currently, three kinds of tests are especially easy to run from
test scripts, and should cover most of your testing needs:

    - test_error()
    - test_p_convergence()
    - test_eigenvalues()

These functions (methods, actually) are part of a ScriptManager
object. Exact documentation for their possible parameters can be
found in the documentation for ScriptManager. 

We will explore them one by one. Here's a really elementary test
script:

    test_error("examples/recwave", "recwave", "recwave4", degree=8,
      max_error=1e-8)

The first parameter gives the directory where code is supposed to be
run, the second one is the stem (i.e. filename w/o extension) of the
.usr file, the third one is the stem of the .rea file.

This call checks that the maximum error returned by the recwave 
examples is less than 1e-8.

Observe how the order of the parameters does not matter as long as you
specify the "name=" part of them.

Here's a slightly more involved example:

    flux_types = [
      {"ifcentral":True, "ifupwind":False},
      {"ifcentral":False, "ifupwind":True}
      ]

    for ft in flux_types:
        test_p_convergence("examples/recwave", "recwave", "recwave4", 
            min_eoc=15, rea_switches=ft)
    
This example checks for both central and upwinding flux if the code
converges with an empirical order of convergence (EOC) of at least 15.
It does so by looping over a list of dictionaries that contain
settings that are passed as rea_switches. These switches are then
used in the runs, where they are automatically written to the .rea files.

As a final example, consider the following, once again simple script:

    test_eigenvalues("examples/2dbox3mix", "2dbox3mix", "2dbox3mix")

This checks the timestepping matrix for "unstable" eigenvalues with
large real parts. The cutoff value is the default here. Note that if
no rea_switches are given, then whatever is currently found in the
.rea file will be used.

All these functions allow many more settings and parameters. Please see the
documentation of the class ScriptManager for complete details.  Also see
examples/testscript.sis for a more elaborate test script.

Note that paths specified in test scripts are always relative
to the root of the nekton tree, not relative to where sisyphus is
run.

There is also a dictionary called PARAMETERS available in the test
script that results from the parsing of the `-s' option to sisyphus.
An option `-s nocentral,nodealias,steps=17' gets turned into the
dictionary

  {'nodealias': None, 'steps': '17', 'nocentral': None}

The meaning of these parameters is entirely up to the test script.
Past the comma-separated key-value format, sisyphus makes no assumption
and takes no notice of these parameters.

The following symbols are also available in script files,
if you need to do more advanced things:

- TestCase
- RetryingTestCase
- NekCompileTestCase
- NekRunTestCase, 
- NekTimeRunTestCase
- NekEigenvalueTestCase
- NekPConvergenceTestCase
- NekPMLTestCase

(If you are looking at the HTML documentation, the symbols above
should be links taking you to the appropriate documentation.)

Additionally, of course, the entire power of Python is at your
disposal to write the test script.

II. Running and output
----------------------

If you seem to have forgotten how to run sisyphus and just need a quick
reminder, run `sisyphus -h'. It will print a brief summary of the 
available options.

Generally, sisyphus is run with a single argument, the test script.
Let's try that:

   $ sisyphus examples/testscript.sis 

We might get:

    Traceback (most recent call last):
      File "bin/sisyphus", line 1304, in ?
        run()
      File "bin/sisyphus", line 1267, in run
        options.results_dir)
      File "bin/sisyphus", line 344, in __init__
        raise RuntimeError, "The results directory `%s' may not exist "
    RuntimeError: The results directory `test-results' may not exist prior 
    to running sisyphus.

This tells us that the default output directory, called `test-results'
exists, and sisyphus is refusing to overwrite it. By specifying the 
option `-p' or `--purge-old', we can force it to delete the existing
directory:

   $ sisyphus -p examples/testscript.sis 

We get:

    compile-recwave-lxi=9... OK
    timerun-recwave4-lxi=9-iocomm=5_dt=0.001_nsteps=100-ifupwind=T_ifcentral=F...
    ... and so on ...

That is, sisyphus is now busy running the tests, and upon completing
each test case, it will print a brief summary of the result. 

Instead of telling sisyphus to purge the previous results, we could also
have specified a different output directory, like so:

   $ sisyphus -o my-new-test-results examples/testscript.sis 

Note that `my-new-test-results' may also not exist, or sisyphus will
complain.

As long as you are running sisyphus inside a Nekton tree, it does not really
matter which directory you are in. It will automatically walk up the directory
tree until it finds the root of the Nekton tree and then interpret all paths in
the script file relative to that.

III. Restricting the set of tests that are run
----------------------------------------------

If you run sisyphus with the option `-n', sisyphus will not actually
execute tests or write results, instead it will only print the list
of tests that *would* be run to standard output:

   $ sisyphus -n -p examples/testscript.sis 

We get something like:

    Would run the following test cases:
    timerun-recwave4-lxi=3-iocomm=5_dt=0.001_nsteps=100-ifupwind=F_ifcentral=T
    pconvergence-2dbox3mix-iocomm=5_dt=0.001_nsteps=100-ifupwind=F_ifcentral=T
    timerun-recwave4-lxi=6-iocomm=5_dt=0.001_nsteps=100-ifupwind=F_ifcentral=T
    compile-2dbox3mix-lxi=4
    compile-2dbox3mix-lxi=5
    compile-2dbox3mix-lxi=6
    ...

The option `-r' lets you specify a shell filename pattern, and only such
test cases whose name matches the pattern will be run. (Additionally,
each test case may *depend* on a few other test cases, such as a 
`run' test case may depend on a `compile' test case. (See TestSupervisor
for details.) The test cases selected by the `-r' option are always
run along with the ones they depend on.

Let's try this:

    $ sisyphus -r 'timerun-2dbox*lxi=3*' -n -p examples/testscript.sis 

    Would run the following test cases:
    timerun-2dbox3mix-lxi=3-iocomm=5_dt=0.001_nsteps=100-ifupwind=T_ifcentral=F
    compile-2dbox3mix-lxi=3
    timerun-2dbox3mix-lxi=3-iocomm=5_dt=0.001_nsteps=100-ifupwind=F_ifcentral=T

You see that the `compile' is run even though it's not been specified.
Be careful to quote the pattern to prevent the shell from expanding
it, leading to potentially confusing results.

IV. If this document seems out of date
--------------------------------------

To update the HTML documentation that you might be looking at, go to
the bin/ subdirectory and run

    $ pydoc -w ./sisyphus.py

This will write a new HTML documentation file sisyphus.html in bin/, 
which you should then move to the doc/ directory.

(Yes, the "./" in front of the file name is definitely necessary.)
"""



import sys
import nektools
from data_tags import *



try:
    import Gnuplot
    HAVE_GNUPLOT = True
except ImportError:
    print "Sisyphus recommends that you have Gnuplot.py for graphing."
    print "You may install it from http://gnuplot-py.sourceforge.net/"
    HAVE_GNUPLOT = False

# generic tool functions -------------------------------------------------
parse_fortran_float = nektools.parse_fortran_float
MySet = nektools.MySet
TaggedDataTable = nektools.TaggedDataTable
my_enumerate = nektools.my_enumerate
sort_by = nektools.sort_by




def get_settings_str(settings):
    """Convert a mapping of settings to a string representation.
    """
    def mystr(val):
        if val == True:
                return "T"
        elif val == False:
                return "F"
        else:
            return str(val)

    return "_".join([
            "%s=%s" % (key, mystr(value))
            for key, value in settings.iteritems()])




def run_command(dir, cmd, env_settings={}):
    """Run a command and capture standard output and standard error.
    Return a tuple `ok, output', where `ok' is a boolean that tells
    whether the execution finished unsignaled and with exit status 0,
    and output is the combined stdout and stderr output of the program.

    `cmd' may either be a string, in which case it is fed to a shell,
    or a list of arguments, in which case it is run directly. `dir'
    specifies the directory in which the command is run. `env' is a
    dictionary that specifies changes to the subprocess's environment.

    Current directory and environment are restored by this function.
    """
    import os
    import popen2

    curdir = os.getcwd()
    env_set = EnvironmentSetAndRestore(env_settings)

    try:
        os.chdir(dir)

        output = "RUNNING "
        if isinstance(cmd, str):
            output += "%s\n" % cmd
        else:
            output += "%s\n" % " ".join(cmd)
        output += "in directory: %s\n" % os.getcwd()

        child_process = popen2.Popen4(cmd)
        output += child_process.fromchild.read()
        wait_status = child_process.wait()

        signaled = os.WIFSIGNALED(wait_status)
        status = os.WEXITSTATUS(wait_status)

        ok = (status == 0) and not signaled

        if status:
            output += "EXIT STATUS %d\n" % status
        if signaled:
            output += "KILLED BY SIGNAL %d\n" % os.WTERMSIG(wait_status)
            if os.WCOREDUMP(wait_status):
                output += "CORE DUMPED\n"
        return ok, output

    finally:
        # restore directory and environment no matter what
        os.chdir(curdir)
        env_set.restore()




if HAVE_GNUPLOT:
    class PNGPlotter(object):
        """A simple wrapper around Gnuplot.py that makes it easy
        to grab the final plot as a PNG file.

        After the PNG file has been obtained, further plotting 
        operations will result in an error.
        """
        def __init__(self):
            self.Filename = "gp-tmp.png"
            self.GP = Gnuplot.Gnuplot()
            self.GP("set terminal png medium")
            self.GP("set output \"%s\"" % self.Filename)

        def title(self, *args):
            self.GP.title(*args)
        def xlabel(self, *args):
            self.GP.xlabel(*args)
        def ylabel(self, *args):
            self.GP.ylabel(*args)
        def plot(self, *args):
            self.GP.plot(*args)
        def __call__(self, *args):
            self.GP(*args)

        def get_png(self):
            del self.GP

            import os
            data = file(self.Filename, "rb").read()
            os.unlink(self.Filename)
            return data
else:
    class PNGPlotter(object):
        """A dummy plotter. We use it if Gnuplot.py is not available.
        """
        def __init__(self):
            pass

        def title(self, *args):
            pass

        def xlabel(self, *args):
            pass
        def ylabel(self, *args):
            pass
        def plot(self, *args):
            pass
        def __call__(self, *args):
            pass

        def get_png(self):
            # This is a 5x5 PNG image of a black X
            return '\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x05'
            '\x00\x00\x00\x05\x08\x00\x00\x00\x00\xa8\x04y9\x00\x00\x00'
            '\x1cIDAT\x08\xd7-\xc6\xb1\x01\x00\x00\x08\xc3 \xfa\xff\xd1q'
            '\x91\x89\x89IZ\x96\xbf\xd8\x01\x98\xe4\t\xfd\x16&\x1e\xa7\x00'
            '\x00\x00\x00IEND\xaeB`\x82'




def all_same(iterable, f=lambda x: x):
    """Return True iff all elements of an iterable are the same. 
    If `f' is specified, return True iff f(x) for all elements
    x is the same.
    """
    iterator = iterable.__iter__()
    try:
        first = iterator.next()
    except StopIteration:
        return True

    for i in iterator:
        if f(i) != f(first):
            return False
    return True




def vandermonde(vector, degree=None):
    """Return a Vandermonde matrix, i.e.

    1 v1 v1^2 ...
    1 v2 v2^2 ...
    ...
    """
    import Numeric as num

    if degree is None:
        degree = len(vector)

    mat = num.zeros((len(vector), degree+1), num.Float64)
    for i, v in my_enumerate(vector):
        for power in range(degree+1):
            mat[i,power] = v**power
    return mat




def fit_polynomial(x_vector, data_vector, degree):
    """Fit a polynomial of degree `degree' at the points 
    in `x_vector' to the `data_vector' by least-squares.

    Return a vector of polynomial coefficients in 
    order of increasing monomial degree.
    """

    import Numeric as num
    import LinearAlgebra as la

    vdm = vandermonde(x_vector, degree)
    vdm_h = num.transpose(num.conjugate(vdm))
    vdm2 = num.matrixmultiply(vdm_h, vdm)
    rhs = num.matrixmultiply(vdm_h, data_vector)
    return la.solve_linear_equations(vdm2, rhs)




def estimate_order_of_convergence(abscissae, errors):
    """Assuming that abscissae and errors are connected by a law of the form

    error = constant * abscissa ^ (-order),

    this function finds, in a least-squares sense, the best approximation of
    constant and order for the given data set. It returns a tuple (constant, order).
    Both inputs must be sequences.
    """

    import Numeric as num

    assert len(abscissae) == len(errors)
    if len(abscissae) <= 1:
        raise RuntimeError, "Need more than one value to guess order of convergence."

    coefficients = fit_polynomial(num.log10(abscissae), 
            num.log10(errors), 1)
    return 10**coefficients[0], -coefficients[1]
  




class Reference(object):
    """A simple value container that makes an rvalue modifiable.
    """
    def __init__(self, value):
        self.Value = value

    def get(self):
        return self.Value

    def set(self, value):
        self.Value = value




class Stopwatch(object):
    """A generic high-precision timekeeping object.
    May be stopped and restarted multiple times.
    Records wall time, subprocesses' userland and kernelspace
    CPU times. 
    """
    def __init__(self):
        self.reset()
        self.Running = False

    def copy(self):
        result = Stopwatch()
        result.Running = self.Running
        result.WallTime = self.WallTime
        result.SubProcUserTime = self.SubProcUserTime
        result.SubProcSystemTime = self.SubProcSystemTime
        if self.Running:
            result.WallStart = self.WallStart
            result.SubProcUserStart = self.SubProcUserStart
            result.SubProcSystemStart = self.SubProcSystemStart
        return result

    def reset(self):
        self.WallTime = 0
        self.SubProcUserTime = 0
        self.SubProcSystemTime = 0

    def start(self):
        import os
        dummy, dummy, self.SubProcUserStart, self.SubProcSystemStart, \
                self.WallStart = os.times()
        self.Running = True

    def stop(self):
        import os
        dummy, dummy, user_end, system_end, wall_end = os.times()

        self.SubProcUserTime += user_end - self.SubProcUserStart
        self.SubProcSystemTime += system_end - self.SubProcSystemStart
        self.WallTime += wall_end - self.WallStart

        del self.SubProcUserStart
        del self.SubProcSystemStart
        del self.WallStart

        self.Running = False

    def __str__(self):
        if self.Running:
            copy = self.copy()
            copy.stop()
            return str(copy)

        if self.SubProcUserTime or self.SubProcSystemTime:
            return "%.1fs wall %.1fs user" % \
                    (self.WallTime, self.SubProcUserTime)
        else:
            return "%.1fs wall" % self.WallTime

    def detailed_str(self):
        if self.Running:
            copy = self.copy()
            copy.stop()
            return copy.detailed_str()

        if self.SubProcUserTime or self.SubProcSystemTime:
            return "%.1fs wall %.1fs user %.1fs system %.1fs cpu" % \
                    (self.WallTime, self.SubProcUserTime, 
                            self.SubProcSystemTime, 
                            self.SubProcUserTime+self.SubProcSystemTime)
        else:
            return "%.1fs wall" % self.WallTime

    def _running(self):
        return self.Running
    running = property(_running, 
            doc="Whether the stopwatch is currently running.")

    def _wall_time(self):
        return self.WallTime
    wall_time = property(_wall_time, 
            doc="How much real time has passed between start and stop.")

    def _subproc_user_time(self):
        return self.SubProcUserTime
    subproc_user_time = property(_subproc_user_time, 
            doc="How much userland CPU time has passed in subprocesses "
            "between start and stop.") 
    def _subproc_system_time(self):
        return self.SubProcSystemTime
    subproc_system_time = property(_subproc_system_time, 
            doc="How much system CPU time has passed in subprocesses "
            "between start and stop.")




def histogram(iterable, f):
    """Return a dictionary that indicates the count of how
    often each value of `f' occurred when letting the argument
    of `f' run over the `iterable'.

    (`f' is a function of one argument.)
    """
    result = {}
    for i in iterable:
        f_i = f(i)
        result[f_i] = result.setdefault(f_i, 0) + 1
    return result




def hist2str(hist):
    """Return a human-readable representation of the histogram()-generated
    `hist'.
    """
    return ", ".join(["%d %s" % (value, key) 
            for key, value in hist.iteritems()])




class EnvironmentSetAndRestore:
    """A class that sets environment variables passed to it,
    and resets them when restore() is called.
    """

    def __init__(self, env_settings):
        import os
        self.PreviousSettings = {}
        for evar, value in env_settings.iteritems():
            if evar in os.environ:
                self.PreviousSettings[evar] = os.environ[evar]
            else:
                self.PreviousSettings[evar] = None

            os.environ[evar] = value

    def restore(self):
        import os
        for evar, value in self.PreviousSettings.iteritems():
            if value is not None:
                os.environ[evar] = value
            else:
                del os.environ[evar]




# abstract testing framework ---------------------------------------------
class TestSupervisor(object):
    """This class keeps a list of TestCase objects, manages their
    dependencies, runs them in the right order, stores and passes 
    results from one TestCase to the next, and outputs the
    HTML summary.

    If you imagine the dependencies as a directed acyclic graph,
    then tests are run in a depth-first order. This is important.
    Consider this case:

       Compile A  Compile B
          |           |
          |           |
        Run A      Run B

    The order (Compile A, Run A, Compile B, Run B) is valid,
    while (Comiple A, Compile B, Run A, Run B) is not. This is
    because the compilation puts the tree in a certain state
    on which "Run A" will depend. If we run "Compile B" in between,
    then this state is no longer present--"Compile B" might be a
    compilation in the same directory, just with different settings.
    """

    def __init__(self, config, results_dir):
        """`config' is a dictionary that carries configuration
        information from the caller to the test cases. Currently,
        there is only one useful key that is evaluated by 
        the test cases, and that is "tree_root", pointing to the
        root of the Nekton tree.

        `results_dir' is the result where the HTML summary is written.
        It may not exist.
        """

        import time

        self.Config = config
        self.TestCases = []
        self.ResultsTable = {}
        self.ResultsDir = results_dir
        self.RunStopwatch = Stopwatch()
        self.RunStart = time.time()
        self.FilesWrittenFor = MySet()
        self.LogsWrittenFor = MySet()
        self.Interrupted = False
        self.LastReportWritten = 0

        import os
        try:
            os.mkdir(self.ResultsDir)
        except OSError:
            raise RuntimeError, "The results directory `%s' may not exist "\
                    "prior to running sisyphus." % self.ResultsDir

    def add_test_case(self, tc):
        self.TestCases.append(tc)

    def roots(self):
        """Return a list of TestCase objects with no dependencies.
        """
        return [tc for tc in self.TestCases if not tc.depends_on()]

    def test_cases(self):
        """Return a list of all test cases.
        """
        return self.TestCases

    def children_map(self):
        """Return a mapping from a test case to the ones that
        depend on it.
        """
        result = {}
        for tc in self.TestCases:
            deps = tc.depends_on()
            for dep in deps:
                result.setdefault(dep, []).append(tc)
        return result

    def restrict_to_subset(self, subset):
        """Beginning with `subset', find all TestCases that
        `subset' directly or indirectly depends on. Then
        restrict the set of tests to run to this.
        """
        deps_of_subset = MySet()
        for tc in subset:
            deps_of_subset.update(tc.recursively_depends_on())
        self.TestCases = list(deps_of_subset | MySet(subset))

    def run(self):
        """Run tests in depth-first order. See TestSupervisor.
        """
        kids_table = self.children_map()

        self.Interrupted = False

        def run_with_kids(root):
            import traceback

            anc_data = []
            anc_data_complete = True

            if self.Interrupted: return

            for ancestor in root.depends_on():
                if ancestor not in self.ResultsTable:
                    # Need to run more ancestors first,
                    # don't run this branch just yet.
                    # We will get to this one via
                    # the last ancestor.
                    return
                else:
                    anc_res = self.ResultsTable[ancestor]
                    if anc_res.State != TestResult.FAILED:
                        anc_data.append((ancestor, anc_res.Data))
                    else:
                        result = TestResult(
                                TestResult.FAILED, 
                                oneliner="prerequisite %s failed" % 
                                ancestor.unique_name())
                        anc_data_complete = False

            sys.stdout.write("%s... " % root.unique_name())
            sys.stdout.flush()

            if anc_data_complete or root.accepts_partial_results():
                # we're actually executing something, so also write report saying so.
                try:
                    self.write_intermediate_results(root)
                except KeyboardInterrupt:
                    self.Interrupted = True
                    return

                sw = Stopwatch()
                sw.start()

                try:
                    result = root.run(anc_data, self.Config)
                except KeyboardInterrupt:
                    result = TestResult(TestResult.NOT_DONE, 
                            oneliner="Interrupted")
                    self.Interrupted = True
                except:
                    extype, value, tb = sys.exc_info()
                    tb_str = "".join(traceback.format_exception(extype, value, tb))
                    result = TestResult(TestResult.FAILED, 
                            oneliner="Failed with %s" % str(extype),
                            logs=[tb_str])

                sw.stop()
                result.Stopwatch = sw

                if not anc_data_complete and result.State == TestResult.GOOD:
                    result.State = TestResult.WARNING
                    result.OneLiner = "%s [incomplete prerequisites]" % \
                            result.OneLiner

            sys.stdout.write("%s\n" % result.OneLiner)

            self.ResultsTable[root] = result

            if self.Interrupted: return

            if root in kids_table:
                for kid in kids_table[root]:
                    run_with_kids(kid)

        self.RunStopwatch.start()
        for root in self.roots():
            run_with_kids(root)
        self.RunStopwatch.stop()
        self.write_results()

    def get_result(self, tc):
        try:
            return self.ResultsTable[tc]
        except KeyError:
            return TestResult(TestResult.NOT_DONE, "Not yet run")

    def summarize_test_states(self, tc_list):
        def tc2state_str(tc):
            return TestResult.state2str(self.get_result(tc).State)
        return hist2str(histogram(tc_list, tc2state_str))

    def write_intermediate_results(self, current_tc=None):
        from time import time
        if abs(time() - self.LastReportWritten) > 60:
            self.LastReportWritten = time()
            self.write_results(current_tc)

    def write_results(self, current_tc=None):
        import os

        def test_dir_rel(tc):
            return tc.unique_name()

        def test_dir(tc):
            import os.path
            dirname = test_dir_rel(tc)
            absname = os.path.join(results_dir, dirname)
            try:
                os.mkdir(absname)
            except OSError:
                pass
            return absname

        def test_link(tc, link_text=None):
            if link_text is None:
                link_text = str(test_indices[tc])
            return "<a href=\"#sumtab-entry-%s\" title=\"%s\">%s</a>" % \
                    (tc.unique_name(), tc.unique_name(), link_text)

        def log_link(tc):
            import os.path

            tres = self.get_result(tc)
            if not tres.Logs or tres.State == TestResult.NOT_DONE: 
                return ""

            if not tc in self.LogsWrittenFor:
                file(os.path.join(test_dir(tc), "log"), "w")\
                        .write((75 * "-" + "\n").join(tres.Logs))
                self.LogsWrittenFor.add(tc)

            return "<a href=\"%s/log\">[log]</a>" % test_dir_rel(tc)

        def category_list(cats):
            result = "<h2>Categories</h2>\n"
            result += "<ul>\n"
            result += "\n".join([
                    "<li><a href=\"#sumtab-%s\"><tt>%s</tt></a></li>\n" % (cat,cat) 
                    for cat in cats
                    ])

            result += "</ul>\n"
            return result

        def summary_table(cat):
            def get_result_local(tc):
                try:
                    return self.ResultsTable[tc]
                except KeyError:
                    if tc is current_tc:
                        return TestResult(TestResult.RUNNING, "Currently running")
                    else:
                        return TestResult(TestResult.NOT_DONE, "Not run yet")

            result = "<h3><a name=\"sumtab-%s\"/>" \
                    "Summary table <tt>%s</tt></h2>\n" % (cat, cat)
            result += """
              <table class="summary_table"><tr>
              <th>#</th>
              <th>Test ID</th>
              <th>Log</th>
              <th>Detail</th>
              <th>Summary</th>
              <th>Times</th>
              <th>Dependencies</th>
              </tr>
              """

            tests = [tc for tc in self.TestCases if tc.category() == cat]
            sort_by(tests, lambda tc: (get_result_local(tc).State, 
                tc.unique_name()))

            for tc in tests:
                tres = get_result_local(tc)
                state = TestResult.state2str(tres.State)

                if tres.Stopwatch is not None:
                    times = str(tres.Stopwatch)
                else:
                    times = ""

                result += """
                  <tr class="%s_test">
                  <td><a name="sumtab-entry-%s"/>%d</td>
                  <td><tt>%s</tt></td>
                  <td>%s</td>
                  <td>%s</td>
                  <td>%s</td>
                  <td>%s</td>
                  <td>%s</td>
                  </tr>
                  """ % (state, tc.unique_name(),
                          test_indices[tc], tc.unique_name(), 
                          log_link(tc), detail_link(tc, "[detail]"),
                          tres.OneLiner, times,
                          ", ".join([test_link(dep) 
                              for dep in tc.depends_on()])
                          )

            result += "</table>"
            result += "%s.<br/>\n" % self.summarize_test_states(tests)
            return result

        def summary_tables(cats):
            result = "<h2>Summary tables</h2>\n"
            for cat in cats:
                result += summary_table(cat)
            return result

        def detail_link(tc, link_text=None):
            import os.path

            tres = self.get_result(tc)
            if not tres.HTML: 
                return ""

            if link_text is None:
                link_text = "<tt>%s</tt>" % tc.unique_name()

            return "<a href=\"#detail-%s\">%s</a>" % (tc.unique_name(), link_text)

        def detail():
            import re
            file_re = re.compile(r"FILE\{([a-z.A-Z0-9]+)\}")
            result = "<h2>Test Details</h2>\n"
            for tc in self.TestCases:
                tcres = self.get_result(tc)
                tc_dir = test_dir_rel(tc)
                if not tcres.HTML:
                    continue
                result += "<a name=\"detail-%s\"/>\n" % tc.unique_name()
                result += "<h3>Detail for <tt>%s</tt></h3>\n" % tc.unique_name()
                html = tcres.HTML
                while True:
                    match = file_re.search(html)
                    if match is None: break
                    html = html[:match.start()] + \
                            "%s/%s" % (tc_dir, match.group(1)) + \
                            html[match.end():]

                result += html

                dep_links = []
                for dep in tc.depends_on():
                    dep_link = detail_link(dep)
                    if dep_link:
                        dep_links.append(dep_link)

                if dep_links:
                    result += """
                    These results are based on:\n
                    <ul>
                    %s
                    </ul>
                    """ % "\n".join(["<li>%s</li>" % dep_link for dep_link in dep_links])

                result += "%s<br/>" % test_link(tc, "Back to summary")

            return result


        def write_files():
            import os.path
            for tc in self.TestCases:
                tcres = self.get_result(tc)

                if tcres.State == TestResult.NOT_DONE:
                    continue
                if tc in self.FilesWrittenFor:
                    continue

                tcdir = test_dir(tc)
                for name, contents in tcres.Files.iteritems():
                    file(os.path.join(tcdir, name), "wb")\
                            .write(contents)

                self.FilesWrittenFor.add(tc)

        def vital_stats():
            import socket
            import time

            def format_time(t):
                return time.strftime("%a, %d %b %Y %H:%M:%S", 
                        time.localtime(t))

            if self.RunStopwatch.running:
                status = "Still running"
            else:
                if self.Interrupted:
                    status = "Stopped by user"
                else:
                    status = "Finished"

            if current_tc is None:
                currently_running = "nothing"
            else:
                currently_running = test_link(current_tc)

            return """
            <h2>Vital statistics</h2>
            <table>
              <tr><td>Status</td><td>%s</td></tr>
              <tr><td>Host</td><td>%s</td></tr>
              <tr><td>Started at</td><td>%s</td></tr>
              <tr><td>Report generated at</td><td>%s</td></tr>
              <tr><td>Time consumed</td><td>%s</td></tr>
              <tr><td>Grand total</td><td>%s</td></tr>
              <tr><td>Currently running</td><td>%s</td></tr>
            </table>
            """ % (status, socket.gethostname(),
                    format_time(self.RunStart),
                    format_time(time.time()),
                    self.RunStopwatch.detailed_str(),
                    self.summarize_test_states(self.TestCases),
                    currently_running)

        test_indices = dict([(tc, idx) for idx, tc in my_enumerate(self.TestCases)])

        results_dir = self.ResultsDir
        html_file = os.path.join(results_dir, "index.html")

        categories = list(MySet([tc.category() for tc in self.TestCases]))
        categories.sort()

        body = \
                vital_stats() + \
                category_list(categories) + \
                summary_tables(categories) + \
                detail()
        
        write_files()

        html = """
        <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
            "http://www.w3.org/TR/html4/strict.dtd">
        <html>
          <head>
            <title>Sisyphus&trade; summary</title>
            <style type="text/css">
              table.summary_table {
                border-style:outset;
                border-width:thin;
              }
              table.summary_table * th {
                border-style:inset;
                border-width:thin;
              }
              table.summary_table * td {
                border-style:inset;
                border-width:thin;
              }
              tr.running_test td {
                background-color:rgb(128,128,255);
              }
              tr.unfinished_test td {
                background-color:rgb(192,192,192);
              }
              tr.failed_test td {
                background-color:rgb(255,128,128);
              }
              tr.warning_test td {
                background-color:rgb(255,255,0);
              }
              tr.good_test td {
                background-color:rgb(0,255,0);
              }
              tr.invalid_test td {
                background-color:rgb(128,128,128);
              }
            </style>
          </head>
          <body>
            <h1>Sisyphus&trade; summary</h1>
            %s
          </body>
        </html>
        """ % body

        file(html_file, "w").write(html)




class TestResult(object):
    """This class is a simple structure that holds the results of running
    a TestCase.

    It has the following data members:

      State: one of the state constants saying what the result of the
        result was, if any.
      OneLiner: a one-line summary of what happened
      Data: the result data to be passed on to depending TestCase instances
      Logs: a list of plain-text (string) logs
      Files: A dictionary mapping files to str objects (which contain the
        file contents). These will end up in the test case's subdirectory
        in the report.
      HTML: Some HTML that will be included in the summary report.
        Any references of the type FILE{name} will be replaced with a
        the complete path to the file (found in Files) of the given name.
      Stopwatch: A Stopwatch instance recording times for the test case.
        A test case is not responsible for providing this--this member
        is automatically managed by the TestSupervisor.
    """

    FAILED = 1
    WARNING = 2
    GOOD = 3
    RUNNING = 4
    NOT_DONE = 5

    def __init__(self, state, oneliner, data=None, logs=[], files={}, html=None):
        self.State = state
        self.Data = data
        self.OneLiner = oneliner
        self.Logs = logs[:]
        self.Files = files.copy()
        self.HTML = html
        self.Stopwatch = None

    def state2str(state):
        if state == TestResult.FAILED:
            return "failed"
        elif state == TestResult.WARNING:
            return "warning"
        elif state == TestResult.GOOD:
            return "good"
        elif state == TestResult.RUNNING:
            return "running"
        elif state == TestResult.NOT_DONE:
            return "unfinished"
        else:
            return "invalid"
    state2str = staticmethod(state2str)




class TestCase(object):
    """Abstract interface for a test case.
    """

    def accepts_partial_results(self):
        """Return True if this TestCase is to be run even if some of
        its prerequisites failed.

        If your derived test case wants partial results, you need to
        override this.
        """
        return False

    def depends_on(self):
        """Return a list of test cases that the given one depends on.
        These test cases will be run before this one, in depth-first
        order. See TestSupervisor.

        If your derived test case has dependencies, you need to override
        this method.
        """
        return []

    def recursively_depends_on(self):
        """Return a list of all test cases that this one depends on,
        directly and indirectly.
        """
        deps = self.depends_on()
        recdeps = []
        for tc in self.depends_on():
            recdeps += tc.recursively_depends_on()
        return deps + recdeps

    def category(self):
        """Return a brief string describing the kind of test
        that is to be conducted.

        You need to override this method in a derived test case.
        """
        raise NotImplementedError

    def unique_id(self):
        """Return a unique id within the given category. 

        You need to override this method in a derived test case.
        """
        raise NotImplementedError

    def unique_name(self):
        """Returns a concatenation of category() and unique_id().
        """
        return "%s-%s" % (self.category(), self.unique_id())

    def run(self, ancestor_results, config):
        """Run the test case, return a TestResult instance.

        `ancestor_results' is a list of tuples (anc, result),
        where `anc' is the ancestor object (as found in the result
        of `depends_on()') and `result' is the `Data' member of
        the TestResult instance it returned.

        You need to override this method in a derived test case.
        """
        raise NotImplementedError




class SuccessfulTestCase(TestCase):
    """A test case that always succeeds.
    """
    def category(self):
        return "dummy"

    def get_settings_str(self):
        return "dummy"

    def unique_id(self):
        return "0"

    def run(self, ancestor_results, config):
        return TestResult(TestResult.GOOD, "OK")




class RetryingTestCase(TestCase):
    """A test case that retries the underlying operation
    a given number of times before actually failing.
    """

    def __init__(self, trycount):
        self.TryCount = trycount

    def combine_logs(self, results):
        logs = []
        for i, result in my_enumerate(results):
            logs.append("*** BEGIN ATTEMPT %d\n" % (i+1))
            logs += result.Logs
            logs.append("*** END ATTEMPT %d\n" % (i+1))
        return logs

    def run_try(self, ancestor_results, config):
        """This method has the same calling sequence and
        results as TestCase.run(). It is the underlying
        operation that is retried.

        You need to override this method in a derived test case.
        """
        raise NotImplementedError

    def run(self, ancestor_results, config):
        results = []
        while len(results) < self.TryCount:
            result = self.run_try(ancestor_results, config)
            results.append(result)
            if result.State == TestResult.GOOD:
                break

        if len(results) > 1:
            if result.State == TestResult.GOOD:
                result.State = TestResult.WARNING
            result.OneLiner += " [%d attempts]" % len(results)
        result.Logs = self.combine_logs(results)
        return result



# nek-specific tests -----------------------------------------------------
class NekCompileTestCase(RetryingTestCase):
    """This test case compiles a given Nekton example and checks
    that no errors occur. It retries the compilation three times
    to deal with the stupid Portland Group license manager.
    """

    def __init__(self, dir, usrstem, polynomial_degree, sizeu_settings={},
            compiler_defines={}):
        RetryingTestCase.__init__(self, 3)
        self.Directory = dir
        self.UsrStem = usrstem
        self.PolynomialDegree = polynomial_degree
        self.SizeuSettings = sizeu_settings.copy()
        self.SizeuSettings["lxi"] = polynomial_degree
        self.CompilerDefines = compiler_defines

    def category(self):
        return "compile"

    def get_settings_str(self):
        return "%s-%s" % (
                get_settings_str(self.SizeuSettings),
                get_settings_str(self.CompilerDefines),
                )

    def unique_id(self):
        return "%s-%s" % (self.UsrStem, self.get_settings_str())

    def _polynomial_degree(self):
        return self.PolynomialDegree
    polynomial_degree = property(_polynomial_degree)

    def run_try(self, ancestor_results, config):
        from os.path import join
        import re

        root = config["tree_root"]
        bindir = join(root, "bin")
        mydir = join(root, self.Directory)
        cleanall = join(bindir, "cleanall")
        makenek = join(bindir, "makenek")

        sizeufile = "SIZEu"
        sizeupath = join(mydir, sizeufile)
        usrfile = "%s.usr" % self.UsrStem
        usrpath = join(mydir, usrfile)

        ok, clean_log = run_command(mydir, cleanall)
        logs = [clean_log]
        if not ok:
            return TestResult(TestResult.FAILED, "cleanall failed", logs=logs)

        sizeu_log = "Changing SIZEu:\n"
        sizeu_lines = file(sizeupath, "r").readlines()
        for key, value in self.SizeuSettings.iteritems():
            sizeu_log += "\nApplying %s=%s...\n" % (key,value)
            pattern = re.compile(r"\(\s*%s\s*\=\s*[-+0-9.ea-zA-Z]+\s*\)" 
                    % re.escape(key), re.IGNORECASE)
            match_count = 0
            for i, line in my_enumerate(sizeu_lines):
                if pattern.search(line) is not None:
                    new_line = pattern.sub("(%s=%s)" % (key,value),
                            line)
                    sizeu_log += "Changing line %d:\n" % (i+1)
                    sizeu_log += "BEFORE: %s" % line
                    sizeu_log += "AFTER : %s" % new_line
                    sizeu_lines[i] = new_line
                    match_count += 1
            if match_count == 0:
                return TestResult(TestResult.FAILED, "No SIZEu match for %s" % key, 
                        logs=logs+[sizeu_log])
            if match_count > 1:
                return TestResult(TestResult.FAILED, "More than one SIZEu match for %s" % key, 
                        logs=logs+[sizeu_log])

        outf = file(sizeupath, "w")
        sizeu_contents = "".join(sizeu_lines)
        outf.write(sizeu_contents)
        outf.close()
        logs.append(sizeu_log)

        env_settings = {"FORTRAN_DEFINES":
            " ".join(["-D%s=%s" % (key, value) 
                    for key, value in self.CompilerDefines.iteritems()])
            }

        ok, make_log = run_command(mydir, [makenek, self.UsrStem],
                env_settings)
        logs.append(make_log)

        files = {}
        files[usrfile] = file(usrpath, "r").read()
        files[sizeufile] = sizeu_contents 

        html = """
        <p>
          <a href=\"FILE{%s}\">.usr file used in compilation</a>
          &middot;
          <a href=\"FILE{%s}\">SIZEu file used in compilation</a>
        </p>
        """ % (usrfile, sizeufile)

        import time
        # wait for NFS cache to settle
        time.sleep(1) 

        if ok:
            return TestResult(TestResult.GOOD, "OK", 
                    logs=logs, files=files, html=html)
        else:
            return TestResult(TestResult.FAILED, "makenek failed", 
                    logs=logs, files=files, html=html)





class NekRunTestCase(TestCase):
    """A base class for all TestCase derivatives that need to run
    Nekton code as part of their checks. This class does not check
    anything by itself.
    """
    def __init__(self, compile_tc, dir, rea, rea_parameters={}, rea_switches={},
            boxres=None):
        self.CompileTC = compile_tc
        self.Directory = dir
        self.ReaStem = rea
        self.ReaParameters = rea_parameters
        self.ReaSwitches = rea_switches
        self.BoxResolution = boxres

    def unique_id(self):
        result = "%s-%s-%s-%s" % (self.ReaStem, 
                self.compile_tc.get_settings_str(),
                get_settings_str(self.ReaParameters),
                get_settings_str(self.ReaSwitches),
                )
        if self.BoxResolution is not None:
            result += "-box%s" % ("x".join(["%d"%i for i in self.BoxResolution]))
        return result

    def _compile_tc(self):
        return self.CompileTC
    compile_tc = property(_compile_tc)

    def _directory(self):
        return self.Directory
    directory = property(_directory)

    def _rea_stem(self):
        return self.ReaStem
    rea_stem = property(_rea_stem)

    def _rea_parameters(self):
        return self.ReaParameters
    rea_parameters = property(_rea_parameters)

    def _rea_switches(self):
        return self.ReaSwitches
    rea_switches = property(_rea_switches)

    def depends_on(self):
        return [self.CompileTC]

    def setup_rea(self, rea):
        for key, value in self.ReaParameters.iteritems():
            rea.set_parameter(key, value)
        for key, value in self.ReaSwitches.iteritems():
            rea.set_switch(key, value)

    def setup_box(self, boxfname):
        import boxfile

        boxf = boxfile.BoxFile(boxfname)
        for i, res in my_enumerate(self.BoxResolution):
            if res is not None:
                boxf.NElements[i] = res
        boxf.write(boxfname)

        return TestResult(TestResult.GOOD, "OK")

    def run(self, ancestor_results, config):
        from os.path import join
        import reafile

        root = config["tree_root"]
        bindir = join(root, "bin")
        mydir = join(root, self.Directory)
        nek = join(bindir, "nek")
        rea_file = "%s.rea" % self.ReaStem
        rea_path = join(mydir, rea_file)
        map_file = "%s.map" % self.ReaStem
        map_path = join(mydir, map_file)
        box_file = "%s.box" % self.ReaStem
        box_path = join(mydir, box_file)

        rea = reafile.REAFile(rea_path)
        rea_original = rea.copy()

        try:
            self.setup_rea(rea)
            rea.write(rea_path)

            if self.BoxResolution is not None:
                self.setup_box(box_path)
                box_contents = file(box_path, "r").read()
            else:
                box_contents = None

            rea_contents = file(rea_path, "r").read()
            try:
                map_contents = file(map_path, "r").read()
            except IOError:
                map_contents = None

            ok, nek_log = run_command(mydir, [nek, "-w", self.ReaStem])

        finally:
            rea_original.write(rea_path)

        if not ok:
            return TestResult(TestResult.FAILED, "nek failed", logs=[nek_log])

        if rea.Dimensions == 2:
            if rea.get_switch("iftm"):
                html = "<p>Used TM calculation in 2D.</p>"
            elif rea.get_switch("ifte"):
                html = "<p>Used TE calculation in 2D.</p>"
            else:
                html = "<p>Used neither TM nor TE in 2D... Oops.</p>"
        else:
            html = "<p>Used full 3D calculation.</p>"

        filelinks = [
          "<a href=\"FILE{%s}\">.rea file used for computation</a>" % rea_file
          ]
        files = {rea_file: rea_contents}

        if box_contents is not None:
            files[box_file] = box_contents
            filelinks.append(
                    "<a href=\"FILE{%s}\">.box file used for computation</a>" % box_file
                    )
        if map_contents is not None:
            files[map_file] = map_contents
            filelinks.append(
                    "<a href=\"FILE{%s}\">.map file used for computation</a>" % map_file
                    )

        html += "<p>%s</p>" % " &middot; ".join(filelinks)
        return TestResult(TestResult.GOOD, "OK", logs=[nek_log],
                html=html, files=files, data=rea)




class NekTimeRunTestCase(NekRunTestCase):
    """This test case runs Nekton in timestepping mode and
    (optionally) checks that the maximum error lies below
    a certain value.
    """

    def __init__(self, compile_tc, dir, reastem, rea_parameters={}, 
            rea_switches={}, max_error=None, max_error_type="L2",
            boxres=None):
        NekRunTestCase.__init__(self, compile_tc, dir, reastem, 
                rea_parameters, rea_switches, boxres)
        self.MaxError = max_error
        self.MaxErrorType = max_error_type

    def category(self):
        return "timerun"

    def setup_rea(self, rea):
        rea.set_switch("ifeig", False)
        rea.set_switch("ifexp", False)
        rea.set_switch("ifrk4", True)
        rea.set_parameter("ifblockoutput", 1)

        NekRunTestCase.setup_rea(self, rea)

    def run(self, ancestor_results, config):
        from os.path import join
        root = config["tree_root"]
        mydir = join(root, self.Directory)

        result = NekRunTestCase.run(self, ancestor_results, config)
        if result.State == TestResult.FAILED:
            return result

        rea = result.Data

        import time

        attempts = 5
        runlog = None
        while attempts and runlog is None:
            try:
                runlog = file(join(mydir, "runlog.dat"), "r").read()
            except IOError:
                pass
            if runlog is None:
                print "WAITING..."
                time.sleep(1)
                attempts -= 1

        if runlog is None:
            result.State = TestResult.FAILED
            result.OneLiner = "Run log data file not found"
            return result

        result.Files["runlog.dat"] = runlog
        result.Data = nektools.parse_run_log(
                runlog, 
                self.CompileTC.UsrStem,
                self.CompileTC.polynomial_degree,
                self.ReaStem,
                self.BoxResolution,
                rea.get_parameter("dt"),
                self.ReaSwitches,
                self.ReaParameters)

        l2errors = result.Data[0].filter(TAG_TYPE, "L2")
        
        steps = l2errors.tag_values(TAG_STEP)
        steps.sort()
        components = l2errors.tag_values(TAG_COMPONENT)
        components.sort()
        pp = PNGPlotter()

        if HAVE_GNUPLOT:
            pp.title("Error development in time (Degree %d)" % 
                    self.CompileTC.polynomial_degree)
            pp.xlabel("Timestep")
            pp.ylabel("L2 Error")

            datasets = [Gnuplot.Data(steps, 
                    l2errors
                    .filter(TAG_COMPONENT, comp)
                    .values_sorted_by(TAG_STEP), 
                    title=comp, with="linespoints")
                for comp in components]
            pp.plot(*datasets)

        result.Files["errors.png"] = pp.get_png()

        result.HTML="<p><img src=\"FILE{errors.png}\"/><br/>\n" + \
                "<a href=\"FILE{runlog.dat}\">Original run log</a></p>\n" + \
                result.HTML

        max_error_here = max(result.Data[0]
                .filter(TAG_TYPE, self.MaxErrorType).
                values())

        if self.MaxError is not None and max_error_here > self.MaxError:
            result.State = TestResult.FAILED
            result.OneLiner = "Maximum %s error (%g) exceeded allowable value (%g)" %\
                    (self.MaxErrorType, max_error_here, self.MaxError)

        return result




class NekEigenvalueTestCase(NekRunTestCase):
    """This TestCase runs Nekton in timestepping-matrix eigenvalue
    computation mode and checks a stability criterion, namely,
    that the maximum real part of all the eigenvlaues is below
    a given value.
    """

    def __init__(self, compile_tc, dir, reastem, rea_parameters={}, 
            rea_switches={}, use_arpack=True, max_real_eigenvalue=1e-13,
            boxres=None):
        NekRunTestCase.__init__(self, compile_tc, dir, reastem, 
                rea_parameters, rea_switches, boxres)
        self.UseArpack = use_arpack
        self.MaxRealEvalue = max_real_eigenvalue

    def category(self):
        return "eig"

    def setup_rea(self, rea):
        rea.set_switch("ifeig", True)
        rea.set_switch("ifexp", False)
        rea.set_switch("ifrk4", False)

        if self.UseArpack:
            rea.set_parameter("ifarpack", 1)
        else:
            rea.set_parameter("ifarpack", 0)

        NekRunTestCase.setup_rea(self, rea)

    def parse_eigenvalues(self, edat):
        lines = edat.split("\n")[0:]
        real_evalues = []
        imag_evalues = []
        for line in lines:
            if len(line.strip()) == 0:
                continue
            words = line.split()
            assert len(words) == 2
            real_evalues.append(float(words[0]))
            imag_evalues.append(float(words[1]))
        return real_evalues, imag_evalues

    def run(self, ancestor_results, config):
        from os.path import join
        root = config["tree_root"]
        mydir = join(root, self.Directory)

        result = NekRunTestCase.run(self, ancestor_results, config)
        if result.State == TestResult.FAILED:
            return result

        try:
            eigenvalues = file(join(mydir, "fort.59"), "r").read()
        except IOError:
            result.State = TestResult.FAILED
            result.OneLiner = "Eigenvalue file not found"
            return result

        result.Files["eigenvalues.dat"] = eigenvalues
        real_evalues, imag_evalues = self.parse_eigenvalues(eigenvalues)
        
        max_re_ev = max(real_evalues)

        if max_re_ev < self.MaxRealEvalue:
            result.State = TestResult.GOOD
            result.OneLiner = "OK"
        else:
            result.State = TestResult.FAILED
            result.OneLiner = "Max real part of eigenvalue too large (%g)" \
                    % max_re_ev

        pp = PNGPlotter()
        if HAVE_GNUPLOT:
            pp.title("Eigenvalue diagram")
            pp.xlabel("Real part")
            pp.ylabel("Imaginary part")

            pp.plot(Gnuplot.Data(real_evalues, imag_evalues))
        result.Files["eigenvalues.png"] = pp.get_png()

        if self.UseArpack:
            package = "ARPACK"
        else:
            package = "LAPACK"

        result.HTML = """
          <img src=\"FILE{eigenvalues.png}\"/><br/>
        <p>Largest occurring real part of any eigenvalue: %g</p>
        <p>Computational package: %s</p>
        <p><a href=\"FILE{eigenvalues.dat}\">Raw eigenvalue data</a></p>
        %s
        """ % (max_re_ev, package, result.HTML)

        return result




class NekPMLTestCase(TestCase):
    """This test case verifies that the PML attenuates a certain 
    time-limited source within a certain time limit, i.e. that the
    maximum Linf-norm drops by at least a certain factor, say 1e-3,
    within the runtime of the underlying TimeRunTestCase

    The NekTimeRunTestCase instances form, of course, dependencies of this
    TestCase.
    """

    def __init__(self, time_test, at_timestep=3000, norm_type="L2",
            min_decay=1e-2, non_decaying_components=[]):
        self.TimeTest = time_test
        self.AtTimestep = at_timestep
        self.NormType = norm_type
        self.MinDecay = min_decay
        self.NonDecayingComponents = non_decaying_components
  
    def depends_on(self):
        return [self.TimeTest]

    def category(self):
        return "pml"

    def unique_id(self):
        return self.TimeTest.unique_id()

    def run(self, ancestor_results, config):
        (ancestor, result), = ancestor_results
        norms = result[0].filter(TAG_TYPE, self.NormType)
        components = norms.tag_values(TAG_COMPONENT)
        components.sort()
        
        ok = True
        oneliner = "OK"
        pml_rows = ""
        count_decayed = 0

        for comp in components:
            comp_norms = norms.filter(TAG_COMPONENT, comp)
            end_norm, = comp_norms\
                    .filter(TAG_STEP, self.AtTimestep) \
                    .values()
            max_norm = max(comp_norms.values())

            if max_norm == 0:
                pml_rows += "<tr><td>%s</td><td>%s</td><td>%s</td></tr>" % \
                        (comp, "None", "Max error is zero - no decay")
            else:
                decay = end_norm / max_norm 
                if decay > self.MinDecay:
                    if comp in self.NonDecayingComponents:
                        remark = "Non-decay OK"
                    else:
                        ok = False
                        oneliner = "Component %s failed decay test" % comp
                        remark = "&gt; %g - Failure" % self.MinDecay
                else:
                    remark = "&lt; %g - OK" % self.MinDecay
                    count_decayed += 1

                pml_rows += "<tr><td>%s</td><td>%g</td><td>%s</td></tr>" % \
                        (comp, decay, remark)

        html = """
        <table>
          <tr><th>Component</th><th>Decay</th><th>Remarks</th></tr>
          %s
        </table>
        """ % pml_rows

        if ok:
            if count_decayed == 0:
                state = TestResult.WARNING
                oneliner = "No component actually decayed"
            else:
                state = TestResult.GOOD
        else:
            state = TestResult.FAILED

        return TestResult(state, oneliner, html=html)




class NekPConvergenceTestCase(TestCase):
    """This test case computes the `empirical order of convergence'
    (EOC) of the error reported by NekTimeRunTestCase instances,
    and compares it to a given minimum value.

    The NekTimeRunTestCase instances form, of course, dependencies of this
    TestCase.
    """

    def __init__(self, time_tests, at_timestep=100, error_type="L2",
            min_eoc=1, non_converging_components=[]):
        self.TimeTests = time_tests
        assert all_same(time_tests, f=lambda tc: tc.directory)
        assert all_same(time_tests, f=lambda tc: tc.rea_stem)
        assert all_same(time_tests, f=lambda tc: tc.rea_parameters)
        assert all_same(time_tests, f=lambda tc: tc.rea_switches)
        self.AtTimestep = at_timestep
        self.ErrorType = error_type
        self.MinEOC = min_eoc
        self.NonConvergingComponents = non_converging_components
  
    def accepts_partial_results(self):
        return True

    def depends_on(self):
        return self.TimeTests

    def category(self):
        return "pconvergence"

    def unique_id(self):
        tc = self.TimeTests[0]
        return "%s-%s-%s" % (tc.rea_stem, get_settings_str(tc.rea_parameters), 
                get_settings_str(tc.rea_switches))

    def run(self, ancestor_results, config):
        if len(ancestor_results) < 2:
            return TestResult(TestResult.FAILED,
                    "can't compute EOC from less than two results")

        errors = reduce(
                lambda x,y: x.union(y),
                [anc_result[0] for anc, anc_result in ancestor_results])
        errors = errors.filter(TAG_TYPE, self.ErrorType)
        errors = errors.filter(TAG_STEP, self.AtTimestep)

        components = errors.tag_values(TAG_COMPONENT)
        components.sort()
        degrees = errors.tag_values(TAG_PDEGREE)
        degrees.sort()
        eoc_rows = ""
        
        ok = True
        oneliner = "OK"

        datasets = []
        count_converging = 0
        for comp in components:
            comp_errors = errors \
                    .filter(TAG_COMPONENT, comp) \
                    .values_sorted_by(TAG_PDEGREE)
            if HAVE_GNUPLOT:
                datasets.append(Gnuplot.Data(degrees, comp_errors, 
                    title=comp, with="linespoints"))

            if 0 in comp_errors:
                eoc_rows += "<tr><td>%s</td><td>%s</td><td>%s</td></tr>" % \
                        (comp, "None", "Zero in dataset")
            else:
                absc, eoc = estimate_order_of_convergence(degrees, comp_errors)

                if eoc < self.MinEOC:
                    if comp in self.NonConvergingComponents:
                        remark = "Non-convergence OK"
                    else:
                        ok = False
                        oneliner = "Component %s failed EOC test" % comp
                        remark = "&lt; %g - Failure" % self.MinEOC
                else:
                    remark = "&gt; %g - OK" % self.MinEOC
                    count_converging += 1

                eoc_rows += "<tr><td>%s</td><td>%g</td><td>%s</td></tr>" % \
                        (comp, eoc, remark)

        eoc_html = """
        <table>
          <tr><th>Component</th><th>EOC</th><th>Remarks</th></tr>
          %s
        </table>
        """ % eoc_rows

        pp = PNGPlotter()
        pp.title("Error vs. polynomial degree after %d timesteps" % 
                self.AtTimestep)
        pp.xlabel("Polynomial degree")
        pp.ylabel("%s Error" % self.ErrorType.upper())
        pp("set logscale xy")
        pp.plot(*datasets)
        errors_png = pp.get_png()

        html = "<img src=\"FILE{errors.png}\"/><br/>\n" + eoc_html

        if ok:
            if count_converging == 0:
                state = TestResult.WARNING
                oneliner = "No component actually converged"
            else:
                state = TestResult.GOOD
        else:
            state = TestResult.FAILED

        return TestResult(state, oneliner, files={"errors.png": errors_png},
                html=html)




class NekRunDataSaver(TestCase):
    """This test case saves a pickled TaggedDataTable of all the 
    run data encountered.
    """

    def __init__(self, time_tests):
        self.TimeTests = time_tests

    def accepts_partial_results(self):
        return True

    def depends_on(self):
        return self.TimeTests

    def category(self):
        return "gather"

    def unique_id(self):
        return "data"
        
    def run(self, ancestor_results, config):
        if len(ancestor_results):
            errors = reduce(
                    lambda x,y: x.union(y),
                    [anc_result[0] for anc, anc_result in ancestor_results])
            times = reduce(
                    lambda x,y: x.union(y),
                    [anc_result[1] for anc, anc_result in ancestor_results])
        else:
            errors = TaggedDataTable()
            times = TaggedDataTable()

        import cPickle as pickle

        return TestResult(TestResult.GOOD, 
                "OK", 
                files={
                "errors.pickle": pickle.dumps(errors),
                "times.pickle": pickle.dumps(times),
                })




class ScriptManager(object):
    """This instance provides the framework within which a
    sisyphus script is executed.
    
    All members of ScriptManager are available without 
    qualification in test scripts, so the method `add_test_case()' 
    is available as `add_test_case()', and `supervisor' is 
    available simply as `supervisor'. When reading the method
    signatures, be aware that you should *not* specify the `self'
    argument--it has been specified for you. (The symbols available
    in the script scope are bound methods.)
    """

    def __init__(self, supervisor):
        """The constructor of the instance. There is no need to call
        this from test scripts.
        """
        self.CompileTCs = {}
        self.Supervisor = supervisor

    def _supervisor(self):
        return self.Supervisor
    supervisor = property(_supervisor, doc="The TestSupervisor "
        "controlling the script.")

    def add_test_case(self, tc):
        """Add the given test case to the TestSupervisor that controls
        this script.
        """
        self.Supervisor.add_test_case(tc)

    def get_compile_tc(self, dir, usr_stem, polynomial_degree, 
            compiler_defines={}):
        """Get a NekCompileTestCase that compiles the code with
        the .usr file specified in `usr_stem', found in the directory
        `dir'. `usr_stem' should not include the .usr extension.
        `polynomial_degree' specifies the parameter `lxi' with which
        the code will be compiled. `compiler_defines' specifies a
        dictionary of compiler macros with values.

        The resulting NekCompileTestCase is automatically added to the 
        supervisor. 
        
        If the same NekCompileTestCase is requested twice, this
        routine will return the same one both times, eliminating
        the creation of a redundant compile task.
        """
        cdef_key = list(compiler_defines)
        cdef_key.sort()
        cdef_key = tuple(cdef_key)
        key = (dir, usr_stem, polynomial_degree, cdef_key)
        try:
            return self.CompileTCs[key]
        except KeyError:
            new_tc = NekCompileTestCase(dir, usr_stem, polynomial_degree,
                    compiler_defines=compiler_defines)
            self.add_test_case(new_tc)
            self.CompileTCs[key] = new_tc
            return new_tc

    def test_error(self, dir, usr_stem, rea_stem, max_error, degree,
            max_error_type="L2",
            at_timestep=100, dt=1e-3, rea_parameters={}, rea_switches={},
            boxres=None):
        """Create a NekTimeRunTestCase that will check if the error
        is less than `expected_error'. `max_error' is compared 
        to the maximum error over all timesteps at which it is computed,
        for all components, including `sum'. `max_error_type' gives
        the type of error, typically "L2" ("Linf" is also possible.)

        `dir', `usr_stem' and `rea_stem' specify the directory, the .usr
        and the .rea file to run. The latter two should not include the
        .usr and .rea extensions.

        `degree' specifies the polynomial degree at which the code is run.
        `at_timestep' specifes for how many time steps the code should run.
        `dt' specifies the length of one timestep.

        By passing a mapping type for `rea_parameters' and/or `rea_switches',
        you may specify additional options to set in the .rea file.
        Example: 

          test_error(..., rea_switches={"ifcentral": True, "ifupwind": False})

        `boxres' finally specifies a tuple of the same length as the 
        dimension of the problem, giving the number of elements in each
        direction that is to be written into the corresponding .box file.
        """

        iocomm = at_timestep/20
        if iocomm == 0:
            iocomm = 1

        my_rea_parameters = rea_parameters.copy()
        my_rea_switches = rea_switches.copy()
        my_rea_parameters["nsteps"] = at_timestep
        if not "iocomm" in rea_parameters:
            my_rea_parameters["iocomm"] = iocomm
        my_rea_parameters["dt"] = dt

        ct = self.get_compile_tc(dir, usr_stem, degree)
        rt = NekTimeRunTestCase(ct, dir, rea_stem,
            rea_parameters=my_rea_parameters, rea_switches=my_rea_switches,
            max_error=max_error, max_error_type=max_error_type,
            boxres=boxres)
        self.add_test_case(rt)
        return rt

    def test_p_convergence(self, dir, usr_stem, rea_stem, min_eoc, 
            non_converging_components=[],
            degrees=[3,6,9], at_timestep=100, dt=1e-3, 
            rea_parameters={}, rea_switches={},
            boxres=None):
        """Create a series of NekTimeRunTestCase instances and a final 
        NekPConvergenceTestCase that will check if the error
        decreases at least by a polynomial order of `min_eoc'.

        Failure is signalled if this is not the case and the 
        component for which it occurs is not part of 
        the list `non_converging_components'.

        `dir', `usr_stem' and `rea_stem' specify the directory, the .usr
        and the .rea file to run. The latter two should not include the
        .usr and .rea extensions.

        `degrees' specifies the polynomial degrees at which the code is run.
        `at_timestep' specifes for how many time steps the code should run.
        `dt' specifies the length of one timestep.

        By passing a mapping type for `rea_parameters' and/or `rea_switches',
        you may specify additional options to set in the .rea file.
        See test_error() for an example.

        See test_error() for the meaning of `boxres'.
        """

        iocomm = at_timestep/20
        if iocomm == 0:
            iocomm = 1

        my_rea_parameters = rea_parameters.copy()
        my_rea_switches = rea_switches.copy()
        my_rea_parameters["nsteps"] = at_timestep
        if not "iocomm" in my_rea_parameters:
            my_rea_parameters["iocomm"] = iocomm
        my_rea_parameters["dt"] = dt

        time_tests = []
        for pd in degrees:
            ct = self.get_compile_tc(dir, usr_stem, pd)
            rt = NekTimeRunTestCase(ct, dir, rea_stem,
                rea_parameters=my_rea_parameters, rea_switches=my_rea_switches,
                boxres=boxres)
            self.add_test_case(rt)
            time_tests.append(rt)
        pt = NekPConvergenceTestCase(time_tests, at_timestep=at_timestep,
            min_eoc=min_eoc, 
            non_converging_components=non_converging_components)
        self.add_test_case(pt)
        return pt

    def test_pml(self, dir, usr_stem, rea_stem, min_decay=3e-1, 
            non_decaying_components=[],
            degrees=[3], at_timestep=1000, dt=1e-3, 
            rea_parameters={}, rea_switches={},
            norm_type="Linf",
            boxres=None):
        """Create a NekTimeRunTestCase instance and a final 
        NekPMLTestCase that will check if the norm of type
        `norm_type' reaches at least `min_decay' times its
        maximum value by timestep `at_timestep'.
        
        Failure is signalled if this is not the case and the 
        component for which it occurs is not part of 
        the list `non_decaying_components'.

        The point of this is to test whether some absorbing boundary
        condition manages to attenuate a solution sufficiently to
        satisfy the above condition.

        `dir', `usr_stem' and `rea_stem' specify the directory, the .usr
        and the .rea file to run. The latter two should not include the
        .usr and .rea extensions.

        `degrees' specifies the polynomial degrees at which the code is run.

        By passing a mapping type for `rea_parameters' and/or `rea_switches',
        you may specify additional options to set in the .rea file.
        See test_error() for an example.

        See test_error() for the meaning of `boxres'.
        """

        iocomm = at_timestep/20
        if iocomm == 0:
            iocomm = 1

        my_rea_parameters = rea_parameters.copy()
        my_rea_switches = rea_switches.copy()
        my_rea_parameters["nsteps"] = at_timestep
        if not "iocomm" in rea_parameters:
            my_rea_parameters["iocomm"] = iocomm
        my_rea_parameters["dt"] = dt

        # set up time-limited gaussian source at the origin
        tstop = dt*at_timestep*0.2
        my_rea_parameters["tmodtyp"] = 3 # Cosine
        my_rea_parameters["smodtyp"] = 2 # Gaussian pulse
        my_rea_parameters["srcscal"] = 1
        my_rea_parameters["fldchs"] = 6 # Ez
        my_rea_parameters["tstart"] = 0
        my_rea_parameters["tstop"] = tstop
        # frequency: one complete cycle till tstop
        # (if this is not a complete cycle, we'll create "charges",
        # i.e. divergences in the field that don't go away,
        # defeating our test.
        my_rea_parameters[64] = 1/float(tstop) 
        my_rea_parameters[66] = 0 # at the origin
        my_rea_parameters[67] = 0
        my_rea_parameters[68] = 0
        my_rea_parameters[69] = 0 # not moving
        my_rea_parameters[70] = 0
        my_rea_parameters[71] = 0
        my_rea_parameters[69] = 0.1 # with small sigma
        my_rea_parameters[70] = 0.1
        my_rea_parameters[71] = 0.1

        tests = []
        for pd in degrees:
            ct = self.get_compile_tc(dir, usr_stem, pd)
            rt = NekTimeRunTestCase(ct, dir, rea_stem,
                rea_parameters=my_rea_parameters, 
                rea_switches=my_rea_switches,
                boxres=boxres)
            self.add_test_case(rt)
            pt = NekPMLTestCase(rt, at_timestep=at_timestep,
                    min_decay=min_decay, 
                    non_decaying_components=non_decaying_components,
                    norm_type=norm_type)
            self.add_test_case(pt)
            tests.append(pt)

        return tests

    def test_eigenvalues(self, dir, usr_stem, rea_stem, degree=3,
            rea_parameters={}, rea_switches={}, max_real_eigenvalue=1e-13,
            use_arpack=True,
            boxres=None):
        """Create a NekEigenvalueTestCase that will check if the largest
        real part of any eigenvalue is less than `max_real_eigenvalue'.

        `dir', `usr_stem' and `rea_stem' specify the directory, the .usr
        and the .rea file to run. The latter two should not include the
        .usr and .rea extensions.

        `degree' specifies the polynomial degree at which the code is run.
        `use_arpack' says whether to compute only a few relevant eigenvectors
        using ARPACK. If False, compute all of them using LAPACK, which may
        be slow and/or memory-intensive.

        By passing a mapping type for `rea_parameters' and/or `rea_switches',
        you may specify additional options to set in the .rea file.
        See test_error() for an example.

        See test_error() for the meaning of `boxres'.
        """

        my_rea_parameters = rea_parameters.copy()
        my_rea_switches = rea_switches.copy()

        cdef = {}
        if use_arpack:
            cdef["USE_ARPACK"] = 1
        else:
            cdef["CEM_EIG_MEMORY_SMALL"] = 1

        ct = self.get_compile_tc(dir, usr_stem, degree, cdef)
        et = NekEigenvalueTestCase(ct, dir, rea_stem, use_arpack=use_arpack,
                rea_parameters=my_rea_parameters, rea_switches=my_rea_switches,
                max_real_eigenvalue=max_real_eigenvalue,
                boxres=boxres)
        self.add_test_case(et)
        return et



def run():
    """The main program."""
    import os
    import os.path
    from optparse import OptionParser
    from fnmatch import fnmatchcase

    description = "Runs test cases in a directory tree, according to a given script."
    usage = "usage: %prog [options] SCRIPTFILE"
    parser = OptionParser(description=description, usage=usage)
    parser.add_option(
	    "-o", "--output", dest="results_dir", default="test-results",
	    help="Name of directory to which reports are written", 
            metavar="OUTPUT_DIR")
    parser.add_option(
	    "-p", "--purge-old", dest="purge_old", action="store_true",
	    help="Purge results directory if it exists.")
    parser.add_option(
	    "-n", "--no-run", dest="no_run", action="store_true",
	    help="Don't run, only show what test cases would be run")
    parser.add_option(
	    "-r", "--restrict", dest="restrict", default="*",
	    help="Restrict the tests to be run by using a shell pattern", 
            metavar="SHELL-PATTERN")
    parser.add_option(
	    "-s", "--script-parameterss", dest="script_params", default="",
	    help="Pass the comma separated list SCRIPT-PARAMS as parameters"
              "to the test script", 
            metavar="SCRIPT-PARAMS")

    options, args = parser.parse_args()

    if len(args) != 1:
        print>> sys.stderr, "%s takes one argument, a script file." % sys.argv[0]
        print>> sys.stderr, "Try `%s -h' for help." % sys.argv[0]
        sys.exit(1)

    if options.purge_old:
        if sys.version_info < (2,3,0):
            print "The -p option is not available if we're running on"
            print "anything less than Python 2.3. This Python is of version"
            print sys.version
            sys.exit(1)
        
        for root, dirs, files in os.walk(options.results_dir, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
                #print "rm", os.path.join(root, name)
            for name in dirs:
                os.rmdir(os.path.join(root, name))
                #print "rmdir", os.path.join(root, name)
        try:
            os.rmdir(options.results_dir)
        except OSError:
            # it's ok if there is nothing to purge.
            pass

    sup = TestSupervisor({
            "tree_root": nektools.find_root(),
            },
            options.results_dir)

    script_mgr = ScriptManager(sup)

    script_params = {}
    for param in options.script_params.split(","):
        key_val = param.split("=")

        if len(key_val) == 1:
            script_params[key_val[0]] = None
        else:
            script_params[key_val[0]] = key_val[1]

    script_scope = {
      "test_p_convergence": script_mgr.test_p_convergence,
      "test_eigenvalues": script_mgr.test_eigenvalues,
      "test_error": script_mgr.test_error,
      "test_pml": script_mgr.test_pml,

      "add_test_case": script_mgr.add_test_case,
      "get_comiple_tc": script_mgr.get_compile_tc,
      "supervisor": script_mgr.supervisor,

      "TestCase": TestCase,
      "RetryingTestCase": RetryingTestCase,
      "NekCompileTestCase": NekCompileTestCase,
      "NekRunTestCase": NekRunTestCase,
      "NekTimeRunTestCase": NekTimeRunTestCase,
      "NekEigenvalueTestCase": NekEigenvalueTestCase,
      "NekPConvergenceTestCase": NekPConvergenceTestCase,
      "NekPMLTestCase": NekPMLTestCase,

      "PARAMETERS": script_params
      }

    execfile(args[0], script_scope)

    run_tcs = [tc for tc in sup.test_cases() \
            if fnmatchcase(tc.unique_name(), options.restrict)]
    sup.restrict_to_subset(run_tcs)

    time_tests = [tc for tc in sup.test_cases() 
        if isinstance(tc, NekTimeRunTestCase)]
    sup.add_test_case(NekRunDataSaver(time_tests))

    if options.no_run:
        print>>sys.stderr, "Would run the following test cases:"
        for tc in sup.test_cases():
            print tc.unique_name()
    else:
        sup.run()

    
if __name__ == "__main__":
    import pdb
    pdb.runcall(run)

