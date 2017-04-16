import os
import subprocess
import fileinput

import pytest

DIR = os.path.dirname(os.path.realpath(__file__))
TOPDIR = os.path.abspath(os.path.join(DIR, '..'))
LOGFILE = os.path.join(TOPDIR, 'tests', 'build.log')
if os.path.isfile(LOGFILE):
    os.remove(LOGFILE)


build_command = pytest.config.getoption('build_command')
arch = pytest.config.getoption('arch')
np = pytest.config.getoption('np')
clean = pytest.config.getoption('clean')
environment = (build_command, arch, np, clean)


testdata = [
    ('2dboxper-te', '2dboxper', {4: 1}),
    ('2dboxper-tm', '2dboxper', {4: 2}),
    ('2dboxpec-te', '2dboxpec', {4: 1}),
    ('2dboxpec-tm', '2dboxpec', {4: 2}),
    ('2dboxpml-te', '2dboxpml', {4: 1}),
    ('2dboxpml-tm', '2dboxpml', {4: 2}),
    ('2ddielectric-te-1mat', '2ddielectric', {4: 1, 70: 0}),
    ('2ddielectric-te-2mat', '2ddielectric', {4: 1, 70: 1}),
    ('2ddielectric-tm-1mat', '2ddielectric', {4: 2, 70: 0}),
    ('2ddielectric-tm-2mat', '2ddielectric', {4: 2, 70: 1}),
    ('3dboxper', '3dboxper', {}),
    ('3dboxpec', '3dboxpec', {}),
    ('3dboxpml', '3dboxpml', {}),
    ('3ddielectric-1mat', '3ddielectric', {70: 0}),
    ('3ddielectric-2mat', '3ddielectric', {70: 1}),
    ('cylwave', 'cylwave', {}),
    ('drude', 'drude', {}),
    ('lorentz', 'lorentz', {}),
    ('graphene', 'graphene', {})
]
testnames = [data[0] for data in testdata]
testdata = [data[1:] + environment for data in testdata]


class BuildError(Exception):
    pass


class ReaState():
    """Context manager for changing an `.rea` file upon entering and
    restoring it upon exiting.

    """
    def __init__(self, name, values):
        self.rea = name + '.rea'
        self.values = values
        self.old_values = {}
        self.offset = 4

    def __enter__(self):
        for i, line in enumerate(fileinput.input(self.rea, inplace=True)):
            try:
                print('  {}'.format(self.values[i-self.offset+1]))
                self.old_values[i] = line
            except KeyError:
                print(line, end='')

    def __exit__(self, exc_type, exc_val, exc_tb):
        for i, line in enumerate(fileinput.input(self.rea, inplace=True)):
            try:
                print(self.old_values[i], end='')
            except KeyError:
                print(line, end='')


def build_test(name, build_command, arch, clean):
    os.chdir(os.path.join(TOPDIR, 'tests', name))
    if clean:
        cleanpath = os.path.join(TOPDIR, 'bin', 'cleanall')
        subprocess.run(cleanpath)
    buildpath = os.path.join(TOPDIR, 'bin', build_command)
    # Use buffering=1 (line buffered) so that everything is written in
    # the correct order.
    with open(LOGFILE, 'a', buffering=1) as logfile:
        logfile.write('***** build log for {} *****\n'.format(name))
        res = subprocess.run([buildpath, '-a', arch], stdout=logfile,
                             stderr=logfile)
    os.chdir(TOPDIR)
    if res.returncode != 0:
        raise BuildError("Build failed; see tests/build.log for details")


def run_test(name, rea_values, np):
    os.chdir(os.path.join(TOPDIR, 'tests', name))
    with ReaState(name, rea_values):
        execpath = os.path.join(TOPDIR, 'bin', 'nek')
        res = subprocess.run([execpath, name, np])
    os.chdir(TOPDIR)
    if res.returncode != 0:
        raise AssertionError("Error too large")


paramstring = 'name, rea_values, build_command, arch, np, clean'
@pytest.mark.parametrize(paramstring, testdata, ids=testnames)
def test_nekcem(name, rea_values, build_command, arch, np, clean):
    build_test(name, build_command, arch, clean)
    run_test(name, rea_values, np)
