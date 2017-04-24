import os
import subprocess
import fileinput
import json

import pytest

DIR = os.path.dirname(os.path.realpath(__file__))
TOPDIR = os.path.abspath(os.path.join(DIR, '..'))

np = pytest.config.getoption('np')
if not np:
    raise ValueError('must specify number of processors')
clean = pytest.config.getoption('clean')
config = pytest.config.getoption('config')
environment = (np, clean, config)

with open(os.path.join(DIR, 'tests.json'), 'r') as json_data:
    testdata = json.load(json_data)
testnames = [key for key in testdata]
print([key for key in testdata])
testdata = [(value["dir"], value["usr"], value["rea"],
             value["params"]) + environment for (key, value) in
            testdata.items()]


class BuildError(Exception):
    pass


class ReaState():
    """Context manager for changing an `.rea` file upon entering and
    restoring it upon exiting.

    """
    def __init__(self, name, params):
        self.rea = name + '.rea'
        self.params = params
        self.old_params = {}
        self.offset = 4

    def __enter__(self):
        for i, line in enumerate(fileinput.input(self.rea, inplace=True)):
            try:
                print('  {}'.format(self.params[i-self.offset+1]))
                self.old_params[i] = line
            except KeyError:
                print(line, end='')

    def __exit__(self, exc_type, exc_val, exc_tb):
        for i, line in enumerate(fileinput.input(self.rea, inplace=True)):
            try:
                print(self.old_params[i], end='')
            except KeyError:
                print(line, end='')


def build_test(directory, usr, np, clean, config):
    os.chdir(os.path.join(TOPDIR, 'tests', directory))
    if config or not os.path.isfile('Makefile'):
        command = os.path.join('..', '..', 'bin', 'configurenek')
        subprocess.call([command, usr])
    if clean:
        subprocess.call(['make', 'clean'])
    args = '-j{}'.format(np)
    with open('compiler.out', 'w') as log:
        code = subprocess.call(['make', args], stdout=log, stderr=log)
    os.chdir(TOPDIR)
    if code != 0:
        raise BuildError('Build failed; see the log for details')


def run_test(directory, rea, params, np):
    os.chdir(os.path.join(TOPDIR, 'tests', directory))
    with ReaState(rea, params):
        nek = os.path.join('..', '..', 'bin', 'nek')
        code = subprocess.call([nek, rea, np])
    os.chdir(TOPDIR)
    if code != 0:
        raise AssertionError('Error too large')


paramstring = 'directory, usr, rea, params, np, clean, config'
@pytest.mark.parametrize(paramstring, testdata, ids=testnames)
def test_nekcem(directory, usr, rea, params, np, clean, config):
    build_test(directory, usr, np, clean, config)
    run_test(directory, rea, params, np)
