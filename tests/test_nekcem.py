import os
import subprocess
import fileinput

import pytest

DIR = os.path.dirname(os.path.realpath(__file__))
TOPDIR = os.path.abspath(os.path.join(DIR, '..'))

np = pytest.config.getoption('np')
if not np:
    raise ValueError('must specify number of processors')
clean = pytest.config.getoption('clean')
config = pytest.config.getoption('config')
environment = (np, clean, config)


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


def build_test(name, np, clean, config):
    os.chdir(os.path.join(TOPDIR, 'tests', name))
    if config or not os.path.isfile('Makefile'):
        command = os.path.join('..', '..', 'bin', 'configurenek')
        subprocess.call([command])
    if clean:
        subprocess.call(['make', 'clean'])
    args = '-j{}'.format(np)
    with open('compiler.out', 'w') as log:
        code = subprocess.call(['make', args], stdout=log, stderr=log)
    os.chdir(TOPDIR)
    if code != 0:
        raise BuildError('Build failed; see the log for details')


def run_test(name, rea_values, np):
    os.chdir(os.path.join(TOPDIR, 'tests', name))
    with ReaState(name, rea_values):
        nek = os.path.join('..', '..', 'bin', 'nek')
        code = subprocess.call([nek, name, np])
    os.chdir(TOPDIR)
    if code != 0:
        raise AssertionError('Error too large')


paramstring = 'name, rea_values, np, clean, config'
@pytest.mark.parametrize(paramstring, testdata, ids=testnames)
def test_nekcem(name, rea_values, np, clean, config):
    build_test(name, np, clean, config)
    run_test(name, rea_values, np)
