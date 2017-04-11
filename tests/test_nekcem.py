import pytest
import os
import subprocess

PWD = os.path.dirname(os.path.realpath(__file__))
TOPDIR = os.path.abspath(os.path.join(PWD, '..'))
LOGFILE = os.path.join(TOPDIR, 'tests', 'build.log')

build_command = pytest.config.getoption('build_command')
arch = pytest.config.getoption('arch')
np = pytest.config.getoption('np')
rebuild = pytest.config.getoption('rebuild')

names = ['2dboxper', '2dboxpec', '2dboxpml', '3dboxper', '3dboxpec',
         'drude']
testdata = [(name, build_command, arch, np, rebuild) for name in names]

if os.path.isfile(LOGFILE):
    os.remove(LOGFILE)


class BuildError(Exception):
    pass


def build_example(name, build_command, arch, rebuild):
    examplepath = os.path.join(TOPDIR, 'tests', name)
    os.chdir(examplepath)
    if rebuild:
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


def run_example(name, np):
    examplepath = os.path.join(TOPDIR, 'tests', name)
    os.chdir(examplepath)
    execpath = os.path.join(TOPDIR, 'bin', 'nek')
    res = subprocess.run([execpath, name, np])
    os.chdir(TOPDIR)
    if res.returncode != 0:
        raise AssertionError("Error too large")


@pytest.mark.parametrize('name, build_command, arch, np, rebuild',
                         testdata, ids=names)
def test_nekcem(name, build_command, arch, np, rebuild):
    build_example(name, build_command, arch, rebuild)
    run_example(name, np)
