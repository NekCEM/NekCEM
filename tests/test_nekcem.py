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

testdata = [(build_command, arch, np, rebuild)]

if os.path.isfile(LOGFILE):
    os.remove(LOGFILE)


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
        try:
            subprocess.run([buildpath, '-a', arch], check=True,
                           stdout=logfile, stderr=logfile)
        finally:
            os.chdir(TOPDIR)


def run_example(name, np):
    examplepath = os.path.join(TOPDIR, 'tests', name)
    os.chdir(examplepath)
    execpath = os.path.join(TOPDIR, 'bin', 'nek')
    try:
        subprocess.run([execpath, name, np], check=True)
    finally:
        os.chdir(TOPDIR)


@pytest.mark.parametrize('build_command, arch, np, rebuild', testdata)
def test_2dboxper(build_command, arch, np, rebuild):
    name = '2dboxper'
    build_example(name, build_command, arch, rebuild)
    run_example(name, np)
