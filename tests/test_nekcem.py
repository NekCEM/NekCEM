import pytest
import os
import subprocess

PWD = os.path.dirname(os.path.realpath(__file__))
TOPDIR = os.path.abspath(os.path.join(PWD, '..'))

build_command = pytest.config.getoption('build_command')
arch = pytest.config.getoption('arch')
np = pytest.config.getoption('np')
rebuild = pytest.config.getoption('rebuild')

testdata = [(build_command, arch, np, rebuild)]


def build_example(name, build_command, arch, rebuild):
    examplepath = os.path.join(TOPDIR, 'tests', name)
    os.chdir(examplepath)
    if rebuild:
        cleanpath = os.path.join('..', '..', 'bin', 'cleanall')
        subprocess.run(cleanpath)
    buildpath = os.path.join('..', '..', 'bin', build_command)
    logfile = open('compiler.log', 'w')
    try:
        subprocess.run([buildpath, '-a', arch], check=True,
                       stdout=logfile, stderr=logfile)
    finally:
        logfile.close()
        os.chdir(TOPDIR)


def run_example(name, np):
    examplepath = os.path.join(TOPDIR, 'tests', name)
    os.chdir(examplepath)
    execpath = os.path.join('..', '..', 'bin', 'nek')
    try:
        subprocess.run([execpath, name, np], check=True)
    finally:
        os.chdir(TOPDIR)


@pytest.mark.parametrize("build_command, arch, np, rebuild", testdata)
def test_2dboxper(build_command, arch, np, rebuild):
    name = '2dboxper'
    build_example(name, build_command, arch, rebuild)
    run_example(name, np)
