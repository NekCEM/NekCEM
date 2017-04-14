def pytest_addoption(parser):
    parser.addoption('--build_command',
                     help=('which build command to use'
                           ' (e.g makenekmpi)'))
    parser.addoption('--arch',
                     help=('which architecture to build for'
                           '(e.g. linux-gnu-mpi)'))
    parser.addoption('--np',
                     help=('how many processors to run with'))
    parser.addoption('--clean', type=bool, default=False,
                     help=('whether to clean before building;'
                           ' default is False'))
