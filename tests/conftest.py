def pytest_addoption(parser):
    parser.addoption('--build_command',
                     help=('which build command to use'
                           ' (e.g makenekmpi)'))
    parser.addoption('--arch',
                     help=('which architecture to build for'
                           '(e.g. linux-gnu-mpi)'))
    parser.addoption('--np',
                     help=('how many processors to run with'))
    parser.addoption('--rebuild', type=bool, default=False,
                     help=('whether to clean and rebuild tests'
                           ' before running; default is False'))
