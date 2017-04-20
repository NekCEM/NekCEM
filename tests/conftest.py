def pytest_addoption(parser):
    parser.addoption('--np',
                     help=('how many processors to run with'))
    parser.addoption('--clean', action="store_true",
                     help=('whether to clean before building;'
                           ' default is False'))
    parser.addoption('--config', action="store_true",
                     help=('whether to run nekconfigure'))
