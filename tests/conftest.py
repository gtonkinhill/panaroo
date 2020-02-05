def pytest_addoption(parser):
    parser.addoption(
        "--datafolder",
        action="append",
        default=[],
        help="data folder to pass to test functions",
    )


def pytest_generate_tests(metafunc):
    if "datafolder" in metafunc.fixturenames:
        metafunc.parametrize("datafolder",
                             metafunc.config.getoption("datafolder"))
