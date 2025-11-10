import os
import pytest

@pytest.fixture(scope="session")
def data_dir():
    """Provides absolute path to the test data directory."""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "data"))

@pytest.fixture(scope="session")
def output_dir():
    """Provides absolute path to the test output directory."""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "output"))

def pytest_addoption(parser):
    parser.addoption(
        "--update-vcf-hashes",
        action="store_true",
        default=False,
        help="Recalculate expected hashes for VCF outputs"
    )