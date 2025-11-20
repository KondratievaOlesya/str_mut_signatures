from __future__ import annotations

import zipfile
from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def data_dir() -> str:
    """Absolute path to the test data directory."""
    return str(Path(__file__).resolve().parent / "data")


@pytest.fixture(scope="session")
def output_dir() -> str:
    """
    Session-wide output directory under tests/output.

    All tests that need to write files in a persistent place should
    use this directory (or its subdirectories).
    """
    root = Path(__file__).resolve().parent
    out_dir = root / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    return str(out_dir)


@pytest.fixture
def temp_output_dir(output_dir: str, request) -> str:
    """
    Per-test temporary output directory under tests/output.

    Creates a subdirectory named after the test node, e.g.:

        tests/output/test_module__test_func/
    """
    base = Path(output_dir)
    # sanitize node name a bit (pytest uses "::" in nodeid)
    name = request.node.name.replace("/", "_").replace("\\", "_").replace("::", "__")
    d = base / name
    d.mkdir(parents=True, exist_ok=True)
    return str(d)


def pytest_addoption(parser):
    parser.addoption(
        "--update-vcf-hashes",
        action="store_true",
        default=False,
        help="Recalculate expected hashes for VCF outputs",
    )


@pytest.fixture(scope="session")
def vcf_dir(data_dir: str, output_dir: str) -> str:
    """
    Unpack tests/data/vcfs/test_input.zip into tests/output/vcfs
    and return the directory containing VCF files.

    If the directory already exists (from a previous run), it is reused.
    """
    data_path = Path(data_dir)
    zip_path = data_path / "vcfs" / "test_input.zip"
    assert zip_path.is_file(), f"Missing test input zip: {zip_path}"

    vcf_root = Path(output_dir) / "vcfs"
    vcf_root.mkdir(parents=True, exist_ok=True)

    # If directory is empty, extract; otherwise assume it's already populated
    if not any(vcf_root.iterdir()):
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(vcf_root)

    inner = list(vcf_root.iterdir())
    if len(inner) == 1 and inner[0].is_dir():
        # Zip contained a single directory; use that
        return str(inner[0])

    # Otherwise, use the top-level extraction dir
    return str(vcf_root)
