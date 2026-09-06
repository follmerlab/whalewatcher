import pathlib
import pytest

DATA = pathlib.Path(__file__).parent / "data"


@pytest.fixture
def data_dir():
    return DATA
