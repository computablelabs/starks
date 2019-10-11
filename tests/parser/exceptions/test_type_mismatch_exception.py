import pytest
from pytest import (
    raises,
)

from starks import (
    compiler,
)
from starks.exceptions import (
    TypeMismatchException,
)

fail_list = [
    """
@public
def test_func() -> int128:
    return (1, 2)
    """,
]


@pytest.mark.parametrize('bad_code', fail_list)
def test_type_mismatch_exception(bad_code):
    with raises(TypeMismatchException):
        compiler.compile_code(bad_code)
