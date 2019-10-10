import pytest

from starks import (
    compiler,
)
from starks.exceptions import (
    StructureException,
)

FAILING_CONTRACTS = [
    """
@public
@constant
@nonreentrant('lock')
def nonreentrant_foo() -> uint256:
    return 1
    """,
]


@pytest.mark.parametrize('failing_contract_code', FAILING_CONTRACTS)
def test_invalid_function_decorators(failing_contract_code):
    with pytest.raises(StructureException):
        compiler.compile_code(failing_contract_code)
