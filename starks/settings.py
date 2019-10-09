import os
from typing import (
    Optional,
)

STARKS_COLOR_OUTPUT = os.environ.get('STARKS_COLOR_OUTPUT', '0') == '1'
STARKS_ERROR_CONTEXT_LINES = int(os.environ.get('STARKS_ERROR_CONTEXT_LINES', '1'))
STARKS_ERROR_LINE_NUMBERS = os.environ.get('STARKS_ERROR_LINE_NUMBERS', '1') == '1'

STARKS_TRACEBACK_LIMIT: Optional[int]

_tb_limit_str = os.environ.get('STARKS_TRACEBACK_LIMIT')
if _tb_limit_str is not None:
    STARKS_TRACEBACK_LIMIT = int(_tb_limit_str)
else:
    STARKS_TRACEBACK_LIMIT = None
