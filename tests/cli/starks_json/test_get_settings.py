#!/usr/bin/env python3

import pytest

from starks.cli.starks_json import (
    get_input_dict_settings,
)
from starks.exceptions import (
    JSONError,
)


def test_unknown_evm():
    with pytest.raises(JSONError):
        get_input_dict_settings({'settings': {'evmVersion': "foo"}})


@pytest.mark.parametrize('evm_version', ['homestead', 'tangerineWhistle', 'spuriousDragon'])
def test_early_evm(evm_version):
    with pytest.raises(JSONError):
        get_input_dict_settings({'settings': {'evmVersion': evm_version}})


@pytest.mark.parametrize('evm_version', ['byzantium', 'constantinople', 'petersburg'])
def test_valid_evm(evm_version):
    get_input_dict_settings({'settings': {'evmVersion': evm_version}})
