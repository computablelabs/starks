import copy
from typing import Dict

from eth_utils.toolz import merge

from starks.eth import constants
from starks.eth.abc import OpcodeAPI
from starks.eth.vm import mnemonics
from starks.eth.vm import opcode_values
from starks.eth.vm.logic import (
    call,
)

from starks.eth.vm.forks.frontier.opcodes import FRONTIER_OPCODES


NEW_OPCODES = {
    opcode_values.DELEGATECALL: call.DelegateCall.configure(
        __name__='opcode:DELEGATECALL',
        mnemonic=mnemonics.DELEGATECALL,
        gas_cost=constants.GAS_CALL,
    )(),
}


HOMESTEAD_OPCODES: Dict[int, OpcodeAPI] = merge(
    copy.deepcopy(FRONTIER_OPCODES),
    NEW_OPCODES
)
