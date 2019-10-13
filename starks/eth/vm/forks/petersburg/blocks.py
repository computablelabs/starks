from rlp.sedes import (
    CountableList,
)
from starks.eth.rlp.headers import (
    BlockHeader,
)
from starks.eth.vm.forks.byzantium.blocks import (
    ByzantiumBlock,
)

from .transactions import (
    PetersburgTransaction,
)


class PetersburgBlock(ByzantiumBlock):
    transaction_class = PetersburgTransaction
    fields = [
        ('header', BlockHeader),
        ('transactions', CountableList(transaction_class)),
        ('uncles', CountableList(BlockHeader))
    ]
