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
    ConstantinopleTransaction,
)


class ConstantinopleBlock(ByzantiumBlock):
    transaction_class = ConstantinopleTransaction
    fields = [
        ('header', BlockHeader),
        ('transactions', CountableList(transaction_class)),
        ('uncles', CountableList(BlockHeader))
    ]
