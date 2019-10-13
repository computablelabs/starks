from rlp.sedes import (
    CountableList,
)
from starks.eth.rlp.headers import (
    BlockHeader,
)
from starks.eth.vm.forks.petersburg.blocks import (
    PetersburgBlock,
)

from .transactions import (
    IstanbulTransaction,
)


class IstanbulBlock(PetersburgBlock):
    transaction_class = IstanbulTransaction
    fields = [
        ('header', BlockHeader),
        ('transactions', CountableList(transaction_class)),
        ('uncles', CountableList(BlockHeader))
    ]
