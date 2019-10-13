from rlp.sedes import (
    CountableList,
)
from starks.eth.rlp.headers import (
    BlockHeader,
)
from starks.eth.vm.forks.frontier.blocks import (
    FrontierBlock,
)
from .transactions import (
    HomesteadTransaction,
)


class HomesteadBlock(FrontierBlock):
    transaction_class = HomesteadTransaction
    fields = [
        ('header', BlockHeader),
        ('transactions', CountableList(transaction_class)),
        ('uncles', CountableList(BlockHeader))
    ]
