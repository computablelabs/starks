from eth_utils import (
    ValidationError,
)

from starks.eth.constants import (
    SECPK1_N,
)

from starks.eth.abc import (
    SignedTransactionAPI,
    StateAPI,
)
from starks.eth.vm.forks.frontier.validation import (
    validate_frontier_transaction,
)


def validate_homestead_transaction(state: StateAPI,
                                   transaction: SignedTransactionAPI) -> None:
    if transaction.s > SECPK1_N // 2 or transaction.s == 0:
        raise ValidationError("Invalid signature S value")

    validate_frontier_transaction(state, transaction)
