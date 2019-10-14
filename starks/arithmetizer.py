from typing import Dict
from starks.poly_utils import multivariates_over
from starks.utils import generate_Xi_s

opcodes_to_polys: Dict[str, str] = {
    'STOP': 'x',
    'ADD': 'x + y',
    'MUL': 'xy',
    'SUB': 'x - y', 
    'DIV': 'x / y', # Is this OK?
    }

def arithmetize(asm, field, width):
    """Compiles Assembly to a polynomial."""
    polysOver = multivariates_over(field, width).factory
    Xs  = generate_Xi_s(field, width)
    return "x"
