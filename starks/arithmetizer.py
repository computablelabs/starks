
opcodes_to_polys: Dict[str, str] = {
    'STOP': 'x',
    'ADD': 'x + y',
    'MUL': 'xy',
    'SUB': 'x - y', 
    'DIV': 'x / y', # Is this OK?
    }

def arithmetize(asm):
    """Compiles Assembly to a polynomial."""
    return "x"
