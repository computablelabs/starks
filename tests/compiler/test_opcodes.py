import starks 


def test_opcodes():
    code = """
@public
def a() -> bool:
    return True
    """

    out = starks.compile_code(code, ['opcodes_runtime', 'opcodes'])

    assert len(out['opcodes']) > len(out['opcodes_runtime'])
    assert out['opcodes_runtime'] in out['opcodes']
