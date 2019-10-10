import starks 


def test_bytecode_runtime():
    code = """
@public
def a() -> bool:
    return True
    """

    out = starks.compile_code(code, ['bytecode_runtime', 'bytecode'])

    assert len(out['bytecode']) > len(out['bytecode_runtime'])
    assert out['bytecode_runtime'][2:] in out['bytecode'][2:]
