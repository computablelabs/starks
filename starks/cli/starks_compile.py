import argparse

# TODO fix this up for the starks library
format_options_help = """Format to print, one or more of:
bytecode (default) - Deployable bytecode
bytecode_runtime   - Bytecode at runtime
abi                - ABI in JSON format
abi_python         - ABI in python format
ast                - AST in JSON format
source_map         - Vyper source map
method_identifiers - Dictionary of method signature to method identifier.
combined_json      - All of the above format options combined as single JSON output.
interface          - Print Vyper interface of a contract
external_interface - Print the External interface of a contract, used for outside contract calls.
opcodes            - List of opcodes as a string
opcodes_runtime    - List of runtime opcodes as a string
ir                 - Print Intermediate Representation in LLL
"""

def _parse_args(argv):

    warnings.simplefilter('always')

    parser = argparse.ArgumentParser(
        description='Starks programming language for Datatrust',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'input_files',
        help='Starks sourcecode to compile',
        nargs='+',
    )
    parser.add_argument(
        '--version',
        action='version',
        version=f'{starks.__version__}+commit.{starks.__commit__}',
    )
    parser.add_argument(
        '-f',
        help=format_options_help,
        default='bytecode', dest='format',
    )
    parser.add_argument(
        '--traceback-limit',
        help='Set the traceback limit for error messages reported by the compiler',
        type=int,
    )
    parser.add_argument(
        '-p',
        help='Set the root path for contract imports',
        default='.', dest='root_folder'
    )

    args = parser.parse_args(argv)

    if args.traceback_limit is not None:
        sys.tracebacklimit = args.traceback_limit
    elif VYPER_TRACEBACK_LIMIT is not None:
        sys.tracebacklimit = VYPER_TRACEBACK_LIMIT
    else:
        # Python usually defaults sys.tracebacklimit to 1000.  We use a default
        # setting of zero so error printouts only include information about where
        # an error occurred in a Vyper source file.
        sys.tracebacklimit = 0

    output_formats = tuple(uniq(args.format.split(',')))

    translate_map = {
        'abi_python': 'abi',
        'json': 'abi',
        'ast': 'ast_dict'
    }
    final_formats = []

    for f in output_formats:
        final_formats.append(translate_map.get(f, f))

    compiled = compile_files(
        args.input_files,
        final_formats,
        args.root_folder,
        args.show_gas_estimates
    )

    if output_formats == ('combined_json',):
        print(json.dumps(compiled))
        return

    for contract_data in list(compiled.values()):
        for f in output_formats:
            o = contract_data[translate_map.get(f, f)]
            if f in ('abi', 'json', 'ast', 'source_map'):
                print(json.dumps(o))
            else:
                print(o)


def uniq(seq: Iterable[T]) -> Iterator[T]:
    """
    Yield unique items in ``seq`` in order.
    """
    seen: Set[T] = set()

    for x in seq:
        if x in seen:
            continue

        seen.add(x)
        yield x

def exc_handler(contract_path: ContractPath, exception: Exception) -> None:
    print(f'Error compiling: {contract_path}')
    raise exception

def compile_files(input_files: Iterable[str],
                  output_formats: OutputFormats,
                  root_folder: str = '.',
                  show_gas_estimates: bool = False) -> OrderedDict:

    if show_gas_estimates:
        parser_utils.LLLnode.repr_show_gas = True

    root_path = Path(root_folder).resolve()
    if not root_path.exists():
        raise FileNotFoundError(f"Invalid root path - '{root_path.as_posix()}' does not exist")

    contract_sources: ContractCodes = OrderedDict()
    for file_name in input_files:
        file_path = Path(file_name)
        try:
            file_str = file_path.resolve().relative_to(root_path).as_posix()
        except ValueError:
            file_str = file_path.as_posix()
        with file_path.open() as fh:
            contract_sources[file_str] = fh.read()

    show_version = False
    if 'combined_json' in output_formats:
        if len(output_formats) > 1:
            raise ValueError("If using combined_json it must be the only output format requested")
        output_formats = ['bytecode', 'bytecode_runtime', 'abi', 'source_map', 'method_identifiers']
        show_version = True

    compiler_data = starks.compile_codes(
        contract_sources,
        output_formats,
        exc_handler=exc_handler,
        interface_codes=get_interface_codes(root_path, contract_sources)
    )
    if show_version:
        compiler_data['version'] = starks.__version__

    return compiler_data
