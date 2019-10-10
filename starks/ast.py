from itertools import (
    chain,
)
import typing

from starks.exceptions import (
    CompilerPanic,
)
from starks.settings import (
    STARKS_ERROR_CONTEXT_LINES,
    STARKS_ERROR_LINE_NUMBERS,
)
from starks.utils import (
    annotate_source_code,
)

BASE_NODE_ATTRIBUTES = (
    'node_id',
    'source_code',
    'col_offset',
    'lineno',
    'end_col_offset',
    'end_lineno',
    'src'
)


class StarksNode:
    __slots__ = BASE_NODE_ATTRIBUTES
    ignored_fields: typing.Tuple = ('ctx', )
    only_empty_fields: typing.Tuple = ()

    @classmethod
    def get_slots(cls):
        return set(chain.from_iterable(
            getattr(klass, '__slots__', [])
            for klass in cls.__class__.mro(cls)
        ))

    def __init__(self, **kwargs):

        for field_name, value in kwargs.items():
            if field_name in self.get_slots():
                setattr(self, field_name, value)
            elif value:
                raise CompilerPanic(
                    f'Unsupported non-empty value field_name: {field_name}, '
                    f' class: {type(self)} value: {value}'
                )

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        for field_name in (i for i in self.get_slots() if i not in BASE_NODE_ATTRIBUTES):
            if getattr(self, field_name, None) != getattr(other, field_name, None):
                return False
        return True

    def __repr__(self):
        cls = type(self)
        class_repr = f'{cls.__module__}.{cls.__qualname__}'

        source_annotation = annotate_source_code(
            self.source_code,
            self.lineno,
            self.col_offset,
            context_lines=STARKS_ERROR_CONTEXT_LINES,
            line_numbers=STARKS_ERROR_LINE_NUMBERS,
        )

        return f'{class_repr}:\n{source_annotation}'


class Module(StarksNode):
    __slots__ = ('body', )


class Name(StarksNode):
    __slots__ = ('id', )


class Subscript(StarksNode):
    __slots__ = ('slice', 'value')


class Index(StarksNode):
    __slots__ = ('value', )


class arg(StarksNode):
    __slots__ = ('arg', 'annotation')


class Tuple(StarksNode):
    __slots__ = ('elts', )


class FunctionDef(StarksNode):
    __slots__ = ('args', 'body', 'returns', 'name', 'decorator_list', 'pos')


class arguments(StarksNode):
    __slots__ = ('args', 'defaults', 'default')
    only_empty_fields = ('vararg', 'kwonlyargs', 'kwarg', 'kw_defaults')


class Import(StarksNode):
    __slots__ = ('names', )


class Call(StarksNode):
    __slots__ = ('func', 'args', 'keywords', 'keyword')


class keyword(StarksNode):
    __slots__ = ('arg', 'value')


class Str(StarksNode):
    __slots__ = ('s', )


class Compare(StarksNode):
    __slots__ = ('comparators', 'ops', 'left', 'right')


class Num(StarksNode):
    __slots__ = ('n', )


class NameConstant(StarksNode):
    __slots__ = ('value', )


class Attribute(StarksNode):
    __slots__ = ('attr', 'value',)


class Op(StarksNode):
    __slots__ = ('op', 'left', 'right')


class BoolOp(Op):
    __slots__ = ('values', )


class BinOp(Op):
    __slots__ = ()


class UnaryOp(Op):
    __slots__ = ('operand', )


class List(StarksNode):
    __slots__ = ('elts', )


class Dict(StarksNode):
    __slots__ = ('keys', 'values')


class Bytes(StarksNode):
    __slots__ = ('s', )


class Add(StarksNode):
    __slots__ = ()


class Sub(StarksNode):
    __slots__ = ()


class Mult(StarksNode):
    __slots__ = ()


class Div(StarksNode):
    __slots__ = ()


class Mod(StarksNode):
    __slots__ = ()


class Pow(StarksNode):
    __slots__ = ()


class In(StarksNode):
    __slots__ = ()


class Gt(StarksNode):
    __slots__ = ()


class GtE(StarksNode):
    __slots__ = ()


class LtE(StarksNode):
    __slots__ = ()


class Lt(StarksNode):
    __slots__ = ()


class Eq(StarksNode):
    __slots__ = ()


class NotEq(StarksNode):
    __slots__ = ()


class And(StarksNode):
    __slots__ = ()


class Or(StarksNode):
    __slots__ = ()


class Not(StarksNode):
    __slots__ = ()


class USub(StarksNode):
    __slots__ = ()


class Expr(StarksNode):
    __slots__ = ('value', )


class Pass(StarksNode):
    __slots__ = ()


class AnnAssign(StarksNode):
    __slots__ = ('target', 'annotation', 'value', 'simple')


class Assign(StarksNode):
    __slots__ = ('targets', 'value')


class If(StarksNode):
    __slots__ = ('test', 'body', 'orelse')


class Assert(StarksNode):
    __slots__ = ('test', 'msg')


class For(StarksNode):
    __slots__ = ('iter', 'target', 'orelse', 'body')


class AugAssign(StarksNode):
    __slots__ = ('op', 'target', 'value')


class Break(StarksNode):
    __slots__ = ()


class Continue(StarksNode):
    __slots__ = ()


class Return(StarksNode):
    __slots__ = ('value', )


class Delete(StarksNode):
    __slots__ = ('targets', )


class stmt(StarksNode):
    __slots__ = ()


class ClassDef(StarksNode):
    __slots__ = ('class_type', 'name', 'body')


class Raise(StarksNode):
    __slots__ = ('exc', )


class Slice(StarksNode):
    only_empty_fields = ('lower', )


class alias(StarksNode):
    __slots__ = ('name', 'asname')


class ImportFrom(StarksNode):
    __slots__ = ('level', 'module', 'names')
