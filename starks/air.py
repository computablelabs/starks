"""This file holds the classes necessary to provide the "algebraic intermediate
representation" (AIR) of the computation. Informally, the key step needed here is the
"arithmetization" of the computation. That is, the computation must be
configured to work within a given finite field. That is, the states of the
computation will be vectors each of whose elements are members of a given
finite field.

In addition, the "transition relation" of how a step of the computation
transforms a given state vector to the next state vector must be represented as
a polynomial relation P such that if x_n and x_{n+1} are two states P(x_n,
x_{n+1}) == 0 if and only if the transition is valid. In addition, boundary
constraint B encode constraints placed on the input and output of the
computation and are also polynomials.

A witness for an AIR is a valid execution trace for the computation. This is a
sequence of states (each of which is a vector of field elements). Note that for
a computation of length T with a state of size w, the trace is O(wT) which
could possibly be very large.
"""

from starks.utils import is_a_power_of_2

def get_computational_trace(inp, steps, constants, step_fn):
  """Get the computational trace for the algebraic intermediate representation.

  Parameters
  ----------
  f: Field
    TODO(rbharath): This shouldn't have to be passed explicitly
  inp: list 
    The input state for the computation
  steps: Int
    The number of steps in the computation
  constants: List
    List of constants defining the computation in question
  step_fn: Function
    A function which maps one state to the next state.
  """
  computational_trace = [inp]
  for i in range(steps - 1):
    poly_constants = constants[i]
    # TODO(rbharath): Is there off-by-one error on round_contants?
    next_state = step_fn(computational_trace[-1], poly_constants)
    computational_trace.append(next_state)
  output = computational_trace[-1]
  print('Done generating computational trace')
  return computational_trace, output

class Computation(object):
  """A simple class defining a computation.
  
  More formally, this class holds an instance of the AIR problem. Recall that
  an instance of the AIR problem is a tuple.

    x = (F, T, w, Ps, C, B)

  Let's define each of these terms in sequence.

  - F is a finite field.
  - T is an integer representing a bound on running time.
  - w is the width of the computation (the dimensionality of the state space)
  - Polys is a list of polynomial constraints. Each polynomial is an element of
    F[X_1,..,X_w,Y_1,...,Y_w] and represents a transition constraint. There are s constraints in Polys
  - C is a monotone boolean circuit. (TODO(rbharath): What does this circuit mean and where is it used?). This encodes valid transition conditions among polynomials. For now, this is simply the AND function, requiring that all constraints in Polys must be met.
  - B is a list of boundary constraints. Each boundary constraint is a tuple
    (i, j, alpha), where i \in [T], j \in [w], alpha \in F.

  Fields
  ------
  dims: Int
    The dimensionality of the state space for the computation.
  inp: Int or List
    Either a single int or a list of integers of length dims
  steps: Int
    An int holding the number of steps of this computation.
  output: Int of List
    Either a single int or a list of integers of length dims
  constants: List
    A list of constants. Each element of constants must be a
    list of length steps
  step_fn: Function
    A function that maps a computation state to the next
    state. A state here is either an int of a list of ints of
    length dims.
  constraint_degree: int
    The degree of the constraint being considered
  extension_factor: Int
    TODO(rbharath): Can this be removed?
  """
  def __init__(self, field, dims, inp, steps, constants, step_fn,
      constraint_degree, extension_factor):
    # TODO(rbharath): Is it necessary to store field explicitly
    self.field = field
    self.dims = dims
    # Handle 1-d case
    if isinstance(inp, int):
      inp = [inp]

    # TODO(rbharath): Can this be removed?
    # Some constraints to make our job easier
    assert steps <= 2**32 // extension_factor
    assert is_a_power_of_2(steps)

    self.inp = inp
    self.steps = steps
    # Check that the constants have right length
    for poly_constants in constants:
      assert len(poly_constants) <= steps
    self.constants = constants

    self.step_fn = step_fn
    self.constraint_degree = constraint_degree
    self.computational_trace, self.output = get_computational_trace(
        inp, steps, constants, step_fn)
    self.extension_factor = extension_factor

    # The AIR variables. TODO(rbharath): Swap the STARK library to use these
    # fields for consistency.
    self.F = field
    self.T = steps
    self.w = dims
    self.Polys = self.generate_constraint_polynomials(self, step_fn)
    self.C = self.generate_monotone_circuit(self, self.Polys)
    self.B = self.generate_boundary_constraint(inp)
  
  def get_witness(self):
    """Returns the witness (computational trace) for this computation."""
    return self.computational_trace

  def generate_boundary_constraint(self, inp):
    pass

  def generate_monotone_circuit(self, Polys):
    pass

  def generate_constraint_polynomials(self, step_fn):
    """Constructs the constraint polynomials for this AIR instance.

    A constraint polynomial is in F[X_1,..,X_w, Y_1,.., Y_w]. An AIR instance
    can hold a set of constraint polynomials each of which enforces a
    transition constraint. Intuitively, X_1,...,X_1 is the vector of current
    states and Y_1,...,Y_w is the vector of the next state. It can be useful at
    times to have more than one constraint for enforcing various transition
    properties.
    """
    pass
