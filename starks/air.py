"""This file holds the classes necessary to provide the
"algebraic intermediate representation" (AIR) of the
computation. Informally, the key step needed here is the
"arithmetization" of the computation. That is, the computation
must be configured to work within a given finite field. That
is, the states of the computation will be vectors each of
whose elements are members of a given finite field.

In addition, the "transition relation" of how a step of the
computation transforms a given state vector to the next state
vector must be represented as a polynomial relation P such
that if x_n and x_{n+1} are two states P(x_n, x_{n+1}) == 0 if
and only if the transition is valid. In addition, boundary
constraint B encode constraints placed on the input and output
of the computation and are also polynomials.

A witness for an AIR is a valid execution trace for the
computation. This is a sequence of states (each of which is a
vector of field elements). Note that for a computation of
length T with a state of size w, the trace is O(wT) which
could possibly be very large.
"""

from typing import List
from typing import Tuple
from starks.utils import is_a_power_of_2
from starks.utils import generate_Xi_s
from starks.poly_utils import multivariates_over
from sympy import *

def get_computational_trace(inp, steps, width, step_polys):
  """Get the computational trace for the algebraic intermediate representation.

  Parameters
  ----------
  f: Field
    TODO(rbharath): This shouldn't have to be passed explicitly
  inp: list 
    The input state for the computation
  steps: Int
    The number of steps in the computation
  step_polys: Function
    A function which maps one state to the next state.
  """
  computational_trace = [inp]
  for i in range(steps - 1):
    # TODO(rbharath): Is there off-by-one error on round_contants?
    next_state = [step_polys[i](computational_trace[-1]) for i in range(width)]
    computational_trace.append(next_state)
  output = computational_trace[-1]
  print('Done generating computational trace')
  return computational_trace, output

class AIR(object):
  """A simple class defining the algebraic intermediate representation of a computation.
  
  More formally, this class holds an instance of the AIR problem. Recall that
  an instance of the AIR problem is a tuple.

    x = (F, T, w, Ps, C, B)

  Let's define each of these terms in sequence.

  - F is a finite field.
  - T is an integer representing a bound on running time.
  - w is the width of the computation (the dimensionality of the state space)
  - Polys is a list of polynomial constraints. Each polynomial is an element of
    F[X_1,..,X_w,Y_1,...,Y_w] and represents a transition constraint. There are s constraints in Polys
  - C is a monotone boolean circuit. This encodes valid transition conditions among polynomials. For now, this is simply the AND function, requiring that all constraints in Polys must be met.
  - B is a list of boundary constraints. Each boundary constraint is a tuple
    (i, j, alpha), where i \in [T], j \in [w], alpha \in F.

  Fields
  ------
  width: Int
    The dimensionality of the state space for the computation.
  inp: Int or List
    Either a single int or a list of integers of length width 
  steps: Int
    An int holding the number of steps of this computation.
  output: Int of List
    Either a single int or a list of integers of length width 
  step_polys: Poly 
    A function that maps a computation state to the next
    state. A state here is either an int of a list of ints of
    length width.
  extension_factor: Int
    TODO(rbharath): Can this be removed?
  """
  def __init__(self, field, width, inp, steps, step_polys, extension_factor):
    self.field = field
    self.width = width
    # Handle 1-d case
    if isinstance(inp, int):
      inp = [inp]

    # Some constraints to make our job easier
    self.t = 9 # this is only a sample value, we need to change it based on a computation
    assert steps == 2**self.t-1
    #assert is_a_power_of_2(steps)

    self.inp = inp
    self.steps = steps

    self.step_polys = step_polys
    self.computational_trace, self.output = get_computational_trace(
        inp, steps, width, step_polys)
    self.extension_factor = extension_factor

    # The AIR variables. TODO(rbharath): Swap the STARK library to use these
    # fields for consistency.
    self.F = field
    self.T = steps
    self.w = width 
    self.Polys = self.generate_constraint_polynomials()
    self.C = self.generate_monotone_circuit(self.Polys)
    self.d = 10 # this is only a sample value, we need to change it based on a computation
    self.CDegree = self.C_degree(self.Polys)
    assert self.CDegree <= 2**self.d
    self.B = self.generate_boundary_constraints()
  
  def generate_witness(self):
    """Returns the witness (computational trace) for this computation."""
    return [[self.computational_trace[i][j] for i in range(self.steps)] for j in range(self.w)]

  def generate_boundary_constraints(self) -> List[Tuple]:
    boundary_constraints = []
    for ind in range(self.w):
      # (i, j, alpha) = (0, ind, inp[ind])
      # This contraints specifies that the input must be fixed
      boundary_constraints.append((0, ind, self.inp[ind]))
    return boundary_constraints

  def generate_monotone_circuit(self, Polys):
    """
    building a monotone boolean circuit over variables Z1, . . . , Zs with multi-input AND and OR gates
    """
    output_eval = []
    for i in range(self.steps - 1):
      Wt = self.computational_trace[i]
      Wt1 = self.computational_trace[i+1]
      temp = []
      for j in range(self.width):
        temp.append(Polys[j](Wt+Wt1))
      output_eval.append(temp)

    output = []
    for i in range(self.steps - 1):
      t_temp = True;
      for j in range(self.width):
        t_temp = t_temp and output_eval[i][j] == 0
      output.append(t_temp)

    for i in range(self.steps - 1):
      if output[i] == False:
        return False

    return True

  def C_degree(self, Polys):
    """
    the degree of the monotone boolean circuit
    """
    return max([Polys[j].degree() for j in range(self.width)])


  def get_degree(self):
    return max([poly.degree() for poly in self.step_polys])

  def generate_constraint_polynomials(self):
    """Constructs the constraint polynomials for this AIR instance.

    A constraint polynomial is in F[X_1,..,X_w, Y_1,.., Y_w].
    An AIR instance can hold a set of constraint polynomials
    each of which enforces a transition constraint.
    Intuitively, X_1,...,X_1 is the vector of current states
    and Y_1,...,Y_w is the vector of the next state. It can be
    useful at times to have more than one constraint for
    enforcing various transition properties.
    """
    # Let's make polynomials in F[x_1,..,x_w,y_1,...,y_w]]
    multi = multivariates_over(self.field, 2*self.width).factory
    multivars = generate_Xi_s(self.field, 2*self.width)
    Xs = multivars[:self.width]
    Ys = multivars[self.width:]
    constraints = []
    for i in range(self.width):
      step_poly = self.step_polys[i]
      Y_i = Ys[i]
      # The constraint_poly
      poly = Y_i - step_poly(Xs)
      constraints.append(poly)
    return constraints
