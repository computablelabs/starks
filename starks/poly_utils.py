"""This file contains a number of polynomial utility functions."""
import random
from typing import List
from typing import Dict
from typing import Tuple
from typing import Callable
from primefac import factorint
from starks.polynomial import Poly
from starks.modp import IntegersModP
from starks.polynomial import polynomials_over
from starks.euclidean import gcd
from starks.numbertype import Field
from starks.numbertype import FieldElement
from starks.numbertype import MultiVarPoly 
from starks.multivariate_polynomial import multivariates_over

def is_irreducible(polynomial: Poly, p: int) -> bool:
  """is_irreducible: Polynomial, int -> bool

  Determine if the given monic polynomial with coefficients in Z/p is
  irreducible over Z/p where p is the given integer
  Algorithm 4.69 in the Handbook of Applied Cryptography
  """
  ZmodP = IntegersModP(p)
  if polynomial.field is not ZmodP:
    raise TypeError("Given a polynomial that's not over %s, but instead %r" %
                    (ZmodP.__name__, polynomial.field.__name__))

  poly = polynomials_over(ZmodP).factory
  x = poly([0, 1])
  power_term = x
  is_unit = lambda p: p.degree() == 0

  for _ in range(int(polynomial.degree() / 2)):
    power_term = power_term.powmod(p, polynomial)
    gcd_over_Zmodp = gcd(polynomial, power_term - x)
    if not is_unit(gcd_over_Zmodp):
      return False

  return True

def generate_irreducible_polynomial(modulus: int, degree: int) -> Poly:
  """ 
  Generate a random irreducible polynomial of a given degree over Z/p, where p
  is given by the integer 'modulus'. This algorithm is expected to terminate
  after 'degree' many irreducibility tests. By Chernoff bounds the probability
  it deviates from this by very much is exponentially small.
  """
  Zp = IntegersModP(modulus)
  Polynomial = polynomials_over(Zp)

  while True:
    coefficients = [Zp(random.randint(0, modulus - 1)) for _ in range(degree)]
    random_monic_polynomial = Polynomial(coefficients + [Zp(1)])

    if is_irreducible(random_monic_polynomial, modulus):
      return random_monic_polynomial

def generate_primitive_polynomial(modulus: int, degree: int) -> Poly:
  """Generates a primitive polynomial over Z/modulus.
  
  Follows algorithm 4.78 in the Handbook of Applied Cryptography
  (http://math.fau.edu/bkhadka/Syllabi/A%20handbook%20of%20applied%20cryptography.pdf).
  Generates a random irreducible polynomial and then checks if it's prime.
  
  """
  Zp = IntegersModP(modulus)
  Polynomial = polynomials_over(Zp)
  while True:
    irred_poly = generate_irreducible_polynomial(modulus, degree)
    if is_primitive(irred_poly, modulus, degree):
      return irred_poly

def is_primitive(irred_poly: Poly, modulus: int, degree: int) -> bool:
  """Returns true if given polynomial is primitve.
  
  Follows algorithm 4.78 in the Handbook of Applied Cryptography
  (http://math.fau.edu/bkhadka/Syllabi/A%20handbook%20of%20applied%20cryptography.pdf).
  """
  # All primitive polynomials are irreducible
  if not is_irreducible(irred_poly, modulus):
    return False
  # factorize p^m - 1
  prime_factors = factorint(modulus**degree - 1)
  # This is returned as dictionary with multiplicities. Turn into list
  prime_factors = [int(factor) for factor in prime_factors.keys()]
  Zp = IntegersModP(modulus)
  polysOver = polynomials_over(Zp)
  # TODO(rbharath): This is x right?
  x = polysOver([0, 1])
  # This is 1 right?
  one = polysOver([1])
  for i, factor in enumerate(prime_factors):
    power = (modulus**degree - 1) // factor
    # TODO(rbharath): This might need a smarter power implementation. In
    # particular, might need to take the modulus at intermediate powers.
    l_x = (x**power) % irred_poly
    if l_x == one:
      return False
  return True

def construct_multivariate_dirac_delta(field: Field, values: List[FieldElement]) -> MultiVarPoly:
  """Constructs the multivariate dirac delta polynomial at 0.

  1_0(x) = \prod_{i=1}^n (1 - x_i^{q-1})

  This can be generalized into the the dirac polynomial at y as follows.

  1_y(x)\prod_{i=1}^n (1 - (x_i - y_i)^{q-1})
  """
  n = len(values)
  multi = multivariates_over(field, n).factory
  q = field.field_size
  base = field(1)
  for i, val in enumerate(values):
    # ith_term = (0,...1,...0) with the 1 in the ith-term
    ith_term = [0] * n
    ith_term[i] = 1 
    term = multi({tuple(ith_term): 1})
    term = field(1) - term**(q-1)
    base = base * term
  return base


def construct_multivariate_coefficients(step_fn: Callable) -> Dict[Tuple[int, ...], FieldElement]:
  """Transforms a function over vector of finite fields into a polynomial.

  Every function f: F_q^n -> F_q is a polynomial if F is a finite field of size
  q. (See Lemma 7 of http://math.uga.edu/~pete/4400ChevalleyWarning.pdf). The
  key trick used in this transformation is the creation of a "dirac-delta"
  multivariate polynomial which is 1 iff all n of its inputs are 0.

  1_0(x) = \prod_{i=1}^n (1 - x_i^{q-1})

  Why does this make sense? For any non-zero element x in F_q, x^{q-1} = 1. How
  can we convert an aribtrary function using these dirac-delta polynomials?

  P_f(x) = \sum_{y \in F_q^n} f(y) \prod_{i=1}^n (1 - (x_i - y_i)^{q-1})

  The idea is that we construct the polynomial term-wise.
  """


def multi_inv(field, values):
  """Use one field inversion to invert many values simultaneously.
  
  TODO(rbharath): Find a reference for this algorithm.
  """
  partials = [field(1)]
  for val in values:
    if val == 0:
      mul_value = 1
    else:
      mul_value = val
    partials.append(partials[-1] * mul_value)
  assert len(partials) == len(values) + 1
  inv = 1 / partials[-1]
  outputs = [0] * len(values)
  for i in range(len(values), 0, -1):
    outputs[i - 1] = partials[i - 1] * inv if values[i - 1] else 0
    if values[i-1] != 0:
      inv = inv * values[i - 1]
  return outputs

def zpoly(modulus, roots):
  """Build a polynomial with the specified roots over Z/modulus.
  
  TODO(rbharath): Find a reference for this implementation. 
  """
  mod = IntegersModP(modulus)
  polysOverMod = polynomials_over(mod).factory
  root = [mod(1)]
  for x in roots:
    root.insert(0, mod(0))
    for j in range(len(root) - 1):
      root[j] -= root[j + 1] * x
  return polysOverMod(root)

def lagrange_interp(modulus, xs, ys):
  """
  Given p+1 y values and x values with no errors, recovers the original
  p+1 degree polynomial. Lagrange interpolation works roughly in the following way.

  1. Suppose you have a set of points, eg. x = [1, 2, 3], y = [2, 5, 10]
  2. For each x, generate a polynomial which equals its corresponding
     y coordinate at that point and 0 at all other points provided.
  3. Add these polynomials together.
  """
  # Generate master numerator polynomial, eg. (x - x1) * (x - x2) * ... * (x - xn)
  root = zpoly(modulus, xs)
  mod = IntegersModP(modulus)
  polysOverMod = polynomials_over(mod).factory
  assert len(root) == len(ys) + 1
  # print(root)
  # Generate per-value numerator polynomials, eg. for x=x2,
  # (x - x1) * (x - x3) * ... * (x - xn), by dividing the master
  # polynomial back by each x coordinate
  nums = [root / polysOverMod([-x, 1]) for x in xs]
  # Generate denominators by evaluating numerator polys at each x
  denoms = [nums[i](xs[i]) for i in range(len(xs))]
  invdenoms = multi_inv(mod, denoms)
  # Generate output polynomial, which is the sum of the per-value numerator
  # polynomials rescaled to have the right y values
  b = [0 for y in ys]
  for i in range(len(xs)):
    yslice = ys[i] * invdenoms[i]
    num_coefficients = nums[i].coefficients
    for j in range(len(ys)):
      if num_coefficients[j] and ys[i]:
        b[j] += num_coefficients[j] * yslice
  return polysOverMod(b)

# Optimized version of the above restricted to deg-4 polynomials
#def lagrange_interp_4(self, xs, ys):
def lagrange_interp_4(modulus, xs, ys):
  mod = IntegersModP(modulus)
  polysOverMod = polynomials_over(mod).factory
  x01, x02, x03, x12, x13, x23 = \
      xs[0] * xs[1], xs[0] * xs[2], xs[0] * xs[3], xs[1] * xs[2], xs[1] * xs[3], xs[2] * xs[3]
  m = modulus
  eq0 = polysOverMod([-x12 * xs[3], (x12 + x13 + x23), -xs[1] - xs[2] - xs[3], 1])
  eq1 = polysOverMod([-x02 * xs[3], (x02 + x03 + x23), -xs[0] - xs[2] - xs[3], 1])
  eq2 = polysOverMod([-x01 * xs[3], (x01 + x03 + x13), -xs[0] - xs[1] - xs[3], 1])
  eq3 = polysOverMod([-x01 * xs[2], (x01 + x02 + x12), -xs[0] - xs[1] - xs[2], 1])
  e0 = eq0(xs[0])
  e1 = eq1(xs[1])
  e2 = eq2(xs[2])
  e3 = eq3(xs[3])
  e01 = e0 * e1
  e23 = e2 * e3
  invall = 1 / (e01 * e23)
  inv_y0 = ys[0] * invall * e1 * e23 
  inv_y1 = ys[1] * invall * e0 * e23 
  inv_y2 = ys[2] * invall * e01 * e3
  inv_y3 = ys[3] * invall * e01 * e2
  return polysOverMod([
      (eq0.coefficients[i] * inv_y0 + eq1.coefficients[i] * inv_y1 + eq2.coefficients[i] * inv_y2 + eq3.coefficients[i] * inv_y3)
      for i in range(4)
  ])


# TODO(rbharath): Does this make a noticeable speed difference? If so add back in later.
## Optimized poly evaluation for degree 4
#def eval_quartic(self, p, x):
#  xsq = x * x % self.modulus
#  xcb = xsq * x
#  return (p[0] + p[1] * x + p[2] * xsq + p[3] * xcb) % self.modulus

# Optimized version of the above restricted to deg-2 polynomials
#def lagrange_interp_2(self, xs, ys):
def lagrange_interp_2(modulus, xs, ys):
  mod = IntegersModP(modulus)
  polysOverMod = polynomials_over(mod).factory
  ###############################################
  if not isinstance(xs, list):
    xs = xs.coefficients
  if not isinstance(ys, list):
    ys = ys.coefficients
  ###############################################
  m = modulus
  eq0 = polysOverMod([-xs[1], 1])
  eq1 = polysOverMod([-xs[0], 1])
  e0 = eq0(xs[0])
  e1 = eq1(xs[1])
  invall = 1/(e0 * e1)
  inv_y0 = ys[0] * invall * e1
  inv_y1 = ys[1] * invall * e0
  return polysOverMod([(eq0.coefficients[i] * inv_y0 + eq1.coefficients[i] * inv_y1) for i in range(2)])

def multi_interp_4(modulus, xsets, ysets):
  """Optimized version of the above restricted to deg-4 polynomials"""
  mod = IntegersModP(modulus)
  polysOverMod = polynomials_over(mod).factory
  data = []
  invtargets = []
  for xs, ys in zip(xsets, ysets):
    x01, x02, x03, x12, x13, x23 = \
        xs[0] * xs[1], xs[0] * xs[2], xs[0] * xs[3], xs[1] * xs[2], xs[1] * xs[3], xs[2] * xs[3]
    m = modulus
    eq0 = polysOverMod([-x12 * xs[3], (x12 + x13 + x23), -xs[1] - xs[2] - xs[3], 1])
    eq1 = polysOverMod([-x02 * xs[3], (x02 + x03 + x23), -xs[0] - xs[2] - xs[3], 1])
    eq2 = polysOverMod([-x01 * xs[3], (x01 + x03 + x13), -xs[0] - xs[1] - xs[3], 1])
    eq3 = polysOverMod([-x01 * xs[2], (x01 + x02 + x12), -xs[0] - xs[1] - xs[2], 1])
    e0 = eq0(xs[0])
    e1 = eq1(xs[1])
    e2 = eq2(xs[2])
    e3 = eq3(xs[3])
    data.append([ys, eq0, eq1, eq2, eq3])
    invtargets.extend([e0, e1, e2, e3])
  invalls = multi_inv(mod, invtargets)
  o = []
  for (i, (ys, eq0, eq1, eq2, eq3)) in enumerate(data):
    invallz = invalls[i * 4:i * 4 + 4]
    inv_y0 = ys[0] * invallz[0]
    inv_y1 = ys[1] * invallz[1]
    inv_y2 = ys[2] * invallz[2]
    inv_y3 = ys[3] * invallz[3]
    o.append(polysOverMod([(eq0.coefficients[i] * inv_y0 + eq1.coefficients[i] * inv_y1 + eq2.coefficients[i] * inv_y2 +
               eq3.coefficients[i] * inv_y3) for i in range(4)]))
  return o
