"""This file contains a number of polynomial utility functions."""
from starks.modp import IntegersModP
from starks.polynomial import polynomials_over

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
  #return [x % self.modulus for x in root]
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
  #nums = [self.div_polys(root, [-x, 1]) for x in xs]
  nums = [root / polysOverMod([-x, 1]) for x in xs]
  # Generate denominators by evaluating numerator polys at each x
  denoms = [nums[i](xs[i]) for i in range(len(xs))]
  invdenoms = multi_inv(mod, denoms)
  # Generate output polynomial, which is the sum of the per-value numerator
  # polynomials rescaled to have the right y values
  b = [0 for y in ys]
  for i in range(len(xs)):
    #yslice = self.mul(ys[i], invdenoms[i])
    yslice = ys[i] * invdenoms[i]
    num_coefficients = nums[i].coefficients
    for j in range(len(ys)):
      if num_coefficients[j] and ys[i]:
        b[j] += num_coefficients[j] * yslice
  #return [x % self.modulus for x in b]
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
  #eq0 = [-xs[1] % m, 1]
  eq0 = polysOverMod([-xs[1], 1])
  #eq1 = [-xs[0] % m, 1]
  eq1 = polysOverMod([-xs[0], 1])
  #e0 = self.eval_poly_at(eq0, xs[0])
  e0 = eq0(xs[0])
  #e1 = self.eval_poly_at(eq1, xs[1])
  e1 = eq1(xs[1])
  #invall = self.inv(e0 * e1)
  invall = 1/(e0 * e1)
  inv_y0 = ys[0] * invall * e1
  inv_y1 = ys[1] * invall * e0
  #return [(eq0[i] * inv_y0 + eq1[i] * inv_y1) % m for i in range(2)]
  #return [(eq0.coefficients[i] * inv_y0 + eq1.coefficients[i] * inv_y1) for i in range(2)]
  return polysOverMod([(eq0.coefficients[i] * inv_y0 + eq1.coefficients[i] * inv_y1) for i in range(2)])

#def multi_interp_4(self, xsets, ysets):
def multi_interp_4(modulus, xsets, ysets):
  """Optimized version of the above restricted to deg-4 polynomials"""
  mod = IntegersModP(modulus)
  polysOverMod = polynomials_over(mod).factory
  data = []
  invtargets = []
  for xs, ys in zip(xsets, ysets):
    x01, x02, x03, x12, x13, x23 = \
        xs[0] * xs[1], xs[0] * xs[2], xs[0] * xs[3], xs[1] * xs[2], xs[1] * xs[3], xs[2] * xs[3]
    #m = self.modulus
    m = modulus
    #eq0 = [-x12 * xs[3] % m, (x12 + x13 + x23), -xs[1] - xs[2] - xs[3], 1]
    eq0 = polysOverMod([-x12 * xs[3], (x12 + x13 + x23), -xs[1] - xs[2] - xs[3], 1])
    #eq1 = [-x02 * xs[3] % m, (x02 + x03 + x23), -xs[0] - xs[2] - xs[3], 1]
    eq1 = polysOverMod([-x02 * xs[3], (x02 + x03 + x23), -xs[0] - xs[2] - xs[3], 1])
    #eq2 = [-x01 * xs[3] % m, (x01 + x03 + x13), -xs[0] - xs[1] - xs[3], 1]
    eq2 = polysOverMod([-x01 * xs[3], (x01 + x03 + x13), -xs[0] - xs[1] - xs[3], 1])
    #eq3 = [-x01 * xs[2] % m, (x01 + x02 + x12), -xs[0] - xs[1] - xs[2], 1]
    eq3 = polysOverMod([-x01 * xs[2], (x01 + x02 + x12), -xs[0] - xs[1] - xs[2], 1])
    #e0 = self.eval_quartic(eq0, xs[0])
    e0 = eq0(xs[0])
    #e1 = self.eval_quartic(eq1, xs[1])
    e1 = eq1(xs[1])
    #e2 = self.eval_quartic(eq2, xs[2])
    e2 = eq2(xs[2])
    #e3 = self.eval_quartic(eq3, xs[3])
    e3 = eq3(xs[3])
    data.append([ys, eq0, eq1, eq2, eq3])
    invtargets.extend([e0, e1, e2, e3])
  #invalls = self.multi_inv(invtargets)
  invalls = multi_inv(mod, invtargets)
  o = []
  for (i, (ys, eq0, eq1, eq2, eq3)) in enumerate(data):
    invallz = invalls[i * 4:i * 4 + 4]
    #inv_y0 = ys[0] * invallz[0] % m
    inv_y0 = ys[0] * invallz[0]
    #inv_y1 = ys[1] * invallz[1] % m
    inv_y1 = ys[1] * invallz[1]
    #inv_y2 = ys[2] * invallz[2] % m
    inv_y2 = ys[2] * invallz[2]
    #inv_y3 = ys[3] * invallz[3] % m
    inv_y3 = ys[3] * invallz[3]
    o.append(polysOverMod([(eq0.coefficients[i] * inv_y0 + eq1.coefficients[i] * inv_y1 + eq2.coefficients[i] * inv_y2 +
               #eq3[i] * inv_y3) % m for i in range(4)])
               eq3.coefficients[i] * inv_y3) for i in range(4)]))
  # assert o == [self.lagrange_interp_4(xs, ys) for xs, ys in zip(xsets, ysets)]
  return o
