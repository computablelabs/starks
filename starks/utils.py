import time
from starks.merkle_tree import blake
from starks.modp import IntegersModP
from typing import Dict
from typing import List
from starks.numbertype import Field
from starks.numbertype import FieldElement
from starks.poly_utils import multivariates_over


def plus_one(num: int) -> int:
    return num + 1

# TODO(rbharath): Wait, does Vitalik's blog post claim that
# the verifier complexity is linear; Nah looks like t*log(t)
# is optimal. Verifier complexity is O(log**2(t)) which should
# be pretty small even for very large computations.
# NOTE(rbharath): These starks here are not zero-knowledge I
# think. Will need to be added onto library later.
def mimc(inp: int, steps: int, round_constants: List[int]):
  """Compute a MIMC permutation for some number of steps"""
  modulus = 2**256 - 2**32 * 351 + 1
  start_time = time.time()
  for i in range(steps - 1):
    inp = (inp**3 + round_constants[i % len(round_constants)]) % modulus
  print("MIMC computed in %.4f sec" % (time.time() - start_time))
  return inp


def get_power_cycle(r: FieldElement, field: Field):
  """
  Get the set of powers of R, until but not including when the
  powers loop back to 1
  """
  o = [field(1), r]
  while o[-1] != field(1):
    o.append(o[-1] * r)
  return o[:-1]

def generate_Xi_s(field: Field, width: int):
  """
  Constructs polynomials X_1,..,X_width

  This helper method generates "index" polynomials. Using
  these polynomials makes it easier to write more complex
  polynomials programmatically.
  """
  polysOver = multivariates_over(field, width).factory
  Xi_s = []
  for i in range(width):
    pre = (0,) * i
    post = (0,) * (width - (i+1))
    index = pre + (1,) + post
    X_i = polysOver({index: field(1)})
    Xi_s.append(X_i)
  return Xi_s


# TODO(rbharath): This function needs to be reworkd to be more general
def get_pseudorandom_indices(entropy, modulus, count, exclude_multiples_of=0):
#def get_pseudorandom_field_elements(entropy, field, count, exclude_multiples_of=0):
  """Extract pseudorandom indices from entropy

  Draws pseudorandom numbers from a given range while avoiding
  a subset (multiples of a forbidden value). For example, to
  sample random indices from a list of length 512 while
  avoiding indices that are multiples of 32.
  """
  assert modulus < 2**24
  data = entropy 
  # Note that we must have len(data) >= 4 * count. This code #
  # expands data to have necessary length. Think of this as an
  # entropy expansion step.
  while len(data) < 4 * count:
    data += blake(data[-32:])
  if exclude_multiples_of == 0:
    return [
        int.from_bytes(data[i:i + 4], 'big') % modulus
        #field(data[i:i + 4])
        for i in range(0, count * 4, 4)
    ]
  else:
    # TODO(rbharath): This is horribly ugly. Figure out how to generalize this...
    real_modulus = modulus * (exclude_multiples_of - 1) // exclude_multiples_of
    o = [
        int.from_bytes(data[i:i + 4], 'big') % real_modulus
        #field(data[i:i + 4])
        for i in range(0, count * 4, 4)
    ]
    return [x + 1 + x // (exclude_multiples_of - 1) for x in o]


def is_a_power_of_2(x):
  return True if x == 1 else False if x % 2 else is_a_power_of_2(x // 2)
