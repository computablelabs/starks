from starks.merkle_tree import blake
import time

# TODO(rbharath): Wait, does Vitalik's blog post claim that
# the verifier complexity is linear; Nah looks like t*log(t)
# is optimal. Verifier complexity is O(log**2(t)) which should
# be pretty small even for very large computations.
# NOTE(rbharath): These starks here are not zero-knowledge I
# think. Will need to be added onto library later.
def mimc(inp, steps, round_constants):
  """Compute a MIMC permutation for some number of steps"""
  modulus = 2**256 - 2**32 * 351 + 1
  start_time = time.time()
  for i in range(steps - 1):
    inp = (inp**3 + round_constants[i % len(round_constants)]) % modulus
  print("MIMC computed in %.4f sec" % (time.time() - start_time))
  return inp


def get_power_cycle(r, modulus):
  """
  Get the set of powers of R, until but not including when the
  powers loop back to 1
  """
  o = [1, r]
  while o[-1] != 1:
    o.append((o[-1] * r) % modulus)
  return o[:-1]


def get_pseudorandom_indices(seed, modulus, count, exclude_multiples_of=0):
  """Extract pseudorandom indices from entropy"""
  assert modulus < 2**24
  data = seed
  while len(data) < 4 * count:
    data += blake(data[-32:])
  if exclude_multiples_of == 0:
    return [
        int.from_bytes(data[i:i + 4], 'big') % modulus
        for i in range(0, count * 4, 4)
    ]
  else:
    real_modulus = modulus * (exclude_multiples_of - 1) // exclude_multiples_of
    o = [
        int.from_bytes(data[i:i + 4], 'big') % real_modulus
        for i in range(0, count * 4, 4)
    ]
    return [x + 1 + x // (exclude_multiples_of - 1) for x in o]


def is_a_power_of_2(x):
  return True if x == 1 else False if x % 2 else is_a_power_of_2(x // 2)
