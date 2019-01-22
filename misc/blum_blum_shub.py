"""Blum-Blum-Shub is a pseudorandom generator.

x_{n+1} = x_n^2 mod M

M = pq is the modulus and is the product of two large primes.

Note that the security of this construction depends on the security of integer
factoring. That means BBS is *not* quantum-secure.
"""

from starks.primality import probablyPrime
import random


def good_prime(p):
  """A good prime is one which is congruent to 3 mod 4.

  This condition guarantees that each quadratic residue for p has one square
  root which is also a quadratic residue.
  """
  return p % 4 == 3 and probablyPrime(p, accuracy=100)


def find_good_prime(num_bits=512):
  """Randomly sample to find a good prime."""
  candidate = 1

  while not good_prime(candidate):
    candidate = random.getrandbits(num_bits)

  return candidate


def make_modulus(num_bits=512):
  """M = pq is the modulus and is the product of two large primes."""
  return find_good_prime(num_bits) * find_good_prime(num_bits)


def parity(n):
  """Computes the sum of bits in n mod 2."""
  # bin(n) returns 0b.... 
  # bin(n)[2:] trims "ob"
  return sum(int(x) for x in bin(n)[2:]) % 2


def blum_blum_shub(modulus_length=512):
  """Returns a function implementing BBS pseudorandom generator.

  x_{n+1} = x_n^2 mod M
  """
  modulus = make_modulus(num_bits=modulus_length)

  def f(inputInt):
    return pow(inputInt, 2, modulus)

  return f


