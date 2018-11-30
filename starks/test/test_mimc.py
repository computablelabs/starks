import unittest
from starks.mimc_stark import mimc
from starks.mimc_stark import mk_mimc_proof

class TestMiMC(unittest.TestCase):
  """
  Basic tests for MiMC and MiMC-stark implementation. 
  """

  def test_mimc(self):
    """
    Basic tests of MiMC.
    """
    inp = 5
    steps = 3
    round_constants = [2, 7]
    val = mimc(inp, steps, round_constants)

  def test_mimc_stark(self):
    """
    Basic tests of MiMC Stark generation
    """
    inp = 5
    steps = 4
    round_constants = [2, 7]
    proof = mk_mimc_proof(inp, steps, round_constants)
    print(len(proof))
    assert 0 == 1
