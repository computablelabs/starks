"""Test the Algebraic linking IOP (ALI) protocol."""

from starks.air import Computation
from starks.apr import APR
from starks.ali import ALIProver
from starks.finitefield import FiniteField 
from starks.poly_utils import multivariates_over

#class TestALI(unittest.TestCase):
#  """Basic tests for ALI classes."""
#
#  def test_ali_prover_constructor(construction):
#    """Tests that ALIProver can be initialized."""
#    width = 2
#    # Set the field small in tests since primitive polynomial generation is slow.
#    p = 7 
#    m = 4
#    steps = 512
#    extension_factor = 8
#    field = FiniteField(p, m)
#    inp = [field(0), field(1)]
#    polysOver = multivariates_over(field, width).factory
#    X_1 = polysOver({(1,0): field(1)})
#    X_2 = polysOver({(0,1): field(1)})
#    step_polys = [X_2, X_1 + X_2] 
#    comp = Computation(field, width, inp, steps, step_polys,
#        extension_factor)
#    apr = APR(comp)
#    ali_prover = ALIProver(apr)
