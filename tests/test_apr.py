"""Test the algebraic placement and routing (APR) transformation."""
import pytest
from starks.finitefield import IntegersModP
from starks.air import AIR 
from starks.apr import construct_neighbors
from starks.apr import APR
from starks.utils import generate_Xi_s
from starks.polynomial import polynomials_over
from starks.poly_utils import multivariates_over
from starks.finitefield import FiniteField 
from starks.poly_utils import generate_primitive_polynomial

def test_construct_neighbors():
    """Test that neighbors are constructed correctly."""
    p = 2
    m = 17
    width = 2
    Tau = list(range(width))
    zeta = generate_primitive_polynomial(p, width)
    Zp = IntegersModP(p)
    basePolys = polynomials_over(Zp)
    ##x^17 + x^3 + 1 is primitive 
    coefficients = [Zp(0)] * 18
    coefficients[0] = Zp(1)
    coefficients[3] = Zp(1)
    coefficients[17] = Zp(1)
    poly = basePolys(coefficients)
    field = FiniteField(p, m, polynomialModulus=poly)
    # g
    g = basePolys([0, 1])
    polysOver = polynomials_over(field).factory
    Nbrs = construct_neighbors(Tau, zeta, g, polysOver)
    assert len(Nbrs) == 3 * len(Tau)

#def test_tilde_expansion():
#    """"Test that that tilde expansion is performed correctly."""
#    assert 0 == 1

#def test_apr_constructor():
#    """Test that the APR class can be initialized."""
#    width = 2
#    # Set the field small in tests since primitive polynomial generation is slow.
#    p = 2
#    m = 17
#    Zp = IntegersModP(p)
#    polysOver = polynomials_over(Zp)
#    #field = FiniteField(p, m)
#    ##x^17 + x^3 + 1 is primitive 
#    coefficients = [Zp(0)] * 18
#    coefficients[0] = Zp(1)
#    coefficients[3] = Zp(1)
#    coefficients[17] = Zp(1)
#    poly = polysOver(coefficients)
#    field = FiniteField(p, m, polynomialModulus=poly)
#    t = 3
#    inp = [field(0), field(1)]
#    polysOver = multivariates_over(field, width).factory
#    [X_1, X_2] = generate_Xi_s(field, width)
#    step_polys = [X_2, X_1 + X_2] 
#    air = AIR(field, width, inp, t, step_polys)
#    apr = APR(air)

#
#  def test_get_witness(self):
#    """Test that a witness can be extracted from the APR instance.
#    
#    # TODO(rbharath): Make this test non-trivial to check the witness has
#    # required properties.
#    """
#    width = 2
#    # Set the field small in tests since primitive polynomial generation is slow.
#    p = 2
#    m = 17
#    #m = 8
#    Zp = IntegersModP(p)
#    polysOver = polynomials_over(Zp)
#    #x^17 + x^3 + 1 is primitive 
#    coefficients = [Zp(0)] * 18
#    coefficients[0] = Zp(1)
#    coefficients[3] = Zp(1)
#    coefficients[17] = Zp(1)
#    poly = polysOver(coefficients)
#    field = FiniteField(p, m, polynomialModulus=poly)
#    steps = 4
#    extension_factor = 8
#    inp = [field(0), field(1)]
#    polysOver = multivariates_over(field, width).factory
#    X_1 = polysOver({(1,0): field(1)})
#    X_2 = polysOver({(0,1): field(1)})
#    step_polys = [X_2, X_1 + X_2] 
#    air = AIR(field, width, inp, steps, step_polys,
#        extension_factor)
#    apr = APR(air)
#    # TODO(rbharath): Uncomment this and reactivate it
#    witness = apr.generate_witness()
