import sympy as sp
from geocusp.residue import blowup_intersection_yw, strict_cusp_blowup_zero_residue, chart_transition_source_to_U12
from geocusp.negative_geometry import symbolic_two_loop_integrand, collective_scaling


def test_exceptional_cross_ratio():
    res, chi = blowup_intersection_yw()
    x,z0,rho,eta = sp.symbols('x z_0 rho eta', nonzero=True)
    assert sp.simplify(chi-eta*z0/x)==0
    assert res!=0


def test_strict_cusp_not_simple_radial_divisor():
    pull, residue = strict_cusp_blowup_zero_residue()
    assert sp.simplify(residue)==0


def test_chart_transition_jacobian_nonzero():
    mapping, inverse, J=chart_transition_source_to_U12()
    assert sp.simplify(J)!=0


def test_two_loop_integrand_exists():
    _, I, X, n=symbolic_two_loop_integrand()
    assert I!=0 and X!=0 and n!=0


def test_collective_tree_scaling_log():
    for L in range(1,6):
        r=collective_scaling(L,[(i,i+1) for i in range(L-1)])
        assert r['scalar_integrand_power']==-2*L
        assert r['radial_form_power']==-1
