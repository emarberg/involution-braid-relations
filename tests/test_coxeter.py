import pytest
import numpy as np

from project.utils import (
    InvalidInputException,
    OperatorException
)

from project.algebra import (
    QuadraticNumber,
    RationalNumber,
    Polynomial
)

from project.coxeter import (
    CoxeterGraph,
    CoxeterVector,
    PartialTransform,
    CoxeterTransform,
    CoxeterWord
)


class TestCoxeterGraph:
    def test_constructor(self):
        g = CoxeterGraph([])
        assert g.generators == []

        g = CoxeterGraph([], generators=[3, 2, 1])
        assert g.generators == [1, 2, 3]
        assert g.get_order(1, 2) == g.get_order(2, 1) == 2

        g = CoxeterGraph([(1, 2, 4), (2, 3, 5), (3, 4, 4)], star=[(1, 4), (2, 3)])
        assert g.generators == [1, 2, 3, 4]
        assert g.get_order(2, 1) == 4
        assert g.get_order(3, 2) == 5
        assert g.get_order(4, 3) == 4
        assert g.star(1) == 4 and g.star(2) == 3 and g.star(3) == 2 and g.star(4) == 1

        # inputs i, j to g.get_order() must both be in generators
        e = None
        try:
            g.get_order(2, 10)
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        e = None
        try:
            g.get_order(10, 2)
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        # input i to g.star() must be in generators
        e = None
        try:
            g.star(10)
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

    def test_constructor_edges_exceptions(self):
        """Test handling of errors in input `edges` for CoxeterGraph constructor."""
        # each i and j for (i, j, m) in `edges` must appear in `generators`
        e = None
        try:
            CoxeterGraph(edges=[(1, 2, 1)], generators=[3, 4, 5])
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        # if (i, j, m) in `edges` has m == 1 then must have i == j
        e = None
        try:
            CoxeterGraph(edges=[(1, 2, 1)])
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        # if (i, j, m) in `edges` has i == j then must have m == 1
        e = None
        try:
            CoxeterGraph(edges=[(1, 1, 2)])
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        # if (i, j, m) is in `edges` then m must be an integer or infinity
        e = None
        try:
            CoxeterGraph(edges=[(1, 2, 3.0)])
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        # `edges` cannot contain (i, j, m) and (j, i, n) if m != n
        e = None
        try:
            CoxeterGraph(edges=[(1, 2, 3), (2, 1, 2)])
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

    def test_constructor_star_exceptions(self):
        """Test handling of errors in input `star` for CoxeterGraph constructor."""
        # if (i, j) is in `star` then i and j must be generators
        e = None
        try:
            CoxeterGraph(edges=[(1, 2, 3)], star=[(3, 4)])
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        # `star` cannot contain (i, j) and (j, k) if i != k
        e = None
        try:
            CoxeterGraph(edges=[(1, 2, 3)], star=[(1, 2), (2, 2)])
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        # `star` must define a group automorphism
        e = None
        try:
            CoxeterGraph(edges=[(1, 2, 3), (2, 3, 4)], star=[(1, 3)])
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

    def test_eq(self):
        """Test == operator for CoxeterGraphs, as well as hash() implementation."""
        g = CoxeterGraph(edges=[(1, 2, 3), (2, 1, 3)],
                         generators=[1, 2, 3, 4],
                         star=[(1, 2), (2, 1), (3, 4), (4, 3)])
        h = CoxeterGraph(edges=[(1, 2, 3), (3, 3, 1), (4, 4, 1)],
                         star=[(1, 2), (3, 4)])
        i = CoxeterGraph(edges=[(1, 2, 3), (3, 3, 1), (4, 4, 1)])
        assert g == h and h == g
        assert g != i and i != g
        assert h != i and i != h
        assert hash(g) == hash(h)

    def test_type_a_small_rank(self):
        """Check that constructions of type A0 and A1 are distinct."""
        g = CoxeterGraph.A(0)
        h = CoxeterGraph.A(1)
        assert g.generators == []
        assert h.generators == [1]
        assert g != h

    def test_convenience_methods(self):
        """
        Test methods for creating common finite Coxeter systems and the method
        returning boolean properties of these systems.
        """
        a5 = CoxeterGraph.A(5)
        b6 = CoxeterGraph.B(6)
        d7 = CoxeterGraph.D(7)
        e6 = CoxeterGraph.E(6)
        e7 = CoxeterGraph.E(7)
        e8 = CoxeterGraph.E(8)
        f4 = CoxeterGraph.F(4)
        g2 = CoxeterGraph.G(2)
        h3 = CoxeterGraph.H(3)
        h4 = CoxeterGraph.H(4)
        i = CoxeterGraph([(1, 2, np.infty)])
        j = CoxeterGraph([(1, 2, 7)])

        A5 = CoxeterGraph.A_tilde(5)  # noqa
        B6 = CoxeterGraph.B_tilde(6)  # noqa
        C6 = CoxeterGraph.C_tilde(6)  # noqa
        D6 = CoxeterGraph.D_tilde(6)  # noqa
        E6 = CoxeterGraph.E_tilde(6)  # noqa
        E7 = CoxeterGraph.E_tilde(7)  # noqa
        E8 = CoxeterGraph.E_tilde(8)  # noqa
        F4 = CoxeterGraph.F_tilde(4)  # noqa
        G2 = CoxeterGraph.G_tilde(2)  # noqa

        a5_twist = CoxeterGraph.A_twist(5)
        b2_twist = CoxeterGraph.B_twist(2)
        d7_twist = CoxeterGraph.D_twist(7)
        e6_twist = CoxeterGraph.E_twist(6)
        f4_twist = CoxeterGraph.F_twist(4)
        g2_twist = CoxeterGraph.G_twist(2)

        A5_twist = CoxeterGraph.A_tilde_twist(5)  # noqa
        B6_twist = CoxeterGraph.B_tilde_twist(6)  # noqa
        C6_twist = CoxeterGraph.C_tilde_twist(6)  # noqa
        D6_twist = CoxeterGraph.D_tilde_twist(6)  # noqa
        D6_half_twist = CoxeterGraph.D_tilde_half_twist(6)  # noqa
        D6_small_twist = CoxeterGraph.D_tilde_small_twist(6)  # noqa
        E6_twist = CoxeterGraph.E_tilde_twist(6)  # noqa
        E7_twist = CoxeterGraph.E_tilde_twist(7)  # noqa

        rA5 = CoxeterGraph.A_tilde_rotate(5)  # noqa
        fA5 = CoxeterGraph.A_tilde_flip(5)  # noqa

        simply_laced = [
            a5, a5_twist, A5, A5_twist, d7, d7_twist, e6, e6_twist, e7, e8,
            E6, E6_twist, E7, E7_twist, E8, D6, rA5, fA5, D6_twist, D6_half_twist, D6_small_twist
        ]
        crystallographic = \
            simply_laced + [
                b2_twist, b6, f4, f4_twist, g2, g2_twist, B6, B6_twist, C6, C6_twist, F4, G2
            ]
        quadratic = crystallographic + [h3, h4, i]
        combined = quadratic + [j]
        for g in combined:
            assert g.is_quadratic() == (g in quadratic)
            assert g.is_simply_laced() == (g in simply_laced)
            assert g.is_crystallographic() == (g in crystallographic)

            # test that can convert to str without errors
            str(g)

        assert h4.get_braid(3, 4) == (3, 4, 3, 4, 3)
        assert f4.get_braid(2, 3) == (2, 3, 2, 3)
        e = None
        try:
            i.get_braid(1, 2)
        except Exception as exception:
            e = exception
        assert type(e) == CoxeterGraph.GetBraidException

    @pytest.mark.parametrize("ctor, base_ctor", [
        (CoxeterGraph.AxA_twist, CoxeterGraph.A),
        (CoxeterGraph.BxB_twist, CoxeterGraph.B),
        (CoxeterGraph.DxD_twist, CoxeterGraph.D)
    ])
    def test_product_systems(self, ctor, base_ctor):
        """
        Test that product Coxeter graphs returned by static methods AxA_twist, BxB_twist, and
        DxD_twist are expected disjoint unions of graphs of type A, B, and D.
        """
        k = 5
        n = 2 * k
        g = ctor(n)
        h = base_ctor(k)
        assert g.generators == list(range(1, n + 1))
        assert all(
            g.get_order(i, j) == g.get_order(i + k, j + k) == h.get_order(i, j)
            for i in range(1, k + 1)
            for j in range(1, k + 1)
        )

    def test_eval_bilinear(self):
        a5 = CoxeterGraph.A(5)
        b6 = CoxeterGraph.B(6)
        g2 = CoxeterGraph.G(2)
        h4 = CoxeterGraph.H(4)
        i = CoxeterGraph([(1, 2, np.infty)])
        j = CoxeterGraph([(1, 2, 7)])
        k = CoxeterGraph([(1, 2, 12)])

        assert a5.eval_bilinear(1, 1) == 1
        assert a5.eval_bilinear(1, 3) == 0
        assert a5.eval_bilinear(1, 2) == -RationalNumber(1) / 2
        assert b6.eval_bilinear(5, 6) == -QuadraticNumber.sqrt(2) / 2
        assert h4.eval_bilinear(3, 4) == -(QuadraticNumber.sqrt(5) + 1) / 4
        assert g2.eval_bilinear(1, 2) == -QuadraticNumber.sqrt(3) / 2
        assert k.eval_bilinear(1, 2) == -(QuadraticNumber.sqrt(6) + QuadraticNumber.sqrt(2)) / 4
        assert i.eval_bilinear(1, 2) == -1

        e = None
        try:
            j.eval_bilinear(1, 2)
        except Exception as exception:
            e = exception
        assert type(e) == CoxeterGraph.EvalBilinearException


class TestCoxeterVector:
    def test_constructor(self):
        g = CoxeterGraph.A(5)

        zero = CoxeterVector(g)
        assert zero.coefficients == {}

        one = CoxeterVector(g, 1)
        assert one.coefficients == {1: 1}
        assert 1 in one and 2 not in one

        nonzero = CoxeterVector(g, 1, Polynomial('x'))
        assert nonzero.coefficients == {1: Polynomial('x')}
        assert 1 in one and 2 not in one

        # second argument to root constructor must belong of g.generators
        e = None
        try:
            CoxeterVector(g, 0)
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

    def test_eq(self):
        g = CoxeterGraph.B(5)
        h = CoxeterGraph.A(5)

        # can only compare roots with other roots with same graph, or the int 0
        r = CoxeterVector(g)
        assert r == 0 and 0 == r
        for zero in [RationalNumber(0), QuadraticNumber(0), Polynomial(0), CoxeterVector(h)]:
            e = None
            try:
                r == zero
            except Exception as exception:
                e = exception
            assert type(e) == OperatorException

            e = None
            try:
                zero == r
            except Exception as exception:
                e = exception
            assert type(e) == OperatorException

        # r and s are linearly independent vectors, so clearly not equal
        r = CoxeterVector(g, 1)
        s = CoxeterVector(g, 2)
        assert r != s

        # roots should be considered equal if their coefficients represent the same numbers
        a = CoxeterVector(g, 1, 2)
        b = CoxeterVector(g, 1, RationalNumber(2))
        c = CoxeterVector(g, 1, QuadraticNumber(2))
        d = CoxeterVector(g, 1, Polynomial(2))

        assert a == b == c == d and d == c == b == a
        assert len({a, b, c, d}) == 1

    def test_eval_bilinear(self):
        g = CoxeterGraph.F(4)
        r = CoxeterVector(g, 1)
        s = CoxeterVector(g, 2)
        t = CoxeterVector(g, 3)

        # <alpha, alpha> == 1 for any root alpha
        assert r.eval_bilinear(r) == r.eval_bilinear(1) == 1
        assert s.eval_bilinear(s) == s.eval_bilinear(2) == 1
        assert t.eval_bilinear(t) == t.eval_bilinear(3) == 1

        # <alpha_i, alpha_j> == -cos(pi / m_ij)
        assert r.eval_bilinear(s) == r.eval_bilinear(2) == -RationalNumber(1) / 2
        assert r.eval_bilinear(t) == r.eval_bilinear(3) == 0
        assert s.eval_bilinear(r) == s.eval_bilinear(1) == -RationalNumber(1) / 2
        assert s.eval_bilinear(t) == s.eval_bilinear(3) == -QuadraticNumber.sqrt(2) / 2
        assert t.eval_bilinear(r) == t.eval_bilinear(1) == 0
        assert t.eval_bilinear(s) == t.eval_bilinear(2) == -QuadraticNumber.sqrt(2) / 2

        # check that <-,-> is bilinear
        assert (r + s).eval_bilinear(s + t) == RationalNumber(1) / 2 - QuadraticNumber.sqrt(2) / 2
        assert (r + s).reflect(3) == r + s + QuadraticNumber.sqrt(2) * t

    def test_eval_bilinear_errors(self):
        """Test error handling for various invalid inputs to CoxeterVector.eval_bilinear()."""
        g = CoxeterGraph.F(4)
        r = CoxeterVector(g, 1)

        # int inputs to r.eval_bilinear or r.reflect must belong to r.graph.generators
        e = None
        try:
            r.eval_bilinear(0)
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        # input to r.reflect() must belong to r.graph.generators
        e = None
        try:
            r.reflect(0)
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        h = CoxeterGraph.H(4)
        u = CoxeterVector(h, 1)

        # cannot compute r.eval_bilinear() if input CoxeterVector has different CoxeterGraph
        e = None
        try:
            r.eval_bilinear(u)
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

    def test_set_variables(self):
        """Tests for CoxeterVector.set_variable method."""
        x = Polynomial({0: 1})
        y = Polynomial({1: 1})

        g = CoxeterGraph.F(4)
        a = CoxeterVector(g, 1, x)
        b = CoxeterVector(g, 2, y)
        c = CoxeterVector(g, 3, x + y)

        assert a.set_variable(0, 3) == CoxeterVector(g, 1, 3)
        assert a.set_variable(1, 3) == a
        assert a.set_variable(0, 0) == CoxeterVector(g)

        assert (a + b + c).set_variable(0, 1) == \
            CoxeterVector(g, 1) + b + (1 + y) * CoxeterVector(g, 3)

    def test_add(self):
        g = CoxeterGraph.F(4)
        a = CoxeterVector(g, 1, 1) - CoxeterVector(g, 2, 1)
        b = CoxeterVector(g, 2, 1) - CoxeterVector(g, 3, 1)
        assert (a + b).coefficients == {1: 1, 3: -1}

        assert a + 0 == 0 + a == a
        e = None
        try:
            a + CoxeterVector(CoxeterGraph.H(3), 3)
        except Exception as exception:
            e = exception
        assert type(e) == OperatorException

    def test_power_error(self):
        """Test that ** operator is not implemented for CoxeterVector objects."""
        g = CoxeterGraph.E(8)
        r = CoxeterVector(g, 1) + CoxeterVector(g, 8)
        try:
            r**2
        except NotImplementedError:
            pass

    def test_convenience_methods(self):
        """Tests for methods is_constant, is_positive, is_negative, and is_valid."""
        g = CoxeterGraph.E(8)

        r = CoxeterVector(g, 1, Polynomial(8))
        s = CoxeterVector(g, 2, Polynomial('x') + 1)
        assert r.is_constant()
        assert not s.is_constant()

        t = -r - s
        u = r - s
        v = CoxeterVector(g, 3, Polynomial('x') - 1)
        assert r.is_positive() and not r.is_negative()
        assert s.is_positive() and not r.is_negative()
        assert not t.is_positive() and t.is_negative()
        assert not u.is_positive() and not u.is_negative()
        assert not v.is_positive() and not v.is_negative()

        z = CoxeterVector(g)
        assert z.is_zero()

        # z is not 'valid' since it is zero
        assert not z.is_valid()
        assert r.is_valid()
        assert s.is_valid()
        assert t.is_valid()
        # u is not 'valid' because its coeff of alpha_1 is > 0 while its coeff of alpha_2 is < 0
        assert not u.is_valid()
        assert v.is_valid()

        # v +/- s is 'valid' since its coeff of alpha_3 is x - 1, which is neither < 0 or > 0
        assert (v + s).is_valid() and (v - s).is_valid()

    def test_repr(self):
        g = CoxeterGraph.A(5)
        r = CoxeterVector(g)
        assert str(r) == '0'


class TestPartialTransform:
    def test_constructor(self):
        g = CoxeterGraph.A(5)
        t = PartialTransform(g)
        assert t.graph == g
        assert t.sigma == {}

        sigma = {1: -CoxeterVector(g, 2), 2: CoxeterVector(g, 1)}
        t = PartialTransform(g, sigma)
        assert t.sigma == sigma
        # check that changes to sigma do not affect T after construction
        del sigma[1]
        assert t.sigma == {1: -CoxeterVector(g, 2), 2: CoxeterVector(g, 1)}

        t[1] *= Polynomial('x')
        # t now has no unconditional descents, but 1 is a conditional descent
        assert t.sigma == {1: -CoxeterVector(g, 2, Polynomial('x')), 2: CoxeterVector(g, 1)}

        t[1] += CoxeterVector(g, 2)
        # t now has no strong conditional descents, but 1 is a weak conditional descent
        assert t.sigma == {1: CoxeterVector(g, 2, 1 - Polynomial('x')), 2: CoxeterVector(g, 1)}

        t[1] += CoxeterVector(g, 2, 2 * Polynomial('x'))
        # t now has no descents, although T is not trivial
        assert t.sigma == {1: CoxeterVector(g, 2, 1 + Polynomial('x')), 2: CoxeterVector(g, 1)}

    def test_constructor_errors(self):
        # Check error handling for invalid constructor inputs.
        g = CoxeterGraph.A(5)
        e = None
        try:
            PartialTransform(g, {0: 0})
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        e = None
        try:
            PartialTransform(g, {1: CoxeterVector(CoxeterGraph.B(5), 1)})
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        t = PartialTransform(g, {1: -CoxeterVector(g, 2), 2: CoxeterVector(g, 1)})
        e = None
        try:
            t[1] = CoxeterVector(CoxeterGraph.B(5), 1)
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

    def test_eq(self):
        g = CoxeterGraph.A(5)
        t = PartialTransform(g, {
            1: CoxeterVector(g, 2, 1 + Polynomial('x')),
            2: CoxeterVector(g, 1)
        })
        u = t
        v = t.copy()

        t[1] = CoxeterVector(g)
        assert u[1] == 0 and u == t
        assert v[1] == CoxeterVector(g, 2, 1 + Polynomial('x')) and t != v
        assert len({t, u, v}) == len({t, v}) == 2

    def test_multiply(self):
        g = CoxeterGraph.A(5)
        t = PartialTransform(g, {1: CoxeterVector(g, 2, 1 + Polynomial('x')), 2: CoxeterVector(g, 1)})

        u = t * 1
        assert u.sigma == {
            1: -CoxeterVector(g, 2, 1 + Polynomial('x')),
            2: CoxeterVector(g, 1) + CoxeterVector(g, 2, 1 + Polynomial('x'))
        }

        v = 1 * t
        assert v.sigma == {
            1: (CoxeterVector(g, 1) + CoxeterVector(g, 2)) * (1 + Polynomial('x')),
            2: -CoxeterVector(g, 1)
        }

        w = 1 * t * 1
        assert w.sigma == {
            1: -(CoxeterVector(g, 1) + CoxeterVector(g, 2)) * (1 + Polynomial('x')),
            2: CoxeterVector(g, 1, Polynomial('x')) + CoxeterVector(g, 2, 1 + Polynomial('x'))
        }

    def test_multiplication_errors(self):
        g = CoxeterGraph.A(5)
        t = PartialTransform(g, {
            1: CoxeterVector(g, 2, 1 + Polynomial('x')),
            2: CoxeterVector(g, 1)
        })

        # cannot multiply by 6 since this is not a generator for the Coxeter system A5
        e = None
        try:
            t * 6
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        e = None
        try:
            6 * t
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        # cannot multiply by 3 since (alpha_2).reflect(3) = alpha_2 + alpha_3 and T[3] is undefined
        e = None
        try:
            t * 3
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

    def test_convenience_methods(self):
        """Tests for various convenience methods of RooTransform objects."""
        g = CoxeterGraph.A(5)
        sigma = {1: CoxeterVector(g, 2, 1 + Polynomial('x')), 2: CoxeterVector(g, 1)}
        t = PartialTransform(g, sigma)

        assert not t.is_constant()
        assert not t.is_complete()
        assert t.is_positive()
        assert not t.is_identity()

        assert len(t) == 2
        assert 1 in t and 2 in t and 3 not in t
        assert {i for i in t} == {1, 2}
        assert set(t.values()) == set(sigma.values())

        u = PartialTransform.identity(g)
        assert u.sigma == {i: CoxeterVector(g, i) for i in g.generators}
        assert u.is_constant()
        assert u.is_complete()
        assert u.is_positive()
        assert u.is_identity()

        # try converting PartialTransforms to str, check that no errors occur
        str(t)
        str(u)

    def test_determinant(self):
        X = Polynomial('X')  # noqa
        Y = Polynomial('Y')  # noqa
        g = CoxeterGraph.A(2)
        t = PartialTransform(g, {1: CoxeterVector(g, 2, 1 + X)})
        assert t.determinant() is None

        t = PartialTransform(g, {1: CoxeterVector(g, 2, 1 + X), 2: CoxeterVector(g, 2, 1 + X**2)})
        assert t.determinant() is None

        t = PartialTransform(g, {1: CoxeterVector(g, 2, 1 + X), 2: CoxeterVector(g, 2, 1 + Y)})
        assert t.determinant() is None

        t = PartialTransform(g, {1: CoxeterVector(g, 2), 2: CoxeterVector(g, 1, -5)})
        assert t.determinant() == 5

        t = PartialTransform(g, {1: CoxeterVector(g, 2, Y), 2: CoxeterVector(g, 1, Y - 5)})
        assert t.determinant() == 5 * Y - Y**2

    def test_get_relations(self):
        g = CoxeterGraph.B(3)
        t = PartialTransform(g, {2: CoxeterVector(g, 2), 3: CoxeterVector(g, 3)})
        assert t.get_relations() == {((2, 3, 2), (3, 2, 3))}

        t = PartialTransform(g, {2: CoxeterVector(g, 3), 3: CoxeterVector(g, 2)})
        assert t.get_relations() == {((2, 3), (3, 2))}

        t = PartialTransform(g, {3: CoxeterVector(g, 2), 2: CoxeterVector(g, 1)})
        assert t.get_relations() == set()


class TestCoxeterTransform:
    def test_constructor(self):
        g = CoxeterGraph.B(3)
        t = CoxeterTransform(g)
        assert t.graph == g
        assert t.sigma == {i: CoxeterVector(g, i) for i in g.generators}
        assert t.right_descents == set()
        assert t.minimal_right_descent is None
        assert t.minimal_reduced_word == ()

        sigma = {i: -CoxeterVector(g, i) for i in g.generators}
        t = CoxeterTransform(g, sigma)
        assert t.graph == g
        assert t.sigma == sigma
        sigma.clear()
        assert t.sigma == {i: -CoxeterVector(g, i) for i in g.generators}
        assert t.right_descents == {1, 2, 3}
        assert t.minimal_right_descent is 1
        assert t.minimal_reduced_word == (3, 2, 3, 1, 2, 3, 1, 2, 1)

    def test_constructor_errors(self):
        # Check error handling for invalid constructor inputs.
        g = CoxeterGraph.B(3)
        e = None
        try:
            CoxeterTransform(g, {0: CoxeterVector(g, 0)})
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        h = CoxeterGraph.A(3)
        e = None
        try:
            CoxeterTransform(g, {
                1: CoxeterVector(g, 1),
                2: CoxeterVector(g, 2),
                3: CoxeterVector(h, 3)
            })
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        e = None
        try:
            CoxeterTransform(g, {
                0: CoxeterVector(g, 1),
                1: CoxeterVector(g, 2),
                2: CoxeterVector(g, Polynomial('x'))
            })
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

    def test_get_inverse(self):
        g = CoxeterGraph.B(3)
        t = CoxeterTransform(g, {i: -CoxeterVector(g, i) for i in g.generators})
        assert t == t.get_inverse()

        u = CoxeterTransform(g)
        assert (u * 3 * 2).get_inverse() == u * 2 * 3

    def test_eq(self):
        g = CoxeterGraph.B(2)
        t = PartialTransform(g, {1: CoxeterVector(g, 1), 2: CoxeterVector(g, 2)})
        u = CoxeterTransform(g, {1: CoxeterVector(g, 1), 2: CoxeterVector(g, 2)})
        assert t == u

        v = u.copy()
        v[1] = -CoxeterVector(g, 1)
        v[2] = -CoxeterVector(g, 2)

        assert t == u and u != v
        assert not hasattr(t, 'minimal_reduced_word')
        assert u.minimal_reduced_word == ()
        assert v.minimal_reduced_word == (2, 1, 2, 1)
        assert len({t, u, v}) == 2

    def test_setitem_errors(self):
        g = CoxeterGraph.B(2)
        v = CoxeterTransform(g, {1: CoxeterVector(g, 1), 2: CoxeterVector(g, 2)})

        # cannot assign value which is non-constant CoxeterVector
        e = None
        try:
            v[1] = -CoxeterVector(g, 1, Polynomial('x'))
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        # cannot assign to index which is not in g.generators
        e = None
        try:
            v[0] = CoxeterVector(g, 1)
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

    def test_span_by_right_relations(self):
        g = CoxeterGraph.A(3)
        relations = {((1, 2), (2, 1)), ((2, 3), (3, 2))}

        a = CoxeterTransform.from_word(g, (1, 2, 3, 1))
        b = CoxeterTransform.from_word(g, (2, 3, 2, 1))
        c = CoxeterTransform.from_word(g, (2, 3, 1, 2))

        assert a.multiply_down((3, 2)) is not None
        assert a.multiply_down((3, 2)).multiply_up((3, 2)) is not None
        assert a.multiply_up((2, 3)) is None
        assert a.span_by_right_relations(relations) == {a, b, c}

    def test_demazure_conjugate(self):
        g = CoxeterGraph.A(5)
        e = CoxeterTransform(g)

        e = e.demazure_conjugate(1)
        assert e.minimal_reduced_word == (1,)

        e = e.demazure_conjugate(2)
        assert e.minimal_reduced_word == (1, 2, 1)

        g = CoxeterGraph.A_twist(5)
        e = CoxeterTransform(g)

        e = e.demazure_conjugate(1)
        assert e.minimal_reduced_word == (5, 1)

        exc = None
        try:
            e.demazure_conjugate(0)
        except Exception as exception:
            exc = exception
        assert type(exc) == InvalidInputException


class TestCoxeterWord:
    def test_init(self):
        g = CoxeterGraph.B(3)
        w = CoxeterWord(g)
        assert w.word == ()
        assert w.left_action == CoxeterTransform.identity(g)
        assert w.right_action == CoxeterTransform.identity(g)
        assert w.is_reduced
        assert w.left_descents == set()
        assert w.right_descents == set()
        assert w.minimal_reduced_word == ()

        w = CoxeterWord(g, (2, 3, 2, 3))
        assert w.word == (2, 3, 2, 3)
        assert w.left_action == CoxeterTransform.from_word(g, (2, 3, 2, 3))
        assert w.right_action == w.left_action.get_inverse()
        assert w.is_reduced
        assert w.left_descents == {2, 3}
        assert w.right_descents == {2, 3}
        assert w.minimal_reduced_word == (2, 3, 2, 3)

        v = CoxeterWord(g, (3, 2, 3, 2))
        assert v.minimal_reduced_word == (2, 3, 2, 3)
        assert v != w and v.left_action == w.left_action and v.right_action == v.right_action
        assert hash(v) != hash(w)

        v = CoxeterWord(g, (3, 2, 3, 2, 2, 2))
        assert not v.is_reduced
        assert v != w and v.left_action == w.left_action and v.right_action == v.right_action
        assert hash(v) != hash(w)

        assert len(w) == 4 and len(v) == 6
        assert str(v) == 'unreduced word (3, 2, 3, 2, 2, 2)'
        assert str(w) == 'reduced word (2, 3, 2, 3)'

    @pytest.mark.parametrize("word, involution_word", [
        ((2, 3), (3, 2, 3)),
        ((3, 1), (3, 1)),
        ((2, 3, 2, 3), (3, 2, 3, 2, 3))
    ])
    def test_to_involution(self, word, involution_word):
        g = CoxeterGraph.B(3)
        w = CoxeterWord(g, word)
        assert w.to_involution().word == involution_word

    def test_get_reduced_words(self):
        g = CoxeterGraph.B(3)
        w = CoxeterWord(g, (3, 2, 3, 2))
        assert w.get_reduced_words() == {(3, 2, 3, 2), (2, 3, 2, 3)}

    def test_copy(self):
        g = CoxeterGraph.B(3)
        w = CoxeterWord(g)
        v = w.copy()
        assert v == w

        w.extend_left(1)
        assert w.word == (1,) and v.word == ()

    def test_multiply(self):
        g = CoxeterGraph.B(3)
        assert CoxeterWord(g) * CoxeterWord(g) == CoxeterWord(g)

        w = ((CoxeterWord(g) * 3) * 2) * 1
        assert w.word == (3, 2, 1)
        assert (w * w).word == (3, 2, 1, 3, 2, 1)

        w = 3 * (2 * (1 * CoxeterWord(g)))
        assert w.word == (3, 2, 1)
        assert (w * w).word == (3, 2, 1, 3, 2, 1)
        assert w.__rmul__(w).word == (3, 2, 1, 3, 2, 1)

        w = 3 * (2 * (1 * ((CoxeterWord(g) * 3) * 2))) * 1
        assert w.word == (3, 2, 1, 3, 2, 1)

    def test_multiply_errors(self):
        g = CoxeterGraph.B(3)
        h = CoxeterGraph.A(3)

        # cannot multiply by non-generator 0, as g.generators = [1, 2, 3]
        e = None
        try:
            CoxeterWord(g) * 0
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        e = None
        try:
            0 * CoxeterWord(g)
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        e = None
        try:
            CoxeterWord(g).extend_right(0)
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        e = None
        try:
            CoxeterWord(g).extend_left(0)
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        # cannot multiply CoxeterWords with different CoxeterGraphs
        e = None
        try:
            CoxeterWord(g).__mul__(CoxeterWord(h))
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

        e = None
        try:
            CoxeterWord(g).__rmul__(CoxeterWord(h))
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException
