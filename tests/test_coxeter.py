import pytest
import numpy as np

from utils import (
    InvalidInputException
)

from algebra import (
    QuadraticNumber,
    RationalNumber,
    Polynomial,
    CoxeterGraph,
    Root,
    RootTransform,
    CoxeterTransform
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

        try:
            g.get_order(2, 10)
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        try:
            g.get_order(10, 2)
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        try:
            g.star(10)
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

    def test_constructor_edges_exceptions(self):
        """Test handling of errors in input `edges` for CoxeterGraph constructor."""
        # each i and j for (i, j, m) in `edges` must appear in `generators`
        try:
            CoxeterGraph(edges=[(1, 2, 1)], generators=[3, 4, 5])
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        # if (i, j, m) in `edges` has m == 1 then must have i == j
        try:
            CoxeterGraph(edges=[(1, 2, 1)])
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        # if (i, j, m) in `edges` has i == j then must have m == 1
        try:
            CoxeterGraph(edges=[(1, 1, 2)])
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        # if (i, j, m) is in `edges` then m must be an integer or infinity
        try:
            CoxeterGraph(edges=[(1, 2, 3.0)])
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        # `edges` cannot contain (i, j, m) and (j, i, n) if m != n
        try:
            CoxeterGraph(edges=[(1, 2, 3), (2, 1, 2)])
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

    def test_constructor_star_exceptions(self):
        """Test handling of errors in input `star` for CoxeterGraph constructor."""
        # if (i, j) is in `star` then i and j must be generators
        try:
            CoxeterGraph(edges=[(1, 2, 3)], star=[(3, 4)])
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        # `star` cannot contain (i, j) and (j, k) if i != k
        try:
            CoxeterGraph(edges=[(1, 2, 3)], star=[(1, 2), (2, 2)])
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        # `star` must define a group automorphism
        try:
            CoxeterGraph(edges=[(1, 2, 3), (2, 3, 4)], star=[(1, 3)])
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

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

    def test_convenience_methods(self):
        """
        Test methods for creating common finite Coxeter systems and the method
        returning boolean properties of these systems.
        """
        a5 = CoxeterGraph.A(5)
        A5 = CoxeterGraph.A_tilde(5)
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

        a5_twist = CoxeterGraph.A2(5)
        d7_twist = CoxeterGraph.D2(7)
        e6_twist = CoxeterGraph.E2(6)
        f4_twist = CoxeterGraph.F2(4)
        g2_twist = CoxeterGraph.G2(2)

        simply_laced = [a5, a5_twist, A5, d7, d7_twist, e6, e6_twist, e7, e8]
        crystallographic = simply_laced + [b6, f4, f4_twist, g2, g2_twist]
        quadratic = crystallographic + [h3, h4, i]
        combined = quadratic + [j]
        for g in combined:
            assert g.is_quadratic() == (g in quadratic)
            assert g.is_simply_laced() == (g in simply_laced)
            assert g.is_crystallographic() == (g in crystallographic)

            # test that can convert to str without errors
            str(g)

        assert h4.get_braid(1, 2) == (1, 2, 1, 2, 1)
        assert f4.get_braid(2, 3) == (2, 3, 2, 3)
        try:
            i.get_braid(1, 2)
        except Exception as e:
            assert type(e) == CoxeterGraph.GetBraidException
        else:
            assert False

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
        assert a5.eval_bilinear(1, 2) == -RationalNumber(1)/2
        assert b6.eval_bilinear(0, 1) == -QuadraticNumber.sqrt(2)/2
        assert h4.eval_bilinear(1, 2) == -(QuadraticNumber.sqrt(5) + 1)/4
        assert g2.eval_bilinear(1, 2) == -QuadraticNumber.sqrt(3)/2
        assert k.eval_bilinear(1, 2) == -(QuadraticNumber.sqrt(6) + QuadraticNumber.sqrt(2))/4
        assert i.eval_bilinear(1, 2) == -1

        try:
            j.eval_bilinear(1, 2)
        except Exception as e:
            assert type(e) == CoxeterGraph.EvalBilinearException
        else:
            assert False


class TestRoot:
    def test_constructor(self):
        g = CoxeterGraph.A(5)

        zero = Root(g)
        assert zero.coefficients == {}

        one = Root(g, 1)
        assert one.coefficients == {1: 1}

        nonzero = Root(g, 1, Polynomial('x'))
        assert nonzero.coefficients == {1: Polynomial('x')}

        try:
            Root(g, 0)
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

    def test_eq(self):
        g = CoxeterGraph.B(5)
        h = CoxeterGraph.A(5)

        # can only compare roots with other roots with same graph, or the int 0
        r = Root(g)
        assert r == 0 and 0 == r
        for zero in [RationalNumber(0), QuadraticNumber(0), Polynomial(0), Root(h)]:
            try:
                r == zero
            except Exception as e:
                assert type(e) == Root.OperatorException
            else:
                assert False

            try:
                zero == r
            except Exception as e:
                assert type(e) == type(zero).OperatorException
            else:
                assert False

        r = Root(g, 1)
        s = Root(g, 2)
        assert r != s

        a = Root(g, 1, 2)
        b = Root(g, 1, RationalNumber(2))
        c = Root(g, 1, QuadraticNumber(2))
        d = Root(g, 1, Polynomial(2))

        assert a == b == c == d and d == c == b == a
        assert len({a, b, c, d}) == 1

    def test_eval_bilinear(self):
        g = CoxeterGraph.F(4)
        r = Root(g, 1)
        s = Root(g, 2)
        t = Root(g, 3)

        assert r.eval_bilinear(r) == r.eval_bilinear(1) == 1
        assert s.eval_bilinear(s) == s.eval_bilinear(2) == 1
        assert t.eval_bilinear(t) == t.eval_bilinear(3) == 1

        assert r.eval_bilinear(s) == r.eval_bilinear(2) == -RationalNumber(1)/2
        assert r.eval_bilinear(t) == r.eval_bilinear(3) == 0
        assert s.eval_bilinear(r) == s.eval_bilinear(1) == -RationalNumber(1)/2
        assert s.eval_bilinear(t) == s.eval_bilinear(3) == -QuadraticNumber.sqrt(2)/2
        assert t.eval_bilinear(r) == t.eval_bilinear(1) == 0
        assert t.eval_bilinear(s) == t.eval_bilinear(2) == -QuadraticNumber.sqrt(2)/2

        assert (r + s).eval_bilinear(s + t) == RationalNumber(1)/2 - QuadraticNumber.sqrt(2)/2
        assert (r + s).reflect(3) == r + s + QuadraticNumber.sqrt(2)*t

        # int inputs to r.eval_bilinear or r.reflect must belong to r.graph.generators
        try:
            r.eval_bilinear(0)
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        try:
            r.reflect(0)
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        h = CoxeterGraph.H(4)
        u = Root(h, 1)

        # cannot compute r.eval_bilinear() if input Root has different CoxeterGraph
        try:
            r.eval_bilinear(u)
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

    def test_set_variables(self):
        """Tests for Root.set_variable and Root.set_variables_to_zero methods."""
        x = Polynomial({0: 1})
        y = Polynomial({1: 1})

        g = CoxeterGraph.F(4)
        a = Root(g, 1, x)
        b = Root(g, 2, y)
        c = Root(g, 3, x + y)

        assert a.set_variable(0, 3) == Root(g, 1, 3)
        assert a.set_variable(1, 3) == a
        assert a.set_variable(0, 0) == a.set_variables_to_zero({0}) == Root(g)

        assert (a + b + c).set_variable(0, 1) == Root(g, 1) + b + (1 + y)*Root(g, 3)
        assert (a + b + c).set_variables_to_zero({0}) == (Root(g, 2) + Root(g, 3))*y
        assert (a + b + c).set_variables_to_zero({1}) == (Root(g, 1) + Root(g, 3))*x
        assert (a + b + c).set_variables_to_zero({0, 1}) == 0

    def test_add(self):
        g = CoxeterGraph.F(4)
        a = Root(g, 1, 1) - Root(g, 2, 1)
        b = Root(g, 2, 1) - Root(g, 3, 1)
        assert (a + b).coefficients == {1: 1, 3: -1}

        assert a + 0 == 0 + a == a
        try:
            a + Root(CoxeterGraph.H(3), 3)
        except Exception as e:
            assert type(e) == Root.OperatorException
        else:
            assert False

    def test_power_error(self):
        """Test that ** operator is not implemented for Root objects."""
        g = CoxeterGraph.E(8)
        r = Root(g, 1) + Root(g, 8)
        try:
            r**2
        except NotImplementedError:
            pass

    def test_convenience_methods(self):
        """Tests for methods is_constant, is_positive, is_negative, and is_valid."""
        g = CoxeterGraph.E(8)

        r = Root(g, 1, Polynomial(8))
        s = Root(g, 2, Polynomial('x') + 1)
        assert r.is_constant()
        assert not s.is_constant()

        t = -r - s
        u = r - s
        v = Root(g, 3, Polynomial('x') - 1)
        assert r.is_positive() and not r.is_negative()
        assert s.is_positive() and not r.is_negative()
        assert not t.is_positive() and t.is_negative()
        assert not u.is_positive() and not u.is_negative()
        assert not v.is_positive() and not v.is_negative()

        z = Root(g)
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


class TestRootTransform:
    pass


class TestCoxeterTransform:
    pass


class TestCoxeterWord:
    pass
