import pytest
import numpy as np

from project.utils import (
    InvalidInputException,
    OperatorException
)

from project.algebra import (
    QuadraticNumber,
    RationalNumber,
    Polynomial,
    CoxeterGraph,
    Root,
    RootTransform,
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

        # input i to g.star() must be in generators
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

        A5 = CoxeterGraph.A_tilde(5)
        B6 = CoxeterGraph.B_tilde(6)
        C6 = CoxeterGraph.C_tilde(6)
        D6 = CoxeterGraph.D_tilde(6)
        E6 = CoxeterGraph.E_tilde(6)
        E7 = CoxeterGraph.E_tilde(7)
        E8 = CoxeterGraph.E_tilde(8)
        F4 = CoxeterGraph.F_tilde(4)
        G2 = CoxeterGraph.G_tilde(2)

        a5_twist = CoxeterGraph.A2(5)
        b2_twist = CoxeterGraph.B2(2)
        d7_twist = CoxeterGraph.D2(7)
        e6_twist = CoxeterGraph.E2(6)
        f4_twist = CoxeterGraph.F2(4)
        g2_twist = CoxeterGraph.G2(2)

        simply_laced = [a5, a5_twist, A5, d7, d7_twist, e6, e6_twist, e7, e8, E6, E7, E8, D6]
        crystallographic = \
            simply_laced + [b2_twist, b6, f4, f4_twist, g2, g2_twist, B6, C6, F4, G2]
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
        assert 1 in one and 2 not in one

        nonzero = Root(g, 1, Polynomial('x'))
        assert nonzero.coefficients == {1: Polynomial('x')}
        assert 1 in one and 2 not in one

        # second argument to root constructor must belong of g.generators
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
                assert type(e) == OperatorException
            else:
                assert False

            try:
                zero == r
            except Exception as e:
                assert type(e) == OperatorException
            else:
                assert False

        # r and s are linearly independent vectors, so clearly not equal
        r = Root(g, 1)
        s = Root(g, 2)
        assert r != s

        # roots should be considered equal if their coefficients represent the same numbers
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

        # <alpha, alpha> == 1 for any root alpha
        assert r.eval_bilinear(r) == r.eval_bilinear(1) == 1
        assert s.eval_bilinear(s) == s.eval_bilinear(2) == 1
        assert t.eval_bilinear(t) == t.eval_bilinear(3) == 1

        # <alpha_i, alpha_j> == -cos(pi / m_ij)
        assert r.eval_bilinear(s) == r.eval_bilinear(2) == -RationalNumber(1)/2
        assert r.eval_bilinear(t) == r.eval_bilinear(3) == 0
        assert s.eval_bilinear(r) == s.eval_bilinear(1) == -RationalNumber(1)/2
        assert s.eval_bilinear(t) == s.eval_bilinear(3) == -QuadraticNumber.sqrt(2)/2
        assert t.eval_bilinear(r) == t.eval_bilinear(1) == 0
        assert t.eval_bilinear(s) == t.eval_bilinear(2) == -QuadraticNumber.sqrt(2)/2

        # check that <-,-> is bilinear
        assert (r + s).eval_bilinear(s + t) == RationalNumber(1)/2 - QuadraticNumber.sqrt(2)/2
        assert (r + s).reflect(3) == r + s + QuadraticNumber.sqrt(2)*t

    def test_eval_bilinear_errors(self):
        """Test error handling for various invalid inputs to Root.eval_bilinear()."""
        g = CoxeterGraph.F(4)
        r = Root(g, 1)

        # int inputs to r.eval_bilinear or r.reflect must belong to r.graph.generators
        try:
            r.eval_bilinear(0)
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        # input to r.reflect() must belong to r.graph.generators
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
        """Tests for Root.set_variable method."""
        x = Polynomial({0: 1})
        y = Polynomial({1: 1})

        g = CoxeterGraph.F(4)
        a = Root(g, 1, x)
        b = Root(g, 2, y)
        c = Root(g, 3, x + y)

        assert a.set_variable(0, 3) == Root(g, 1, 3)
        assert a.set_variable(1, 3) == a
        assert a.set_variable(0, 0) == Root(g)

        assert (a + b + c).set_variable(0, 1) == Root(g, 1) + b + (1 + y)*Root(g, 3)
        # assert (a + b + c).set_variables_to_zero({0}) == (Root(g, 2) + Root(g, 3))*y
        # assert (a + b + c).set_variables_to_zero({1}) == (Root(g, 1) + Root(g, 3))*x
        # assert (a + b + c).set_variables_to_zero({0, 1}) == 0

    def test_add(self):
        g = CoxeterGraph.F(4)
        a = Root(g, 1, 1) - Root(g, 2, 1)
        b = Root(g, 2, 1) - Root(g, 3, 1)
        assert (a + b).coefficients == {1: 1, 3: -1}

        assert a + 0 == 0 + a == a
        try:
            a + Root(CoxeterGraph.H(3), 3)
        except Exception as e:
            assert type(e) == OperatorException
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

    def test_repr(self):
        g = CoxeterGraph.A(5)
        r = Root(g)
        assert str(r) == '0'


class TestRootTransform:
    def test_constructor(self):
        g = CoxeterGraph.A(5)
        T = RootTransform(g)
        assert T.graph == g
        assert T.sigma == {}
        assert T.unconditional_descents == set()
        assert T.strong_conditional_descents == set()
        assert T.weak_conditional_descents == set()

        sigma = {1: -Root(g, 2), 2: Root(g, 1)}
        T = RootTransform(g, sigma)
        assert T.sigma == sigma
        # check that changes to sigma do not affect T after construction
        del sigma[1]
        assert T.sigma == {1: -Root(g, 2), 2: Root(g, 1)}
        assert T.unconditional_descents == {1}
        assert T.strong_conditional_descents == {1}
        assert T.weak_conditional_descents == {1}

        T[1] *= Polynomial('x')
        # T now has no unconditional descents, but 1 is a conditional descent
        assert T.sigma == {1: -Root(g, 2, Polynomial('x')), 2: Root(g, 1)}
        assert T.unconditional_descents == set()
        assert T.strong_conditional_descents == {1}
        assert T.weak_conditional_descents == {1}

        T[1] += Root(g, 2)
        # T now has no strong conditional descents, but 1 is a weak conditional descent
        assert T.sigma == {1: Root(g, 2, 1 - Polynomial('x')), 2: Root(g, 1)}
        assert T.unconditional_descents == set()
        assert T.strong_conditional_descents == set()
        assert T.weak_conditional_descents == {1}

        T[1] += Root(g, 2, 2*Polynomial('x'))
        # T now has no descents, although T is not trivial
        assert T.sigma == {1: Root(g, 2, 1 + Polynomial('x')), 2: Root(g, 1)}
        assert T.unconditional_descents == set()
        assert T.strong_conditional_descents == set()
        assert T.weak_conditional_descents == set()

    def test_constructor_errors(self):
        # Check error handling for invalid constructor inputs.
        g = CoxeterGraph.A(5)
        try:
            RootTransform(g, {0: 0})
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        try:
            RootTransform(g, {1: Root(CoxeterGraph.B(5), 1)})
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        T = RootTransform(g, {1: -Root(g, 2), 2: Root(g, 1)})
        try:
            T[1] = Root(CoxeterGraph.B(5), 1)
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

    def test_eq(self):
        g = CoxeterGraph.A(5)
        T = RootTransform(g, {1: Root(g, 2, 1 + Polynomial('x')), 2: Root(g, 1)})
        U = T
        V = T.copy()

        T[1] = Root(g)
        assert U[1] == 0 and U == T
        assert V[1] == Root(g, 2, 1 + Polynomial('x')) and T != V
        assert len({T, U, V}) == len({T, V}) == 2

    def test_multiply(self):
        g = CoxeterGraph.A(5)
        T = RootTransform(g, {1: Root(g, 2, 1 + Polynomial('x')), 2: Root(g, 1)})

        U = T * 1
        assert U.sigma == {
            1: -Root(g, 2, 1 + Polynomial('x')),
            2: Root(g, 1) + Root(g, 2, 1 + Polynomial('x'))
        }

        V = 1 * T
        assert V.sigma == {
            1: (Root(g, 1) + Root(g, 2)) * (1 + Polynomial('x')),
            2: -Root(g, 1)
        }

        W = 1 * T * 1
        assert W.sigma == {
            1: -(Root(g, 1) + Root(g, 2))*(1 + Polynomial('x')),
            2: Root(g, 1, Polynomial('x')) + Root(g, 2, 1 + Polynomial('x'))
        }

    def test_multiplication_errors(self):
        g = CoxeterGraph.A(5)
        T = RootTransform(g, {1: Root(g, 2, 1 + Polynomial('x')), 2: Root(g, 1)})

        # cannot multiply by 6 since this is not a generator for the Coxeter system A5
        try:
            T * 6
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        try:
            6 * T
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        # cannot multiply by 3 since (alpha_2).reflect(3) = alpha_2 + alpha_3 and T[3] is undefined
        try:
            T * 3
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

    def test_convenience_methods(self):
        """Tests for various convenience methods of RooTransform objects."""
        g = CoxeterGraph.A(5)
        sigma = {1: Root(g, 2, 1 + Polynomial('x')), 2: Root(g, 1)}
        T = RootTransform(g, sigma)

        assert not T.is_constant()
        assert not T.is_complete()
        assert T.is_positive()
        assert not T.is_identity()

        assert len(T) == 2
        assert 1 in T and 2 in T and 3 not in T
        assert {i for i in T} == {1, 2}
        assert set(T.values()) == set(sigma.values())

        U = RootTransform.identity(g)
        assert U.sigma == {i: Root(g, i) for i in g.generators}
        assert U.is_constant()
        assert U.is_complete()
        assert U.is_positive()
        assert U.is_identity()

        # check (non-)conversion of RootTransform to CoxeterTransform
        try:
            T.to_coxeter_transform()
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        V = U * 1 * 2 * 3
        assert V.to_coxeter_transform().minimal_reduced_word == (1, 2, 3)

        # try converting RootTransforms to str, check that no errors occur
        str(U)
        str(V)


class TestCoxeterTransform:
    def test_constructor(self):
        g = CoxeterGraph.B(3)
        T = CoxeterTransform(g)
        assert T.graph == g
        assert T.sigma == {i: Root(g, i) for i in g.generators}
        assert T.right_descents == set()
        assert T.minimal_right_descent is None
        assert T.minimal_reduced_word == ()

        sigma = {i: -Root(g, i) for i in g.generators}
        T = CoxeterTransform(g, sigma)
        assert T.graph == g
        assert T.sigma == sigma
        sigma.clear()
        assert T.sigma == {i: -Root(g, i) for i in g.generators}
        assert T.right_descents == {0, 1, 2}
        assert T.minimal_right_descent is 0
        assert T.minimal_reduced_word == (2, 1, 0, 1, 2, 1, 0, 1, 0)

    def test_constructor_errors(self):
        # Check error handling for invalid constructor inputs.
        g = CoxeterGraph.B(3)
        try:
            CoxeterTransform(g, {0: Root(g, 0)})
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        h = CoxeterGraph.A(3)
        try:
            CoxeterTransform(g, {0: Root(g, 0), 1: Root(g, 1), 2: Root(h, 2)})
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        try:
            CoxeterTransform(g, {0: Root(g, 0), 1: Root(g, 1), 2: Root(g, Polynomial('x'))})
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

    def test_get_inverse(self):
        g = CoxeterGraph.B(3)
        T = CoxeterTransform(g, {i: -Root(g, i) for i in g.generators})
        assert T == T.get_inverse()

        U = CoxeterTransform(g)
        assert (U * 0 * 1).get_inverse() == U * 1 * 0

    def test_eq(self):
        g = CoxeterGraph.B(2)
        T = RootTransform(g, {0: Root(g, 0), 1: Root(g, 1)})
        U = CoxeterTransform(g, {0: Root(g, 0), 1: Root(g, 1)})
        assert T == U

        V = U.copy()
        V[0] = -Root(g, 0)
        V[1] = -Root(g, 1)

        assert T == U and U != V
        assert not hasattr(T, 'minimal_reduced_word')
        assert U.minimal_reduced_word == ()
        assert V.minimal_reduced_word == (1, 0, 1, 0)
        assert len({T, U, V}) == 2

    def test_setitem_errors(self):
        g = CoxeterGraph.B(2)
        V = CoxeterTransform(g, {0: Root(g, 0), 1: Root(g, 1)})

        # cannot assign value which is non-constant Root
        try:
            V[0] = -Root(g, 0, Polynomial('x'))
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        # cannot assign to index which is not in g.generators
        try:
            V[2] = Root(g, 1)
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

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

        g = CoxeterGraph.A2(5)
        e = CoxeterTransform(g)

        e = e.demazure_conjugate(1)
        assert e.minimal_reduced_word == (5, 1)

        try:
            e.demazure_conjugate(0)
        except Exception as exc:
            assert type(exc) == InvalidInputException
        else:
            assert False


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

        w = CoxeterWord(g, (1, 0, 1, 0))
        assert w.word == (1, 0, 1, 0)
        assert w.left_action == CoxeterTransform.from_word(g, (0, 1, 0, 1))
        assert w.right_action == w.left_action.get_inverse()
        assert w.is_reduced
        assert w.left_descents == {0, 1}
        assert w.right_descents == {0, 1}
        assert w.minimal_reduced_word == (0, 1, 0, 1)

        v = CoxeterWord(g, (0, 1, 0, 1))
        assert v != w and v.left_action == w.left_action and v.right_action == v.right_action
        assert hash(v) != hash(w)

        v = CoxeterWord(g, (0, 1, 0, 1, 1, 1))
        assert not v.is_reduced
        assert v != w and v.left_action == w.left_action and v.right_action == v.right_action
        assert hash(v) != hash(w)

        assert len(w) == 4 and len(v) == 6
        assert str(v) == 'unreduced word (0, 1, 0, 1, 1, 1)'
        assert str(w) == 'reduced word (1, 0, 1, 0)'

    @pytest.mark.parametrize("word, involution_word", [
        ((1, 0), (0, 1, 0)),
        ((0, 2), (0, 2)),
        ((1, 0, 1, 0), (0, 1, 0, 1, 0))
    ])
    def test_to_involution(self, word, involution_word):
        g = CoxeterGraph.B(3)
        w = CoxeterWord(g, word)
        assert w.to_involution().word == involution_word

    def test_get_reduced_words(self):
        g = CoxeterGraph.B(3)
        w = CoxeterWord(g, (1, 0, 1, 0))
        assert w.get_reduced_words() == {(1, 0, 1, 0), (0, 1, 0, 1)}

    def test_copy(self):
        g = CoxeterGraph.B(3)
        w = CoxeterWord(g)
        v = w.copy()
        assert v == w

        w.extend_left(0)
        assert w.word == (0,) and v.word == ()

    def test_multiply(self):
        g = CoxeterGraph.B(3)
        assert CoxeterWord(g) * CoxeterWord(g) == CoxeterWord(g)

        w = ((CoxeterWord(g) * 0) * 1) * 2
        assert w.word == (0, 1, 2)
        assert (w * w).word == (0, 1, 2, 0, 1, 2)

        w = 0 * (1 * (2 * CoxeterWord(g)))
        assert w.word == (0, 1, 2)
        assert (w * w).word == (0, 1, 2, 0, 1, 2)
        assert w.__rmul__(w).word == (0, 1, 2, 0, 1, 2)

        w = 0 * (1 * (2 * ((CoxeterWord(g) * 0) * 1))) * 2
        assert w.word == (0, 1, 2, 0, 1, 2)

    def test_multiply_errors(self):
        g = CoxeterGraph.B(3)
        h = CoxeterGraph.A(3)

        # cannot multiply by non-generator 3, as g.generators = [0, 1, 2]
        try:
            CoxeterWord(g) * 3
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        try:
            3 * CoxeterWord(g)
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        try:
            CoxeterWord(g).extend_right(3)
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        try:
            CoxeterWord(g).extend_left(3)
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        # cannot multiply CoxeterWords with different CoxeterGraphs
        try:
            CoxeterWord(g).__mul__(CoxeterWord(h))
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

        try:
            CoxeterWord(g).__rmul__(CoxeterWord(h))
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False
