import pytest
import numpy as np

from utils import (
    InvalidInputException
)

from algebra import (
    QuadraticNumber,
    RationalNumber,
    Monomial,
    Polynomial,
    CoxeterGraph,
    Root
)

from solver import (
    ConstraintsManager
)


class TestConstraintsManager:
    def test_add_zero_constraint(self):
        manager = ConstraintsManager()
        a = 12
        b = RationalNumber(7, 8)
        c = QuadraticNumber.sqrt(11)
        d = Polynomial('x')
        e = Polynomial('x') * Polynomial('y')
        f = Polynomial('x')**3

        g = CoxeterGraph.A(5)
        r = Root(g, 1, a) + Root(g, 2, b) + Root(g, 3, c) + Root(g, 4, d) + Root(g, 5, e)

        manager.add_zero_constraint(a)
        assert a in manager.linear_constraints

        manager.add_zero_constraint(b)
        assert b in manager.linear_constraints

        manager.add_zero_constraint(c)
        assert c in manager.linear_constraints

        manager.add_zero_constraint(d)
        assert d in manager.linear_constraints

        manager.add_zero_constraint(e)
        assert e in manager.quadratic_constraints

        try:
            manager.add_zero_constraint(f)
        except Exception as exc:
            assert type(exc) == InvalidInputException
        else:
            assert False

        try:
            manager.add_zero_constraint(None)
        except Exception as exc:
            assert type(exc) == InvalidInputException
        else:
            assert False

        manager.linear_constraints = set()
        manager.quadratic_constraints = set()
        assert len(manager.linear_constraints) == len(manager.quadratic_constraints) == 0

        manager.add_zero_constraint(r)
        assert manager.linear_constraints == {a, b, c, d}
        assert manager.quadratic_constraints == {e}

    def test_add_nonpositive_constraint(self):
        manager = ConstraintsManager()
        a = -3
        b = RationalNumber(-2)
        c = QuadraticNumber(-1)
        d = -Polynomial('x')
        e = -Polynomial('x') * Polynomial('y')
        f = -Polynomial('x')**3

        g = CoxeterGraph.A(5)
        r = Root(g, 1, a) + Root(g, 4, d) + Root(g, 5, e)

        manager.add_nonpositive_constraint(a)
        manager.add_nonpositive_constraint(b)
        manager.add_nonpositive_constraint(c)
        manager.add_nonpositive_constraint(d)
        manager.add_nonpositive_constraint(e)
        manager.add_nonpositive_constraint(f)
        assert manager.nonpositive_constraints == {a, b, c, d, e, f}

        try:
            manager.add_nonpositive_constraint(None)
        except Exception as exc:
            assert type(exc) == InvalidInputException
        else:
            assert False

        manager.nonpositive_constraints = set()
        assert len(manager.nonpositive_constraints) == 0

        manager.add_nonpositive_constraint(r)
        assert manager.nonpositive_constraints == {a, d, e}

        manager.add_nonpositive_constraint(-e)
        assert -e in manager.quadratic_constraints

        manager.nonpositive_constraints = set()
        assert len(manager.nonpositive_constraints) == 0

        manager.add_nonpositive_constraint(c)
        manager.add_nonpositive_constraint(b)
        manager.add_nonpositive_constraint(a)
        manager.add_nonpositive_constraint(d)
        manager.add_nonpositive_constraint(d + e)
        assert manager.nonpositive_constraints == {c, d}

    def test_add_nonzero_constraint(self):
        manager = ConstraintsManager()
        a = -3
        b = -Polynomial('x')
        c = -Polynomial('x') * Polynomial('y')

        g = CoxeterGraph.A(5)
        r = Root(g, 1, a) + Root(g, 4, b) + Root(g, 5, c)

        manager.add_nonzero_constraint(r)
        assert manager.nonzero_constraints == {r}

        try:
            manager.add_nonzero_constraint(c)
        except Exception as exc:
            assert type(exc) == InvalidInputException
        else:
            assert False

    def test_repr(self):
        manager = ConstraintsManager()
        assert str(manager) == ''

        manager.add_zero_constraint(Polynomial('x'))
        manager.add_zero_constraint(Polynomial('x') * Polynomial('y'))
        manager.add_nonpositive_constraint(Polynomial('x') - 1)
        manager.add_nonzero_constraint(Root(CoxeterGraph.A(5), 1, Polynomial('y')))
        assert str(manager) != ''

    def test_copy(self):
        manager = ConstraintsManager()
        a = 12
        b = RationalNumber(7, 8)
        c = QuadraticNumber.sqrt(11)
        d = Polynomial('x')
        e = Polynomial('x') * Polynomial('y')

        manager.add_zero_constraint(a)
        manager.add_zero_constraint(b)
        manager.add_zero_constraint(c)

        other = manager.copy()
        assert manager == other

        other.add_zero_constraint(d)
        other.add_zero_constraint(e)

        assert manager.linear_constraints == {a, b, c}
        assert other.linear_constraints == {a, b, c, d}
        assert manager.quadratic_constraints == set()
        assert other.quadratic_constraints == {e}
        assert manager != other

    def test_remove_vacuous_constraints(self):
        manager = ConstraintsManager()

        # check that method removes zero polynomials from sets of linear and quadratic constraints
        manager.add_zero_constraint(Polynomial(0))
        assert manager.linear_constraints == {0}
        manager.remove_vacuous_constraints()
        assert manager.linear_constraints == set()

        manager.quadratic_constraints = {Polynomial(0)}
        assert manager.quadratic_constraints == {0}
        manager.remove_vacuous_constraints()
        assert manager.quadratic_constraints == set()

        # check that method removes from set of nonpositive constraints any polynomials which have
        # all nonpositive coefficients or which are bounded above by other nonpositive constraints
        x = Polynomial('x')
        manager.add_nonpositive_constraint(-x - 1)
        manager.add_nonpositive_constraint(x - 2)
        manager.add_nonpositive_constraint(x - 1)
        assert manager.nonpositive_constraints == {-x - 1, x - 2, x - 1}
        manager.remove_vacuous_constraints()
        assert manager.nonpositive_constraints == {x - 1}

        # check that method removes nonzero roots from set of nonzero constraints
        g = CoxeterGraph.A(5)
        r = Root(g, 1)
        manager.add_nonzero_constraint(r)
        assert manager.nonzero_constraints == {r}
        manager.remove_vacuous_constraints()
        assert manager.nonzero_constraints == set()

    def test_simplify(self):
        g = CoxeterGraph.A(5)
        x = Polynomial('x')
        y = Polynomial('y')
        manager = ConstraintsManager()

        manager.add_zero_constraint(x + y - 2)
        manager.add_zero_constraint(x - y + 2)
        manager.add_zero_constraint(y)
        manager.add_nonpositive_constraint(3 - y)
        manager.add_zero_constraint(x + y**2)
        manager.add_nonzero_constraint(Root(g, 1, x))

        variable_substitutions = set(manager.simplify())
        assert variable_substitutions == {(Monomial('x'), 0), (Monomial('y'), 2)}
        assert manager.linear_constraints == {2}
        assert manager.quadratic_constraints == {4}
        assert manager.nonpositive_constraints == {1}
        assert manager.nonzero_constraints == {Root(g)}
        assert not manager.is_valid()
