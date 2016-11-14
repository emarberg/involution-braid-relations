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
