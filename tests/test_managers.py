import pytest

from project.utils import (
    reverse_tuple,
    InvalidInputException
)

from project.algebra import (
    QuadraticNumber,
    RationalNumber,
    Monomial,
    Polynomial,
    CoxeterGraph,
    CoxeterTransform,
    Root
)

from project.managers import (
    ConstraintsManager,
    PartialBraid,
    BraidQueue
)


@pytest.mark.parametrize("tup, expected", [
    ((), ()),
    ((1,), (1,)),
    ((1, 2), (2, 1)),
    ((1, 2, 3), (3, 2, 1)),
])
def test_reverse_tuple(tup, expected):
    assert reverse_tuple(tup) == expected


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

        # zero contraints must be linear or quadratic Polynomials
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

        # check that adding Root as nonpositive constraint introduces constraint for each coeff
        manager.add_nonpositive_constraint(r)
        assert manager.nonpositive_constraints == {a, d, e}

        # check that if nonpositive constraints contains bothx and -x,
        # then x gets added as a zero constraints
        manager.add_nonpositive_constraint(-e)
        assert -e in manager.quadratic_constraints

        manager.nonpositive_constraints = set()
        assert len(manager.nonpositive_constraints) == 0

        manager.add_nonpositive_constraint(c)
        manager.add_nonpositive_constraint(b)
        manager.add_nonpositive_constraint(a)
        manager.add_nonpositive_constraint(d)
        manager.add_nonpositive_constraint(d + e)
        # since a <= b < c and d + e <= d, we should observe only the following constraints:
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

        # input to manager.add_nonzero_constraint must be a Root
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

    def test_eq(self):
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
        """Tests for method of ConstraintsManager which reduces various constraints."""
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


class TestPartialBraid:
    def test_partial_braid_constructor_errors(self):
        g = CoxeterGraph.A(3)

        # input (s, t) must both be in g.generators
        try:
            PartialBraid(g, s=0, t=1)
        except Exception as e:
            assert type(e) == InvalidInputException
        else:
            assert False

    @pytest.mark.parametrize("m, is_fixer, expected", [
        (5, True, (1, 2, 1)),
        (5, False, (1, 2, 1)),
        (6, True, (2, 1, 2, 1)),
        (6, False, (1, 2, 1)),
    ])
    def test_extend_words(self, m, is_fixer, expected):
        g = CoxeterGraph([(1, 2, m)])

        state = PartialBraid(g, 1, 2)
        assert state.word_s.word == () and state.word_s.word == ()

        state._extend_words(is_fixer)
        assert state.word_s.word == expected

        other = tuple(map(lambda i: 3-i, expected))
        assert state.word_t.word == other


class TestBraidQueue:
    def test_A3(self):
        g = CoxeterGraph.A(3)

        # test algorithm where trivial output is expected
        q = BraidQueue(g, 1, 3, verbose_level=BraidQueue.VERBOSE_LEVEL_HIGH)
        q.go()
        assert q.sufficient_relations == set() and q.minimal_relations == []

        # test algorithm in small case where nontrivial output is expected
        q = BraidQueue(g, verbose_level=BraidQueue.VERBOSE_LEVEL_HIGH)
        assert {(state.s, state.t) for state in q.queue} == {(1, 2), (1, 3), (2, 3)}
        assert q.sufficient_relations == set()
        assert q.minimal_relations == []

        q.go()
        assert q.sufficient_relations == {
            ((1, 2), (2, 1)),
            ((2, 3), (3, 2))
        }
        assert q.minimal_relations == [
            ((1, 2), (2, 1)),
            ((2, 3), (3, 2))
        ]

    def test_B3(self):
        g = CoxeterGraph.B(3)
        q = BraidQueue(g, verbose_level=BraidQueue.VERBOSE_LEVEL_LOW)
        q.go()
        assert q.sufficient_relations == {
            ((1, 2), (2, 1)),
            ((0, 1, 0), (1, 0, 1)),
            ((0, 1, 2, 0, 1, 0), (0, 1, 2, 1, 0, 1))
        }
        assert q.minimal_relations == [
            ((1, 2), (2, 1)),
            ((0, 1, 0), (1, 0, 1)),
            ((0, 1, 2, 0, 1, 0), (0, 1, 2, 1, 0, 1))
        ]

    def test_2A3(self):
        # Test algorithm in small twisted case
        g = CoxeterGraph.A2(3)
        q = BraidQueue(g, verbose_level=BraidQueue.VERBOSE_LEVEL_LOW)
        q.go(do_sanity_check=True)
        assert q.sufficient_relations == {
            ((1,), (3,)),
            ((2, 1, 2, 3), (2, 1, 3, 2)),
            ((2, 3, 1, 2), (2, 3, 2, 1))
        }
        assert q.minimal_relations == [
            ((1,), (3,)),
            ((2, 1, 2, 3), (2, 1, 3, 2))
        ]

        # check that the following does not cause any errors
        q.next()

    def test_get_next_level_of_involutions_to_atoms(self):
        g = CoxeterGraph.A(2)

        level = BraidQueue._get_next_level_of_involutions_to_atoms(g)
        assert level == {CoxeterTransform(g): {CoxeterTransform(g)}}

        level = BraidQueue._get_next_level_of_involutions_to_atoms(g, level)
        assert level == {
            CoxeterTransform.from_word(g, (1,)): {CoxeterTransform.from_word(g, (1,))},
            CoxeterTransform.from_word(g, (2,)): {CoxeterTransform.from_word(g, (2,))}
        }

        level = BraidQueue._get_next_level_of_involutions_to_atoms(g, level)
        assert level == {
            CoxeterTransform.from_word(g, (1, 2, 1)): {
                CoxeterTransform.from_word(g, (1, 2)),
                CoxeterTransform.from_word(g, (2, 1))
            },
        }

        level = BraidQueue._get_next_level_of_involutions_to_atoms(g, level)
        assert level == {}
