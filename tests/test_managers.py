import pytest

from project.utils import (
    reverse_tuple,
    InvalidInputException
)

from project.algebra import (
    QuadraticNumber,
    RationalNumber,
    Monomial,
    Polynomial
)

from project.coxeter import (
    CoxeterGraph,
    CoxeterTransform,
    CoxeterWord,
    CoxeterVector,
    PartialTransform
)

from project.managers import (
    ConstraintsManager,
    BraidSystem,
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
        e = Polynomial('x')**3

        g = CoxeterGraph.A(5)
        r = CoxeterVector(g, 1, a) + \
            CoxeterVector(g, 2, b) + \
            CoxeterVector(g, 3, c) + \
            CoxeterVector(g, 4, d)

        manager.add_zero_constraint(a)
        assert a in manager.zero_constraints

        manager.add_zero_constraint(b)
        assert b in manager.zero_constraints

        manager.add_zero_constraint(c)
        assert c in manager.zero_constraints

        manager.add_zero_constraint(d)
        assert d in manager.zero_constraints

        # zero contraints must be linear Polynomials
        exception = None
        try:
            manager.add_zero_constraint(e)
        except Exception as exc:
            exception = exc
        assert type(exception) == InvalidInputException

        exception = None
        try:
            manager.add_zero_constraint(None)
        except Exception as exc:
            exception = exc
        assert type(exception) == InvalidInputException

        manager.zero_constraints = set()
        assert len(manager.zero_constraints) == 0

        manager.add_zero_constraint(r)
        assert manager.zero_constraints == {a, b, c, d}

    def test_add_nonpositive_constraint(self):
        manager = ConstraintsManager()
        a = -3
        b = RationalNumber(-2)
        c = QuadraticNumber(-1)
        d = -Polynomial('x')
        e = -Polynomial('x') * Polynomial('y')
        f = -Polynomial('x')**3

        g = CoxeterGraph.A(5)
        r = CoxeterVector(g, 1, a) + CoxeterVector(g, 4, 1 + d) + CoxeterVector(g, 5, 1 + e)

        manager.add_nonpositive_constraint(a)
        manager.add_nonpositive_constraint(b)
        manager.add_nonpositive_constraint(c)
        manager.add_nonpositive_constraint(d)
        manager.add_nonpositive_constraint(e)
        manager.add_nonpositive_constraint(f)
        assert manager.nonpositive_constraints == set()

        exception = None
        try:
            manager.add_nonpositive_constraint(None)
        except Exception as exc:
            exception = exc
        assert type(exception) == InvalidInputException

        manager.nonpositive_constraints = set()
        assert len(manager.nonpositive_constraints) == 0

        # check that adding CoxeterVector as constraint introduces constraint for each coeff
        manager.add_nonpositive_constraint(-r)
        assert manager.nonpositive_constraints == {-1 - d, -1 - e}
        assert manager.zero_constraints == {-a}

    def test_add_nonzero_constraint(self):
        manager = ConstraintsManager()
        a = -3
        b = -Polynomial('x')
        c = -Polynomial('x') * Polynomial('y')

        g = CoxeterGraph.A(5)
        r = CoxeterVector(g, 1, a) + CoxeterVector(g, 4, b) + CoxeterVector(g, 5, c)

        manager.add_nonzero_constraint(r)
        assert manager.nonzero_roots == {r}

        # input to manager.add_nonzero_constraint must be a CoxeterVector
        exception = None
        try:
            manager.add_nonzero_constraint(c)
        except Exception as exc:
            exception = exc
        assert type(exception) == InvalidInputException

    def test_repr(self):
        manager = ConstraintsManager()
        assert str(manager) == ''

        manager.add_zero_constraint(Polynomial('x'))
        manager.add_nonpositive_constraint(Polynomial('x') - 1)
        manager.add_nonzero_constraint(CoxeterVector(CoxeterGraph.A(5), 1, Polynomial('y')))
        assert str(manager) != ''

    def test_remove_vacuous_constraints(self):
        manager = ConstraintsManager()

        # check that method removes zero polynomials from sets of linear constraints
        manager.add_zero_constraint(Polynomial(0))
        assert manager.zero_constraints == {0}
        manager.remove_vacuous_constraints()
        assert manager.zero_constraints == set()

        # check that method removes from set of nonpositive constraints any polynomials which have
        # all nonpositive coefficients or which are bounded above by other nonpositive constraints
        x = Polynomial('x')
        manager.add_nonpositive_constraint(-x - 1)
        manager.add_nonpositive_constraint(x - 2)
        manager.add_nonpositive_constraint(x - 1)
        assert manager.nonpositive_constraints == {x - 2, x - 1}
        manager.remove_vacuous_constraints()
        assert manager.nonpositive_constraints == {x - 1}

        # check that method removes nonzero roots from set of nonzero constraints
        g = CoxeterGraph.A(5)
        r = CoxeterVector(g, 1)
        manager.add_nonzero_constraint(r)
        assert manager.nonzero_roots == {r}
        manager.remove_vacuous_constraints()
        assert manager.nonzero_roots == set()

    def test_simplify(self):
        """Tests for method of ConstraintsManager which reduces various constraints."""
        g = CoxeterGraph.A(5)
        x = Polynomial('x')
        y = Polynomial('y')
        z = Polynomial('z')
        manager = ConstraintsManager()

        manager.add_zero_constraint(x + y - 2)
        manager.add_zero_constraint(x - y + 2)
        manager.add_zero_constraint(y + z)
        manager.add_nonpositive_constraint(3 - y)
        manager.add_nonzero_constraint(CoxeterVector(g, 1, x))

        variable_substitutions = set(manager.simplify())
        assert variable_substitutions == {
            (Monomial('x'), 0),
            (Monomial('y'), 2),
            (Monomial('z'), -2)
        }
        assert manager.zero_constraints == set()
        assert manager.nonpositive_constraints == {1}
        assert manager.nonzero_roots == {CoxeterVector(g)}
        assert not manager.is_valid()

    def test_is_value_nonpositive(self):
        x = Polynomial('x')
        y = Polynomial('y')
        manager = ConstraintsManager()

        assert manager.is_value_nonpositive(-1)
        assert manager.is_value_nonpositive(-x - 1)
        assert not manager.is_value_nonpositive(x - 1)

        manager.nonpositive_constraints = {3 + x}
        assert manager.is_value_nonpositive(5 * x)
        assert manager.is_value_nonpositive(5)
        assert not manager.is_value_nonpositive(y)
        assert not manager.is_value_nonpositive(x**2)

    def test_is_value_negative(self):
        x = Polynomial('x')
        y = Polynomial('y')
        manager = ConstraintsManager()

        assert manager.is_value_negative(-1)
        assert manager.is_value_negative(-x - 1)
        assert not manager.is_value_negative(-x)

        manager.nonpositive_constraints = {x}
        assert not manager.is_value_negative(2 * x)

        manager.nonpositive_constraints = {-3 + x}
        assert not manager.is_value_negative(-6 + 2 * x)
        assert not manager.is_value_negative(-5 + 2 * x)
        assert manager.is_value_negative(-7 + 2 * x)

        assert not manager.is_value_negative(-y)
        assert not manager.is_value_negative(-x**2)


class TestBraidSystem:
    def test_braid_system_constructor_errors(self):
        g = CoxeterGraph.A(3)

        # input (s, t) must both be in g.generators
        e = None
        try:
            BraidSystem(g, s=0, t=1)
        except Exception as exception:
            e = exception
        assert type(e) == InvalidInputException

    @pytest.mark.parametrize("m, is_fixer, expected", [
        (5, True, (1, 2, 1)),
        (5, False, (1, 2, 1)),
        (6, True, (2, 1, 2, 1)),
        (6, False, (1, 2, 1)),
    ])
    def test_extend_words(self, m, is_fixer, expected):
        g = CoxeterGraph([(1, 2, m)])

        state = BraidSystem(g, 1, 2, is_fixer)
        assert state.word_s.word == expected

        other = tuple(map(lambda i: 3 - i, expected))
        assert state.word_t.word == other

    def test_get_unconditional_descent(self):
        g = CoxeterGraph.A(3)
        state = BraidSystem(g, s=1, t=2)
        state._extend_words()

        state.sigma = PartialTransform(g, {1: -CoxeterVector(g, 1), 3: -CoxeterVector(g, 3)})
        assert state.get_unconditional_descent() == 1

        state.sigma = PartialTransform(g, {1: -CoxeterVector(g, 1), 2: -CoxeterVector(g, 3)})
        assert state.get_unconditional_descent() == 1

        state.sigma = PartialTransform(
            g, {1: CoxeterVector(g, 1), 3: -CoxeterVector(g, 3, 1 + Polynomial('x'))}
        )
        assert state.get_unconditional_descent() == 3

        state.sigma = PartialTransform(
            g, {1: CoxeterVector(g, 1), 3: -CoxeterVector(g, 3, Polynomial('x'))}
        )
        assert state.get_unconditional_descent() is None

        state.sigma = PartialTransform(
            g, {1: CoxeterVector(g, 1), 3: CoxeterVector(g, 3, 1 - Polynomial('x'))}
        )
        assert state.get_unconditional_descent() is None

    def test_get_conditional_descent(self):
        x = Polynomial('x')
        g = CoxeterGraph.A(3)
        state = BraidSystem(g, s=1, t=2)

        state.sigma = PartialTransform(g, {1: -CoxeterVector(g, 1), 3: -CoxeterVector(g, 3)})
        assert state.get_conditional_descent() == 1

        state.sigma = PartialTransform(
            g, {1: -CoxeterVector(g, 1) + CoxeterVector(g, 2) - CoxeterVector(g, 3, x)}
        )
        assert state.get_conditional_descent() == 1

        state.sigma = PartialTransform(g, {
            1: (1 - x) * CoxeterVector(g, 2),
            2: (x - 1) * (CoxeterVector(g, 1) - CoxeterVector(g, 2) + CoxeterVector(g, 3))
        })
        assert state.get_conditional_descent() is None

        state.constraints.nonpositive_constraints.add(x - 1)
        assert state.get_conditional_descent() == 2

    def test_get_children_from_conditional_descent(self):
        x = Polynomial('x')
        g = CoxeterGraph.A(3)
        state = BraidSystem(g, s=1, t=3)
        state.word_s = CoxeterWord(g, (1,))
        state.word_t = CoxeterWord(g, (3,))

        alpha = CoxeterVector(g, 3)
        beta = \
            CoxeterVector(g, 1, -x) + \
            CoxeterVector(g, 2, 1 - 2 * x) + \
            CoxeterVector(g, 3, -x)
        gamma = CoxeterVector(g, 1)
        state.sigma = PartialTransform(g, {1: alpha, 2: beta, 3: gamma})

        children = state._get_children_from_conditional_descent()
        assert len(children) == 3
        assert \
            sorted((c.word_s.word, c.word_t.word) for c in children) == \
            [((1,), (3,)), ((2, 1), (2, 3)), ((2, 1), (2, 3))]

    def test_is_sigma_valid(self):
        g = CoxeterGraph.A(3)
        state = BraidSystem(g, s=1, t=2)

        state.sigma = PartialTransform(g, {1: CoxeterVector(g, 1) - CoxeterVector(g, 2)})
        assert not state._is_sigma_valid()

    def test_is_realized(self):
        g = CoxeterGraph.A(3)
        state = BraidSystem(g, s=1, t=2)
        state.sigma = PartialTransform(g, {
            1: CoxeterVector(g, 3),
            2: CoxeterVector(g, 2),
            3: CoxeterVector(g, 1)
        })
        assert not state.is_realized()

        state.sigma = PartialTransform.identity(g)
        state.word_s = CoxeterWord(g, (1, 2))
        state.word_t = CoxeterWord(g, (3, 2))
        assert not state.is_realized()

        state.word_s = CoxeterWord(g, (1, 2))
        state.word_t = CoxeterWord(g, (2, 1))
        assert state.is_realized()

    def test_is_descent_periodic(self):
        g = CoxeterGraph.G_tilde(2)
        state = BraidSystem(g, s=1, t=2)
        state.word_s = CoxeterWord(g, (3, 2, 1, 2, 1, 2, 3, 2, 1, 2, 1, 2, 3, 1, 2, 1))
        state.word_t = CoxeterWord(g, (3, 2, 1, 2, 1, 2, 3, 2, 1, 2, 1, 2, 3, 2, 1, 2))

        alpha = \
            CoxeterVector(g, 1, 5 * QuadraticNumber.sqrt(3)) + \
            CoxeterVector(g, 2, 11) + \
            CoxeterVector(g, 3, 6)
        beta = \
            CoxeterVector(g, 1, -16 - 6 * QuadraticNumber.sqrt(3)) + \
            CoxeterVector(g, 2, -12 - 11 * QuadraticNumber.sqrt(3)) + \
            CoxeterVector(g, 3, -6 - 6 * QuadraticNumber.sqrt(3))
        gamma = \
            CoxeterVector(g, 1, 17 + 11 * QuadraticNumber.sqrt(3)) + \
            CoxeterVector(g, 2, 22 + 11 * QuadraticNumber.sqrt(3)) + \
            CoxeterVector(g, 3, 11 + 6 * QuadraticNumber.sqrt(3))

        state.sigma = PartialTransform(g, {1: alpha, 2: beta, 3: gamma})
        history = [state]
        for _ in range(16):
            new = history[-1]
            descent = new.get_unconditional_descent()
            commutes = new.sigma[descent] == -CoxeterVector(new.graph, new.graph.star(descent))
            print(descent, new)
            new = new._branch_from_descent(descent, commutes=commutes).reduce()
            history.append(new)
        assert new.is_descent_periodic(history)

        q = BraidQueue(g)
        q.descent_periodic_systems = [state]
        q.minimize_relations()

        state.sigma = PartialTransform(g, {1: 2 * alpha, 2: beta, 3: gamma})
        assert not state.is_descent_periodic([state])

        state.sigma = PartialTransform(g, {1: Polynomial('x') * alpha, 2: beta, 3: gamma})
        assert not state.is_descent_periodic([state])

    def test_is_redundant(self):
        g = CoxeterGraph.A_twist(5)
        alpha = {i: CoxeterVector(g, i) for i in g.generators}
        state = BraidSystem(g, 3, 4, True)
        state.word_s = CoxeterWord(g, (1, 2, 4, 3))
        state.word_t = CoxeterWord(g, (1, 2, 3, 4))

        state.sigma = PartialTransform(g, {1: alpha[5], 2: alpha[4], 3: alpha[3], 4: alpha[2]})
        assert state.is_redundant()

    def test_repr(self):
        a = 'initial sigma: alpha_s -> alpha_s^*, alpha_t -> alpha_t^*'
        b = 'initial sigma: alpha_s -> alpha_t^*, alpha_t -> alpha_s^*'
        g = CoxeterGraph.D(5)
        state = BraidSystem(g, 3, 4, True)
        assert a in str(state) and b not in str(state)
        state = BraidSystem(g, 3, 4, False)
        assert a not in str(state) and b in str(state)


class TestBraidQueue:
    def test_A3(self):  # noqa
        g = CoxeterGraph.A(3)

        # test algorithm where trivial output is expected
        q = BraidQueue(g, 1, 3)
        q.go()
        assert q.sufficient_relations == {
            ((1, 2), (2, 1)),
            ((2, 3), (3, 2))
        }
        assert q.minimal_relations == [
            ((1, 2), (2, 1)),
            ((2, 3), (3, 2))
        ]

        # test algorithm in small case where nontrivial output is expected
        q = BraidQueue(g)
        assert {(state.s, state.t) for state in q.queue} == {(1, 2), (2, 3)}
        assert q.sufficient_relations == {((1, 2), (2, 1)), ((2, 3), (3, 2))}
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

    def test_B3(self):  # noqa
        g = CoxeterGraph.B(3)
        q = BraidQueue(g)
        q.go()
        assert q.minimal_relations == [
            ((1, 2), (2, 1)),
            ((2, 3, 2), (3, 2, 3)),
            ((3, 2, 1, 2, 3, 2), (3, 2, 1, 3, 2, 3))
        ]

    def test_2A3(self):  # noqa
        # Test algorithm in small twisted case
        g = CoxeterGraph.A_twist(3)
        q = BraidQueue(g)
        assert {(state.s, state.t) for state in q.queue} == {(1, 2), (2, 3), (1, 3)}
        assert q.sufficient_relations == {((1,), (3,))}

        q.go(verify=True)
        assert q.sufficient_relations == {
            ((1,), (3,)),
            ((2, 1, 2, 3), (2, 1, 3, 2))
        }
        assert q.minimal_relations == [
            ((1,), (3,)),
            ((2, 1, 2, 3), (2, 1, 3, 2))
        ]

        # check that the following does not cause any errors
        q.next()

    def test_D4(self):  # noqa
        # Test algorithm in nontrivial case
        g = CoxeterGraph.D(4)
        q = BraidQueue(g, 2, 3)
        q.go(verify=False)
        assert q.minimal_relations == [
            ((1, 2), (2, 1)),
            ((2, 3), (3, 2)),
            ((2, 4), (4, 2)),
            ((3, 1, 2, 4, 2, 1, 2, 3), (3, 1, 2, 4, 2, 1, 3, 2))
        ]

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

    def test_verify_relations(self):
        g = CoxeterGraph.B(3)
        q = BraidQueue(g)

        # check error handling when too many relations are present
        q.minimal_relations = [((1, 2), (2, 1)), ((2, 3), (3, 2))]
        try:
            q.verify_relations(None)
        except Exception as e:
            assert str(e) == 'Error: minimal relations do not preserve all sets of atoms.'

        # check error handling when too few relations are present
        q.minimal_relations = [((2, 3, 2), (3, 2, 3)), ((1, 2), (2, 1))]
        try:
            q.verify_relations(None)
        except Exception as e:
            assert str(e) == 'Error: minimal relations fail to span all sets of atoms.'

        # does not raise exception if we limit the length of atoms to check
        q.verify_relations(upper_length=3)

        # also raise no exceptions if atom length is unlimited but we use sufficient relations
        q.minimal_relations += [((3, 2, 1, 2, 3, 2), (3, 2, 1, 3, 2, 3))]
        q.verify_relations(None)
