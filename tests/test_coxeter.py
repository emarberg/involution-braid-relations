import pytest

from algebra import (
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

        # each i and j for (i, j, m) in `triples` must appear in `generators`
        try:
            CoxeterGraph(triples=[(1, 2, 1)], generators=[3, 4, 5])
        except Exception as e:
            assert str(e).startswith(
                'Invalid input to CoxeterGraph: `generators` must contain i and j for all'
            )
        else:
            assert False

        # if (i, j, m) in `triples` has m == 1 then must have i == j
        try:
            CoxeterGraph(triples=[(1, 2, 1)])
        except Exception as e:
            assert str(e).startswith('Invalid input to CoxeterGraph: `triples` contains invalid')
        else:
            assert False

        # if (i, j, m) in `triples` has i == j then must have m == 1
        try:
            CoxeterGraph(triples=[(1, 1, 2)])
        except Exception as e:
            assert str(e).startswith('Invalid input to CoxeterGraph: `triples` contains invalid')
        else:
            assert False

        # if (i, j, m) is in `triples` then m must be an integer or infinity
        try:
            CoxeterGraph(triples=[(1, 2, 3.0)])
        except Exception as e:
            assert str(e).startswith('Invalid input to CoxeterGraph: `triples` contains invalid')
        else:
            assert False

        # `triples` cannot contain (i, j, m) and (j, i, n) if m != n
        try:
            CoxeterGraph(triples=[(1, 2, 3), (2, 1, 2)])
        except Exception as e:
            assert str(e).startswith(
                'Invalid input to CoxeterGraph: `triples` contains inconsistent'
            )
        else:
            assert False

        # if (i, j) is in `star` then i and j must be generators
        try:
            CoxeterGraph(triples=[(1, 2, 3)], star=[(3, 4)])
        except Exception as e:
            assert str(e).startswith('Invalid input to CoxeterGraph: `star` contains pairs which')
        else:
            assert False

        # `star` cannot contain (i, j) and (j, k) if i != k
        try:
            CoxeterGraph(triples=[(1, 2, 3)], star=[(1, 2), (2, 2)])
        except Exception as e:
            assert str(e).startswith('Invalid input to CoxeterGraph: `star` contains inconsistent')
        else:
            assert False

        # `star` must define a group automorphism
        try:
            CoxeterGraph(triples=[(1, 2, 3), (2, 3, 4)], star=[(1, 3)])
        except Exception as e:
            assert str(e).startswith(
                'Invalid input to CoxeterGraph: `star = [(1, 3)]` does not define automorphism'
            )
        else:
            assert False


class TestRoot:
    pass


class TestRootTransform:
    pass


class TestCoxeterTransform:
    pass


class TestCoxeterWord:
    pass
