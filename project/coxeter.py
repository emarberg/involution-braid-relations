import numpy as np

from project.utils import (
    reverse_tuple,
    get_braids,
    InvalidInputException,
    OperatorException,
    NumberMixin,
    VectorMixin
)

from project.algebra import (
    RationalNumber,
    QuadraticNumber,
    Polynomial,
    Matrix
)


class CoxeterGraph:

    """
    Class for objects representing Coxeter graphs with a distinguisehd graph involution.
    This data classifies twisted Coxeter systems up to isomorphism.
    """

    class EvalBilinearException(Exception):
        def __init__(self):
            super(CoxeterGraph.EvalBilinearException, self).__init__(
                'Error in CoxeterGraph.eval_bilinear: '
                'method currently not supported for parameters (i, j) '
                'when m_ij is not 1, 2, 3, 4, 5, 6, 12, or infinity')

    class GetBraidException(Exception):
        def __init__(self):
            super(CoxeterGraph.GetBraidException, self).__init__(
                'Error in CoxeterGraph.get_braid: m_ij must be finite for input parameters (i, j)')

    def __init__(self, edges, generators=None, star=None):
        """
        Input `edges` should be list of int tuples (i, j, m) where i, j are indices of simple
        generators and m is order of s_i s_j. Triples with m == 1 or 2 can be omitted.
        If (i, j, m) is included then the reverse triple (j, i, m) also may be omitted.

        Input `generators` should be list of integers indexing simple generators of the group.
        If not provided, this list will be formed from the set of numbers i or j in `edges`.
        If i or j is not in `generators` for some (i, j, m) in `edges`, an error will be raised.

        Input `star` should be list of pairs (i, j) such that the involution * : S -> S
        with i^* = j and j^* = i extends to an automorphism of W. If `star` is not given,
        * is defined to be the identity map.
        """
        self._setup_generators(edges, generators)
        self._setup_orders(edges)
        self._setup_star(star)

    def _setup_generators(self, edges, generators):
        """Assign sorted list of simple generator indices to self.generators."""
        generators_from_edges = {t[0] for t in edges} | {t[1] for t in edges}
        try:
            if generators is not None:
                generators = set(generators)
                assert generators_from_edges.issubset(generators)
                self.generators = sorted(generators)
            else:
                self.generators = sorted(generators_from_edges)
        except:
            raise InvalidInputException(self, (edges, generators))

    def _setup_orders(self, edges):
        """Construct dictionary with orders m_ij of products of simple generators."""
        self.orders = {}
        for i, j, m in edges:
            valid_order = (type(m) == int and 1 <= m) or m == np.infty
            valid_order = valid_order and ((m == 1) == (i == j))
            valid_generators = i in self.generators and j in self.generators

            if not (valid_order and valid_generators) or self.orders.get((i, j), m) != m:
                raise InvalidInputException(self, '%s in `tuples`' % str((i, j, m)))
            elif i != j and m != 2:
                self.orders[(i, j)] = m
                self.orders[(j, i)] = m

    def _setup_star(self, star):
        """Construct dictionary with images of the *-involution 'star'."""
        self._star = {}
        if star:
            generators_from_star = {t[0] for t in star} | {t[1] for t in star}
            if not generators_from_star.issubset(set(self.generators)):
                raise InvalidInputException(self, star)
            for i, j in star:
                if self._star.get(i, j) == j and self._star.get(j, i) == i:
                    self._star[i] = j
                    self._star[j] = i
                else:
                    raise InvalidInputException(self, '%s in `star`' % str((i, j)))

        # validate that input `star` encodes automorphism
        for i in self.generators:
            for j in self.generators:
                if self.get_order(i, j) != self.get_order(self.star(i), self.star(j)):
                    raise InvalidInputException(self, star)

    def __eq__(self, other):
        return \
            type(other) == CoxeterGraph and \
            self.generators == other.generators and \
            all(self.star(i) == other.star(i) for i in self.generators) and \
            all(self.get_order(i, j) == other.get_order(i, j)
                for i in self.generators
                for j in self.generators)

    def __hash__(self):
        gens = tuple(self.generators)
        orders = tuple(tuple(self.get_order(i, j) for j in gens) for i in gens)
        star = tuple(self.star(i) for i in self.generators)
        return hash((gens, orders, star))

    def star(self, i):
        if i not in self.generators:
            raise InvalidInputException(self, i, 'star')
        else:
            return self._star.get(i, i)

    def get_braid(self, i, j):
        """Returns alternating tuple (i, j, i, j, ...) of length m_ij."""
        gens = [i, j]
        m = self.get_order(i, j)
        if m == np.infty:
            raise CoxeterGraph.GetBraidException
        else:
            return tuple(gens[t % 2] for t in range(m))

    def get_order(self, i, j):
        """Return order of s_i * s_j in W."""
        if i not in self.generators:
            raise InvalidInputException(self, i, 'get_order')
        elif j not in self.generators:
            raise InvalidInputException(self, j, 'get_order')
        elif i == j:
            return 1
        else:
            return self.orders.get((i, j), 2)

    def get_semiorder(self, i, j, generators_are_fixed):
        """
        Returns twisted involution length of longest element in dihedral group <s_i, s_j>
        with respect to the automorphism fixing s_i and s_j (if generators_are_fixed is True)
        or transposing s_i and s_j (if generators_are_fixed is False).
        """
        m = self.get_order(i, j)
        if m % 2 != 0:
            return (m + 1) // 2
        elif generators_are_fixed:
            return m // 2 + 1
        else:
            return m // 2

    def eval_bilinear(self, i, j):
        """
        Returns -cos(pi/m_ij) as an int, RationalNumber, or QuadraticNumber.
        For simplicity, an error is raised if m_ij is not 1, 2, 3, 4, 5, 6, 12, or infinity,
        even though -cos(pi/m_ij) can be a QuadraticNumber for other values of m_ij.
        """
        m = self.get_order(i, j)
        if m == 1:
            return 1
        elif m == 2:
            return 0
        elif m == 3:
            return RationalNumber(-1, 2)
        elif m == 4:
            return -QuadraticNumber.sqrt(2) / 2
        elif m == 5:
            return -(QuadraticNumber.sqrt(5) + 1) / 4
        elif m == 6:
            return -QuadraticNumber.sqrt(3) / 2
        elif m == 12:
            return -(QuadraticNumber.sqrt(6) + QuadraticNumber.sqrt(2)) / 4
        elif m == np.infty:
            return -1
        else:
            raise CoxeterGraph.EvalBilinearException

    def is_simply_laced(self):
        return set(self.orders.values()).issubset({1, 2, 3})

    def is_crystallographic(self):
        return set(self.orders.values()).issubset({1, 2, 3, 4, 6})

    def is_quadratic(self):
        return set(self.orders.values()).issubset({1, 2, 3, 4, 5, 6, 12, np.infty})

    def __repr__(self):
        entries = [str(self.get_order(i, j)) for i in self.generators for j in self.generators]
        if entries:
            max_len = max(list(map(len, entries)))
            entries = list(map(lambda x: (max_len - len(x)) * ' ' + x, entries))
        n = len(self.generators)
        s = ''
        while entries:
            s += ' '.join(entries[:n]) + '\n'
            entries = entries[n:]
        if s.endswith('\n'):
            s = s[:-1]
        s += '\n\n' + ', '.join([str(i) + '^* = ' + str(self.star(i)) for i in self.generators])
        return s

    def get_max_involution_word_length(self, limit=1000):
        max_len = 0
        e = CoxeterTransform(self)
        while e.right_descents != set(self.generators) and (limit is None or max_len < limit):
            i = next(iter(set(self.generators) - e.right_descents))
            e = e.demazure_conjugate(i)
            max_len += 1
        return max_len

    def get_inverse_atoms(self, relations, start_word):
        """
        Return set of CoxeterTransforms which represent inverses of the elements spanned by
        the input `relations` plus the orindary braid relations.
        """
        start = CoxeterTransform.from_word(self, reverse_tuple(start_word))
        relations = {(reverse_tuple(a), reverse_tuple(b)) for a, b in relations}
        return start.span_by_right_relations(relations)

    def get_half_braid_relations(self):
        """Returns half-braid relations of twisted Coxeter system."""
        relations = set()
        for i in self.generators:
            for j in self.generators:
                # relations encode symmetric data, so can skip half of them
                if i >= j:
                    continue

                if {i, j} == {self.star(i), self.star(j)}:
                    n = self.get_semiorder(i, j, self.star(i) == i)
                    if n < self.get_order(i, j):
                        relations.add(get_braids(i, j, n))
        return relations

    @staticmethod
    def A(n, star=None):  # noqa
        """
        Dynkin diagram labeling is:

            1--2--...--n

        """
        edges = [(i, i + 1, 3) for i in range(1, n)]
        return CoxeterGraph(edges, star=star)

    @staticmethod
    def A_twist(n):  # noqa
        star = [(i, n + 1 - i) for i in range(1, n + 1)]
        return CoxeterGraph.A(n, star=star)

    @staticmethod
    def B_twist(n=2):  # noqa
        assert n == 2
        return CoxeterGraph([(0, 1, 4)], star=[(0, 1)])

    @staticmethod
    def B(n):  # noqa
        """
        Dynkin diagram labeling is:

            0==1--2--...--(n-1)

        """
        assert 2 <= n
        edges = [(i, i + 1, 3) for i in range(1, n - 1)] + [(0, 1, 4)]
        return CoxeterGraph(edges)

    @staticmethod
    def D(n, star=None):  # noqa
        """
        Dynkin diagram labeling is:

               0
               |
            1--2--3--...--(n-1)

        """
        assert 4 <= n
        edges = [(i, i + 1, 3) for i in range(1, n - 1)] + [(0, 2, 3)]
        return CoxeterGraph(edges, star=star)

    @staticmethod
    def D_twist(n):  # noqa
        assert 4 <= n
        star = [(0, 1)] + [(i, i) for i in range(2, n)]
        return CoxeterGraph.D(n, star=star)

    @staticmethod
    def E(n, star=None):  # noqa
        """
        Dynkin diagram labeling is:

                  3
                  |
            1--2--4--5--...--n

        """
        assert n in [6, 7, 8]
        edges = [(1, 2, 3), (2, 4, 3)] + [(i, i + 1, 3) for i in range(3, n)]
        return CoxeterGraph(edges, star=star)

    @staticmethod
    def E_twist(n=6):  # noqa
        assert n == 6
        star = [(1, 6), (2, 5), (3, 3), (4, 4)]
        return CoxeterGraph.E(n, star=star)

    @staticmethod
    def F(n=4, star=None):  # noqa
        assert n == 4
        edges = [(1, 2, 3), (2, 3, 4), (3, 4, 3)]
        return CoxeterGraph(edges, star=star)

    @staticmethod
    def F_twist(n=4):  # noqa
        assert n == 4
        star = [(1, 4), (2, 3)]
        return CoxeterGraph.F(n, star=star)

    @staticmethod
    def G(n=2, star=None):  # noqa
        assert n == 2
        return CoxeterGraph([(1, 2, 6)], star=star)

    @staticmethod
    def G_twist(n=2):  # noqa
        assert n == 2
        return CoxeterGraph.G(n, star=[(1, 2)])

    @staticmethod
    def H(n):  # noqa
        """
        Dynkin diagram labeling is:

            1==2--...--n

        where edge == has label 5.
        """
        assert n in [3, 4]
        if n == 3:
            edges = [(1, 2, 5), (2, 3, 3)]
        elif n == 4:
            edges = [(1, 2, 5), (2, 3, 3), (3, 4, 3)]
        return CoxeterGraph(edges)

    @staticmethod
    def A_tilde(n, star=None):  # noqa
        assert n >= 2
        edges = [(i, i + 1, 3) for i in range(1, n)] + [(n, n + 1, 3), (n + 1, 1, 3)]
        return CoxeterGraph(edges, star=star)

    @staticmethod
    def A_tilde_twist(n):  # noqa
        star = [(i, n + 2 - i) for i in range(1, n + 2)]
        return CoxeterGraph.A_tilde(n, star)

    @staticmethod
    def A_tilde_flip(n):  # noqa
        assert n % 2 != 0
        star = [(i, n + 1 - i) for i in range(1, n + 1)]
        return CoxeterGraph.A_tilde(n, star)

    @staticmethod
    def A_tilde_rotate(n):  # noqa
        assert n % 2 != 0
        r = (n + 1) // 2
        star = [(i, (i - 1 + r) % (n + 1) + 1) for i in range(1, n + 2)]
        return CoxeterGraph.A_tilde(n, star)

    @staticmethod
    def B_tilde(n, star=None):  # noqa
        """
        Dynkin diagram labeling is:

                           n
                           |
            0==1--2--...--(n-2)--(n-1)

        """
        assert n >= 3
        edges = [(i, i + 1, 3) for i in range(1, n - 1)] + [(0, 1, 4)] + [(n - 2, n, 3)]
        return CoxeterGraph(edges, star=star)

    @staticmethod
    def B_tilde_twist(n):  # noqa
        star = [(n - 1, n)]
        return CoxeterGraph.B_tilde(n, star)

    @staticmethod
    def C_tilde(n, star=None):  # noqa
        """
        Dynkin diagram labeling is:

            0==1--2--...--(n-1)==n

        """
        assert n >= 2
        edges = [(i, i + 1, 3) for i in range(1, n - 1)] + [(0, 1, 4)] + [(n - 1, n, 4)]
        return CoxeterGraph(edges, star=star)

    @staticmethod
    def C_tilde_twist(n):  # noqa
        star = [(i, n - i) for i in range(n + 1)]
        return CoxeterGraph.C_tilde(n, star)

    @staticmethod
    def D_tilde(n, star=None):  # noqa
        """
        Dynkin diagram labeling is:

               0           n
               |           |
            1--2--3--...--(n-2)--(n-1)

        """
        assert n >= 4
        edges = [(i, i + 1, 3) for i in range(1, n - 1)] + [(0, 2, 3)] + [(n - 2, n, 3)]
        return CoxeterGraph(edges, star=star)

    @staticmethod
    def D_tilde_twist(n):  # noqa
        star = [(i, n - i) for i in range(n + 1)]
        return CoxeterGraph.D_tilde(n, star)

    @staticmethod
    def D_tilde_small_twist(n):  # noqa
        return CoxeterGraph.D_tilde(n, star=[(0, 1), (n - 1, n)])

    @staticmethod
    def D_tilde_half_twist(n):  # noqa
        return CoxeterGraph.D_tilde(n, star=[(0, 1)])

    @staticmethod
    def E_tilde(n, star=None):  # noqa
        """
        Dynkin diagram labelings are:

                  *
                  |
                  3
                  |
            1--2--4--5--6

                     3
                     |
            *--1--2--4--5--6--7

                  3
                  |
            1--2--4--5--6--7--8--*

        """
        assert n in [6, 7, 8]
        edges = [('1', '2', 3), ('2', '4', 3)] + [(str(i), str(i + 1), 3) for i in range(3, n)]
        if n == 6:
            edges += [('3', '*', 3)]
        elif n == 7:
            edges += [('*', '1', 3)]
        elif n == 8:
            edges += [('8', '*', 3)]
        return CoxeterGraph(edges, star=star)

    @staticmethod
    def E_tilde_twist(n):  # noqa
        assert n in [6, 7]
        if n == 6:
            star = [('1', '6'), ('2', '5')]
        elif n == 7:
            star = [('*', '7'), ('1', '6'), ('2', '5')]
        return CoxeterGraph.E_tilde(n, star)

    @staticmethod
    def F_tilde(n=4):  # noqa
        """
        Dynkin diagram labeling is:

            1--2==3--4--5

        """
        assert n == 4
        edges = [(1, 2, 3), (2, 3, 4), (3, 4, 3), (4, 5, 3)]
        return CoxeterGraph(edges)

    @staticmethod
    def G_tilde(n=2):  # noqa
        """
        Dynkin diagram labeling is:

            1≡≡2--3

        """
        assert n == 2
        return CoxeterGraph([(1, 2, 6), (2, 3, 3)])


class CoxeterVector(VectorMixin, NumberMixin):

    """
    Class for objects encoding individual vectors in the geometric representation
    of a provided Coxeter system, i.e., linear combinations of simple roots in the
    root system associated to the input Coxeter group. Coefficients in these vectors
    should be ints, RationalNumbers, QuadraticNumbers, or Polynomials.
    """

    @property
    def coefficients(self):
        return self._coefficients

    @coefficients.setter
    def coefficients(self, value):
        self._coefficients = value

    def __init__(self, coxeter_graph, index=None, coeff=1):
        """
        If index is None, returns zero vector. Otherwise, index should be an element of
        coxeter_graph.generators.
        """
        self.graph = coxeter_graph
        if index is None or coeff == 0:
            self._coefficients = {}
        elif index in coxeter_graph.generators:
            self._coefficients = {index: coeff}
        else:
            raise InvalidInputException(self, index)

    def is_comparable(self, other):
        if type(other) == int and other == 0:
            return True
        elif type(other) == CoxeterVector and other.graph == self.graph:
            return True
        else:
            return False

    def eval_bilinear(self, other):
        """
        Returns value of standard bilinear form evaluated at self and other, where other is another
        CoxeterVector object (or a generator index which will be converted to a CoxeterVector).
        """
        try:
            is_generator = (other in self.graph.generators)
        except:
            is_generator = False

        if is_generator:
            other = CoxeterVector(self.graph, other)

        if type(other) == CoxeterVector and other.graph == self.graph:
            ans = 0
            for i, u in self:
                for j, v in other:
                    ans += u * v * self.graph.eval_bilinear(i, j)
            return ans
        else:
            raise InvalidInputException(self, other, 'eval_bilinear')

    def reflect(self, i):
        """Given index i, returns (self - 2 <self, alpha_i> alpha_i)."""
        if i in self.graph.generators:
            v = 2 * self.eval_bilinear(i)
            return self - CoxeterVector(self.graph, i, v)
        else:
            raise InvalidInputException(self, i, 'reflect')

    def is_constant(self):
        """Returns true if all coefficients are constant polynomials."""
        return not any(type(v) == Polynomial and not v.is_constant() for i, v in self)

    def is_positive(self):
        return (not self.is_zero()) and all(0 <= v for _, v in self)

    def is_negative(self):
        return (not self.is_zero()) and all(v <= 0 for _, v in self)

    def is_valid(self):
        """
        A CoxeterVector object, which is a linear combination sum_i c_i alpha_i of simple roots,
        is 'invalid' if it is zero or if for some i and j it holds that -c_i and c_j are
        both polynomials whose coefficients are all nonnegative and whose constant
        coefficients are positive.

        An invalid linear combination does not specialize to a positive or negative root
        for any choice of nonnegative variable inputs.
        """
        if self.is_zero():
            return False
        if any(v < 0 for _, v in self) and any(0 < v for _, v in self):
            return False
        return True

    def set_variable(self, variable, value):
        new = CoxeterVector(self.graph)
        for i, v in self:
            if type(v) == Polynomial:
                v = v.set_variable(variable, value)
            if v != 0:
                new.coefficients[i] = v
        return new

    def __add__(self, other):
        if other == 0:
            other = CoxeterVector(self.graph)
        if type(other) == CoxeterVector and self.graph == other.graph:
            indices = set(self.coefficients.keys()) | set(other.coefficients.keys())
            new = CoxeterVector(self.graph)
            new.coefficients = {i: self[i] + other[i] for i in indices if (self[i] + other[i]) != 0}
            return new
        else:
            raise OperatorException(self, other, '__add__')

    def __mul__(self, other):
        """Implements scalar multiplication."""
        new = CoxeterVector(self.graph)
        if other != 0:
            new.coefficients = {i: other * self[i] for i in self.coefficients}
        return new

    def __hash__(self):
        return hash((self.graph,) + tuple(self[i] for i in self.graph.generators))

    def __pow__(self, exponent):
        """Overrides NumberMixin.__pow__ implementation."""
        raise NotImplementedError  # pragma: no cover

    def __repr__(self):
        return super(CoxeterVector, self).__repr__()[1:-1]  # slicing removes '(' and ')' characters

    @classmethod
    def get_index_repr(cls, index):
        return 'alpha_' + str(index)


class TransformMixin:

    """
    Mixin class for methods shared by the PartialTransform and CoxeterTransform classes below,
    implementing common operations for linear transforms of (subspaces of) the goemetric
    representation of some Coxeter system. Most such methods are pretty trivial, with the
    notable exception of the __mul__ and __rmul__ functions.
    """

    @property
    def graph(self):
        """Returns CoxeterGraph."""
        raise NotImplementedError  # pragma: no cover

    @property
    def sigma(self):
        """Returns dict with keys in self.graph.generators, and whose values are CoxeterVectors."""
        raise NotImplementedError  # pragma: no cover

    @classmethod
    def identity(cls, coxeter_graph):
        """Returns object corresponding to identity map on geometric representation."""
        sigma = {i: CoxeterVector(coxeter_graph, i) for i in coxeter_graph.generators}
        return cls(coxeter_graph, sigma)

    def __eq__(self, other):
        return isinstance(other, TransformMixin) and \
            other.graph == self.graph and other.sigma == self.sigma

    def __hash__(self):
        return hash((self.graph,) + tuple(sorted(self.sigma.items())))

    def __len__(self):
        return len(self.sigma)

    def __getitem__(self, i):
        return self.sigma[i]

    def __contains__(self, i):
        return i in self.sigma

    def __iter__(self):
        return self.sigma.__iter__()

    def values(self):
        return self.sigma.values()

    def is_identity(self):
        if sorted(self.sigma.keys()) == self.graph.generators:
            return all(self.sigma[i] == CoxeterVector(self.graph, i) for i in self.graph.generators)
        return False

    def __mul__(self, j):
        """Returns composition of transform with image of s_j in geometric representation."""
        if j not in self.graph.generators:
            raise InvalidInputException(self, j, '__mul__')
        new = {}
        for i in self.sigma:
            root = CoxeterVector(self.graph, i).reflect(j)
            for k, v in root:
                if k not in self.sigma:
                    raise InvalidInputException(self, j, '__mul__')
                new[i] = new.get(i, 0) + self.sigma[k] * v
        return self.__class__(self.graph, new)

    def __rmul__(self, j):
        """Returns composition of image of s_j in geometric representation with transform."""
        if j not in self.graph.generators:
            raise InvalidInputException(self, j, '__rmul__')
        new = {}
        for i in self.sigma:
            new[i] = self.sigma[i].reflect(j)
        return self.__class__(self.graph, new)

    def __repr__(self):
        s = '{\n'
        for i in self.graph.generators:
            if i in self.sigma:
                s += '  alpha_%s -> %s\n' % (i, self.sigma[i])
            else:
                s += '  alpha_%s -> ?\n' % i
        s += '}\n'
        return s


class PartialTransform(TransformMixin):

    """
    Class for objects which represent linear transformations T : U -> V where V is the
    geometric representation of a Coxeter system and U is a subspace of V spanned by some
    subset of simple roots.
    """

    def __init__(self, coxeter_graph, sigma={}):
        """
        Input `sigma` should be dict whose keys are elements of coxeter_graph.generators and
        whose values are CoxeterVectors whose graph fields match coxeter_graph.
        """
        if all(i in coxeter_graph.generators for i in sigma) and \
           all(type(r) == CoxeterVector and r.graph == coxeter_graph for r in sigma.values()):
            self._graph = coxeter_graph
            self._sigma = sigma.copy()
        else:
            raise InvalidInputException(self, (coxeter_graph, sigma))

    @property
    def sigma(self):
        return self._sigma

    @property
    def graph(self):
        return self._graph

    def copy(self):
        return PartialTransform(self.graph, self.sigma)

    def __setitem__(self, i, val):
        if i in self.graph.generators and type(val) == CoxeterVector and val.graph == self.graph:
            self._sigma[i] = val
        else:
            raise InvalidInputException(self, (i, val), '__setitem__')

    def is_constant(self):
        return all(r.is_constant() for r in self.sigma.values())

    def is_complete(self):
        return set(self.sigma.keys()) == set(self.graph.generators)

    def is_positive(self):
        return all(r.is_positive() for r in self.sigma.values())

    def get_variables(self):
        """Return set of indices for variables occuring in our (polynomials) coefficients."""
        ans = set()
        for i in self.sigma:
            for j, f in self.sigma[i]:
                if type(f) == Polynomial:
                    ans = ans | f.get_variables()
        return ans

    def determinant(self):
        """
        Returns determinant of the linear transformation defined by PartialTransform,
        if its domain is the vector space spanned by all simple roots, and if all
        entries in its matrix (in the basis of simple roots) are linear or constant
        polynomails in the same variable. Returns None if these conditions do not hold.
        """
        if not self.is_complete():
            return None

        matrix = [
            [self.sigma[i][j] for j in self.graph.generators]
            for i in self.graph.generators
        ]
        # check if matrix has entries which are non-linear polynomials
        if any(type(x) == Polynomial and x.degree() > 1 for row in matrix for x in row):
            return None

        variables = self.get_variables()
        # entries of matrix involve more than one variable
        if len(variables) > 1:
            return None
        elif len(variables) == 1:
            variable = next(iter(variables))
        else:
            variable = None

        return Matrix(matrix, variable).determinant()

    def get_relations(self):
        """
        Returns set of pairs of form ((s, t, s, ...), (t, s, t, ...)) where s, t belong
        to self.graph.generators such that (1) the PartialTransform maps {alpha_s, alpha_t}
        to {alpha_s^*, alpha_t^*} and (2) the length of each tuple in the pair is
        self.graph.get_semiorder(s, t, b), where b is True or False according to
        whether alpha_s -> alpha_s^* or alpha_s -> alpha_t^*, and this length is < m(s,t).
        """
        relations = set()
        for i in self.sigma:
            for j in self.sigma:
                # relations encode symmetric data, so can skip half of them
                if i >= j:
                    continue

                a = CoxeterVector(self.graph, self.graph.star(i))
                b = CoxeterVector(self.graph, self.graph.star(j))

                if {self.sigma[i], self.sigma[j]} == {a, b}:
                    n = self.graph.get_semiorder(i, j, self.sigma[i] == a)
                    relations.add(get_braids(i, j, n))
        return relations


class CoxeterTransform(TransformMixin):

    """
    Class for objects which encode the (invertible) linear transformations which arise as images
    of elements of a Coxeter group in its geometric representation.
    """

    def __init__(self, coxeter_graph, sigma=None):
        """
        Input `sigma` should be a dict whose key set is set(oxeter_graph.generators) and
        whose values are CoxeterVectors with constant coefficients. The constructor does
        not check if the provided inputs define a linear transformation from a group element.
        """
        if sigma:
            keys_valid = set(sigma.keys()) == set(coxeter_graph.generators)
            roots_valid = all(
                type(r) == CoxeterVector and r.graph == coxeter_graph
                for r in sigma.values()
            )
            nonvariable = all(r.is_constant() for r in sigma.values())
            if keys_valid and roots_valid and nonvariable:
                self._sigma = sigma.copy()
            else:
                raise InvalidInputException(self, sigma)
        else:
            self._sigma = {i: CoxeterVector(coxeter_graph, i) for i in coxeter_graph.generators}
        self._graph = coxeter_graph
        self._right_descents = None
        self._minimal_right_descent = None
        self._minimal_reduced_word = None

    @property
    def sigma(self):
        return self._sigma

    @property
    def graph(self):
        return self._graph

    @classmethod
    def from_word(cls, coxeter_graph, word):
        g = CoxeterTransform(coxeter_graph)
        for i in word:
            g *= i
        return g

    def copy(self):
        other = CoxeterTransform(self.graph, self.sigma)
        other._right_descents = self._right_descents
        other._minimal_right_descent = self._minimal_right_descent
        other._minimal_reduced_word = self._minimal_reduced_word
        return other

    def _unset_cached_properties(self):
        self._right_descents = None
        self._minimal_right_descent = None
        self._minimal_reduced_word = None

    def __setitem__(self, i, value):
        if not (value.is_constant() and type(value) == CoxeterVector and value.graph == self.graph):
            raise InvalidInputException(self, (i, value), '__setitem__')
        elif i not in self.graph.generators:
            raise InvalidInputException(self, (i, value), '__setitem__')
        else:
            self._sigma[i] = value
            self._unset_cached_properties()

    @property
    def right_descents(self):
        if self._right_descents is None:
            self._right_descents = {i for i in self.sigma if self.sigma[i].is_negative()}
        return self._right_descents

    @property
    def minimal_right_descent(self):
        if self._minimal_right_descent is None and self.right_descents:
            self._minimal_right_descent = min(self.right_descents)
        return self._minimal_right_descent

    @property
    def minimal_reduced_word(self):
        """Returns lexicographically minimal reduced word for self, read right to left."""
        if self._minimal_reduced_word is None:
            i = self.minimal_right_descent
            if i is None:
                self._minimal_reduced_word = ()
            else:
                self._minimal_reduced_word = (self * i).minimal_reduced_word + (i,)
        return self._minimal_reduced_word

    def get_inverse(self):
        return self.from_word(self.graph, reverse_tuple(self.minimal_reduced_word))

    def multiply_up(self, word):
        """
        With u = self and v = CoxeterTransform(word), returns the
        product uv if ell(uv) = ell(u) + ell(v) and None otherwise.
        Input `word` should be list or tuple of integers.
        """
        if word and word[0] not in self.right_descents:
            return (self * word[0]).multiply_up(word[1:])
        elif word and word[0] in self.right_descents:
            return None
        else:
            return self

    def multiply_down(self, word):
        """
        With u = self and v = CoxeterTransform(word), returns the
        product uv if ell(uv) = ell(u) - ell(v) and None otherwise.
        Input `word` should be iterable over list or tuple of integers.
        """
        if word and word[0] not in self.right_descents:
            return None
        elif word and word[0] in self.right_descents:
            return (self * word[0]).multiply_down(word[1:])
        else:
            return self

    def toggle_right_segment(self, from_segment, to_segment):
        y = self.multiply_down(reverse_tuple(from_segment))
        if y is not None:
            y = y.multiply_up(to_segment)
        return y

    def span_by_right_relations(self, relations):
        """
        Given set `relations` consisting of pairs (a, b) where a and b are tuples of
        simple generators, and CoxeterTransform `start`, returns equivalence class of
        CoxeterTransforms containing `start` spanned by the relations defined by setting x ~ y
        when x has a reduced word ending with a, y has a reduced word ending with b, and
        xa = yb, for some (a, b) in `relations`.
        """
        generated = set()
        to_add = {self}
        while to_add:
            to_add.difference_update({None})
            to_add.difference_update(generated)
            generated.update(to_add)
            next_add = set()
            for x in to_add:
                for word_a, word_b in relations:
                    y = x.toggle_right_segment(word_a, word_b)
                    z = x.toggle_right_segment(word_b, word_a)
                    next_add.update({y, z})
            to_add = next_add
        return generated

    def demazure_conjugate(self, i):
        if i not in self.graph.generators:
            raise InvalidInputException(self, i, 'demazure_conjugate')
        j = self.graph.star(i)
        ans = j * self * i
        if ans == self:
            ans = self * i
        return ans


class CoxeterWord:

    """
    Class for objects representing finite sequences of simple generators in a Coxeter group,
    i.e., factorizations of Coxeter group elements into simple generators.
    """

    def __init__(self, coxeter_graph, word=()):
        self.graph = coxeter_graph
        self.word = ()
        self.left_action = CoxeterTransform.identity(coxeter_graph)
        self.right_action = CoxeterTransform.identity(coxeter_graph)
        self.is_reduced = True
        for i in word:
            self.extend_right(i)

    @property
    def left_descents(self):
        """Returns left descent set of group element represented by a CoxeterWord."""
        return self.right_action.right_descents

    @property
    def right_descents(self):
        """Returns right descent set of group element represented by a CoxeterWord."""
        return self.left_action.right_descents

    @property
    def minimal_reduced_word(self):
        """Returns lexicographically minimal reduced word for a CoxeterWord, read left to right."""
        return reverse_tuple(self.right_action.minimal_reduced_word)

    def __eq__(self, other):
        return type(other) == CoxeterWord and self.graph == other.graph and self.word == other.word

    def __hash__(self):
        return hash((self.graph, self.word))

    def __len__(self):
        return len(self.word)

    def to_involution(self):
        new = CoxeterWord(self.graph)
        for i in self.word:
            alpha = CoxeterVector(self.graph, self.graph.star(i))
            if new.left_action[i] not in [alpha, -alpha]:
                new.extend_left(self.graph.star(i))
            new.extend_right(i)
        return new

    def get_reduced_words(self):
        """Computes and returns set of all reduced words for group element of a CoxeterWord."""
        ans = set()
        to_add = {self.word}
        while to_add:
            ans.update(to_add)
            next_to_add = set()
            for word in to_add:
                for i in range(len(word) - 1):
                    s, t = word[i:i + 2]
                    m = self.graph.get_order(s, t)
                    if word[i:i + m] == self.graph.get_braid(s, t):
                        new_word = word[:i] + self.graph.get_braid(t, s) + word[i + m:]
                        next_to_add.add(new_word)
            to_add = next_to_add - ans
        return ans

    def copy(self):
        other = CoxeterWord(self.graph)
        other.word = self.word
        other.left_action = self.left_action.copy()
        other.right_action = self.right_action.copy()
        other.is_reduced = self.is_reduced
        return other

    def extend_right(self, j):
        """Append j to self.word and update other fields."""
        self.word = self.word + (j,)
        self.is_reduced = self.is_reduced and (j not in self.right_descents)
        self.left_action = self.left_action * j
        self.right_action = j * self.right_action

    def extend_left(self, j):
        """Prepend j to self.word and update other fields."""
        self.word = (j,) + self.word
        self.is_reduced = self.is_reduced and (j not in self.left_descents)
        self.left_action = j * self.left_action
        self.right_action = self.right_action * j

    def _is_valid_argument(self, other):
        return \
            (type(other) == CoxeterWord and other.graph == self.graph) or \
            (type(other) == int and other in self.graph.generators)

    def __mul__(self, other):
        """Input `other` must be int or CoxeterWord."""
        if not self._is_valid_argument(other):
            raise InvalidInputException(self, type(other), '__mul__')

        if type(other) == int:
            new = self.copy()
            new.extend_right(other)
        else:
            new = self.copy()
            for i in other.word:
                new.extend_right(i)
        return new

    def __rmul__(self, other):
        """Input `other` must be int or CoxeterWord."""
        if not self._is_valid_argument(other):
            raise InvalidInputException(self, type(other), '__rmul__')

        if type(other) == int:
            new = self.copy()
            new.extend_left(other)
        else:
            new = self.copy()
            for i in reverse_tuple(other.word):
                new.extend_left(i)
        return new

    def __repr__(self):
        letters = map(str, self.word)
        return (not self.is_reduced) * 'un' + 'reduced word (' + ', '.join(letters) + ')'
