import numpy as np

from project.utils import (
    InvalidInputException,
    ZeroDivisionException,
    OperatorException,
    CannotFactorException,
    OperatorMixin,
    NumberMixin,
    VectorMixin
)


class RationalNumber(OperatorMixin, NumberMixin):

    """
    Class for objects representing rational numbers, with
    the usual field operations and standard total ordering.
    """

    @classmethod
    def one(cls):
        return RationalNumber(1)

    def __init__(self, p=0, q=1):
        """Returns RationalNumber with value p/q. Inputs must be ints or RationalNumbers."""
        if type(p) == type(q) == RationalNumber:
            p, q = (p.numerator * q.denominator, q.numerator * p.denominator)
        elif type(p) == int and type(q) == RationalNumber:
            p, q = (p * q.denominator, q.numerator)
        elif type(q) == int and type(p) == RationalNumber:
            p, q = (p.numerator, q * p.denominator)

        if type(p) == int and type(q) == int and q != 0:
            p, q = self.reduce(p, q)
            self.numerator = p
            self.denominator = q
        elif q == 0:
            raise ZeroDivisionException(self)
        else:
            raise InvalidInputException(self, (type(p), type(q)))

    @classmethod
    def reduce(cls, p, q):
        """Given integers p, q, returns pair (p/d, q/d) where d = gcd(p, q)."""
        if q < 0:
            q = q * -1
            p = p * -1
        d = np.math.gcd(p, q)
        return p // d, q // d

    def is_integer(self):
        return self.denominator == 1

    def is_rational(self):
        return True

    def __hash__(self):
        """If self is an integer, its hash is the same as the corresponding int."""
        if self.is_integer():
            return hash(self.numerator)
        else:
            return hash((self.numerator, self.denominator))

    def eq__int(self, other):
        return self.numerator == other and self.denominator == 1

    def eq__rational_number(self, other):
        return self.numerator == other.numerator and self.denominator == other.denominator

    def eq__quadratic_number(self, other):
        return QuadraticNumber(self) == other

    def eq__polynomial(self, other):
        return Polynomial(self) == other

    def lt__int(self, other):
        return self.numerator < other * self.denominator

    def lt__rational_number(self, other):
        return self.numerator * other.denominator < other.numerator * self.denominator

    def lt__quadratic_number(self, other):
        return QuadraticNumber(self) < other

    def lt__polynomial(self, other):
        return Polynomial(self) < other

    def __le__(self, other):
        if type(other) == Polynomial:
            return Polynomial(self) <= other
        else:
            return self == other or self < other

    def add__int(self, other):
        p, q = self.numerator + other * self.denominator, self.denominator
        return RationalNumber(p, q)

    def add__rational_number(self, other):
        p = self.numerator * other.denominator + other.numerator * self.denominator
        q = self.denominator * other.denominator
        return RationalNumber(p, q)

    def add__quadratic_number(self, other):
        return other.add__rational_number(self)

    def add__polynomial(self, other):
        return other.add__rational_number(self)

    def mul__int(self, other):
        p, q = self.numerator * other, self.denominator
        return RationalNumber(p, q)

    def mul__rational_number(self, other):
        p, q = self.numerator * other.numerator, self.denominator * other.denominator
        return RationalNumber(p, q)

    def mul__quadratic_number(self, other):
        return other.mul__rational_number(self)

    def mul__polynomial(self, other):
        return other.mul__rational_number(self)

    def truediv__int(self, other):
        return RationalNumber(self, other)

    def truediv__rational_number(self, other):
        return RationalNumber(self, other)

    def truediv__quadratic_number(self, other):
        return QuadraticNumber(self).truediv__quadratic_number(other)

    def truediv__polynomial(self, other):
        if other.is_constant():
            return self / other.get_constant_part()
        else:
            raise OperatorException(self, other, '__truediv__')

    def __rtruediv__(self, other):
        return other * RationalNumber(1, self)

    def __repr__(self):
        if self.numerator == 0:
            return '0'
        elif self.denominator != 1:
            return str(self.numerator) + '/' + str(self.denominator)
        else:
            return str(self.numerator)


class PrimeFactorization:

    """
    Class for objects to store and access the unique factorization
    of a non-zero integer into prime-power factors.
    """

    _factorization_cache = {}  # updated whenever a new PrimeFactorization is created

    def __init__(self, i=1):
        """Input i must be nonzero int."""
        if type(i) == int and i != 0:
            self.factorization = PrimeFactorization.get_prime_factorization(i)
            self.n = i
        else:
            raise InvalidInputException(self, type(i))

    def __getitem__(self, p):
        return self.factorization.get(p, 0)

    def __eq__(self, other):
        return type(other) == PrimeFactorization and self.n == other.n

    def __lt__(self, other):
        return self.n < other.n

    def __hash__(self):
        return self.n

    def __mul__(self, other):
        if type(other) == PrimeFactorization:
            ans = PrimeFactorization()
            ans.n = self.n * other.n
            factors = set(self.factorization.keys()) | set(other.factorization.keys())
            ans.factorization = {p: self[p] + other[p] for p in factors if p != -1}
            if ans.n < 0:
                ans.factorization[-1] = 1
            return ans
        else:
            raise Exception('Cannot multiply PrimeFactorization with %s' % type(other))

    def __repr__(self):
        return repr(self.factorization)

    def get_square_free_part(self):
        """
        Returns PrimeFactorization of integer n with smallest absolute value
        such that self/n is a square-free positive integer.
        """
        ans = PrimeFactorization()
        ans.factorization = {p: e % 2 for p, e in self.factorization.items() if e % 2 != 0}
        ans.n = 1
        for p in ans.factorization:
            ans.n *= p
        return ans

    def get_truncated_square_root(self):
        """Returns sqrt(self/n) where self.get_square_free_part() == PrimeFactorization(n)."""
        r = 1
        for p, e in self.factorization.items():
            r *= p**(e // 2)
        return r

    @classmethod
    def get_divisor_exponent(cls, n, p):
        """Return maximum positive integer e such that p**e divides n."""
        if n % p != 0:
            return 0
        if abs(n) == abs(p):
            return 1
        else:
            return 1 + cls.get_divisor_exponent(n // p, p)

    @classmethod
    def get_prime_factorization(cls, i):
        """
        Returns dict with prime integer keys (with -1 considered to be prime) and positive
        integer values which gives the prime factorization of the input integer i.
        """
        if i in cls._factorization_cache:
            return cls._factorization_cache[i].copy()
        else:
            i_input = i
            factorization = {}
            possible_factors = range(2, abs(i) + 1)
            while possible_factors:
                p = possible_factors[0]
                e = cls.get_divisor_exponent(i, p)
                if e != 0:
                    factorization[p] = e
                    i = i // p**e
                possible_factors = [a for a in possible_factors[1:] if a <= abs(i) and a % p != 0]
            if i_input < 0:
                factorization[-1] = 1
            cls._factorization_cache[i_input] = factorization
            return factorization.copy()


class QuadraticNumber(VectorMixin, OperatorMixin, NumberMixin):

    """
    Class for objects representing rational linear combinations of square roots of integers,
    with the usual field operations and (for real-valued objects) the usual total ordering
    on real numbers. QuadraticNumbers are implemented as vectors spanned by
    PrimeFactorizations with RationalNumber coefficients.
    """

    class ImaginaryComparisonException(Exception):
        def __init__(self, diff):
            super(QuadraticNumber.ImaginaryComparisonException, self).__init__(
                'Cannot compare quadratic numbers with non-real difference %s' % diff)

    class IndeterminateComparisonException(Exception):
        def __init__(self, a, b):
            super(QuadraticNumber.IndeterminateComparisonException, self).__init__(
                'Cannot determine inequality %s < %s' % (a, b))

    @property
    def coefficients(self):
        return self._coefficients

    @coefficients.setter
    def coefficients(self, value):
        self._coefficients = value

    @classmethod
    def one(cls):
        return QuadraticNumber(1)

    def __init__(self, i=0):
        """Input i must be int or RationalNumber."""
        if type(i) == int:
            i = RationalNumber(i)
        if type(i) == RationalNumber:
            self._coefficients = {}
            if i != 0:
                self._coefficients[PrimeFactorization(1)] = i
        else:
            raise InvalidInputException(self, type(i))

    def is_comparable(self, other):
        return type(other) in [int, RationalNumber, QuadraticNumber, Polynomial]

    def __repr__(self):
        if self.is_rational():
            return str(self.get_rational_part())
        else:
            return super(QuadraticNumber, self).__repr__()

    def __hash__(self):
        """If self is rational, its hash is the same as the corresponding RationalNumber."""
        if self.is_rational():
            return hash(self.get_rational_part())
        else:
            return hash(tuple(sorted(self.coefficients.items())))

    def lt__int(self, other):
        return self.lt__quadratic_number(QuadraticNumber(other))

    def lt__rational_number(self, other):
        return self.lt__quadratic_number(QuadraticNumber(other))

    def lt__quadratic_number(self, other):
        """
        Tries to compare self and QuadraticNumber other. Raises exception if self - other is
        not real-valued or too complicated for us to compute its sign, and otherwise returns
        True/False if the difference is negative/positive.

        Our approach to determining (whether we can compute) the sign of self - other is overly
        simplistic: we raise an exception to indicate failure if the difference is a sum of too
        many distinct radicals. There are efficient algorithms that can solve this problem in
        general (see http://cs.smith.edu/~orourke/TOPP/P33.html), but our simpler implementation
        is good enough for the applications of this program.
        """
        diff = other - self
        if diff.is_rational():
            return 0 < diff[PrimeFactorization(1)]
        elif not diff.is_real():
            raise QuadraticNumber.ImaginaryComparisonException(diff)
        elif all(0 < c for c in diff.coefficients.values()):
            return True
        elif all(c < 0 for c in diff.coefficients.values()):
            return False

        positive_part, negative_part, boolean = diff._get_normalized_parts()
        if boolean is not None:
            return boolean
        elif max(len(positive_part), len(negative_part)) > 2:
            raise QuadraticNumber.IndeterminateComparisonException(negative_part, positive_part)
        else:
            return (negative_part**2).lt__quadratic_number(positive_part**2)

    def _get_normalized_parts(self):
        """
        Decompose QuadraticNumber into positive and negative parts. If either
        or these is rational, renormalize both values and attempt to compare.
        Returns renormalized values and output of comparison, which may be None.
        """
        positive_part, negative_part = self.decompose()
        boolean = None
        if positive_part.is_rational():
            negative_part = negative_part / positive_part
            positive_part = QuadraticNumber(1)
            if any(pf.n * coeff**2 > 1 for pf, coeff in negative_part.coefficients.items()):
                boolean = False
        elif negative_part.is_rational():
            positive_part = positive_part / negative_part
            negative_part = QuadraticNumber(1)
            if any(pf.n * coeff**2 > 1 for pf, coeff in positive_part.coefficients.items()):
                boolean = True
        return positive_part, negative_part, boolean

    def lt__polynomial(self, other):
        return Polynomial(self).lt__polynomial(other)

    def __le__(self, other):
        if type(other) == Polynomial:
            return Polynomial(self) <= other
        else:
            return self == other or self < other

    def add__int(self, other):
        return self.add__quadratic_number(QuadraticNumber(other))

    def add__rational_number(self, other):
        return self.add__quadratic_number(QuadraticNumber(other))

    def add__quadratic_number(self, other):
        keys = set(self.coefficients.keys()) | set(other.coefficients.keys())
        new = QuadraticNumber()
        new.coefficients = {i: self[i] + other[i] for i in keys if (self[i] + other[i]) != 0}
        return new

    def add__polynomial(self, other):
        return other.add__quadratic_number(self)

    def mul__int(self, other):
        return self.mul__quadratic_number(QuadraticNumber(other))

    def mul__rational_number(self, other):
        return self.mul__quadratic_number(QuadraticNumber(other))

    def mul__quadratic_number(self, other):
        new = QuadraticNumber()
        for factors_1, coeff_1 in other:
            for factors_2, coeff_2 in self:
                factors = factors_1 * factors_2
                square_free = factors.get_square_free_part()
                coeff = coeff_1 * coeff_2 * factors.get_truncated_square_root()
                new[square_free] += coeff
        return new

    def mul__polynomial(self, other):
        return other.mul__quadratic_number(self)

    def truediv__int(self, other):
        return self.mul__rational_number(RationalNumber(1, other))

    def truediv__rational_number(self, other):
        return self.mul__rational_number(RationalNumber(1, other))

    def truediv__quadratic_number(self, other):
        if self == 0:
            return QuadraticNumber(0)
        elif other.is_rational():
            return self / other.get_rational_part()
        else:
            conjugate = other.conjugate()
            return (self * conjugate) / (other * conjugate)

    def truediv__polynomial(self, other):
        if other.is_constant():
            return self / other.get_constant_part()
        else:
            raise OperatorException(self, other, '__truediv__')

    def __rtruediv__(self, other):
        conjugate = self.conjugate()
        return (other * conjugate) / (self * conjugate)

    @classmethod
    def sqrt(cls, i):
        """
        Returns QuadraticNumber which is the positive square root of the input, if this exists.
        The input i must be an int, a RationalNumber, or a QuadraticNumber of the form
        a + b sqrt(c) for rational numbers a, b, c (we do not check whether the input has a
        square root if its form is more complicated than this, even though one might exist.).
        """
        if type(i) == QuadraticNumber and i.is_rational():
            i = i.get_rational_part()

        if type(i) == RationalNumber:
            denom = i.denominator
            i = i.numerator * i.denominator
        else:
            denom = 1

        if i == 0 and type(i) in [int, QuadraticNumber]:
            return QuadraticNumber()

        if type(i) == int:
            pf = PrimeFactorization(i)
            square_free = pf.get_square_free_part()
            num = pf.get_truncated_square_root()
            ans = QuadraticNumber()
            ans.coefficients[square_free] = RationalNumber(num, denom)
            return ans

        # try to compute sqrt when i = a + b sqrt(c) for rational numbers a, b, c:
        # in this case we can have sqrt(i) = sqrt(a) * ( sqrt(x) +/- sqrt(1-x) ) for some x
        if type(i) == QuadraticNumber and len(i) == 2 and PrimeFactorization(1) in i:
            u = PrimeFactorization(1)
            v = next(iter(i.coefficients.keys() - {u}))
            a, b = i[u], i[v]

            r = b / a
            q = v.n * r**2
            sign = (-1)**(r < 0)

            c = (1 + QuadraticNumber.sqrt(1 - q)) / 2
            if c.is_rational():
                rescale = QuadraticNumber.sqrt(a)
                return rescale * (QuadraticNumber.sqrt(c) + sign * QuadraticNumber.sqrt(1 - c))

        raise Exception('Cannot compute square root of `%s`' % i)

    def is_rational(self):
        return all(pf.n == 1 for pf in self.coefficients)

    def is_real(self):
        return all(0 < pf.n for pf in self.coefficients)

    def decompose(self):
        """Returns pair (x, y) of QuadraticNumbers with positive coeff such that self = x - y."""
        positive_part, negative_part = QuadraticNumber(), QuadraticNumber()
        positive_part.coefficients = {i: v for i, v in self.coefficients.items() if 0 < v}
        negative_part.coefficients = {i: -v for i, v in self.coefficients.items() if v < 0}
        return positive_part, negative_part

    def get_rational_part(self):
        return self[PrimeFactorization(1)]

    def conjugate(self):
        """
        When self has form (a_0 + a_1 sqrt(n_1) + a_2 sqrt(n_2) + ...), this method
        returns product of all sums (a_0 +/- a_1 sqrt(n_1) +/- a_2 sqrt(n_2) +/- ...).
        """
        a = self.get_rational_part()
        b = self - a
        n = len(b)
        ans = QuadraticNumber(1)
        for i in range(1, 2**n):
            coefficients = b.coefficients.copy()
            signs = [(-1)**((i >> k) % 2) for k in range(n)]
            for j in sorted(coefficients):
                coefficients[j] *= signs[0]
                signs = signs[1:]
            term = QuadraticNumber()
            term.coefficients = coefficients
            ans *= (a + term)
        return ans

    @classmethod
    def get_index_repr(cls, index):
        return (index.n != 1) * ('sqrt(' + str(index.n) + ')')


class Monomial:

    """
    Class for objects representing monomials, that is, products of commuting indeterminates (and
    their inverses). Formally, a monomial is just a dict with integer values, but unlike a dict,
    objects of this class are hashable and orderable, can be multiplied with each other, and so on.
    """

    def __init__(self, exponents=None):
        """
        Input must be int or nonempty string (taken to be name of a single variable with
        exponent 1) or a dict with int values and keys which are all strings or all ints.
        There are some minor conditions on which strings can serve as variables names.
        """
        if exponents is None:
            self.exponents = {}
        else:
            try:
                assert exponents != '' and type(exponents) in [str, int, dict]

                if type(exponents) != dict:
                    self.exponents = {exponents: 1}
                else:
                    self.exponents = {i: e for i, e in exponents.items() if e != 0}

                assert not any(type(e) != int for e in self.exponents.values())

                all_ints = all(type(i) == int for i in self.exponents.keys())
                all_strs = all(type(i) == str for i in self.exponents.keys())
                assert all_ints or all_strs
                if all_strs:
                    # don't allow indices like 'x[...]' since leads to ambiguous str conversion
                    assert not any(s.startswith('x[') and s.endswith(']') for s in self.exponents)
            except:
                raise InvalidInputException(self, exponents)

    def __eq__(self, other):
        if type(other) != Monomial:
            return other == 1 and len(self.exponents) == 0
        else:
            return self.exponents == other.exponents

    def __lt__(self, other):
        return tuple(sorted(self.exponents.items())) < tuple(sorted(other.exponents.items()))

    def __hash__(self):
        return hash(tuple(sorted(self.exponents.items())))

    def __mul__(self, other):
        if type(other) == Monomial:
            keys = set(self.exponents.keys()) | other.exponents.keys()
            exponents = {i: self[i] + other[i] for i in keys if self[i] + other[i] != 0}
            return Monomial(exponents)
        else:
            raise Exception('Cannot multiply Monomial with `%s`' % type(other))

    def __pow__(self, other):
        if type(other) != int:
            raise Exception('Cannot exponentiate Monomial by `%s`' % type(other))
        else:
            return Monomial({i: self[i] * other for i in self.exponents})

    def __getitem__(self, i):
        return self.exponents.get(i, 0)

    def __repr__(self):
        if not self.exponents:
            return '1'
        terms = []
        for i, e in sorted(self.exponents.items()):
            base = self.index_to_string(i)
            if e != 1:
                base = base + '**' + str(e)
            terms += [base]
        return ' * '.join(terms)

    def degree(self):
        deg = 0
        for e in self.exponents.values():
            if e < 0:
                return None
            deg += e
        return deg

    @classmethod
    def index_to_string(cls, i):
        s = str(i)
        try:
            int(i)
        except:
            return s
        else:
            return 'x[%s]' % repr(i)


class Polynomial(VectorMixin, OperatorMixin, NumberMixin):

    """
    Class for objects representing polynomials in commuting variables. Polynomials are
    implemented here as vectors spanned by Monomials with QuadraticNumber coefficients.

    The most important quirk in our implementation is how the comparators <, <=, >, and >= are
    defined. Specifically:

        f < g returns True iff (f - g) has all positive real coefficients and nonzero constant term.
        f <= g returns True iff (f - g) is zero or has all positive real coefficients.

    The reverse operators > and >= are defined symmetrically. These definitions give us a simple
    way to check if f has positive or nonnegative values for all nonnegative inputs (this happens
    if 0 < f or 0 <= f). Note, however, that (f <= g) is not the same as (f < g or f == g).
    """

    @property
    def leading_coefficient(self):
        if self.coefficients:
            monomial = max(self.coefficients.keys())
            return self.coefficients[monomial]
        return 0

    @property
    def coefficients(self):
        return self._coefficients

    @coefficients.setter
    def coefficients(self, value):
        self._coefficients = value

    @classmethod
    def one(cls):
        return Polynomial(1)

    def __init__(self, i=None):
        """
        Input must be int, RationalNumber, QuadraticNumber, dict, or Monomial.
        In the first three cases, returns constant polynomial with given constant term.
        In last two cases, returns polynomial representing the given monomial.
        """
        if i is None or (type(i) in [int, RationalNumber, QuadraticNumber] and i == 0):
            self._coefficients = {}
        elif type(i) == str or type(i) == dict:
            try:
                monomial = Monomial(i)
            except:
                raise InvalidInputException(self, type(i))
            else:
                self._coefficients = {monomial: 1}
        elif type(i) == Monomial:
            self._coefficients = {i: 1}
        elif type(i) in [int, RationalNumber]:
            self._coefficients = {Monomial(): QuadraticNumber(i)}
        elif type(i) == QuadraticNumber:
            self._coefficients = {Monomial(): i}
        else:
            raise InvalidInputException(self, type(i))

    def is_comparable(self, other):
        return type(other) in [int, RationalNumber, QuadraticNumber, Polynomial]

    def lt__int(self, other):
        return self.lt__helper(other)

    def lt__rational_number(self, other):
        return self.lt__helper(other)

    def lt__quadratic_number(self, other):
        return self.lt__helper(other)

    def lt__polynomial(self, other):
        return self.lt__helper(other)

    def lt__helper(self, other):
        diff = other - self
        return 0 < diff.get_constant_part() and all(0 < c for c in diff.coefficients.values())

    def __le__(self, other):
        diff = other - self
        return all(0 < c for c in diff.coefficients.values())

    def __hash__(self):
        if self.is_constant():
            return hash(self.get_constant_part())
        else:
            def key(dict_items):
                index, coeff = dict_items
                return (index, hash(coeff))

            return hash(tuple(sorted(self.coefficients.items(), key=key)))

    def add__int(self, other):
        return self.add__polynomial(Polynomial(other))

    def add__rational_number(self, other):
        return self.add__polynomial(Polynomial(other))

    def add__quadratic_number(self, other):
        return self.add__polynomial(Polynomial(other))

    def add__polynomial(self, other):
        new = Polynomial()
        keys = set(self.coefficients.keys()) | other.coefficients.keys()
        new.coefficients = {i: self[i] + other[i] for i in keys if self[i] + other[i] != 0}
        return new

    def mul__int(self, other):
        return self.mul__polynomial(Polynomial(other))

    def mul__rational_number(self, other):
        return self.mul__polynomial(Polynomial(other))

    def mul__quadratic_number(self, other):
        return self.mul__polynomial(Polynomial(other))

    def mul__polynomial(self, other):
        new = Polynomial()
        for monomial_1, coeff_1 in other:
            for monomial_2, coeff_2 in self:
                new[monomial_1 * monomial_2] += coeff_1 * coeff_2
        return new

    def truediv__int(self, other):
        return self.truediv__helper(QuadraticNumber(other))

    def truediv__rational_number(self, other):
        return self.truediv__helper(QuadraticNumber(other))

    def truediv__quadratic_number(self, other):
        return self.truediv__helper(other)

    def truediv__polynomial(self, other):
        if other.is_constant():
            return self / other[1]
        else:
            raise OperatorException(self, other, '__truediv__')

    def truediv__helper(self, other):
        new = Polynomial()
        for i, v in self:
            new[i] = v / other
        return new

    def __rtruediv__(self, other):
        if self.is_constant() and type(other) in [int, RationalNumber, QuadraticNumber]:
            return Polynomial(other / self.get_constant_part())
        else:
            raise OperatorException(self, other, '__rtruediv__')

    def __getitem__(self, i):
        if i == 1:
            i = Monomial()
        elif type(i) in [dict, str]:
            i = Monomial(i)
        return super(Polynomial, self).__getitem__(i)

    def __repr__(self):
        if len(self) == 0:
            return '0'
        sep = ' * '
        s = ''
        for monomial, coeff in sorted(self.coefficients.items()):
            if type(coeff) == QuadraticNumber and coeff.is_rational():
                coeff = coeff.get_rational_part()
            coeff_is_rational = type(coeff) in [RationalNumber, int]
            if coeff == 1:
                s += ' + ' + str(monomial)
            elif coeff == -1:
                s += ' - ' + str(monomial)
            elif monomial == 1:
                if coeff_is_rational and coeff < 0:
                    s += ' - ' + str(-coeff)
                else:
                    s += ' + ' + str(coeff)
            elif coeff_is_rational and coeff < 0:
                s += ' - ' + str(-coeff) + sep + str(monomial)
            elif coeff_is_rational:
                s += ' + ' + str(coeff) + sep + str(monomial)
            else:
                s += ' + ' + str(coeff) + sep + str(monomial)
        if s.startswith(' - '):
            s = '-' + s[3:]
        elif s.startswith(' + '):
            s = s[3:]
        return s

    def get_variables(self):
        """Return set of integers indexing the indeterminates that appear in this polynomial."""
        return set(i for monomial in self.coefficients for i in monomial.exponents)

    def _get_monic_quadratic_coefficients(self):
        """
        Assuming self is univariate and non-constant, returns real QuadraticNumbers a, b, c and
        (monomial) Polynomial x such that f = ax^2 + bx + c is monic and self is a constant
        multiple of f. If such numbers do not exist, raises CannotFactorException.
        """
        def to_quadratic_number(i):
            if type(i) != QuadraticNumber:
                i = QuadraticNumber(i)
            return i

        x = next(iter(self.get_variables()))
        a = to_quadratic_number(self[{x: 2}])
        b = to_quadratic_number(self[{x: 1}])
        c = to_quadratic_number(self[{}])
        x = Polynomial({x: 1})

        # check that self is linear or quadratic
        if self != a * x**2 + b * x + c:
            raise CannotFactorException(self)

        # normalize coefficients so that a * x**2 + b * x + c is monic
        if a != 0:
            b /= a
            c /= a
            a = QuadraticNumber(1)
        elif b != 0:
            c /= b
            b = QuadraticNumber(1)

        if not (a.is_real() and b.is_real() and c.is_real()):
            raise CannotFactorException(self)

        return a, b, c, x

    def get_real_quadratic_factors(self):
        """
        If self is a nonzero constant multiple of a polynomial in one variable with degree at most
        two and real coefficients, then this returns set of self's real-rooted monic linear factors.
        Returns empty set if self has no real roots, and raises Exception if cannot uniquely factor.
        """
        # check that self is not zero
        if self == 0:
            raise CannotFactorException(self)

        # case: polynomial is nonzero and constant, so has no real roots
        if self.is_constant():
            return set()

        # case: multiple variables
        if len(self.get_variables()) > 1:
            # is polynomial divisible by some variable x?
            factorable_variables = set(next(iter(self.coefficients)).exponents.keys())
            for monomial in self.coefficients:
                factorable_variables &= monomial.exponents.keys()
            if factorable_variables:
                v = next(iter(factorable_variables))
                return {self * Polynomial(Monomial({v: -1})), Polynomial({v: 1})}
            raise CannotFactorException(self)

        # if we get here, self is univariate and not constant
        a, b, c, x = self._get_monic_quadratic_coefficients()

        # case: discriminant is negative, so polynomial has no real roots
        if (b**2 - 4 * a * c) < 0:
            return set()

        # case: polynomial is linear, so has one root
        if a == 0:
            return {x + c}

        # finally, try to apply quadratic formula (note that a == 1)
        try:
            r = (-b + QuadraticNumber.sqrt(b**2 - 4 * c)) / 2
            s = -b - r
            return {x - r, x - s}
        except Exception as e:
            # raise broader exception as failure indicates that roots are not QuadraticNumbers
            raise Exception(
                ('Cannot compute roots of factorable polynomial: %s\n' % str(self)) +
                ('Error message: %s' % str(e)))

    def degree(self):
        """Returns degree of polynomialor None if self is a Laurent polynomial."""
        monomial_degrees = {0}
        for monomial in self.coefficients:
            monomial_degrees.add(monomial.degree())
        if None in monomial_degrees:
            return None
        else:
            return max(monomial_degrees)

    def is_constant(self):
        return list(self.coefficients.keys()) in [[], [Monomial()]]

    def get_constant_part(self):
        return self[Monomial()]

    def is_rational(self):
        if self.is_constant():
            v = self.get_constant_part()
            return type(v) == int or v.is_rational()
        return False

    def set_variable(self, variable, value):
        """
        Returns polynomial given by replacing given variable in self with provided value.
        The input variable must be a linear Monomial, or a valid input to Monomial.__init__
        that constructs a linear Monomial.
        """
        try:
            assert type(value) in [int, RationalNumber, QuadraticNumber, Polynomial]

            if type(variable) != Monomial:
                variable = Monomial(variable)

            assert list(variable.exponents.values()) == [1]

            variable_index = next(iter(variable.exponents.keys()))
        except:
            raise InvalidInputException(self, (variable, value), 'set_variable')

        new = Polynomial()
        for monomial, coeff in self:
            e = monomial[variable_index]
            if e != 0 and value != 0:
                exponents = {i: e for i, e in monomial.exponents.items() if i != variable_index}
                new += Polynomial(exponents) * coeff * (value**e)
            elif e < 0 and value == 0:
                raise Exception('Division by zero when setting variable in `%s`' % str(self))
            elif e == 0:
                new[monomial] += coeff
        return new


class Matrix:

    """
    Custom class implementing simple matrix operations -- essentially just determinants, inverses,
    multplication of matrices by vectors, and Vandermonde matrices.
    """

    def __init__(self, matrix, variable=None):
        """
        Constructs matrix from list of lists of Polynomials in the single provided variable.
        The default value of None for `variable` indicates that the entries of the input matrix
        are constant. It is assumed that all entries of matrix are constant or linear polynomials.
        """
        if len({len(row) for row in matrix}) > 1:
            raise InvalidInputException(self, type(matrix))
        self.rows = matrix
        self.variable = variable
        self.num_rows = len(matrix)
        if matrix:
            self.num_cols = len(matrix[0])
        else:
            self.num_cols = 0

    def __repr__(self):
        max_len = 0
        for row in self.rows:
            for e in row:
                max_len = max([max_len, len(str(e))])

        def pad(x):
            ans = str(x)
            return ans + (max_len - len(ans)) * ' '

        return '\n'.join(['[' + '  '.join([pad(e) for e in row]) + ']' for row in self.rows])

    def copy(self):
        matrix = [[e for e in row] for row in self.rows]
        return Matrix(matrix, self.variable)

    def __eq__(self, other):
        if type(other) != Matrix or self.num_rows != other.num_rows:
            return False
        for i in range(len(self.rows)):
            if self.rows[i] != other.rows[i]:
                return False
        return True

    def __getitem__(self, i):
        return self.rows[i]

    def __mul__(self, other):
        """Other should be a matrix of an iterable representing a vector."""
        if type(other) == Matrix:
            variables = {self.variable, other.variable} - {None}
            if self.num_cols != other.num_rows or len(variables) == 2:
                raise OperatorException(self, other, '__mul__')
            if variables:
                var = variables.pop()
            else:
                var = None
            return Matrix([
                [
                    sum(self[i][j] * other[j][k] for j in range(self.num_cols))
                    for k in range(other.num_cols)
                ]
                for i in range(self.num_rows)
            ], variable=var)
        elif len(other) != self.num_cols:
            raise OperatorException(self, other, '__mul__')
        return [sum(row[i] * other[i] for i in range(len(other))) for row in self.rows]

    @classmethod
    def vandermonde_matrix(cls, n):
        return Matrix([[i**j for j in range(n)] for i in range(n)])

    @classmethod
    def vandermonde_inverse(cls, n):
        return cls.vandermonde_matrix(n).inverse()

    def inverse(self):
        n = self.num_rows
        if n != self.num_cols:
            return None

        # augment matrix by identity matrix on the right
        matrix = self.copy()
        for i in range(n):
            matrix.rows[i] += i * [0] + [1] + (n - 1 - i) * [0]
        matrix.num_cols = 2 * n

        # row reduce
        for column in range(n):
            i = matrix._find_invertible_entry_row(column)
            if i is None:
                return None
            matrix._cancel_column(i, column)
            matrix._swap_rows(i, column)

        # cut off first n columns to get inverse
        for i in range(n):
            matrix.rows[i] = matrix.rows[i][n:]
        matrix.num_cols = n
        return matrix

    def determinant(self):
        matrix = self.copy()
        det = QuadraticNumber(1)
        for column in range(matrix.num_cols):
            i = matrix._find_invertible_entry_row(column)
            j = matrix._find_nonzero_entry_row(column)
            if i is None and j is None:
                return 0
            elif i is None and j is not None:
                det *= matrix[j][column]
                det *= matrix._swap_rows(j, column)
            else:
                det *= matrix._cancel_column(i, column)
                det *= matrix._swap_rows(i, column)
        return det

    def _swap_rows(self, i, j):
        """Swaps rows i and j of given matrix, returns -1 if i != j and otherwise 1."""
        if i == j:
            return 1
        else:
            row_i = self.rows[i]
            self.rows[i] = self.rows[j]
            self.rows[j] = row_i
            return -1

    def _normalize_row(self, i, column):
        """Divides all entries in row i by scalar in given column of that row. Returns scalar."""
        scalar = self[i][column]
        if type(scalar) == int:
            scalar = RationalNumber(scalar)
        return self._scale_row(i, scalar)

    def _scale_row(self, i, scalar):
        """Divides all entries in row i by given scalar, and returns scalar."""
        for j in range(self.num_cols):
            self.rows[i][j] /= scalar
        return scalar

    def _add_row(self, src_row, dest_row, scalar):
        """Add scalar multiple of src_row of matrix to dest_row."""
        assert src_row != dest_row
        for k in range(self.num_cols):
            self.rows[dest_row][k] += scalar * self.rows[src_row][k]

    def _cancel_column(self, i, column):
        """
        Rescales row i so that entry in given column is 1, then adds copies of rescaled row to
        cancel all other entries in given column. Returns original value of self[i][column].
        """
        ans = self._normalize_row(i, column)
        for j in range(self.num_rows):
            scalar = self[j][column]
            if j == i or scalar == 0:
                continue
            self._add_row(i, j, -scalar)
        return ans

    def _find_nonzero_entry_row(self, column):
        """Return first row index i >= column such that self[i][column] is nonzero, or None."""
        for i in range(column, self.num_rows):
            scalar = self[i][column]
            if scalar != 0:
                return i
        return None

    def _find_invertible_entry_row(self, column):
        """
        Tries to find index of row with invertible entry in given column. Returns index if found.
        If all entries in column are (linear) polynomials, adds multiples of one row to the
        others in order to cancel entries in all but this row, then returns this row's index.
        If this fails, returns None.
        """
        for i in range(column, self.num_rows):
            scalar = self[i][column]
            if scalar != 0 and (type(scalar) != Polynomial or scalar.is_constant()):
                return i

        if self.variable is None:
            return None

        for i in range(column, self.num_rows):
            for j in range(i + 1, self.num_rows):
                f = self[i][column]
                g = self[j][column]
                if f == 0 or g == 0:
                    continue

                assert f.degree() <= 1 and g.degree() <= 1

                a = f[{self.variable: 1}]
                b = g[{self.variable: 1}]
                if type(b) == int:
                    b = RationalNumber(b)

                self._add_row(i, j, -b / a)
                if self[j][column] != 0:
                    return j
