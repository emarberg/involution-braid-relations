import numpy as np

from project.utils import (
    reverse_tuple,
    IndeterminatePowerException,
    InvalidInputException,
    ZeroDivisionException,
    NumberMixin,
    VectorMixin
)


class OperatorMixin:
    """
    Mixin class which helps implement polymorphic binary operations ==, <, +, *, and /
    between int, RationalNumber, QuadraticNumber, and Polynomial objects.
    """

    def __eq__(self, other):
        if type(other) == int:
            return self.eq__int(other)
        elif type(other) == RationalNumber:
            return self.eq__rational_number(other)
        elif type(other) == QuadraticNumber:
            return self.eq__quadratic_number(other)
        elif type(other) == Polynomial:
            return self.eq__polynomial(other)
        else:
            raise self.OperatorException(other)

    def eq__int(self, other):
        """Evaluates self == other under assumption that type(other) is int."""
        raise NotImplementedError

    def eq__rational_number(self, other):
        """Evaluates self == other under assumption that type(other) is RationalNumber."""
        raise NotImplementedError

    def eq__quadratic_number(self, other):
        """Evaluates self == other under assumption that type(other) is QuadraticNumber."""
        raise NotImplementedError

    def eq__polynomial(self, other):
        """Evaluates self == other under assumption that type(other) is Polynomial."""
        raise NotImplementedError

    def __lt__(self, other):
        if type(other) == int:
            return self.lt__int(other)
        elif type(other) == RationalNumber:
            return self.lt__rational_number(other)
        elif type(other) == QuadraticNumber:
            return self.lt__quadratic_number(other)
        elif type(other) == Polynomial:
            return self.lt__polynomial(other)
        else:
            raise self.OperatorException(other, '__lt__')

    def lt__int(self, other):
        """Evaluates self < other under assumption that type(other) is int."""
        raise NotImplementedError

    def lt__rational_number(self, other):
        """Evaluates self <other under assumption that type(other) is RationalNumber."""
        raise NotImplementedError

    def lt__quadratic_number(self, other):
        """Evaluates self < other under assumption that type(other) is QuadraticNumber."""
        raise NotImplementedError

    def lt__polynomial(self, other):
        """Evaluates self < other under assumption that type(other) is Polynomial."""
        raise NotImplementedError

    def __add__(self, other):
        if type(other) == int:
            return self.add__int(other)
        elif type(other) == RationalNumber:
            return self.add__rational_number(other)
        elif type(other) == QuadraticNumber:
            return self.add__quadratic_number(other)
        elif type(other) == Polynomial:
            return self.add__polynomial(other)
        else:
            raise self.OperatorException(other, '__add__')

    def add__int(self, other):
        """Evaluates self + other under assumption that type(other) is int."""
        raise NotImplementedError

    def add__rational_number(self, other):
        """Evaluates self + other under assumption that type(other) is RationalNumber."""
        raise NotImplementedError

    def add__quadratic_number(self, other):
        """Evaluates self + other under assumption that type(other) is QuadraticNumber."""
        raise NotImplementedError

    def add__polynomial(self, other):
        """Evaluates self + other under assumption that type(other) is Polynomial."""
        raise NotImplementedError

    def __mul__(self, other):
        if type(other) == int:
            return self.mul__int(other)
        elif type(other) == RationalNumber:
            return self.mul__rational_number(other)
        elif type(other) == QuadraticNumber:
            return self.mul__quadratic_number(other)
        elif type(other) == Polynomial:
            return self.mul__polynomial(other)
        elif type(other) == Root:
            return self.mul__root(other)
        else:
            raise self.OperatorException(other, '__mul__')

    def mul__int(self, other):
        """Evaluates self * other under assumption that type(other) is int."""
        raise NotImplementedError

    def mul__rational_number(self, other):
        """Evaluates self * other under assumption that type(other) is RationalNumber."""
        raise NotImplementedError

    def mul__quadratic_number(self, other):
        """Evaluates self * other under assumption that type(other) is QuadraticNumber."""
        raise NotImplementedError

    def mul__polynomial(self, other):
        """Evaluates self * other under assumption that type(other) is Polynomial."""
        raise NotImplementedError

    def mul__root(self, other):
        """Evaluates self * other under assumption that type(other) is Root."""
        return other * self

    def __truediv__(self, other):
        if other == 0:
            raise ZeroDivisionException(self)
        elif type(other) == int:
            return self.truediv__int(other)
        elif type(other) == RationalNumber:
            return self.truediv__rational_number(other)
        elif type(other) == QuadraticNumber:
            return self.truediv__quadratic_number(other)
        elif type(other) == Polynomial:
            return self.truediv__polynomial(other)
        else:
            raise self.OperatorException(other, '__truediv__')

    def truediv__int(self, other):
        """Evaluates self * other under assumption that type(other) is int."""
        raise NotImplementedError

    def truediv__rational_number(self, other):
        """Evaluates self * other under assumption that type(other) is RationalNumber."""
        raise NotImplementedError

    def truediv__quadratic_number(self, other):
        """Evaluates self * other under assumption that type(other) is QuadraticNumber."""
        raise NotImplementedError

    def truediv__polynomial(self, other):
        """Evaluates self * other under assumption that type(other) is Polynomial."""
        raise NotImplementedError


class RationalNumber(OperatorMixin, NumberMixin):

    class OperatorException(Exception):
        def __init__(self, other, operator='__eq__'):
            method = 'RationalNumber.' + operator
            super(RationalNumber.OperatorException, self).__init__(
                'Cannot evaluate %s with input of type `%s`' % (method, type(other)))

    def __init__(self, p=0, q=1):
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
            raise RationalNumber.OperatorException(other, '__truediv__')

    def __rtruediv__(self, other):
        return other * RationalNumber(1, self)

    def __pow__(self, exponent):
        if type(exponent) != int:
            raise RationalNumber.OperatorException(exponent, '__pow__')
        elif exponent == 0 and self != 0:
            return RationalNumber(1)
        elif exponent == 0 and self == 0:
            raise IndeterminatePowerException

        x = super(RationalNumber, self).__pow__(abs(exponent))
        if exponent < 0:
            return RationalNumber(1, x)
        else:
            return x

    def __repr__(self):
        if self.numerator == 0:
            return '0'
        elif self.denominator != 1:
            return str(self.numerator) + '/' + str(self.denominator)
        else:
            return str(self.numerator)


class PrimeFactorization:
    _factorization_cache = {}

    def __init__(self, i=1):
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
            ans.factorization = {p: self[p] + other[p] for p in factors}
            if ans[-1] != 0 and ans[-1] % 2 == 0:
                del ans.factorization[-1]
            return ans
        else:
            raise Exception('Cannot multiply PrimeFactorization with %s' % type(other))

    def __repr__(self):
        return repr(self.factorization)

    def get_square_free_part(self):
        ans = PrimeFactorization()
        ans.factorization = {p: e % 2 for p, e in self.factorization.items() if e % 2 != 0}
        ans.n = 1
        for p in ans.factorization:
            ans.n *= p
        return ans

    def get_truncated_square_root(self):
        r = 1
        for p, e in self.factorization.items():
            r *= p**(e // 2)
        return r

    @classmethod
    def get_divisor_exponent(cls, n, p):
        """Return maximum positive integer e such that p^e divides n."""
        if n % p != 0:
            return 0
        if abs(n) == abs(p):
            return 1
        else:
            return 1 + cls.get_divisor_exponent(n // p, p)

    @classmethod
    def get_prime_factorization(cls, i):
        if i in cls._factorization_cache:
            return cls._factorization_cache[i].copy()
        else:
            i_input = i
            factorization = {}
            N = range(2, abs(i) + 1)
            while N:
                p = N[0]
                e = cls.get_divisor_exponent(i, p)
                if e != 0:
                    factorization[p] = e
                    i = i // p**e
                N = [a for a in N[1:] if a <= abs(i) and a % p != 0]
            if i_input < 0:
                factorization[-1] = 1
            cls._factorization_cache[i_input] = factorization
            return factorization.copy()


class QuadraticNumber(VectorMixin, OperatorMixin, NumberMixin):

    class ImaginaryComparisonException(Exception):
        def __init__(self, diff):
            super(QuadraticNumber.ImaginaryComparisonException, self).__init__(
                'Cannot compare quadratic numbers with non-real difference %s' % diff)

    class IndeterminateComparisonException(Exception):
        def __init__(self, a, b):
            super(QuadraticNumber.IndeterminateComparisonException, self).__init__(
                'Cannot determine inequality %s < %s' % (a, b))

    def __init__(self, i=0):
        if type(i) == int:
            i = RationalNumber(i)
        if type(i) == RationalNumber:
            self.coefficients = {}
            if i != 0:
                self.coefficients[PrimeFactorization(1)] = i
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
        if self.is_rational():
            return hash(self.get_rational_part())
        else:
            return hash(tuple(sorted(self.coefficients.items())))

    def lt__int(self, other):
        return self.lt__quadratic_number(QuadraticNumber(other))

    def lt__rational_number(self, other):
        return self.lt__quadratic_number(QuadraticNumber(other))

    def lt__quadratic_number(self, other):
        diff = other - self
        if diff.is_rational():
            return 0 < diff[PrimeFactorization(1)]
        elif not diff.is_real():
            raise QuadraticNumber.ImaginaryComparisonException(diff)
        elif all(0 < c for c in diff.coefficients.values()):
            return True
        elif all(c < 0 for c in diff.coefficients.values()):
            return False

        positive_part, negative_part = diff.decompose()
        if max(len(positive_part), len(negative_part)) > 2:
            raise QuadraticNumber.IndeterminateComparisonException(negative_part, positive_part)
        else:
            return (negative_part**2).lt__quadratic_number(positive_part**2)

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

    def __pow__(self, exponent):
        if type(exponent) == int:
            if exponent == 0 and self != 0:
                return QuadraticNumber(1)
            elif exponent == 0 and self == 0:
                raise IndeterminatePowerException
            elif exponent < 0:
                return 1 / super(QuadraticNumber, self).__pow__(-exponent)
        return super(QuadraticNumber, self).__pow__(exponent)

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
        raise QuadraticNumber.OperatorException(self, other, '__truediv__')

    def __rtruediv__(self, other):
        conjugate = self.conjugate()
        return (other * conjugate) / (self * conjugate)

    @classmethod
    def sqrt(cls, i):
        """TODO: clean up this method."""
        if type(i) == QuadraticNumber and i.is_rational():
            i = i.get_rational_part()
        if type(i) == RationalNumber:
            denom = i.denominator
            i = i.numerator * i.denominator
        else:
            denom = 1

        if i == 0 and type(i) in [int, QuadraticNumber]:
            return QuadraticNumber()
        elif type(i) == int:
            pf = PrimeFactorization(i)
            square_free = pf.get_square_free_part()
            ans = QuadraticNumber()
            ans.coefficients[square_free] = RationalNumber(pf.get_truncated_square_root(), denom)
            return ans
        elif type(i) == QuadraticNumber:
            q = i / (3 + cls.sqrt(5))
            if q.is_rational():
                return cls.sqrt(q) * (cls.sqrt(2) + cls.sqrt(10)) / 2
            q = i / (7 + 3 * cls.sqrt(5))
            if q.is_rational():
                return cls.sqrt(q) * (3 * cls.sqrt(2) + cls.sqrt(10)) / 2
        raise Exception('Cannot compute square root of `%s`' % i)

    def is_rational(self):
        return all(pf.n == 1 for pf in self.coefficients)

    def is_real(self):
        return all(0 < pf.n for pf in self.coefficients)

    def decompose(self):
        positive_part, negative_part = QuadraticNumber(), QuadraticNumber()
        positive_part.coefficients = {i: v for i, v in self.coefficients.items() if 0 < v}
        negative_part.coefficients = {i: -v for i, v in self.coefficients.items() if v < 0}
        return positive_part, negative_part

    def get_rational_part(self):
        return self[PrimeFactorization(1)]

    def conjugate(self):
        """
        If self has form (a_0 + a_1 sqrt(n_1) + a_2 sqrt(n_2) + ...) then
        returns product of all conjugates (a_0 +/- a_1 sqrt(n_1) +/- a_2 sqrt(n_2) + ...).
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
    def __init__(self, exponents=None):
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
                base = base + '^' + str(e)
            terms += [base]
        return ' '.join(terms)

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
    def __init__(self, i=None):
        if i is None or (type(i) in [int, RationalNumber, QuadraticNumber] and i == 0):
            self.coefficients = {}
        elif type(i) == str or type(i) == dict:
            try:
                monomial = Monomial(i)
            except:
                raise InvalidInputException(self, type(i))
            else:
                self.coefficients = {monomial: 1}
        elif type(i) == Monomial:
            self.coefficients = {i: 1}
        elif type(i) in [int, RationalNumber]:
            self.coefficients = {Monomial(): QuadraticNumber(i)}
        elif type(i) == QuadraticNumber:
            self.coefficients = {Monomial(): i}
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

    def __pow__(self, other):
        if type(other) == int and other == 0:
            return Polynomial(1)
        else:
            return super(Polynomial, self).__pow__(other)

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
            raise Polynomial.OperatorException(self, other, '__truediv__')

    def truediv__helper(self, other):
        new = Polynomial()
        for i, v in self:
            new[i] = v / other
        return new

    def __rtruediv__(self, other):
        if self.is_constant() and type(other) in [int, RationalNumber, QuadraticNumber]:
            return Polynomial(other / self.get_constant_part())
        else:
            raise Polynomial.OperatorException(self, other, '__rtruediv__')

    def __getitem__(self, i):
        if i == 1:
            i = Monomial()
        elif type(i) in [dict, str]:
            i = Monomial(i)

        return super(Polynomial, self).__getitem__(i)

    def __repr__(self):
        if len(self) == 0:
            return '0'
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
                s += ' - ' + str(-coeff) + str(monomial)
            elif coeff_is_rational:
                s += ' + ' + str(coeff) + str(monomial)
            else:
                s += ' + ' + str(coeff) + '' + str(monomial)
        if s.startswith(' - '):
            s = '-' + s[3:]
        elif s.startswith(' + '):
            s = s[3:]
        return s

    def get_variables(self):
        """Return set of integers indexing the indeterminates that appear in this polynomial."""
        return set(i for monomial in self.coefficients for i in monomial.exponents)

    def get_quadratic_discriminant(self):
        x = next(iter(self.get_variables()))
        a = self[{x: 2}]
        b = self[{x: 1}]
        c = self[{}]
        return b**2 - 4 * a * c

    def is_factorable(self):
        if len(self.get_variables()) != 1:
            return False

        x = next(iter(self.get_variables()))
        a = self[{x: 2}]
        b = self[{x: 1}]
        c = self[{}]

        return self == a * Polynomial({x: 2}) + b * Polynomial({x: 1}) + c * Polynomial({})

    def get_factors(self):
        """TODO: improve this method."""
        if len(self.get_variables()) != 1:
            raise Exception('Cannot factor `%s`' % str(self))

        def to_quadratic_number(i):
            if type(i) != QuadraticNumber:
                i = QuadraticNumber(i)
            return i

        x = next(iter(self.get_variables()))
        a = to_quadratic_number(self[{x: 2}])
        b = to_quadratic_number(self[{x: 1}])
        c = to_quadratic_number(self[{}])
        x = Polynomial({x: 1})

        if self != a * x**2 + b * x + c:
            raise Exception('Cannot factor `%s`' % str(self))

        if a != 0:
            b /= a
            c /= a

        if a == 0 and b != 0:
            return {x + c/b}
        elif c == 0:
            return {x, x + b}
        else:
            try:
                r = (-b + QuadraticNumber.sqrt(b**2 - 4 * c)) / 2
                s = (-b - QuadraticNumber.sqrt(b**2 - 4 * c)) / 2
                return {x - r, x - s}
            except:
                raise Exception('Cannot factor `%s`' % str(self))

    def degree(self):
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
            if type(v) == int or v.is_rational():
                return True
        return False

    def set_variable(self, variable, value):
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


class CoxeterGraph:

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
        Input `edges` should be list of tuples (i, j, m) where i, j are indices of simple
        generators and m is order of s_i s_j. Triples with m == 1 or 2 can be omitted.
        If (i, j, m) is included then the reverse triple (j, i, m) also may be omitted.

        Input `generators` should be list of integers indexing simple generators of the group.
        If not provided, this list will be formed from the set of numbers i or j in `edges`.
        If i or j does not belong to `generators` for some (i, j, m) in `edges`, an error
        will be raised.

        Input `star` should be list of pairs (i, j) such that the involution * : S -> S
        with i^* = j and j^* = i extends to an automorphism of W. If `star` is not given,
        * is defined to be the identity map.
        """
        self._setup_generators(edges, generators)
        self._setup_orders(edges)
        self._setup_star(star)

    def _setup_generators(self, edges, generators):
        # assign sorted list of simple generator indices to self.generators
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
        # construct dictionary with orders m_ij of products of simple generators
        self.orders = {}
        for i, j, m in edges:
            # value of m must be an integer with 1 <= m <= infty
            valid_order = (type(m) == int and 1 <= m) or m == np.infty
            # must have m == 1 iff i == j
            valid_order = valid_order and ((m == 1) == (i == j))
            valid_generators = i in self.generators and j in self.generators

            if not (valid_order and valid_generators) or self.orders.get((i, j), m) != m:
                raise InvalidInputException(self, '%s in `tuples`' % str((i, j, m)))
            elif i != j and m != 2:
                self.orders[(i, j)] = m
                self.orders[(j, i)] = m

    def _setup_star(self, star):
        # construct dictionary with images of the *-involution 'star'
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
        m = self.get_order(i, j)
        if m % 2 != 0:
            return (m + 1) // 2
        elif generators_are_fixed:
            return m // 2 + 1
        else:
            return m // 2

    def eval_bilinear(self, i, j):
        """Returns -cos(pi/m_ij)."""
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

    @staticmethod
    def A(n, star=None):
        edges = [(i, i + 1, 3) for i in range(1, n)]
        return CoxeterGraph(edges, star=star)

    @staticmethod
    def A2(n):
        star = [(i, n + 1 - i) for i in range(1, n + 1)]
        return CoxeterGraph.A(n, star=star)

    @staticmethod
    def B2(n=2):
        assert n == 2
        return CoxeterGraph([(0, 1, 4)], star=[(0, 1)])

    @staticmethod
    def B(n):
        """
        Dynkin diagram labeling is:

            0==1--2--...--(n-1)

        """
        assert 2 <= n
        edges = [(i, i + 1, 3) for i in range(1, n - 1)] + [(0, 1, 4)]
        return CoxeterGraph(edges)

    @staticmethod
    def D(n, star=None):
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
    def D2(n):
        assert 4 <= n
        star = [(0, 1)] + [(i, i) for i in range(2, n)]
        return CoxeterGraph.D(n, star=star)

    @staticmethod
    def E(n, star=None):
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
    def E2(n=6):
        assert n == 6
        star = [(1, 6), (2, 5), (3, 3), (4, 4)]
        return CoxeterGraph.E(n, star=star)

    @staticmethod
    def F(n=4, star=None):
        assert n == 4
        edges = [(1, 2, 3), (2, 3, 4), (3, 4, 3)]
        return CoxeterGraph(edges, star=star)

    @staticmethod
    def F2(n=4):
        assert n == 4
        star = [(1, 4), (2, 3)]
        return CoxeterGraph.F(n, star=star)

    @staticmethod
    def G(n=2, star=None):
        assert n == 2
        return CoxeterGraph([(1, 2, 6)], star=star)

    @staticmethod
    def G2(n=2):
        assert n == 2
        return CoxeterGraph.G(n, star=[(1, 2)])

    @staticmethod
    def H(n):
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
    def A_tilde(n, star=None):
        assert n >= 2
        edges = [(i, i + 1, 3) for i in range(1, n)] + [(n, n + 1, 3), (n + 1, 1, 3)]
        return CoxeterGraph(edges, star=star)

    @staticmethod
    def B_tilde(n, star=None):
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
    def C_tilde(n, star=None):
        """
        Dynkin diagram labeling is:

            0==1--2--...--(n-1)==n

        """
        assert n >= 3
        edges = [(i, i + 1, 3) for i in range(1, n - 1)] + [(0, 1, 4)] + [(n - 1, n, 4)]
        return CoxeterGraph(edges, star=star)

    @staticmethod
    def D_tilde(n, star=None):
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
    def E_tilde(n, star=None):
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
    def F_tilde(n=4):
        """
        Dynkin diagram labeling is:

            1--2==3--4--5

        """
        assert n == 4
        edges = [(1, 2, 3), (2, 3, 4), (3, 4, 3), (4, 5, 3)]
        return CoxeterGraph(edges)

    @staticmethod
    def G_tilde(n=2):
        """
        Dynkin diagram labeling is:

            1≡≡2--3

        """
        assert n == 2
        return CoxeterGraph([(1, 2, 6), (2, 3, 3)])


class Root(VectorMixin, NumberMixin):

    def __init__(self, coxeter_graph, index=None, coeff=1):
        self.graph = coxeter_graph
        if index is None or coeff == 0:
            self.coefficients = {}
        elif index in coxeter_graph.generators:
            self.coefficients = {index: coeff}
        else:
            raise InvalidInputException(self, index)

    def is_comparable(self, other):
        if type(other) == int and other == 0:
            return True
        elif type(other) == Root and other.graph == self.graph:
            return True
        else:
            return False

    def eval_bilinear(self, other):
        try:
            is_generator = (other in self.graph.generators)
        except:
            is_generator = False

        if is_generator:
            other = Root(self.graph, other)

        if type(other) == Root and other.graph == self.graph:
            ans = 0
            for i, u in self:
                for j, v in other:
                    ans += u * v * self.graph.eval_bilinear(i, j)
            return ans
        else:
            raise InvalidInputException(self, other, 'eval_bilinear')

    def reflect(self, index):
        if index in self.graph.generators:
            v = 2 * self.eval_bilinear(index)
            return self - Root(self.graph, index, v)
        else:
            raise InvalidInputException(self, index, 'reflect')

    def is_constant(self):
        return not any(type(v) == Polynomial and not v.is_constant() for i, v in self)

    def is_positive(self):
        return (not self.is_zero()) and all(0 <= v for _, v in self)

    def is_negative(self):
        return (not self.is_zero()) and all(v <= 0 for _, v in self)

    def is_valid(self):
        """
        A Root object, which is a linear combination sum_i c_i alpha_i of simple roots,
        is 'invalid' if it is zero or if for some i and j it holds that -c_i and c_j are
        both polynomials whose coefficients are all nonnegative and whose constant
        coefficients are positive.

        An 'invalid' linear combination does not specialize to a positive or negative root
        for any choice of nonnegative indeterminate values.
        """
        if self.is_zero():
            return False
        if any(v < 0 for _, v in self) and any(0 < v for _, v in self):
            return False
        return True

    def set_variable(self, variable, value):
        new = Root(self.graph)
        for i, v in self:
            if type(v) == Polynomial:
                v = v.set_variable(variable, value)
            if v != 0:
                new.coefficients[i] = v
        return new

    def __add__(self, other):
        if other == 0:
            other = Root(self.graph)
        if type(other) == Root and self.graph == other.graph:
            indices = set(self.coefficients.keys()) | set(other.coefficients.keys())
            new = Root(self.graph)
            new.coefficients = {i: self[i] + other[i] for i in indices if (self[i] + other[i]) != 0}
            return new
        else:
            raise Root.OperatorException(self, other, '__add__')

    def __mul__(self, other):
        new = Root(self.graph)
        if other != 0:
            new.coefficients = {i: other * self[i] for i in self.coefficients}
        return new

    def __hash__(self):
        return hash((self.graph,) + tuple(self[i] for i in self.graph.generators))

    def __pow__(self, exponent):
        raise NotImplementedError

    def __repr__(self):
        return super(Root, self).__repr__()[1:-1]

    @classmethod
    def get_index_repr(cls, index):
        return 'alpha_' + str(index)


class RootTransform:
    def __init__(self, coxeter_graph, sigma={}):
        if all(i in coxeter_graph.generators for i in sigma) and \
           all(type(r) == Root and r.graph == coxeter_graph for r in sigma.values()):
            self.graph = coxeter_graph
            self.sigma = sigma.copy()
            self._unconditional_descents = None
            self._strong_conditional_descents = None
            self._weak_conditional_descents = None
        else:
            raise InvalidInputException(self, (coxeter_graph, sigma))

    def __eq__(self, other):
        return isinstance(other, RootTransform) and \
            other.graph == self.graph and other.sigma == self.sigma

    def __getitem__(self, i):
        return self.sigma[i]

    def __hash__(self):
        return hash((self.graph,) + tuple(sorted(self.sigma.items())))

    def _unset_cached_properties(self):
        self._unconditional_descents = None
        self._strong_conditional_descents = None
        self._weak_conditional_descents = None

    def copy(self):
        other = RootTransform(self.graph, self.sigma)
        other._unconditional_descents = self._unconditional_descents
        other._strong_conditional_descents = self._unconditional_descents
        other._weak_conditional_descents = self._weak_conditional_descents
        return other

    def __setitem__(self, i, value):
        if i in self.graph.generators and type(value) == Root and value.graph == self.graph:
            self.sigma[i] = value
            self._unset_cached_properties()
        else:
            raise InvalidInputException(self, (i, value), '__setitem__')

    @property
    def unconditional_descents(self):
        if self._unconditional_descents is None:
            self._unconditional_descents = {
                i for i in self.sigma
                if any(f < 0 for f in self.sigma[i].coefficients.values())
            }
        return self._unconditional_descents

    @property
    def strong_conditional_descents(self):
        if self._strong_conditional_descents is None:
            self._strong_conditional_descents = {
                i for i in self.sigma
                if any(f <= 0 for f in self.sigma[i].coefficients.values())
            }
        return self._strong_conditional_descents

    @property
    def weak_conditional_descents(self):
        if self._weak_conditional_descents is None:
            self._weak_conditional_descents = {
                i for i in self.sigma
                if not any(0 < f for f in self.sigma[i].coefficients.values())
            }
        return self._weak_conditional_descents

    @classmethod
    def identity(cls, coxeter_graph):
        sigma = {i: Root(coxeter_graph, i) for i in coxeter_graph.generators}
        return cls(coxeter_graph, sigma)

    def __mul__(self, j):
        if j not in self.graph.generators:
            raise InvalidInputException(self, j, '__mul__')
        new = {}
        for i in self.sigma:
            root = Root(self.graph, i).reflect(j)
            for k, v in root:
                if k not in self.sigma:
                    raise InvalidInputException(self, j, '__mul__')
                new[i] = new.get(i, 0) + self.sigma[k] * v
        return self.__class__(self.graph, new)

    def __rmul__(self, j):
        if j not in self.graph.generators:
            raise InvalidInputException(self, j, '__rmul__')
        new = {}
        for i in self.sigma:
            new[i] = self.sigma[i].reflect(j)
        return self.__class__(self.graph, new)

    def to_coxeter_transform(self):
        return CoxeterTransform(self.graph, self.sigma)

    def is_constant(self):
        return all(r.is_constant() for r in self.sigma.values())

    def is_complete(self):
        return set(self.sigma.keys()) == set(self.graph.generators)

    def is_positive(self):
        return all(r.is_positive() for r in self.sigma.values())

    def __len__(self):
        return len(self.sigma)

    def __contains__(self, i):
        return i in self.sigma

    def __iter__(self):
        return self.sigma.__iter__()

    def values(self):
        return self.sigma.values()

    def is_identity(self):
        if sorted(self.sigma.keys()) == self.graph.generators:
            return all(self.sigma[i] == Root(self.graph, i) for i in self.graph.generators)
        return False

    def __repr__(self):
        s = '{\n'
        for i in sorted(self.sigma):
            s += '  alpha_%s -> %s\n' % (i, self.sigma[i])
        s += '}'
        return s


class CoxeterTransform(RootTransform):
    def __init__(self, coxeter_graph, sigma=None):
        if sigma:
            keys_valid = set(sigma.keys()) == set(coxeter_graph.generators)
            roots_valid = all(type(r) == Root and r.graph == coxeter_graph for r in sigma.values())
            nonvariable = all(r.is_constant() for r in sigma.values())
            if keys_valid and roots_valid and nonvariable:
                self.sigma = sigma.copy()
            else:
                raise InvalidInputException(self, sigma)
        else:
            self.sigma = {i: Root(coxeter_graph, i) for i in coxeter_graph.generators}
        self.graph = coxeter_graph
        self._right_descents = None
        self._minimal_right_descent = None
        self._minimal_reduced_word = None

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
        if not (value.is_constant() and type(value) == Root and value.graph == self.graph):
            raise InvalidInputException(self, (i, value), '__setitem__')
        elif i not in self.graph.generators:
            raise InvalidInputException(self, (i, value), '__setitem__')
        else:
            self.sigma[i] = value
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
        return self.right_action.right_descents

    @property
    def right_descents(self):
        return self.left_action.right_descents

    @property
    def minimal_reduced_word(self):
        """Returns lexicographically minimal reduced word for self, read left to right."""
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
            alpha = Root(self.graph, self.graph.star(i))
            if new.left_action[i] not in [alpha, -alpha]:
                new.extend_left(self.graph.star(i))
            new.extend_right(i)
        return new

    def get_reduced_words(self):
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
