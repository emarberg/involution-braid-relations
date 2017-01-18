import numpy as np

from project.utils import (
    IndeterminatePowerException,
    InvalidInputException,
    ZeroDivisionException,
    OperatorException,
    OperatorMixin,
    NumberMixin,
    VectorMixin
)


class RationalNumber(OperatorMixin, NumberMixin):

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
            raise OperatorException(self, other, '__truediv__')

    def __rtruediv__(self, other):
        return other * RationalNumber(1, self)

    def __pow__(self, exponent):
        if type(exponent) != int:
            raise OperatorException(self, exponent, '__pow__')
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

    def __init__(self, i=0):
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
        if other.is_constant():
            return self / other.get_constant_part()
        else:
            raise OperatorException(self, other, '__truediv__')

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

    @property
    def coefficients(self):
        return self._coefficients

    @coefficients.setter
    def coefficients(self, value):
        self._coefficients = value

    def __init__(self, i=None):
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
            return {x + c / b}
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
