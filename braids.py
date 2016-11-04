import numpy as np
import time
import itertools


def reverse_tuple(input_tuple):
    return tuple(reversed(input_tuple))


class VectorMixin:
    def is_zero(self):
        return len(self.coefficients) == 0

    def __setitem__(self, index, value):
        if value != 0:
            self.coefficients[index] = value
        elif index in self.coefficients:
            del self.coefficients[index]

    def __getitem__(self, index):
        return self.coefficients.get(index, 0)

    def __iter__(self):
        return self.coefficients.items().__iter__()

    def __contains__(self, i):
        return i in self.coefficients

    def __len__(self):
        return len(self.coefficients)

    def __eq__(self, other):
        return len(other - self) == 0

    @classmethod
    def get_index_repr(cls, index):
        """Return nice string representation of given `index`."""
        raise NotImplementedError

    def __repr__(self):
        s = ''
        for i, v in sorted(self.coefficients.items()):
            alpha = self.get_index_repr(i)
            if alpha == '' and 0 < v:
                s += ' + ' + str(v)
            elif alpha == '' and v < 0:
                s += ' - ' + str(-v)
            elif v == 1:
                s += ' + ' + alpha
            elif v == -1:
                s += ' - ' + alpha
            elif 0 < v and type(v) in [int, RationalNumber]:
                s += ' + ' + str(v) + '*' + alpha
            elif v < 0 and type(v) in [int, RationalNumber]:
                s += ' - ' + str(-v) + '*' + alpha
            else:
                s += ' + (' + str(v) + ')*' + alpha
        if s == '':
            return '0'
        elif s.startswith(' - '):
            return '-' + s[3:]
        elif s.startswith(' + '):
            return s[3:]
        else:
            return s


class NumberMixin:
    def __lt__(self, other):
        raise NotImplementedError

    def __le__(self, other):
        return self == other or self < other

    def __gt__(self, other):
        return -self < -other

    def __ge__(self, other):
        return -self <= -other

    def __neg__(self):
        return self * -1

    def __add__(self, other):
        raise NotImplementedError

    def __radd__(self, other):
        return self + other

    def __mul__(self, other):
        raise NotImplementedError

    def __rmul__(self, other):
        return self * other

    def __sub__(self, other):
        return self + other * -1

    def __rsub__(self, other):
        return -(self - other)

    def __pow__(self, exponent):
        assert type(exponent) == int and 0 < exponent
        if exponent == 1:
            return self + 0
        elif exponent % 2 == 0:
            x = self**(exponent//2)
            return x * x
        else:
            x = self**((exponent-1)//2)
            return x * x * self


class RationalNumber(NumberMixin):
    def __init__(self, p=0, q=1):
        if type(p) == RationalNumber:
            q = q * p.denominator
            p = p.numerator
        if type(q) == RationalNumber:
            p = p * q.denominator
            q = q.numerator
        if type(p) == int and type(q) == int and q != 0:
            p, q = self.reduce(p, q)
            self.numerator = p
            self.denominator = q
        else:
            raise Exception('Invalid input types to RationalNumber: %s' % str((type(p), type(q))))

    @classmethod
    def reduce(cls, p, q):
        if q < 0:
            q = q * -1
            p = p * -1
        d = np.math.gcd(p, q)
        return p//d, q//d

    @classmethod
    def to_rational_number(cls, other):
        if type(other) == RationalNumber:
            return other
        else:
            return RationalNumber(other)

    def is_integer(self):
        return self.denominator == 1

    def __hash__(self):
        if self.is_integer():
            return hash(self.numerator)
        else:
            return hash((self.numerator, self.denominator))

    def __eq__(self, other):
        if type(other) == QuadraticNumber:
            return QuadraticNumber.to_quadratic_number(self) == other
        if type(other) == Polynomial:
            return Polynomial.to_polynomial(self) == other
        else:
            other = self.to_rational_number(other)
            return self.numerator == other.numerator and self.denominator == other.denominator

    def __lt__(self, other):
        if type(other) == QuadraticNumber:
            return QuadraticNumber.to_quadratic_number(self) < other
        elif type(other) == Polynomial:
            return Polynomial.to_polynomial(self) < other
        else:
            other = self.to_rational_number(other)
            return (self - other).numerator < 0

    def __add__(self, other):
        if type(other) in [QuadraticNumber, Polynomial]:
            return other + self
        else:
            other = self.to_rational_number(other)
            n = self.numerator*other.denominator + other.numerator*self.denominator
            p, q = self.reduce(n, self.denominator*other.denominator)
            return RationalNumber(p, q)

    def __mul__(self, other):
        if type(other) in [Root, QuadraticNumber, Polynomial]:
            return other * self
        else:
            other = self.to_rational_number(other)
            p, q = self.reduce(self.numerator*other.numerator, self.denominator*other.denominator)
            return RationalNumber(p, q)

    def __truediv__(self, other):
        if type(other) == QuadraticNumber:
            return QuadraticNumber(self) / other
        elif type(other) == Polynomial and other.is_constant():
            return self / other.get_constant_part()
        else:
            return RationalNumber(self, other)

    def __rtruediv__(self, other):
        return other * RationalNumber(1, self)

    def __pow__(self, exponent):
        if exponent == 0 and self != 0:
            return RationalNumber(1)
        elif exponent == 0 and self == 0:
            raise Exception('Cannot compute indeterminate power 0^0')

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


PRIME_FACTORIZATION_CACHE = {}


class PrimeFactorization:
    def __init__(self, i=0):
        if type(i) == int:
            self.factorization = PrimeFactorization.get_prime_factorization(i)
            if i < 0:
                self.factorization[-1] = 1
            self.n = i
        else:
            raise Exception('Invalid input type to PrimeFactorization: %s' % type(i))

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
            factors = set(self.factorization.keys()).union(set(other.factorization.keys()))
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
            r *= p**(e//2)
        return r

    @classmethod
    def get_divisor_exponent(cls, n, p):
        """Return maximum positive integer e such that p^e divides n."""
        if n % p != 0:
            return 0
        if abs(n) == abs(p):
            return 1
        else:
            return 1 + cls.get_divisor_exponent(n//p, p)

    @classmethod
    def get_prime_factorization(cls, i):
        if i in PRIME_FACTORIZATION_CACHE:
            return PRIME_FACTORIZATION_CACHE[i]
        else:
            i_input = i
            factorization = {}
            N = range(2, abs(i)+1)
            while N:
                p = N[0]
                e = cls.get_divisor_exponent(i, p)
                if e != 0:
                    factorization[p] = e
                    i = i//p**e
                N = [a for a in N[1:] if a <= i and a % p != 0]
            PRIME_FACTORIZATION_CACHE[i_input] = factorization
            return factorization


class QuadraticNumber(VectorMixin, NumberMixin):
    def __init__(self, i=0):
        if type(i) == int:
            i = RationalNumber(i)
        if type(i) == RationalNumber:
            self.coefficients = {}
            if i != 0:
                self.coefficients[PrimeFactorization(1)] = i
        else:
            raise Exception('Invalid input type to QuadraticNumber: %s' % type(i))

    @classmethod
    def to_quadratic_number(cls, other):
        if type(other) == QuadraticNumber:
            return other
        else:
            return QuadraticNumber(other)

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

        if type(i) == int:
            pf = PrimeFactorization(i)
            square_free = pf.get_square_free_part()
            ans = QuadraticNumber()
            ans.coefficients[square_free] = RationalNumber(pf.get_truncated_square_root(), denom)
            return ans
        elif type(i) == QuadraticNumber:
            q = i / (3 + cls.sqrt(5))
            if q.is_rational():
                return cls.sqrt(q) * (cls.sqrt(2) + cls.sqrt(10))/2
            q = i / (7 + 3*cls.sqrt(5))
            if q.is_rational():
                return cls.sqrt(q) * (3*cls.sqrt(2) + cls.sqrt(10))/2
        else:
            raise Exception('Cannot compute square root of `%s`' % i)

    def is_rational(self):
        return all(pf.n == 1 for pf in self.coefficients)

    def __hash__(self):
        if self.is_rational():
            return hash(self.get_rational_part())
        else:
            return hash(tuple(sorted(self.coefficients.items())))

    def is_real(self):
        return all(0 < pf.n for pf in self.coefficients)

    def to_float(self):
        if self.is_real():
            ans = 0.0
            for pf, coeff in self.coefficients.items():
                ans += np.sqrt(pf.n) * coeff.numerator / coeff.denominator
            return ans
        else:
            raise Exception('Cannot convert non-real QuadraticNumber to float.')

    def decompose(self):
        positive_part, negative_part = QuadraticNumber(), QuadraticNumber()
        positive_part.coefficients = {i: v for i, v in self.coefficients.items() if 0 < v}
        negative_part.coefficients = {i: -v for i, v in self.coefficients.items() if v < 0}
        return positive_part, negative_part

    def __lt__(self, other):
        if type(other) == Polynomial:
            return Polynomial.to_polynomial(self) < other

        other = self.to_quadratic_number(other)
        diff = other - self

        if diff.is_rational():
            return 0 < diff[PrimeFactorization(1)]
        elif not diff.is_real():
            raise Exception('Cannot compare quadratic numbers have non-real difference.')
        elif all(0 < c for c in diff.coefficients.values()):
            return True
        elif all(c < 0 for c in diff.coefficients.values()):
            return False

        positive_part, negative_part = diff.decompose()
        if max(len(positive_part), len(negative_part)) > 3:
            raise Exception('Cannot determine inequality %s < %s.' % (negative_part, positive_part))
        else:
            return negative_part**2 < positive_part**2

    def __add__(self, other):
        if type(other) in [Polynomial]:
            return other + self
        else:
            other = self.to_quadratic_number(other)
            keys = set(self.coefficients.keys()).union(set(other.coefficients.keys()))
            new = QuadraticNumber()
            new.coefficients = {i: self[i] + other[i] for i in keys if (self[i] + other[i]) != 0}
            return new

    def __mul__(self, other):
        if type(other) in [Root, Polynomial]:
            return other * self
        else:
            other = self.to_quadratic_number(other)
            new = QuadraticNumber()
            for factors_1, coeff_1 in other:
                for factors_2, coeff_2 in self:
                    factors = factors_1*factors_2
                    square_free = factors.get_square_free_part()
                    coeff = coeff_1 * coeff_2 * factors.get_truncated_square_root()
                    new[square_free] += coeff
            return new

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

    def __truediv__(self, other):
        if not (type(other) in [int, RationalNumber, QuadraticNumber] and other != 0):
            raise Exception('Cannot divide Quadratic number by `%s`' % str(other))
        if self == 0:
            return QuadraticNumber(0)
        elif type(other) in [int, RationalNumber]:
            return self * RationalNumber(1, other)
        elif other.is_rational():
            return self / other.get_rational_part()
        else:
            conjugate = other.conjugate()
            return (self * conjugate) / (other * conjugate)

    def __rtruediv__(self, other):
        conjugate = self.conjugate()
        return (other * conjugate) / (self * conjugate)

    @classmethod
    def get_index_repr(cls, index):
        return (index.n != 1) * ('sqrt(' + str(index.n) + ')')


class Monomial:
    def __init__(self, exponents=None):
        if exponents is None:
            self.exponents = {}
        elif type(exponents) == str:
            e = Monomial.string_to_index(exponents)
            self.exponents = {e: 1}
        elif type(exponents) == int:
            self.exponents = {exponents: 1}
        elif type(exponents) != dict or not all(type(i) == int for i in exponents):
            raise Exception('Invalid input `%s` to Monomial.' % exponents)
        else:
            self.exponents = exponents

    @classmethod
    def string_to_index(cls, s):
        if s.isalpha() and len(s) == 1:
            return ord(s)
        else:
            raise Exception('Invalid input to Monomial.string_to_index: `%s`' % s)

    @classmethod
    def index_to_string(cls, n):
        if ord('a') <= n <= ord('z') or ord('A') <= n <= ord('Z'):
            return chr(n)
        elif n < 0:
            return 'y_' + str(-n)
        else:
            return 'x_' + str(n)

    def __eq__(self, other):
        if type(other) != Monomial:
            return other == 1 and len(self.exponents) == 0
        else:
            return self.exponents == other.exponents

    def __lt__(self, other):
        return repr(self) < repr(other)

    def __hash__(self):
        return hash(tuple(sorted(self.exponents.items())))

    def __mul__(self, other):
        if type(other) == Monomial:
            keys = set(self.exponents.keys()).union(other.exponents.keys())
            exponents = {i: self[i] + other[i] for i in keys if self[i] + other[i] != 0}
            return Monomial(exponents)
        else:
            raise Exception('Cannot multiply Monomial with `%s`.' % type(other))

    def __pow__(self, other):
        if other == 0:
            return Monomial()
        elif type(other) == int:
            return Monomial({i: self[i]*other for i in self.exponents})
        else:
            raise Exception('Cannot exponentiate Monomial by `%s`.' % type(other))

    def __getitem__(self, i):
        return self.exponents.get(i, 0)

    def __repr__(self):
        if not self.exponents:
            return '1'
        s = ''
        for i, e in sorted(self.exponents.items()):
            base = self.index_to_string(i)
            if e != 1:
                base = base + '^' + str(e)
            s += base
        return s


class Polynomial(VectorMixin, NumberMixin):
    def __init__(self, i=None):
        if i is None or i == 0:
            self.coefficients = {}
        elif (type(i) == str and len(i) == 1 and i.isalpha()) or type(i) == dict:
            self.coefficients = {Monomial(i): 1}
        elif type(i) == Monomial:
            self.coefficients = {i: 1}
        elif type(i) in [int, RationalNumber]:
            self.coefficients = {Monomial(): QuadraticNumber(i)}
        elif type(i) == QuadraticNumber:
            self.coefficients = {Monomial(): i}
        else:
            raise Exception('Invalid input type `%s` for Polynomial' % type(i))

    def __lt__(self, other):
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

    def __add__(self, other):
        other = self.to_polynomial(other)
        new = Polynomial()
        keys = set(self.coefficients.keys()).union(other.coefficients.keys())
        new.coefficients = {i: self[i] + other[i] for i in keys if self[i] + other[i] != 0}
        return new

    def __mul__(self, other):
        if type(other) == Root:
            return other * self
        else:
            other = self.to_polynomial(other)
            new = Polynomial()
            for monomial_1, coeff_1 in other:
                for monomial_2, coeff_2 in self:
                    new[monomial_1*monomial_2] += coeff_1 * coeff_2
            return new

    def __truediv__(self, other):
        if type(other) == Polynomial and other.is_constant():
            other = other[1]

        if not (type(other) in [int, RationalNumber, QuadraticNumber] and other != 0):
            raise Exception('Cannot divide Polynomial by `%s`' % str(other))
        else:
            other = QuadraticNumber.to_quadratic_number(other)
            new = Polynomial()
            for i, v in self:
                new[i] = v/other
            return new

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
                s += ' + (' + str(coeff) + ')' + str(monomial)
        if s.startswith(' - '):
            s = '-' + s[3:]
        elif s.startswith(' + '):
            s = s[3:]
        return s

    def get_variables(self):
        """Return set of integers indexing the indeterminates that appear in this polynomial."""
        return set(i for monomial in self.coefficients for i in monomial.exponents)

    def is_factorable(self):
        if len(self.get_variables()) != 1:
            return False

        x = next(iter(self.get_variables()))
        a = self[{x: 2}]
        b = self[{x: 1}]
        c = self[{}]

        if self != a*Polynomial({x: 2}) + b*Polynomial({x: 1}) + c*Polynomial({}):
            return False
        if a == 0 and b == 0:
            return False
        return True

    def get_factors(self):
        """TODO: improve this method."""
        if len(self.get_variables()) != 1:
            raise Exception('Cannot factor `%s`' % str(self))

        x = next(iter(self.get_variables()))
        a = self[{x: 2}]
        b = self[{x: 1}]
        c = self[{}]
        x = Polynomial({x: 1})

        if self != a*x**2 + b*x + c:
            raise Exception('Cannot factor `%s`' % str(self))
        # print(a, b, c, c == 0, x - b/a)
        if a != 0:
            b /= a
            c /= a
        # print(a, b, c, c == 0, x - b/a)

        if a == 0 and b == 0:
            raise Exception('Constant polynomial has no roots.')
        elif a == 0:
            assert self == x + c/b
            return [x + c/b]
        elif c == 0:
            assert self == a * x * (x + b)
            return [x, x + b]
        else:
            r = (-b + QuadraticNumber.sqrt(b**2 - 4*c)) / 2
            s = (-b - QuadraticNumber.sqrt(b**2 - 4*c)) / 2
            if r == s:
                assert self == a * (x - r) * (x - r)
                return [x - r]
            else:
                assert self == a * (x - r) * (x - s)
                return [x - r, x - s]

    def is_linear(self):
        for monomial in self.coefficients:
            if len(monomial.exponents) > 1:
                return False
            if any(e != 1 for i, e in monomial.exponents.items()):
                return False
        return True

    def is_constant(self):
        return list(self.coefficients.keys()) in [[], [Monomial()]]

    def get_constant_part(self):
        return self[Monomial()]

    def set_variable(self, variable, value):
        if type(variable) == int:
            pass
        elif type(variable) == str:
            variable = Monomial.string_to_index(variable)
        elif type(variable) == Monomial and list(variable.exponents.values()) == [1]:
            variable = next(iter(variable.exponents.keys()))
        else:
            raise Exception('Invalid input `%s` to Polynomial.set_variable.' % variable)

        new = Polynomial()
        for monomial, coeff in self:
            e = monomial[variable]
            if e != 0 and value != 0:
                exponents = {i: e for i, e in monomial.exponents.items() if i != variable}
                new += Polynomial(exponents) * coeff * (value**e)
            elif e < 0 and value == 0:
                raise Exception('Division by zero when setting variable in `%s`' % str(self))
            elif e == 0:
                new[monomial] += coeff
        return new

    def set_variables_to_zero(self, variables):
        """Input `variables` should be set of integers."""
        new = Polynomial()
        for i, v in self:
            if set(i.exponents.keys()).isdisjoint(set(variables)):
                new[i] = v
        return new

    @classmethod
    def to_polynomial(cls, other):
        if type(other) == Polynomial:
            return other
        else:
            return Polynomial(other)


class CoxeterGraph:
    def __init__(self, triples, star=None):
        """
        Input `triples` should be list of tuples (i, j, m) where i, j are indices of simple
        generators and m is order of s_i s_j. Triples with m == 1 or 2 can be omitted.
        If (i, j, m) is included then the reverse tuple (j, i, m) may also be omitted.

        Input `star` should be list of pairs (i, j) such that the involution * : S -> S
        with i^* = j and j^* = i extends to an automorphism of W. If `star` is not given,
        * is defined to be the identity map.
        """
        self.generators = sorted({t[0] for t in triples}.union({t[1] for t in triples}))

        # construct dictionary with orders m_ij of products of simple generators
        self.orders = {}
        for i, j, m in triples:
            if i != j and m != 2:
                # check that m=m_ij is a positive int or infinity
                valid_order = (type(m) == int and 3 <= m) or m == np.infty
                # check that triples does not contain inconsistent entries
                if valid_order and self.orders.get((i, j), m) == m:
                    self.orders[(i, j)] = m
                    self.orders[(j, i)] = m
                else:
                    raise Exception('Invalid input `triples = %s` to CoxeterGraph.' % triples)

        # construct dictionary with images of the *-involution 'star'
        self._star = {}
        if star:
            for i, j in star:
                if self._star.get(i, j) == j and self._star.get(j, i) == i:
                    self._star[i] = j
                    self._star[j] = i
                else:
                    raise Exception('Invalid input `star = %s` to CoxeterGraph.' % star)

        # validate that input `star` encodes automorphism
        for i in self.generators:
            for j in self.generators:
                if self.get_order(i, j) != self.get_order(self.star(i), self.star(j)):
                    raise Exception('Invalid input `star = %s` to CoxeterGraph.' % star)

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
        return self._star.get(i, i)

    def get_braid(self, i, j):
        gens = [i, j]
        m = self.get_order(i, j)
        if m == np.infty:
            raise Exception('Error in CoxeterGraph.get_braid: order m_ij must be finite.')
        else:
            return tuple(gens[t % 2] for t in range(m))

    def get_order(self, i, j):
        """Return order of s_i*s_j in W."""
        if i == j:
            return 1
        else:
            return self.orders.get((i, j), 2)

    def eval_bilinear(self, i, j):
        """Returns -cos(pi/m_ij)."""
        m = self.get_order(i, j)
        if self.is_simply_laced():
            if m == 1:
                return RationalNumber(1)
            elif m == 2:
                return RationalNumber(0)
            else:
                return RationalNumber(-1, 2)
        elif self.is_quadratic():
            if m == 1:
                return QuadraticNumber(1)
            elif m == 2:
                return QuadraticNumber(0)
            elif m == 3:
                return -QuadraticNumber(1)/2
            elif m == 4:
                return -QuadraticNumber.sqrt(2)/2
            elif m == 5:
                return -(QuadraticNumber.sqrt(5) + 1)/4
            elif m == 6:
                return -QuadraticNumber.sqrt(3)/2
            elif m == 12:
                return -(QuadraticNumber.sqrt(6) + QuadraticNumber.sqrt(2))/4
            elif m == np.infty:
                return QuadraticNumber(-1)
        else:
            raise Exception('Error in CoxeterGraph.eval_bilinear: '
                            'currently, m_ij must be 1, 2, ..., 6 or 12 or infinity.')
            #  return -np.cos(np.pi / self.get_order(i, j))

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
            entries = list(map(lambda x: (max_len - len(x))*' ' + x, entries))
        n = len(self.generators)
        s = ''
        while entries:
            s += ' '.join(entries[:n]) + '\n'
            entries = entries[n:]
        if s.endswith('\n'):
            s = s[:-1]
        s += '\n\n' + ', '.join([str(i) + '^* = ' + str(self.star(i)) for i in self.generators])
        return s

    @staticmethod
    def A(n, star=None):
        triples = [(i, i+1, 3) for i in range(1, n)]
        return CoxeterGraph(triples, star)

    @staticmethod
    def A2(n):
        star = [(i, n+1-i) for i in range(1, n+1)]
        return CoxeterGraph.A(n, star)

    @staticmethod
    def A_tilde(n):
        triples = [(i, i+1, 3) for i in range(1, n)] + [(n, 1, 3)]
        return CoxeterGraph(triples)

    @staticmethod
    def B(n):
        """
        Dynkin diagram labeling is:

            0==1--2--...--(n-1)

        """
        assert 2 <= n
        triples = [(i, i+1, 3) for i in range(1, n-1)] + [(0, 1, 4)]
        return CoxeterGraph(triples)

    @staticmethod
    def B_tilde(n):
        # TODO
        pass

    @staticmethod
    def C_tilde(n):
        # TODO
        pass

    @staticmethod
    def D(n, star=None):
        """
        Dynkin diagram labeling is:

               0
               |
            1--2--3--...--(n-1)

        """
        assert 4 <= n
        triples = [(i, i+1, 3) for i in range(1, n-1)] + [(0, 2, 3)]
        return CoxeterGraph(triples, star)

    @staticmethod
    def D2(n):
        assert 4 <= n
        star = [(0, 1)] + [(i, i) for i in range(2, n)]
        return CoxeterGraph.D(n, star)

    @staticmethod
    def D_tilde(n):
        # TODO
        pass

    @staticmethod
    def E(n, star=None):
        """
        Dynkin diagram labeling is:

                  3
                  |
            1--2--4--5--...--n

        """
        assert n in [6, 7, 8]
        triples = [(1, 2, 3), (2, 4, 3), (3, 4, 3)] + [(i, i+1, 3) for i in range(4, n)]
        return CoxeterGraph(triples, star)

    @staticmethod
    def E2(n):
        assert n == 6
        star = [(1, 6), (2, 5), (3, 3), (4, 4)]
        return CoxeterGraph.E(n, star)

    @staticmethod
    def E_tilde(n):
        # TODO
        pass

    @staticmethod
    def F(n=4, star=None):
        assert n == 4
        triples = [(1, 2, 3), (2, 3, 4), (3, 4, 3)]
        return CoxeterGraph(triples, star)

    @staticmethod
    def F2(n=4):
        assert n == 4
        star = [(1, 4), (2, 3)]
        return CoxeterGraph.F(n, star)

    @staticmethod
    def F_tilde(n):
        # TODO
        pass

    @staticmethod
    def G(n=2):
        assert n == 2
        return CoxeterGraph([(1, 2, 6)])

    @staticmethod
    def G2(n=2):
        assert n == 2
        return CoxeterGraph([(1, 2, 6)], [(1, 2)])

    @staticmethod
    def G_tilde(n):
        # TODO
        pass

    @staticmethod
    def H(n):
        """
        Dynkin diagram labeling is:

            1==2--...--n

        where edge == has label 5.
        """
        assert n in [3, 4]
        if n == 3:
            triples = [(1, 2, 5), (2, 3, 3)]
        elif n == 4:
            triples = [(1, 2, 5), (2, 3, 3), (3, 4, 3)]
        return CoxeterGraph(triples)


class Root(VectorMixin, NumberMixin):
    def __init__(self, coxeter_graph, index=None, coeff=1):
        self.graph = coxeter_graph
        if index is None or coeff == 0:
            self.coefficients = {}
        elif index in coxeter_graph.generators:
            self.coefficients = {index: coeff}
        else:
            raise Exception('Invalid `index = %s` in constuctor for Root.' % index)

    def __eq__(self, other):
        if other == 0 or type(other) == Root:
            return len(self - other) == 0
        else:
            return False

    def eval_bilinear(self, other):
        if type(other) == int:
            other = Root(self.graph, other)

        if type(other) == Root and other.graph == self.graph:
            ans = 0
            for i, u in self:
                for j, v in other:
                    ans += u * v * self.graph.eval_bilinear(i, j)
            return ans
        else:
            raise Exception('Cannot evaluate bilinear form with input `other = %s`.' % other)

    def reflect(self, index):
        if index in self.graph.generators:
            v = 2 * self.eval_bilinear(index)
            return self - Root(self.graph, index, v)
        else:
            raise Exception('Cannot reflect by root alpha_i with `i = %s`.' % index)

    def is_constant(self):
        return not any(type(v) == Polynomial and not v.is_constant() for i, v in self)

    def is_positive(self):
        return (not self.is_zero()) and all(0 <= v for _, v in self)

    def is_negative(self):
        return (not self.is_zero()) and all(v <= 0 for _, v in self)

    def is_valid(self):
        if self.is_zero():
            return False
        if any(v < 0 for _, v in self) and any(0 < v for _, v in self):
            return False
        return True

    def set_variables_to_zero(self, variables):
        new = Root(self.graph)
        for i, v in self:
            if type(v) == Polynomial:
                v = v.set_variables_to_zero(variables)
            if v != 0:
                new.coefficients[i] = v
        return new

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
            indices = set(self.coefficients.keys()).union(set(other.coefficients.keys()))
            new = Root(self.graph)
            new.coefficients = {i: self[i] + other[i] for i in indices if (self[i] + other[i]) != 0}
            return new
        else:
            raise Exception('Cannot add `%s` to Root.' % other)

    def __mul__(self, other):
        new = Root(self.graph)
        if other != 0:
            new.coefficients = {i: other*self[i] for i in self.coefficients}
        return new

    def __hash__(self):
        return hash((self.graph,) + tuple(self[i] for i in self.graph.generators))

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
            raise Exception('Invalid inputs (%s, %s) to RootTransform.' % (coxeter_graph, sigma))

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
            raise Exception(
                'Invalid inputs (%s, %s) to %s.__setitem__.' % (i, value, self.__name__)
            )

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
        if j in self.graph.generators:
            new = {}
            for i in self.sigma:
                root = Root(self.graph, i).reflect(j)
                for k, v in root:
                    new[i] = new.get(i, 0) + self.sigma[k] * v
            return self.__class__(self.graph, new)
        else:
            raise Exception('Cannot multiply %s by %s.' % (self.__name__, j))

    def __rmul__(self, j):
        if j in self.graph.generators:
            new = {}
            for i in self.sigma:
                new[i] = self.sigma[i].reflect(j)
            return self.__class__(self.graph, new)
        else:
            raise Exception('Cannot right multiply %s by %s.' % (self.__name__, j))

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
        return all(self.sigma.get(i, None) == Root(self.graph, i) for i in self.graph.generators)

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
                raise Exception('Invalid input `sigma = %s` to CoxeterTransform.' % sigma)
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
            raise Exception('Invalid value `%s` in CoxeterTransform.__setitem__.' % value)
        elif i not in self.graph.generators:
            raise Exception('Invalid index `%s` in CoxeterTransform.__setitem__.' % i)
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


class CoxeterWord:
    def __init__(self, coxeter_graph, word=()):
        self.graph = coxeter_graph
        self.word = []
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

    def to_involution_word(self):
        new = CoxeterWord(self.graph)
        for i in self.word:
            alpha = Root(self.graph, self.graph.star(i))
            if new.left_action[i] not in [alpha, -alpha]:
                new.extend_left(self.graph.star(i))
            new.extend_right(i)
        return new

    def get_reduced_words(self):
        ans = set()
        to_add = set([tuple(self.word)])
        while to_add:
            ans.update(to_add)
            next_to_add = set()
            for word in to_add:
                for i in range(len(word)-1):
                    s, t = word[i:i+2]
                    m = self.graph.get_order(s, t)
                    if word[i:i+m] == self.graph.get_braid(s, t):
                        new_word = word[:i] + self.graph.get_braid(t, s) + word[i+m:]
                        next_to_add.add(new_word)
            to_add = next_to_add.difference(ans)
        return ans

    def copy(self):
        other = CoxeterWord(self.graph)
        other.word = self.word[:]
        other.left_action = self.left_action.copy()
        other.right_action = self.right_action.copy()
        other.is_reduced = self.is_reduced
        return other

    def extend_right(self, j):
        """Replace self.word by self.word + [j] and update other fields."""
        self.word = self.word + [j]
        self.is_reduced = self.is_reduced and (j not in self.right_descents)
        self.left_action = self.left_action * j
        self.right_action = j * self.right_action

    def extend_left(self, j):
        """Replace self.word by [j] + self.word and update other fields."""
        self.word = [j] + self.word
        self.is_reduced = self.is_reduced and (j not in self.left_descents)
        self.left_action = j * self.left_action
        self.right_action = self.right_action * j

    def _check_valid(self, other):
        check_a = type(other) == CoxeterWord and other.graph == self.graph
        check_b = type(other) == int and other in self.graph.generators
        assert check_a or check_b

    def __mul__(self, other):
        self._check_valid(other)

        if type(other) == int:
            new = self.copy()
            new.extend_right(other)
        else:
            new = self.copy()
            for i in other.word:
                new.extend_right(i)
        return new

    def __rmul__(self, other):
        self._check_valid(other)

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
        return (not self.is_reduced)*'un' + 'reduced word [' + ', '.join(letters) + ']'


class ConstraintsManager:
    def __init__(self):
        self.pivots = []
        # list of linear Polynomials which must be <= 0
        self.nonpositive_constraints = set()
        # list of linear Polynomials which must be == 0
        self.zero_constraints = set()
        # list of Roots which must be != 0
        self.nonzero_constraints = set()
        # list of quadratic Polynomials which must be == 0
        self.quadratic_constraints = set()

    def __eq__(self, other):
        return \
            ConstraintsManager == type(other) and \
            self.pivots == other.pivots and \
            self.nonpositive_constraints == other.nonpositive_constraints and \
            self.zero_constraints == other.zero_constraints and \
            self.nonzero_constraints == other.nonzero_constraints and \
            self.quadratic_constraints == other.quadratic_constraints

    def copy(self):
        other = ConstraintsManager()
        other.pivots = self.pivots.copy()
        other.nonpositive_constraints = self.nonpositive_constraints.copy()
        other.zero_constraints = self.zero_constraints.copy()
        other.nonzero_constraints = self.nonzero_constraints.copy()
        other.quadratic_constraints = self.quadratic_constraints.copy()
        return other

    def add_zero_constraint(self, c):
        if type(c) == Root:
            for v in c.coefficients.values():
                self.add_zero_constraint(v)
        else:
            c = Polynomial.to_polynomial(c)
            if c.is_linear():
                self.zero_constraints.add(c)
            else:
                self.quadratic_constraints.add(c)

    def add_nonpositive_constraint(self, c):
        if type(c) == Root:
            for v in c.coefficients.values():
                self.add_nonpositive_constraint(v)
            return
        c = Polynomial.to_polynomial(c)
        if -c in self.nonpositive_constraints:
            self.add_zero_constraint(c)
        elif not any(c <= f for f in self.nonpositive_constraints):
            self.nonpositive_constraints.add(c)

    def add_nonzero_constraint(self, root):
        self.nonzero_constraints.add(root)

    def row_reduce_zero_constraints(self):
        degenerate_constraints = set()
        while self.zero_constraints:
            row = self.zero_constraints.pop()
            variables = row.get_variables()
            if not variables:
                if row != 0:
                    degenerate_constraints.add(row)
                continue

            var = Monomial({variables.pop(): 1})
            substitution = -row / row[var]
            substitution[var] = 0

            for i in range(len(self.pivots)):
                prev_var, prev_subst = self.pivots[i]
                c = prev_subst[var]
                if c != 0:
                    new_subst = prev_subst + c*substitution
                    # note that since new_subst was assigned to newly created in previous line,
                    # the following deletion does not affect any other existing polynomial.
                    new_subst[var] = 0
                    self.pivots[i] = (prev_var, new_subst)

            self.pivots.append((var, substitution))
            self.row_reduce_by_pivot(var, substitution)

        self.zero_constraints = degenerate_constraints

    def row_reduce_by_pivot(self, var, substitution):
        self.zero_constraints = {
            c for c in self.reduced_row_generator(var, substitution, self.zero_constraints)
            if c != 0
        }
        self.nonpositive_constraints = {
            c for c in self.reduced_row_generator(var, substitution, self.nonpositive_constraints)
            if not (c <= 0)
        }

    @classmethod
    def reduced_row_generator(cls, var, substitution, old_constraints):
        for constraint in old_constraints:
            c = constraint[var]
            if c != 0:
                new_constraint = constraint + c*substitution
                # note that since new_subst was assigned to newly created in previous line,
                # the following deletion does not affect any other existing polynomial.
                new_constraint[var] = 0
                yield new_constraint
            else:
                yield constraint

    def remove_vacuous_constraints(self):
        self.zero_constraints = {
            f for f in self.zero_constraints
            if not (f == 0)
        }
        self.nonpositive_constraints = {
            f for f in self.nonpositive_constraints
            if not (f <= 0) and not any(f <= g and f != g for g in self.nonpositive_constraints)
        }
        # if 10 < len(self.nonpositive_constraints):
        # self.nonpositive_constraints = [
        #     f for f in self.nonpositive_constraints
        #     if not any(f <= g and f != g for g in self.nonpositive_constraints)
        # ]
        self.nonzero_constraints = {
            r for r in self.nonzero_constraints
            if not any(v < 0 or 0 < v for v in r.coefficients.values())
        }
        self.quadratic_constraints = {
            f for f in self.quadratic_constraints if not (f == 0)
        }

    def __repr__(self):
        s = '\nconstraints:\n'
        i = 1

        def pad(j):
            return (3 - len(str(j)))*' ' + str(j)

        for var, substitution in self.pivots:
            s += '%s. %s = %s\n' % (pad(i), var, substitution)
            i += 1

        for c in self.nonpositive_constraints:
            s += '%s. 0 >= %s\n' % (pad(i), c)
            i += 1

        for c in self.zero_constraints:
            s += '%s. 0 == %s\n' % (pad(i), c)
            i += 1

        for c in self.quadratic_constraints:
            s += '%s. 0 == %s\n' % (pad(i), c)
            i += 1

        for c in self.nonzero_constraints:
            s += '%s. 0 != %s\n' % (pad(i), c)
            i += 1

        if s == '\nconstraints:\n':
            return ''
        else:
            return s

    def is_valid(self):
        if any(0 < f for _, f in self.pivots) or \
           any(0 < f for f in self.nonpositive_constraints) or \
           any(0 < f or f < 0 for f in self.zero_constraints) or \
           any(0 < f or f < 0 for f in self.quadratic_constraints) or \
           any(r == 0 for r in self.nonzero_constraints):
            return False
        return True


class BraidResolver:
    def __init__(self, coxeter_graph, s, t):
        if s in coxeter_graph.generators and t in coxeter_graph.generators:
            self.graph = coxeter_graph
            self.s = s
            self.t = t
            self.sigma = RootTransform(coxeter_graph)
            self.word_s = CoxeterWord(coxeter_graph)
            self.word_t = CoxeterWord(coxeter_graph)
            self.constraints = ConstraintsManager()
        else:
            raise Exception('Invalid inputs `s = %s`, `t = %s` in BraidResolver.' % (s, t))

    def __eq__(self, other):
        return \
            BraidResolver == type(other) and \
            self.graph == other.graph and \
            self.sigma == other.sigma and \
            self.word_s.left_action == other.word_s.left_action and \
            self.word_t.left_action == other.word_t.left_action and \
            self.constraints == other.constraints

    def __len__(self):
        return len(self.sigma)

    def __repr__(self):
        s = '\n'
        s += 'Resolver State:\n'
        s += '***************\n'
        s += 's = %s, word_s = %s\n' % (self.s, self.word_s)
        s += 't = %s, word_t = %s\n' % (self.t, self.word_t)
        s += '\n'
        s += 'sigma = %s' % self.sigma
        s += '\n'
        s += '      unconditional descent: %s\n' % str(self.get_unconditional_descent())
        s += ' strong conditional descent: %s\n' % str(self.get_strong_conditional_descent())
        s += '  weak conditional descents: %s\n' % list(self.get_weak_conditional_descents())
        s += str(self.constraints)
        s += '***************\n'
        return s

    def copy(self):
        other = BraidResolver(self.graph, self.s, self.t)
        other.sigma = self.sigma.copy()
        other.word_s = self.word_s.copy()
        other.word_t = self.word_t.copy()
        other.constraints = self.constraints.copy()
        return other

    def clear_constraints(self):
        self.constraints = ConstraintsManager()

    def get_semiorder(self, is_fixer=True):
        m = self.graph.get_order(self.s, self.t)
        if m % 2 != 0:
            return (m+1)//2
        elif is_fixer:
            return m//2 + 1
        else:
            return m//2

    def get_unconditional_descent(self):
        descents_to_avoid = self.word_s.left_descents.union(self.word_t.left_descents)
        unconditional = \
            list(self.sigma.unconditional_descents.intersection(descents_to_avoid)) + \
            list(self.sigma.unconditional_descents.difference(descents_to_avoid))
        if unconditional:
            return unconditional[0]
        return None

    def get_strong_conditional_descent(self):
        for i in self.sigma:
            root = Root(self.graph)
            for j, f in self.sigma[i]:
                if f <= 0 or any(f <= g for g in self.constraints.nonpositive_constraints):
                    root += Root(self.graph, j, f)
            if root != 0:
                return i, root
        return None

    def get_weak_conditional_descents(self):
        descents_to_avoid = self.word_s.left_descents.union(self.word_t.left_descents)
        return self.sigma.weak_conditional_descents.difference(descents_to_avoid)

    def is_quadratic_constraint_factorable(self):
        if self.constraints.quadratic_constraints:
            c = next(iter(self.constraints.quadratic_constraints))
            return c.is_factorable()
        else:
            return False

    def branch(self, should_filter=True):
        t0 = time.time()
        children, label = self._get_children()

        t1 = time.time()
        for child in children:
            child.reduce()

        t2 = time.time()
        children = [child for child in children if not should_filter or child.is_valid()]

        t3 = time.time()
        descr = '\n'
        descr += 'BRANCHING: %s\n' % label
        descr += '  CONSTRUCT: %s seconds\n' % (t1-t0)
        descr += '  REDUCTION: %s seconds\n' % (t2-t1)
        descr += '  VALIDITY : %s seconds' % (t3-t2)
        return children, descr

    def _get_children(self):
        children = []

        if len(self.sigma) == 0:
            return self._get_children_first_iteration(), 'first iteration'

        if self.is_quadratic_constraint_factorable():
            return self._get_children_from_quadratic_constraint(), 'reducing quadratic constraint'

        unconditional = self.get_unconditional_descent()
        if unconditional:
            children, iterations = self._get_children_from_unconditional_descents()
            return children, 'unconditional descent, depth=%s' % iterations

        strong = self.get_strong_conditional_descent()
        if strong:
            i, root = strong
            children = self._get_children_from_strong_conditional_descents(i, root)
            return children, 'strong conditional descent'

        weak = self.get_weak_conditional_descents()
        if weak:
            children = self._get_children_from_weak_conditional_descents(weak)
            return children, 'weak conditional descents'

        if self.sigma.is_constant() and not self.sigma.is_complete():
            return self._get_children_from_new_descent(), 'new descent'

        if self.sigma.is_constant() and not self.sigma.is_identity():
            return [], 'empty'

        raise Exception('Bad case: %s' % self)

    def _get_children_first_iteration(self):
        gens = [self.s, self.t]
        alpha = Root(self.graph, self.graph.star(self.s))
        beta = Root(self.graph, self.graph.star(self.t))

        fixer = BraidResolver(self.graph, self.s, self.t)
        fixer.sigma = RootTransform(self.graph, {self.s: alpha, self.t: beta})
        for i in range(self.get_semiorder(is_fixer=True)):
            fixer.word_s.extend_left(gens[i % 2])
            fixer.word_t.extend_left(gens[(i+1) % 2])

        transposer = BraidResolver(self.graph, self.s, self.t)
        transposer.sigma = RootTransform(self.graph, {self.s: beta, self.t: alpha})
        for i in range(self.get_semiorder(is_fixer=False)):
            transposer.word_s.extend_left(gens[i % 2])
            transposer.word_t.extend_left(gens[(i+1) % 2])

        return [fixer, transposer]

    def _get_children_from_quadratic_constraint(self):
        children = []
        constraint = next(iter(self.constraints.quadratic_constraints))
        for factor in constraint.get_factors():
            child = self.copy()
            child.constraints.add_zero_constraint(factor)
            children.append(child)
        return children

    def _get_children_from_unconditional_descents(self):
        children, current, iterations = [], [self], 0
        while current:
            new = []
            for state in current:
                state.reduce()
                if not state.is_valid():
                    continue
                descent = state.get_unconditional_descent()
                if descent:
                    new += state._get_children_from_next_unconditional_descent(descent)
                else:
                    children.append(state)
            current = new
            iterations += 1
        return children, iterations

    def _get_children_from_next_unconditional_descent(self, i):
        children = []
        if self.sigma[i].is_constant():
            commutes = self.sigma[i] == Root(self.graph, self.graph.star(i), -1)
            children.append(self._branch_from_descent(i, translate=commutes))
        else:
            children.append(self._branch_from_descent(i, True))
            children.append(self._branch_from_descent(i, False))
        return children

    def _get_children_from_strong_conditional_descents(self, i, nonzero_root):
        constraints = [(j, f) for (j, f) in nonzero_root]

        child_a = self._branch_from_descent(i, True)
        child_a.constraints.add_nonzero_constraint(nonzero_root)

        child_b = self._branch_from_descent(i, False)
        child_b.constraints.add_nonzero_constraint(nonzero_root)

        child_c = self.copy()
        for _, f in constraints:
            child_c.constraints.add_zero_constraint(f)

        return [child_a, child_b, child_c]

    def _get_children_from_weak_conditional_descents(self, descents):
        children = []
        for i in descents:
            children.append(self._branch_from_descent(i, True))
            children.append(self._branch_from_descent(i, False))
            children.append(self._branch_from_descent(i, descent=False))
        return children

    def _branch_from_descent(self, i, translate=True, descent=True):
        beta = self.sigma[i]
        new = self.copy()
        if not descent:
            new.constraints.add_nonpositive_constraint(-beta)
        else:
            new.constraints.add_nonpositive_constraint(beta)
            new.word_s.extend_left(i)
            new.word_t.extend_left(i)

            alpha = Root(new.graph, new.graph.star(i))
            if translate:
                new.constraints.add_zero_constraint(beta + alpha)
                new.sigma = new.sigma * i
            else:
                new.constraints.add_nonzero_constraint(beta + alpha)
                new.sigma = self.graph.star(i) * new.sigma * i
        return new

    def _get_children_from_new_descent(self):
        g = self.graph
        children = []

        candidates = [
            i for i in g.generators
            if i not in self.sigma and any(g.get_order(i, j) != 2 for j in self.sigma)
        ]
        for i in candidates:
            # add child with new weak descent i
            child = self.copy()
            child.clear_constraints()
            child.sigma[i] = sum([-Polynomial({j: 1})*Root(g, j) for j in g.generators])
            for j in child.sigma:
                child.constraints.add_zero_constraint(
                    child.sigma[i].eval_bilinear(child.sigma[j]) - g.eval_bilinear(i, j)
                )
            children.append(child)

        # add child with no new weak descents
        child = self.copy()
        for i in g.generators:
            if i not in child.sigma:
                child.sigma[i] = Root(g, i)
        children.append(child)
        return children

    def is_valid(self):
        if not self.are_descents_valid():
            return False
        # invalid if word_s or word_t is not reduced
        if not self.word_s.is_reduced or not self.word_t.is_reduced:
            return False
        if not self.is_sigma_valid():
            return False
        # ignore if any cannot be satisfied
        if not self.constraints.is_valid():
            return False

        return True

    def are_descents_valid(self):
        """
        Returns False if word_s and word_t share a right descent, or have (distinct)
        right descents with product of order greater than m_st. These cases may be excluded
        by induction.
        """
        if not self.word_s.right_descents.isdisjoint(self.word_t.right_descents):
            return False
        if any(self.graph.get_order(self.s, self.t) < self.graph.get_order(u, v)
           for u in self.word_s.right_descents
           for v in self.word_t.right_descents):
                return False
        return True

    def is_sigma_valid(self):
        # invalid if sigma sends any root to 0 or non-negative/positive combination of simple roots
        if any(not root.is_valid() for root in self.sigma.values()):
            return False
        # if sigma sends all roots to positive roots, and no variables remain
        if self.sigma.is_constant() and self.sigma.is_complete() and self.sigma.is_positive():
            # invalid if sigma is not trivial
            if not self.sigma.is_identity():
                return False
            # invalid if word_s and word_t are not involution words for the same element
            if not self.is_valid_involution_word_pair():
                return False
        return True

    def is_valid_involution_word_pair(self):
        x = self.word_s.to_involution_word()
        y = self.word_t.to_involution_word()
        return x.is_reduced and y.is_reduced and x.left_action == y.left_action

    def is_final(self):
        return self.is_valid() and self.sigma.is_identity()

    def reduce(self):
        self.constraints.row_reduce_zero_constraints()
        self.simplify_sigma()
        zeros = set()
        for f in self.constraints.nonpositive_constraints:
            if 0 <= f:
                zeros.update(f.get_variables())
        if zeros:
            self.set_variables_to_zero(zeros)
            self.constraints.row_reduce_zero_constraints()
            self.simplify_sigma()
        self.constraints.remove_vacuous_constraints()

    def simplify_sigma(self):
        for var, substitution in self.constraints.pivots:
            for i in self.sigma:
                self.sigma[i] = self.sigma[i].set_variable(var, substitution)

            self.constraints.nonzero_constraints = {
                r.set_variable(var, substitution) for r in self.constraints.nonzero_constraints
            }
            self.constraints.quadratic_constraints = {
                f.set_variable(var, substitution) for f in self.constraints.quadratic_constraints
            }
        self.constraints.pivots = []

    def set_variables_to_zero(self, variables):
        self.constraints.nonpositive_constraints = {
            f.set_variables_to_zero(variables) for f in self.constraints.nonpositive_constraints
        }
        self.constraints.zero_constraints = {
            f.set_variables_to_zero(variables) for f in self.constraints.zero_constraints
        }
        self.constraints.nonzero_constraints = {
            r.set_variables_to_zero(variables) for r in self.constraints.nonzero_constraints
        }
        self.constraints.quadratic_constraints = {
            f.set_variables_to_zero(variables) for f in self.constraints.quadratic_constraints
        }
        for i in self.sigma:
            self.sigma[i] = self.sigma[i].set_variables_to_zero(variables)


class ResolverQueue:

    VERBOSE_LEVEL_NONE = 0
    VERBOSE_LEVEL_MEDIUM = 1
    VERBOSE_LEVEL_HIGH = 2

    def __init__(self, coxeter_graph, s=None, t=None, verbose_level=VERBOSE_LEVEL_MEDIUM):
        self.graph = coxeter_graph
        if not s or not t:
            self.queue = [
                BraidResolver(coxeter_graph, s, t)
                for s in coxeter_graph.generators
                for t in coxeter_graph.generators
                if s < t
            ]
        else:
            self.queue = [BraidResolver(coxeter_graph, s, t)]
        self.final = set()
        self.minimal = None
        self.verbose_level = verbose_level

    def __len__(self):
        return len(self.queue)

    def _update(self, children):
        added = []
        for child in children:
            if child in self.queue:
                continue
            self._print_verbose(child)
            if child.is_final():
                self._add_final(child)
            else:
                self._insert(child)
                added.append(child)
        if not added:
            self._print_verbose('\n(no states added)\n')

    def _insert(self, child):
        i = 0
        while i < len(self.queue) and len(self.queue[i]) < len(child):
            i += 1
        self.queue.insert(i, child)

    def get_progress(self):
        d = {}
        for state in self.queue:
            l = len(state)
            d[l] = d.get(l, 0) + 1
        e = {}
        for state in self.queue:
            l = len(state.word_s)
            e[l] = e.get(l, 0) + 1
        if not d:
            return '0'
        else:
            d_str = ' '.join([str(ell) + '^' + str(mul) for ell, mul in sorted(d.items())])
            e_str = ' '.join([str(ell) + '^' + str(mul) for ell, mul in sorted(e.items())])
            return d_str + ' = ' + e_str

    def _add_final(self, child):
        u = tuple(child.word_s.word)
        v = tuple(child.word_t.word)
        if u < v:
            self.final.add((u, v))
        else:
            self.final.add((v, u))

    def get_inverse_atoms(self, relations, start_word):
        start = CoxeterTransform.from_word(self.graph, reverse_tuple(start_word))
        relations = {(reverse_tuple(a), reverse_tuple(b)) for a, b in relations}
        return start.span_by_right_relations(relations)

    def are_atoms_connected(self, relations, start_word, target_word):
        target = CoxeterTransform.from_word(self.graph, reverse_tuple(target_word))
        return target in self.get_inverse_atoms(relations, start_word)

    def finalize_relation(self, relations):
        finalized = []
        for u, v in relations:
            # find largest common initial prefix of u and v
            i = len(u)
            while 0 < i and u[:i] != v[:i]:
                i -= 1
            start_word = u[:i]
            inverse_atoms = self.get_inverse_atoms(relations, start_word)
            w = min(reverse_tuple(a.minimal_reduced_word) for a in inverse_atoms)
            finalized += [(w + u[i:], w + v[i:])]
        return finalized

    def minimize_relations(self, final):
        necessary, redundant = self._get_necessary_relations(final)
        return self._finalize_necessary_relations(final, necessary, redundant)

    def _get_necessary_relations(self, final):
        self._print("  1. Get relations which cannot be generated by all others combined.")
        necessary, redundant = [], []
        for i in range(len(final)):
            t0 = time.time()
            tag = '     * Relation %s/%s:' % (i+1, len(final))
            u, v = final[i]
            initial = [(a, b) for a, b in final[:i] if (a, b) not in redundant]
            complement = [(a, b) for a, b in (final[:i] + final[i+1:]) if (a, b) not in redundant]
            if self.are_atoms_connected(initial, u, v):
                redundant.append((u, v))
                label = "redundant"
            elif not self.are_atoms_connected(complement, u, v):
                necessary.append((u, v))
                label = "independent"
            else:
                label = "not redundant or independent"
            t1 = time.time()
            self._print("%s (%s), computation took %s seconds" % (tag, label, t1 - t0))
        return necessary, redundant

    def _finalize_necessary_relations(self, final, necessary, redundant):
        self._print("\n  2. Get smallest set of relations which generate all the others.")
        rest = [x for x in final if x not in necessary and x not in redundant]
        for i in range(len(rest)+1):
            for current in itertools.combinations(rest, i):
                candidate = necessary + list(current)
                self._print("     * Trying combination of relations:")
                for rel in candidate:
                    self._print("         %s <---> %s" % rel)
                complement = [x for x in rest if x not in current]
                if all(self.are_atoms_connected(candidate, u, v) for u, v in complement):
                    self._print("       Works!")
                    return self.finalize_relation(candidate)
                else:
                    self._print("       Does not span.")
        raise Exception('Error in `ResolverQueue._finalize_necessary_relations`: returning None.')

    def _print_verbose(self, string):
        self._print(string, self.VERBOSE_LEVEL_HIGH)

    def _print(self, string, level=VERBOSE_LEVEL_MEDIUM):
        if level <= self.verbose_level:
            print(string)

    def next(self):
        if len(self) == 0:
            self._print('Queue is empty.')
        else:
            self._print_verbose('')
            self._print_verbose('-----------')
            self._print_verbose('Next state:')
            self._print_verbose('-----------')
            self._print_verbose(self.queue[0])

            children, descr = self.queue[0].branch()
            self.queue = self.queue[1:]
            self._print(descr)

            self._print_verbose('')
            self._print_verbose('-------------')
            self._print_verbose('Child states:')
            self._print_verbose('-------------')
            self._update(children)
            self._print('States in queue: %s = %s' % (len(self), self.get_progress()))
            self._print('Final states: %s' % len(self.final))

    def go(self):
        t1 = time.time()
        while 0 < len(self):
            self.next()
        t2 = time.time()

        self._print('')
        self._print('Duration: %s seconds' % (t2-t1))

        final = sorted(self.final, key=lambda x: (len(x[0]), x))
        self._print('')
        self._print('-----------------------')
        self._print('Twisted Coxeter system:')
        self._print('-----------------------')
        self._print(self.graph)

        self._print('')
        self._print('---------------------')
        self._print('Sufficient relations:')
        self._print('---------------------')
        for u, v in final:
            self._print('%s <---> %s' % (u, v))

        self._print('')
        self._print('Finding minimal relations.')
        self._print('')
        self.minimal = self.minimize_relations(final)
        self._print('')
        self._print('Duration: %s seconds' % (time.time()-t2))

        print('')
        print('-----------------------------')
        print('Minimal sufficient relations:')
        print('-----------------------------')
        for u, v in self.minimal:
            print('%s <---> %s' % (u, v))
