
def reverse_tuple(input_tuple):
    return tuple(reversed(input_tuple))


class IndeterminatePowerException(Exception):
    def __init__(self):
        super(IndeterminatePowerException, self).__init__('Cannot compute indeterminate power 0**0')


class NegativePowerException(Exception):
    def __init__(self):
        super(NegativePowerException, self).__init__(
            '** not implemented when exponent is non-positive integer')


class InvalidInputException(Exception):
    def __init__(self, obj, inputs, method='__init__'):
        super(InvalidInputException, self).__init__(
            'Invalid inputs to %s.%s: %s' % (obj.__class__.__name__, method, str(inputs)))


class ZeroDivisionException(Exception):
    def __init__(self, n):
        super(ZeroDivisionException, self).__init__('Cannot divide %s by 0' % n.__class__.__name__)


class OperatorException(Exception):
    def __init__(self, a, b, operator='__eq__'):
        method = a.__class__.__name__ + '.' + operator
        super(OperatorException, self).__init__(
            'Cannot evaluate %s with input of type `%s`' % (method, type(b)))


class OperatorMixin:

    """
    Mixin class implementing polymorphic binary operations ==, <, +, *, and /
    between int, RationalNumber, QuadraticNumber, and Polynomial objects.
    """

    def __eq__(self, other):
        if type(other) == int:
            return self.eq__int(other)
        elif type(other).__name__ == 'RationalNumber':
            return self.eq__rational_number(other)
        elif type(other).__name__ == 'QuadraticNumber':
            return self.eq__quadratic_number(other)
        elif type(other).__name__ == 'Polynomial':
            return self.eq__polynomial(other)
        else:
            raise OperatorException(self, other)

    def eq__int(self, other):
        """Evaluates self == other under assumption that type(other) is int."""
        raise NotImplementedError  # pragma: no cover

    def eq__rational_number(self, other):
        """Evaluates self == other under assumption that type(other) is RationalNumber."""
        raise NotImplementedError  # pragma: no cover

    def eq__quadratic_number(self, other):
        """Evaluates self == other under assumption that type(other) is QuadraticNumber."""
        raise NotImplementedError  # pragma: no cover

    def eq__polynomial(self, other):
        """Evaluates self == other under assumption that type(other) is Polynomial."""
        raise NotImplementedError  # pragma: no cover

    def __lt__(self, other):
        if type(other) == int:
            return self.lt__int(other)
        elif type(other).__name__ == 'RationalNumber':
            return self.lt__rational_number(other)
        elif type(other).__name__ == 'QuadraticNumber':
            return self.lt__quadratic_number(other)
        elif type(other).__name__ == 'Polynomial':
            return self.lt__polynomial(other)
        else:
            raise OperatorException(self, other, '__lt__')

    def lt__int(self, other):
        """Evaluates self < other under assumption that type(other) is int."""
        raise NotImplementedError  # pragma: no cover

    def lt__rational_number(self, other):
        """Evaluates self <other under assumption that type(other) is RationalNumber."""
        raise NotImplementedError  # pragma: no cover

    def lt__quadratic_number(self, other):
        """Evaluates self < other under assumption that type(other) is QuadraticNumber."""
        raise NotImplementedError  # pragma: no cover

    def lt__polynomial(self, other):
        """Evaluates self < other under assumption that type(other) is Polynomial."""
        raise NotImplementedError  # pragma: no cover

    def __add__(self, other):
        if type(other) == int:
            return self.add__int(other)
        elif type(other).__name__ == 'RationalNumber':
            return self.add__rational_number(other)
        elif type(other).__name__ == 'QuadraticNumber':
            return self.add__quadratic_number(other)
        elif type(other).__name__ == 'Polynomial':
            return self.add__polynomial(other)
        else:
            raise OperatorException(self, other, '__add__')

    def add__int(self, other):
        """Evaluates self + other under assumption that type(other) is int."""
        raise NotImplementedError  # pragma: no cover

    def add__rational_number(self, other):
        """Evaluates self + other under assumption that type(other) is RationalNumber."""
        raise NotImplementedError  # pragma: no cover

    def add__quadratic_number(self, other):
        """Evaluates self + other under assumption that type(other) is QuadraticNumber."""
        raise NotImplementedError  # pragma: no cover

    def add__polynomial(self, other):
        """Evaluates self + other under assumption that type(other) is Polynomial."""
        raise NotImplementedError  # pragma: no cover

    def __mul__(self, other):
        if type(other) == int:
            return self.mul__int(other)
        elif type(other).__name__ == 'RationalNumber':
            return self.mul__rational_number(other)
        elif type(other).__name__ == 'QuadraticNumber':
            return self.mul__quadratic_number(other)
        elif type(other).__name__ == 'Polynomial':
            return self.mul__polynomial(other)
        elif type(other).__name__ == 'Root':
            return other * self
        else:
            raise OperatorException(self, other, '__mul__')

    def mul__int(self, other):
        """Evaluates self * other under assumption that type(other) is int."""
        raise NotImplementedError  # pragma: no cover

    def mul__rational_number(self, other):
        """Evaluates self * other under assumption that type(other) is RationalNumber."""
        raise NotImplementedError  # pragma: no cover

    def mul__quadratic_number(self, other):
        """Evaluates self * other under assumption that type(other) is QuadraticNumber."""
        raise NotImplementedError  # pragma: no cover

    def mul__polynomial(self, other):
        """Evaluates self * other under assumption that type(other) is Polynomial."""
        raise NotImplementedError  # pragma: no cover

    def __truediv__(self, other):
        if other == 0:
            raise ZeroDivisionException(self)
        elif type(other) == int:
            return self.truediv__int(other)
        elif type(other).__name__ == 'RationalNumber':
            return self.truediv__rational_number(other)
        elif type(other).__name__ == 'QuadraticNumber':
            return self.truediv__quadratic_number(other)
        elif type(other).__name__ == 'Polynomial':
            return self.truediv__polynomial(other)
        else:
            raise OperatorException(self, other, '__truediv__')

    def truediv__int(self, other):
        """Evaluates self * other under assumption that type(other) is int."""
        raise NotImplementedError  # pragma: no cover

    def truediv__rational_number(self, other):
        """Evaluates self * other under assumption that type(other) is RationalNumber."""
        raise NotImplementedError  # pragma: no cover

    def truediv__quadratic_number(self, other):
        """Evaluates self * other under assumption that type(other) is QuadraticNumber."""
        raise NotImplementedError  # pragma: no cover

    def truediv__polynomial(self, other):
        """Evaluates self * other under assumption that type(other) is Polynomial."""
        raise NotImplementedError  # pragma: no cover


class VectorMixin:

    """
    Mixin class implementing common magic functions for dict-like objects
    that represent elements of a vector space.
    """

    @property
    def coefficients(self):
        """This should return a dict whose values are number-like objects."""
        raise NotImplementedError  # pragma: no cover

    @coefficients.setter
    def coefficients(self, value):
        raise NotImplementedError  # pragma: no cover

    def __eq__(self, other):
        if self.is_comparable(other):
            return len(other - self) == 0
        else:
            raise OperatorException(self, other)

    def is_comparable(self, other):
        """Returns True if we can evaluate ==, etc, between self and other."""
        raise NotImplementedError  # pragma: no cover

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

    @classmethod
    def get_index_repr(cls, index):
        """Return nice string representation of given `index`."""
        raise NotImplementedError  # pragma: no cover

    @classmethod
    def is_rational_coeff(cls, v):
        return type(v) == int or (hasattr(v, 'is_rational') and v.is_rational())

    def __repr__(self):
        """Rather complicated convenience method for printing VectorMixin objects nicely."""
        s = ''
        for i, v in sorted(self.coefficients.items()):
            alpha = self.get_index_repr(i)
            if alpha == '' and 0 < v:
                s += ' + ' + str(v)
            elif alpha == '' and v < 0:
                s += ' - ' + str(-v)  # pragma: no cover
            elif v == 1:
                s += ' + ' + alpha
            elif v == -1:
                s += ' - ' + alpha
            elif 0 < v and VectorMixin.is_rational_coeff(v):
                s += ' + ' + str(v) + '*' + alpha
            elif v < 0 and VectorMixin.is_rational_coeff(v):
                s += ' - ' + str(-v) + '*' + alpha
            else:
                s += ' + (' + str(v) + ')*' + alpha
        if s == '':
            s = '0'
        elif s.startswith(' - '):
            s = '-' + s[3:]
        elif s.startswith(' + '):
            s = s[3:]

        return '(' + s + ')'


class NumberMixin:

    """
    Mixin class implementing common magic functions for number-like objects;
    technically, for objects that represent elements of a partially ordered ring.

    The __eq__ operator is not implemented here or required to be implemented,
    in contrast to most other operations (namely, < <= > >= - + * **). This
    convention ensures that when both VectorMixin and NumberMixin are inherited,
    neither class overrides any methods of the other.
    """

    def __lt__(self, other):
        raise NotImplementedError  # pragma: no cover

    def __le__(self, other):
        raise NotImplementedError  # pragma: no cover

    def __gt__(self, other):
        return -self < -other

    def __ge__(self, other):
        return -self <= -other

    def __neg__(self):
        return self * -1

    def __add__(self, other):
        raise NotImplementedError  # pragma: no cover

    def __radd__(self, other):
        return self + other

    def __mul__(self, other):
        raise NotImplementedError  # pragma: no cover

    def __rmul__(self, other):
        return self * other

    def __sub__(self, other):
        return self + (other * -1)

    def __rsub__(self, other):
        return other + (self * -1)

    @classmethod
    def one(cls):
        """Returns multiplicative identity element of ring to which NumberMixin belongs."""
        raise NotImplementedError

    def __pow__(self, exponent):
        if type(exponent) != int:
            raise OperatorException(self, exponent, '__pow__')
        elif exponent == 0 and self == 0:
            raise IndeterminatePowerException
        elif exponent == 0:
            return self.one()
        elif exponent < 0:
            try:
                return 1 / self**(-exponent)
            except:
                raise NegativePowerException
        elif exponent == 1:
            return self
        elif exponent % 2 == 0:
            x = self**(exponent // 2)
            return x * x
        else:
            x = self**((exponent - 1) // 2)
            return x * x * self
