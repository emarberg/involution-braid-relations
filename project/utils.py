
def reverse_tuple(input_tuple):
    return tuple(reversed(input_tuple))


class IndeterminatePowerException(Exception):
    def __init__(self):
        super(IndeterminatePowerException, self).__init__('Cannot compute indeterminate power 0**0')


class InvalidInputException(Exception):
    def __init__(self, obj, inputs, method='__init__'):
        super(InvalidInputException, self).__init__(
            'Invalid inputs to %s.%s: %s' % (obj.__class__.__name__, method, str(inputs)))


class ZeroDivisionException(Exception):
    def __init__(self, n):
        super(ZeroDivisionException, self).__init__('Cannot divide %s by 0' % n.__class__.__name__)


class VectorMixin:

    class OperatorException(Exception):
        def __init__(self, a, b, operator='__eq__'):
            method = a.__class__.__name__ + '.' + operator
            super(VectorMixin.OperatorException, self).__init__(
                'Cannot evaluate %s with input of type `%s`' % (method, type(b)))

    def __eq__(self, other):
        if self.is_comparable(other):
            return len(other - self) == 0
        else:
            raise VectorMixin.OperatorException(self, other)

    def is_comparable(self, other):
        """Returns True if we can evaluate ==, etc, between self and other."""
        raise NotImplementedError

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
        raise NotImplementedError

    @classmethod
    def is_rational_coeff(cls, v):
        try:
            return type(v) == int or v.is_rational()
        except:
            return False

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

    class PowerException(Exception):
        def __init__(self):
            super(NumberMixin.PowerException, self).__init__(
                '** not implemented when exponent is non-positive or non-integer')

    def __lt__(self, other):
        raise NotImplementedError

    def __le__(self, other):
        raise NotImplementedError

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
        return self + (other*-1)

    def __rsub__(self, other):
        return other + (self*-1)

    def __pow__(self, exponent):
        if type(exponent) != int or exponent <= 0:
            raise NumberMixin.PowerException
        if exponent == 1:
            return self + 0
        elif exponent % 2 == 0:
            x = self**(exponent//2)
            return x * x
        else:
            x = self**((exponent-1)//2)
            return x * x * self