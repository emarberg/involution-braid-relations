import mock
import pytest

from ..braids import RationalNumber, QuadraticNumber, Polynomial


class TestRationalNumbers:

    @pytest.mark.parametrize("p, q, numerator, denominator", [
        (1, 2, 1, 2),
        (-1, 2, -1, 2),
        (1, -2, -1, 2),
        (48, 96, 1, 2),
        (-48, -96, 1, 2),
        (RationalNumber(1, 2), 3, 1, 6),
        (RationalNumber(2, 3), RationalNumber(4, 7), 7, 6),
        (-1, RationalNumber(15, 6), -2, 5),
        (-4 * 9 * 5 * 7, -2 * 3 * 25 * 49, 2 * 3, 5 * 7)
    ])
    def test_constructor(self, p, q, numerator, denominator):
        """Test constructor for RationalNumber with a variety of valid inputs."""
        x = RationalNumber(p, q)
        assert x.numerator == numerator and x.denominator == denominator

    @pytest.mark.parametrize("p, q", [
        (1, 0),
        (1.0, 20),
        (QuadraticNumber(), 1),
        (Polynomial(), 1),
        (RationalNumber(1, 2), None),
        (None, RationalNumber(1, 2)),
        (1, RationalNumber(0, 1)),
        (RationalNumber(7, 8), RationalNumber(0, -100))
    ])
    def test_constructor_errors(self, p, q):
        """Test error handling of invalid arguments passed to QuadraticNumber constructor."""
        try:
            # denominator must be nonzero
            RationalNumber(p, q)
        except Exception as e:
            assert str(e).startswith('Invalid input types to RationalNumber: ')
        else:
            assert False

    def test_hash(self):
        """Test that hashes for RationalNumbers which are integers match the hashes of the ints."""
        N = 100
        for p in range(-N, N):
            for q in range(-N, N):
                try:
                    x = RationalNumber(p, q)
                except:
                    assert q == 0
                else:
                    if p % q == 0:
                        assert hash(x) == hash(p//q)
                        assert x.is_integer()

        assert len({1, RationalNumber(1, 1)}) == 1

    def test_addition(self):
        """Simple tests for addition and substraction of RationalNumbers."""
        a = RationalNumber(3, 7)
        b = -a
        c = RationalNumber(4, 9)

        assert a + 0 == 0 + a == a
        assert a + b == b + a == 0 == RationalNumber()
        assert a + 1 == 1 + a == RationalNumber(10, 7)
        assert a - 2 == -2 + a == RationalNumber(-11, 7)
        assert a + c == c + a == RationalNumber(3*9 + 4*7, 7*9)
        assert a - c == -c + a == RationalNumber(3*9 - 4*7, 7*9)

        # cannot add RationalNumber and float
        try:
            a + 1.2
        except:
            pass
        else:
            assert False

    def test_multiplication(self):
        """Simple tests for multiplication and division of RationalNumbers."""
        a = RationalNumber(3, 7)
        b = -a
        c = RationalNumber(4, 9)

        assert a * 1 == 1 * a == a
        assert a * b == b * a == RationalNumber(-9, 49)
        assert a * -2 == -2 * a == RationalNumber(-6, 7)
        assert a / 2 == RationalNumber(3, 14)
        assert a / c == RationalNumber(27, 28)
        assert c / a == RationalNumber(28, 27)
        assert a * c == c * a == RationalNumber(4, 21)

        # cannot multiply RationalNumber and float
        try:
            a * -3.8
        except:
            pass
        else:
            assert False

    @pytest.mark.parametrize("zero", [
        0,
        0.0,
        RationalNumber(0, 1),
        QuadraticNumber(0),
        Polynomial(),
    ])
    def test_division_errors(self, zero):
        """Test error handling for attempted division of RationalNumber by zero."""
        a = RationalNumber(3, 7)
        try:
            a / zero
        except:
            pass
        else:
            assert False

    def test_power(self):
        """Simple tests for exponentiation of RationalNumbers."""
        q = RationalNumber(3, 2)**4
        assert q.numerator == 81 and q.denominator == 16

        q = RationalNumber(3, 2)**-4
        assert q.numerator == 16 and q.denominator == 81

        q = RationalNumber(3, 2)**0
        assert q.numerator == 1 and q.denominator == 1

        # cannot compute 0**0
        try:
            RationalNumber(0, 1)**0
        except:
            pass
        else:
            assert False

        # cannot compute RationlNumber**RationalNumber, even if exponent is integer
        try:
            RationalNumber(5, 2)**RationalNumber(6, 2)
        except:
            pass
        else:
            assert False

    def test_repr(self):
        """Simple tests for conversion of RationalNumbers to str."""
        assert str(RationalNumber(1, 2)) == '1/2'
        assert str(RationalNumber(-1, 2)) == '-1/2'
        assert str(RationalNumber(2, 1)) == '2'
        assert str(RationalNumber(2, -1)) == '-2'


class TestPrimeFactorization:
    pass


class TestQuadraticNumbers:
    pass


class TestMonomial:
    pass


class TestPolynomial:
    pass


class TestComparisons:
    pass
