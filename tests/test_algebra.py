import mock
import pytest

from ..braids import RationalNumber, QuadraticNumber, Polynomial


class TestRationalNumbers:

    def test_constructor(self):
        """Test constructor for RationalNumber with a variety of valid inputs."""
        q = RationalNumber()
        assert q.numerator == 0 and q.denominator == 1

        q = RationalNumber(12)
        assert q.numerator == 12 and q.denominator == 1

        q = RationalNumber(1, 2)
        assert q.numerator == 1 and q.denominator == 2

        q = RationalNumber(-1, 2)
        assert q.numerator == -1 and q.denominator == 2

        q = RationalNumber(1, -2)
        assert q.numerator == -1 and q.denominator == 2

        q = RationalNumber(48, 96)
        assert q.numerator == 1 and q.denominator == 2

        q = RationalNumber(RationalNumber(1, 2), 3)
        assert q.numerator == 1 and q.denominator == 6

        q = RationalNumber(RationalNumber(2, 3), RationalNumber(4, 7))
        assert q.numerator == 7 and q.denominator == 6

        q = RationalNumber(-1, RationalNumber(15, 6))
        assert q.numerator == -2 and q.denominator == 5

        q = RationalNumber(-4 * 9 * 5 * 7, -2 * 3 * 25 * 49)
        assert q.numerator == 2 * 3 and q.denominator == 5 * 7

    @pytest.mark.parametrize("p, q", [
        (1, 0),
        (1.0, 20),
        (QuadraticNumber(), 1),
        (Polynomial(), 1),
        (RationalNumber(1, 2), None),
        (None, RationalNumber(1, 2)),
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

    def test_to_rational_number(self):
        RationalNumber(1, 2)

    def test_is_integer(self):
        RationalNumber(1, 2)

    def test_hash(self):
        RationalNumber(1, 2)

    def test_comparisons(self):
        RationalNumber(1, 2)

    def test_addition(self):
        RationalNumber(1, 2)

    def test_multiplication(self):
        RationalNumber(1, 2)

    def test_division(self):
        RationalNumber(1, 2)

    def test_power(self):
        RationalNumber(1, 2)

    def test_repr(self):
        RationalNumber(1, 2)


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
