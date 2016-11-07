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
